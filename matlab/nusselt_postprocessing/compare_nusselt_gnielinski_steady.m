function out = compare_nusselt_gnielinski_steady(nusselt_file, reynolds_file, varargin)
%COMPARE_NUSSELT_GNIELINSKI_STEADY
% Compara el Nusselt de simulación con el Nusselt de Gnielinski
% usando SOLO la parte steady-state y calculando Gnielinski con una
% correlación estacionaria de fricción de Darcy en función de Re.
%
% Gnielinski:
%   Nu = [(f/8)*(Re-1000)*Pr] / [1 + 12.7*sqrt(f/8)*(Pr^(2/3)-1)]
%
% con
%   f = (0.79*log(Re) - 1.64)^(-2)    [Petukhov, Darcy]
%
% USO:
%   out = compare_nusselt_gnielinski_steady("Nusselt","Reynolds", "Ncycles_remove", 1, "restart_mode", "overwrite", "Pr", 0.7, "plot", true);
%
% PARÁMETROS:
%   "nusselt_t_col"         : columna tiempo en Nusselt   (default 1)
%   "nusselt_val_col"       : columna valor Nusselt       (default 2)
%   "re_t_col"              : columna tiempo en Reynolds  (default 1)
%   "re_val_col"            : columna valor Reynolds      (default 2)
%   "Pr"                    : Prandtl                     (default 0.7)
%   "Ncycles_remove"        : nº ciclos a eliminar        (default 1)
%   "dt_uniform"            : dt uniforme común           (default [])
%   "min_points_cycle"      : mínimo puntos/ciclo         (default 50)
%   "interp_method"         : "pchip" o "spline"          (default "pchip")
%   "restart_mode"          : "merge","overwrite","last"  (default "merge")
%   "plot"                  : true/false                  (default true)
%   "re_min_valid"          : Re mínimo validez Gnielinski(default 3000)
%   "re_max_valid"          : Re máximo validez Gnielinski(default 5e6)
%   "line_density"          : densidad rayado             (default 60)
%
% SALIDA:
%   out.t_u, out.Nu_u, out.Re_u
%   out.fDarcy_u, out.Nu_gn_u
%   out.T_est, out.f_est, out.t_cut
%   out.t_start_ss, out.t_end_ss, out.N_full_cycles
%   out.t_ss, out.Nu_sim_ss, out.Re_ss, out.Nu_gn_ss
%   out.deltaNu_ss, out.relDiff_ss
%   out.meanNu_sim_ss, out.meanNu_gn_ss
%   out.meanDelta_ss, out.meanRelDiff_ss

close all

% -------------------- Parse inputs --------------------
p = inputParser;

p.addRequired("nusselt_file",  @(s)ischar(s) || isstring(s));
p.addRequired("reynolds_file", @(s)ischar(s) || isstring(s));

p.addParameter("nusselt_t_col",   1, @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter("nusselt_val_col", 2, @(x)isnumeric(x) && isscalar(x) && x>=1);

p.addParameter("re_t_col",   1, @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter("re_val_col", 2, @(x)isnumeric(x) && isscalar(x) && x>=1);

p.addParameter("Pr", 0.7, @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter("Ncycles_remove", 1, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter("dt_uniform", [], @(x)isnumeric(x) && (isempty(x) || (isscalar(x) && x>0)));
p.addParameter("min_points_cycle", 50, @(x)isnumeric(x) && isscalar(x) && x>=10);
p.addParameter("interp_method", "pchip", @(s)any(strcmpi(string(s),["pchip","spline"])));
p.addParameter("restart_mode", "merge", @(s)any(strcmpi(string(s),["merge","overwrite","last"])));
p.addParameter("plot", true, @(x)islogical(x) && isscalar(x));

p.addParameter("re_min_valid", 3000, @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter("re_max_valid", 5e6,  @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter("line_density", 60, @(x)isnumeric(x) && isscalar(x) && x>5);

p.parse(nusselt_file, reynolds_file, varargin{:});

nusselt_file  = string(p.Results.nusselt_file);
reynolds_file = string(p.Results.reynolds_file);

Pr             = p.Results.Pr;
Ncycles_remove = p.Results.Ncycles_remove;
dt_user        = p.Results.dt_uniform;
minPtsCycle    = p.Results.min_points_cycle;
interpMethod   = lower(string(p.Results.interp_method));
restartMode    = lower(string(p.Results.restart_mode));
doPlot         = p.Results.plot;

reMinValid  = p.Results.re_min_valid;
reMaxValid  = p.Results.re_max_valid;
lineDensity = p.Results.line_density;

% -------------------- Read and clean each signal --------------------
[t_nu, Nu_clean] = read_and_clean_series( ...
    nusselt_file, p.Results.nusselt_t_col, p.Results.nusselt_val_col, restartMode);

[t_re, Re_clean] = read_and_clean_series( ...
    reynolds_file, p.Results.re_t_col, p.Results.re_val_col, restartMode);

% -------------------- Common time base --------------------
t_start = max([t_nu(1), t_re(1)]);
t_end   = min([t_nu(end), t_re(end)]);

if t_end <= t_start
    error("No hay solape temporal común entre Nusselt y Reynolds tras la limpieza.");
end

if isempty(dt_user)
    dt_nu = median(diff(t_nu));
    dt_re = median(diff(t_re));
    dt_u = min([dt_nu, dt_re]);
else
    dt_u = dt_user;
end

% -------------------- Estimate period from Nusselt --------------------
t_tmp = (t_start:dt_u:t_end).';
if numel(t_tmp) < 5
    error("La malla temporal común provisional tiene muy pocos puntos.");
end

Nu_tmp = interp1_local(t_nu, Nu_clean, t_tmp, interpMethod);

Nu_detr = Nu_tmp - mean(Nu_tmp, "omitnan");
Nfft = numel(Nu_detr);

if Nfft < 8
    error("No hay suficientes puntos para estimar el periodo con FFT.");
end

w = hann(Nfft);
Y = fft(Nu_detr .* w);
P2 = abs(Y/Nfft);
P1 = P2(1:floor(Nfft/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (0:floor(Nfft/2))' / (Nfft*dt_u);

if numel(f) < 3
    error("No hay suficiente resolución para estimar el periodo.");
end

P1(1) = 0;
[~, kmax] = max(P1);
f_est = f(kmax);

if f_est <= 0
    warning("No se pudo estimar un periodo dominante. Se usará T_est = NaN.");
    T_est = NaN;
else
    T_est = 1 / f_est;
end

% -------------------- Adjust dt to ensure enough points/cycle --------------------
if isfinite(T_est) && T_est > 0
    dt_max = T_est / minPtsCycle;
    if dt_u > dt_max
        dt_u = dt_max;
    end
end

% -------------------- Final common uniform grid --------------------
t_u = (t_start:dt_u:t_end).';
if numel(t_u) < 5
    error("La malla temporal común final tiene muy pocos puntos.");
end

Nu_u = interp1_local(t_nu, Nu_clean, t_u, interpMethod);
Re_u = interp1_local(t_re, Re_clean, t_u, interpMethod);

% -------------------- Gnielinski standard with Darcy(Petukhov) --------------------
fDarcy_u = (0.79*log(Re_u) - 1.64).^(-2);

Nu_gn_u = ((fDarcy_u./8) .* (Re_u - 1000) .* Pr) ./ ...
          (1 + 12.7 .* sqrt(fDarcy_u./8) .* (Pr^(2/3) - 1));

mask_valid_u = isfinite(Nu_u) & isfinite(Re_u) & isfinite(fDarcy_u) & ...
               (Re_u >= reMinValid) & (Re_u <= reMaxValid) & ...
               (fDarcy_u > 0);

Nu_gn_u(~mask_valid_u) = NaN;

% -------------------- Steady-state cut: min-to-min, robust --------------------
if isfinite(T_est) && T_est > 0
    t_cut = t_u(1) + Ncycles_remove * T_est;

    winSmooth = max(5, round(0.05*T_est/dt_u));
    Nu_det = smoothdata(Nu_u, "movmean", winSmooth);

    minSep_samples = max(5, round(0.6*T_est/dt_u));
    isMin = islocalmin(Nu_det, 'MinSeparation', minSep_samples);

    t_mins = t_u(isMin);

    idx_min_after_cut = find(t_mins >= t_cut, 1, 'first');

    if isempty(idx_min_after_cut)
        t_start_ss = t_cut;
        t_end_ss   = t_u(end);
        idx_ss = t_u >= t_start_ss;
        N_full_cycles = NaN;
    else
        t_start_ss = t_mins(idx_min_after_cut);

        idx_last_min = find(t_mins > t_start_ss, 1, 'last');

        if isempty(idx_last_min)
            t_end_ss = t_u(end);
            idx_ss = t_u >= t_start_ss;
            N_full_cycles = 0;
        else
            t_end_ss = t_mins(idx_last_min);
            idx_ss = (t_u >= t_start_ss) & (t_u < t_end_ss);
            N_full_cycles = idx_last_min - idx_min_after_cut;
        end
    end
else
    t_cut = t_u(1);
    t_start_ss = t_u(1);
    t_end_ss = t_u(end);
    idx_ss = true(size(t_u));
    N_full_cycles = NaN;
end

t_ss      = t_u(idx_ss);
Nu_sim_ss = Nu_u(idx_ss);
Re_ss     = Re_u(idx_ss);
fDarcy_ss = fDarcy_u(idx_ss);
Nu_gn_ss  = Nu_gn_u(idx_ss);
mask_valid_ss = mask_valid_u(idx_ss);

if numel(t_ss) < 5
    warning("La parte steady-state tiene pocos puntos.");
end

% -------------------- Comparison on steady-state only --------------------
deltaNu_ss = Nu_sim_ss - Nu_gn_ss;
relDiff_ss = 100 * deltaNu_ss ./ Nu_gn_ss;

mask_better_ss = deltaNu_ss > 0 & isfinite(deltaNu_ss);
mask_worse_ss  = deltaNu_ss < 0 & isfinite(deltaNu_ss);
mask_equal_ss  = abs(deltaNu_ss) <= 1e-12 & isfinite(deltaNu_ss);

meanNu_sim_ss = mean(Nu_sim_ss, "omitnan");
meanNu_gn_ss  = mean(Nu_gn_ss, "omitnan");

meanDelta_ss   = mean(deltaNu_ss, "omitnan");
meanRelDiff_ss = mean(relDiff_ss, "omitnan");

areaPositive_ss = trapz(t_ss(mask_better_ss), deltaNu_ss(mask_better_ss));
areaNegative_ss = trapz(t_ss(mask_worse_ss),  -deltaNu_ss(mask_worse_ss));

% -------------------- Output struct --------------------
out = struct();

out.nusselt_file  = nusselt_file;
out.reynolds_file = reynolds_file;

out.restart_mode = restartMode;
out.interp_method = interpMethod;
out.Pr = Pr;
out.Ncycles_remove = Ncycles_remove;
out.dt_uniform = dt_u;

out.t_nu = t_nu;
out.Nu_clean = Nu_clean;

out.t_re = t_re;
out.Re_clean = Re_clean;

out.f_est = f_est;
out.T_est = T_est;
out.t_cut = t_cut;
out.t_start_ss = t_start_ss;
out.t_end_ss = t_end_ss;
out.N_full_cycles = N_full_cycles;

out.t_u = t_u;
out.Nu_u = Nu_u;
out.Re_u = Re_u;
out.fDarcy_u = fDarcy_u;
out.Nu_gn_u = Nu_gn_u;
out.mask_valid_u = mask_valid_u;

out.idx_ss = idx_ss;
out.t_ss = t_ss;
out.Nu_sim_ss = Nu_sim_ss;
out.Re_ss = Re_ss;
out.fDarcy_ss = fDarcy_ss;
out.Nu_gn_ss = Nu_gn_ss;
out.mask_valid_ss = mask_valid_ss;

out.deltaNu_ss = deltaNu_ss;
out.relDiff_ss = relDiff_ss;

out.mask_better_ss = mask_better_ss;
out.mask_worse_ss = mask_worse_ss;
out.mask_equal_ss = mask_equal_ss;

out.meanNu_sim_ss = meanNu_sim_ss;
out.meanNu_gn_ss = meanNu_gn_ss;
out.meanDelta_ss = meanDelta_ss;
out.meanRelDiff_ss = meanRelDiff_ss;
out.areaPositive_ss = areaPositive_ss;
out.areaNegative_ss = areaNegative_ss;

% -------------------- Print results --------------------
fprintf("\n--- RESULTADOS COMPARACIÓN NUSSELT vs GNIELINSKI (STEADY-STATE) ---\n");
fprintf("Modo restart: %s\n", restartMode);
fprintf("Prandtl: %.6g\n", Pr);
fprintf("Periodo estimado: T = %.6g (f = %.6g)\n", T_est, f_est);
fprintf("dt uniforme común: %.6g\n", dt_u);
fprintf("Corte teórico inicial: t_cut = %.6g\n", t_cut);
fprintf("Inicio steady-state real: t_start_ss = %.6g\n", t_start_ss);
fprintf("Fin steady-state real:    t_end_ss   = %.6g\n", t_end_ss);
fprintf("Nº ciclos completos en steady-state: %g\n", N_full_cycles);
fprintf("Intervalo steady-state: t in [%.6g, %.6g] (N=%d)\n", t_ss(1), t_ss(end), numel(t_ss));
fprintf("Mean Nu_sim_ss: %.8g\n", meanNu_sim_ss);
fprintf("Mean Nu_gn_ss:  %.8g\n", meanNu_gn_ss);
fprintf("Mean(Nu_sim_ss - Nu_gn_ss): %.8g\n", meanDelta_ss);
fprintf("Mean(100*(Nu_sim_ss-Nu_gn_ss)/Nu_gn_ss): %.8g %%\n", meanRelDiff_ss);
fprintf("Área ventaja Nu_sim sobre Gnielinski: %.8g\n", areaPositive_ss);
fprintf("Área ventaja Gnielinski sobre Nu_sim: %.8g\n", areaNegative_ss);
fprintf("Fracción de tiempo Nu_sim_ss > Nu_gn_ss: %.4f\n", mean(mask_better_ss, "omitnan"));
fprintf("Fracción de tiempo Nu_sim_ss < Nu_gn_ss: %.4f\n\n", mean(mask_worse_ss, "omitnan"));

% -------------------- Plots --------------------
if doPlot
    figure;
    hold on; box on; grid on;

    hatch_between_intervals(t_ss, Nu_sim_ss, Nu_gn_ss, mask_better_ss, [0 0.5 0], lineDensity);
    hatch_between_intervals(t_ss, Nu_sim_ss, Nu_gn_ss, mask_worse_ss,  [0.8 0 0], lineDensity);

    plot(t_ss, Nu_sim_ss, 'b-', 'LineWidth', 1.4, 'DisplayName', 'Nu simulación (steady-state)');
    plot(t_ss, Nu_gn_ss, '-', 'Color', [0.85 0.4 0], 'LineWidth', 1.4, 'DisplayName', 'Nu Gnielinski');

    xlabel('t');
    ylabel('Nu');
    title('Comparación Nusselt simulación vs Gnielinski (steady-state), DOBLE HOLD (XX), \alpha = 0.0033 \beta = 1');
    legend('Location','best');

    figure;
    hold on; box on; grid on;
    yline(0,'k--');
    plot(t_ss, deltaNu_ss, 'm-', 'LineWidth', 1.2, 'DisplayName', '\DeltaNu = Nu_{sim} - Nu_{Gn}');
    xlabel('t');
    ylabel('\DeltaNu');
    title('Diferencia absoluta en steady-state');
    legend('Location','best');

    figure;
    hold on; box on; grid on;
    yline(0,'k--');
    plot(t_ss, relDiff_ss, 'c-', 'LineWidth', 1.2, ...
        'DisplayName', '100 (Nu_{sim}-Nu_{Gn}) / Nu_{Gn}');
    xlabel('t');
    ylabel('Diferencia relativa [%]');
    title('Diferencia relativa respecto a Gnielinski (steady-state)');
    legend('Location','best');
end

end

% =====================================================================
% =========================== HELPERS ==================================
% =====================================================================

function [t_clean, y_clean] = read_and_clean_series(filename, t_col, y_col, restartMode)

M = readmatrix(filename);

ncol = size(M,2);
if ncol < max(t_col, y_col)
    error("El archivo %s no tiene suficientes columnas.", filename);
end

t = M(:,t_col);
y = M(:,y_col);

mask = isfinite(t) & isfinite(y);
t = t(mask);
y = y(mask);

if numel(t) < 5
    error("Serie demasiado corta en archivo %s tras limpiar NaNs/Infs.", filename);
end

switch lower(string(restartMode))
    case "last"
        startIdx = 1;
        for i = 2:numel(t)
            if t(i) <= t(i-1)
                startIdx = i;
            end
        end

        t_clean = t(startIdx:end);
        y_clean = y(startIdx:end);

        keep = [true; diff(t_clean) > 0];
        t_clean = t_clean(keep);
        y_clean = y_clean(keep);

    case {"merge","overwrite"}
        restartIdx = find(diff(t) <= 0) + 1;
        segStarts = [1; restartIdx(:)];
        segEnds   = [restartIdx(:)-1; numel(t)];

        t_clean = [];
        y_clean = [];

        timeTolFactor = 1e-12;

        for k = 1:numel(segStarts)
            idx1 = segStarts(k);
            idx2 = segEnds(k);

            tk = t(idx1:idx2);
            yk = y(idx1:idx2);

            if numel(tk) >= 2
                keepk = [true; diff(tk) > 0];
                tk = tk(keepk);
                yk = yk(keepk);
            end

            if isempty(tk)
                continue;
            end

            if isempty(t_clean)
                t_clean = tk;
                y_clean = yk;
                continue;
            end

            tol = timeTolFactor * max(1, abs(t_clean(end)));

            switch lower(string(restartMode))
                case "merge"
                    mask_new = tk > (t_clean(end) + tol);
                    t_clean = [t_clean; tk(mask_new)];
                    y_clean = [y_clean; yk(mask_new)];

                case "overwrite"
                    firstNewTime = tk(1);
                    tol2 = timeTolFactor * max(1, abs(firstNewTime));

                    mask_keep_old = t_clean < (firstNewTime - tol2);
                    t_clean = [t_clean(mask_keep_old); tk];
                    y_clean = [y_clean(mask_keep_old); yk];
            end
        end

        if numel(t_clean) >= 2
            keep = [true; diff(t_clean) > 0];
            t_clean = t_clean(keep);
            y_clean = y_clean(keep);
        end

    otherwise
        error("restart_mode no reconocido.");
end

if numel(t_clean) < 5
    error("Tras gestionar reinicios en %s, la serie es demasiado corta.", filename);
end

end

function yq = interp1_local(t, y, tq, method)
switch lower(string(method))
    case "pchip"
        yq = interp1(t, y, tq, "pchip");
    case "spline"
        yq = interp1(t, y, tq, "spline");
    otherwise
        error("Método de interpolación no reconocido.");
end
end

function hatch_between_intervals(t, y1, y2, mask, hatchColor, lineDensity)

if ~any(mask)
    return;
end

segments = logical_segments(mask);

for s = 1:size(segments,1)
    i1 = segments(s,1);
    i2 = segments(s,2);

    tt  = t(i1:i2);
    yy1 = y1(i1:i2);
    yy2 = y2(i1:i2);

    finiteMask = isfinite(tt) & isfinite(yy1) & isfinite(yy2);
    tt  = tt(finiteMask);
    yy1 = yy1(finiteMask);
    yy2 = yy2(finiteMask);

    if numel(tt) < 2
        continue;
    end

    xp = [tt; flipud(tt)];
    yp = [yy1; flipud(yy2)];

    patch('XData', xp, 'YData', yp, ...
          'FaceColor', hatchColor, ...
          'FaceAlpha', 0.12, ...
          'EdgeColor', 'none', ...
          'HandleVisibility', 'off');

    xmin = tt(1);
    xmax = tt(end);
    ymin = min([yy1; yy2]);
    ymax = max([yy1; yy2]);

    if xmax <= xmin || ymax <= ymin
        continue;
    end

    nLines = max(8, round((xmax - xmin) / (t(end)-t(1)) * lineDensity));
    offsets = linspace(ymin-(ymax-ymin), ymax, nLines);

    for c = 1:numel(offsets)
        xline_ = [xmin, xmax];
        yline_ = [offsets(c), offsets(c) + (ymax-ymin)];

        [xc, yc] = clip_line_to_band(xline_, yline_, tt, yy1, yy2);

        if numel(xc) >= 2
            plot(xc, yc, '-', 'Color', hatchColor, 'LineWidth', 0.5, ...
                 'HandleVisibility', 'off');
        end
    end
end

end

function segs = logical_segments(mask)
d = diff([false; mask(:); false]);
starts = find(d == 1);
ends   = find(d == -1) - 1;
segs = [starts, ends];
end

function [xc, yc] = clip_line_to_band(xline_, yline_, tt, y1, y2)

ns = 300;
xs = linspace(xline_(1), xline_(2), ns).';
ys = linspace(yline_(1), yline_(2), ns).';

if xs(1) < tt(1) || xs(end) > tt(end)
    validx = xs >= tt(1) & xs <= tt(end);
    xs = xs(validx);
    ys = ys(validx);
    if numel(xs) < 2
        xc = [];
        yc = [];
        return;
    end
end

yu = interp1(tt, y1, xs, 'linear');
yl = interp1(tt, y2, xs, 'linear');

ymax = max(yu, yl);
ymin = min(yu, yl);

inside = ys >= ymin & ys <= ymax;

if ~any(inside)
    xc = [];
    yc = [];
    return;
end

segs = logical_segments(inside);
[~, idx] = max(segs(:,2) - segs(:,1));
ii = segs(idx,1):segs(idx,2);

xc = xs(ii);
yc = ys(ii);
end