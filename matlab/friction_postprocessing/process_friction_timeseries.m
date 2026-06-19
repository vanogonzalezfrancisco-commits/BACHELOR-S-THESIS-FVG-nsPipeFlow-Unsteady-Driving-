function out = process_friction_timeseries(filename, varargin)
%PROCESS_FRICTION_TIMESERIES
%
% Reads the nsPipeFlow friction output, reconstructs restarts, interpolates
% to a uniform time grid, detects the deceleration branch from Ub, reflects
% Cf to negative values during that branch, removes initial cycles and
% computes the cycle-averaged friction coefficient.
%
% Expected friction file columns:
%   1 -> time
%   2 -> Ub
%   3 -> Uc
%   4 -> Ufr
%   5 -> Cf
%
% Backflow/friction-sign correction criterion:
%   1. Define global low and high levels from the clean Ub signal.
%   2. The correction window starts when Ub leaves the high branch.
%   3. The correction window ends when Ub reaches the low branch.
%
% Main output:
%   out.meanCf_signed_uniform_ss
%
% Example:
%   out = process_friction_timeseries("friction", ...
%       "Ncycles_remove", 1, ...
%       "restart_mode", "overwrite", ...
%       "high_fraction", 0.85, ...
%       "low_fraction", 0.001, ...
%       "smooth_fraction_T", 0.005, ...
%       "plot", true);

close all

%% -------------------- Parse inputs --------------------

p = inputParser;
p.addRequired("filename", @(s)ischar(s) || isstring(s));

p.addParameter("Ncycles_remove", 1, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter("dt_uniform", [], @(x)isnumeric(x) && (isempty(x) || (isscalar(x) && x>0)));
p.addParameter("plot", true, @(x)islogical(x) && isscalar(x));
p.addParameter("min_points_cycle", 50, @(x)isnumeric(x) && isscalar(x) && x>=10);
p.addParameter("interp_method", "spline", @(s)any(strcmpi(string(s), ["spline","pchip"])));
p.addParameter("restart_mode", "overwrite", @(s)any(strcmpi(string(s), ["merge","overwrite","last"])));

p.addParameter("high_fraction", 0.95, @(x)isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter("low_fraction", 0.05, @(x)isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter("range_percentiles", [1 99], @(x)isnumeric(x) && numel(x)==2 && x(1)>=0 && x(2)<=100 && x(1)<x(2));
p.addParameter("smooth_fraction_T", 0.02, @(x)isnumeric(x) && isscalar(x) && x>=0);

p.parse(filename, varargin{:});

Ncycles_remove   = p.Results.Ncycles_remove;
dt_uniform_user  = p.Results.dt_uniform;
doPlot           = p.Results.plot;
minPtsCycle      = p.Results.min_points_cycle;
interpMethod     = lower(string(p.Results.interp_method));
restartMode      = lower(string(p.Results.restart_mode));
highFraction     = p.Results.high_fraction;
lowFraction      = p.Results.low_fraction;
rangePercentiles = p.Results.range_percentiles;
smoothFractionT  = p.Results.smooth_fraction_T;

if lowFraction >= highFraction
    error("low_fraction must be smaller than high_fraction.");
end

%% -------------------- Load data --------------------

M = readmatrix(filename);

if size(M,2) < 5
    error("The friction file must have at least 5 columns: time, Ub, Uc, Ufr, Cf.");
end

t_raw = M(:,1);
X_raw = M(:,2:5);

mask = all(isfinite([t_raw, X_raw]), 2);
t_raw = t_raw(mask);
X_raw = X_raw(mask,:);

if numel(t_raw) < 5
    error("The friction signal is too short after removing NaNs/Infs.");
end

% The Cf stored by nsPipeFlow is positive because the Fortran routine uses abs().
X_raw(:,4) = abs(X_raw(:,4));

%% -------------------- Step 1: Handle restarts --------------------

[t_clean, X_clean] = clean_restarts_matrix(t_raw, X_raw, restartMode);

Ub_clean     = X_clean(:,1);
Uc_clean     = X_clean(:,2);
Ufr_clean    = X_clean(:,3);
Cf_abs_clean = abs(X_clean(:,4));

if numel(t_clean) < 5
    error("The friction signal is too short after restart reconstruction.");
end

%% -------------------- Step 2: Estimate dominant period from Ub --------------------

dt_med = median(diff(t_clean));
[T_est, f_est] = estimate_dominant_period(t_clean, Ub_clean, dt_med);

%% -------------------- Step 3: Uniform time grid --------------------

dt_u = choose_uniform_dt(dt_uniform_user, dt_med, T_est, minPtsCycle);
t_u = (t_clean(1):dt_u:t_clean(end)).';

X_u = interp1(t_clean, [Ub_clean, Uc_clean, Ufr_clean, Cf_abs_clean], t_u, interpMethod);

Ub_u     = X_u(:,1);
Uc_u     = X_u(:,2);
Ufr_u    = X_u(:,3);
Cf_abs_u = abs(X_u(:,4));

%% -------------------- Step 4: Smooth Ub for robust branch detection --------------------

if smoothFractionT > 0 && isfinite(T_est) && T_est > 0
    smooth_window = max(5, round(smoothFractionT*T_est/dt_u));
    Ub_smooth = smoothdata(Ub_u, "movmean", smooth_window);
else
    smooth_window = 1;
    Ub_smooth = Ub_u;
end

%% -------------------- Step 5: Define global thresholds --------------------

Ub_min_global = percentile_local(Ub_smooth, rangePercentiles(1));
Ub_max_global = percentile_local(Ub_smooth, rangePercentiles(2));
Ub_range_global = Ub_max_global - Ub_min_global;

if Ub_range_global <= 0 || ~isfinite(Ub_range_global)
    error("Could not define a valid global Ub range.");
end

Ub_low_global  = Ub_min_global + lowFraction  * Ub_range_global;
Ub_high_global = Ub_min_global + highFraction * Ub_range_global;

%% -------------------- Step 6: Detect deceleration/backflow windows --------------------

[detected_windows, negative_mask] = detect_deceleration_windows(t_u, Ub_smooth, Ub_low_global, Ub_high_global);

Cf_signed_u = Cf_abs_u;
Cf_signed_u(negative_mask) = -abs(Cf_signed_u(negative_mask));

%% -------------------- Step 7: Detect cycle minima for averaging window --------------------

if isfinite(T_est) && T_est > 0
    minSep_samples = max(5, round(0.6*T_est/dt_u));
else
    minSep_samples = max(5, round(0.1*numel(t_u)));
end

idx_mins = find(islocalmin(Ub_smooth, "MinSeparation", minSep_samples));

if numel(idx_mins) < 2 && isfinite(T_est) && T_est > 0
    idx_mins = detect_minima_by_period(t_u, Ub_smooth, T_est);
end

if numel(idx_mins) < 2
    error("Could not detect enough cycle minima from Ub.");
end

%% -------------------- Step 8: Remove transient cycles and keep full cycles --------------------

[t_cut, t_start_ss, t_end_ss, idx_ss, N_full_cycles] = build_steady_state_window(t_u, idx_mins, T_est, Ncycles_remove);

t_ss = t_u(idx_ss);

Cf_abs_ss = Cf_abs_u(idx_ss);
Cf_signed_ss = Cf_signed_u(idx_ss);

%% -------------------- Step 9: Averages --------------------

meanCf_abs_clean_weighted = trapz(t_clean, Cf_abs_clean) / (t_clean(end) - t_clean(1));
meanCf_abs_uniform_full = mean(Cf_abs_u, "omitnan");
meanCf_signed_uniform_full = mean(Cf_signed_u, "omitnan");

if numel(t_ss) >= 2
    meanCf_abs_uniform_ss = mean(Cf_abs_ss, "omitnan");
    meanCf_signed_uniform_ss = mean(Cf_signed_ss, "omitnan");
    meanCf_abs_uniform_ss_trapz = trapz(t_ss, Cf_abs_ss) / (t_ss(end) - t_ss(1));
    meanCf_signed_uniform_ss_trapz = trapz(t_ss, Cf_signed_ss) / (t_ss(end) - t_ss(1));
else
    meanCf_abs_uniform_ss = NaN;
    meanCf_signed_uniform_ss = NaN;
    meanCf_abs_uniform_ss_trapz = NaN;
    meanCf_signed_uniform_ss_trapz = NaN;
end

%% -------------------- Output struct --------------------

out = struct();

out.filename = string(filename);
out.restart_mode = restartMode;
out.interp_method = interpMethod;

out.t_clean = t_clean;
out.Ub_clean = Ub_clean;
out.Uc_clean = Uc_clean;
out.Ufr_clean = Ufr_clean;
out.Cf_abs_clean = Cf_abs_clean;

out.dt_uniform = dt_u;
out.t_u = t_u;
out.Ub_u = Ub_u;
out.Uc_u = Uc_u;
out.Ufr_u = Ufr_u;
out.Cf_abs_u = Cf_abs_u;
out.Cf_signed_u = Cf_signed_u;

out.Ub_smooth = Ub_smooth;
out.smooth_window = smooth_window;

out.Ub_min_global = Ub_min_global;
out.Ub_max_global = Ub_max_global;
out.Ub_low_global = Ub_low_global;
out.Ub_high_global = Ub_high_global;

out.idx_mins = idx_mins;
out.t_mins = t_u(idx_mins);

out.detected_windows = detected_windows;
out.negative_mask = negative_mask;

out.high_fraction = highFraction;
out.low_fraction = lowFraction;
out.range_percentiles = rangePercentiles;

out.f_est = f_est;
out.T_est = T_est;
out.Ncycles_remove = Ncycles_remove;

out.t_cut = t_cut;
out.t_start_ss = t_start_ss;
out.t_end_ss = t_end_ss;
out.N_full_cycles = N_full_cycles;

out.t_ss = t_ss;
out.Cf_abs_ss = Cf_abs_ss;
out.Cf_signed_ss = Cf_signed_ss;

out.meanCf_abs_clean_weighted = meanCf_abs_clean_weighted;
out.meanCf_abs_uniform_full = meanCf_abs_uniform_full;
out.meanCf_signed_uniform_full = meanCf_signed_uniform_full;

out.meanCf_abs_uniform_ss = meanCf_abs_uniform_ss;
out.meanCf_signed_uniform_ss = meanCf_signed_uniform_ss;

out.meanCf_abs_uniform_ss_trapz = meanCf_abs_uniform_ss_trapz;
out.meanCf_signed_uniform_ss_trapz = meanCf_signed_uniform_ss_trapz;

%% -------------------- Plots --------------------

if doPlot
    FS = 28;
    FSsmall = 22;

    %% Figure 1: corrected signed friction coefficient

    figure('Color','w');

    plot(t_u, Cf_abs_u, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Stored $C_f$'); hold on;
    plot(t_u, Cf_signed_u, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Corrected signed $C_f$');

    yl = ylim;

    for k = 1:size(detected_windows,1)
        patch([detected_windows(k,1) detected_windows(k,2) detected_windows(k,2) detected_windows(k,1)], ...
              [yl(1) yl(1) yl(2) yl(2)], ...
              [0.85 0.85 0.85], ...
              'FaceAlpha', 0.35, ...
              'EdgeColor', 'none', ...
              'HandleVisibility', 'off');
    end

    uistack(findobj(gca,'Type','line'), 'top');

    xlabel('$t$', 'Interpreter','latex', 'FontSize', FS);
    ylabel('$C_f$', 'Interpreter','latex', 'FontSize', FS);

    format_axes_latex(FS);

    lgd = legend('Interpreter','latex', 'Location','best');
    lgd.FontSize = FS;

    set(gcf, 'Units', 'normalized', 'OuterPosition', [0.05 0.05 0.90 0.85]);

    %% Figure 2: comparison between bulk velocity and friction correction

    scaleCf = 20;

    figure('Color','w');

    plot(t_u, Ub_u, 'k-', 'LineWidth', 1.8, 'DisplayName', '$U_b$'); hold on;
    plot(t_u, Cf_abs_u * scaleCf, 'b-', 'LineWidth', 1.8, 'DisplayName', '$C_{f,\mathrm{DNS}} \times 20$');
    plot(t_u, Cf_signed_u * scaleCf, 'r-', 'LineWidth', 1.8, 'DisplayName', '$C_{f,\mathrm{processed}} \times 20$');

    yline(Ub_high_global, 'k--', 'LineWidth', 1.4, 'HandleVisibility', 'off');
    yline(Ub_low_global, 'k--', 'LineWidth', 1.4, 'HandleVisibility', 'off');

    xlabel('$t$', 'Interpreter','latex', 'FontSize', FS);

    format_axes_latex(FS);

    xText = 5250;
    yTextOffset = -0.005;

    text(xText, Ub_high_global + yTextOffset, '$U_{b,\mathrm{high}}$', ...
        'Interpreter','latex', ...
        'FontSize', FSsmall, ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','left');

    text(xText, Ub_low_global + yTextOffset, '$U_{b,\mathrm{low}}$', ...
        'Interpreter','latex', ...
        'FontSize', FSsmall, ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','left');

    lgd = legend('Interpreter','latex', 'Location','best');
    lgd.FontSize = FS;

    set(gcf, 'Units', 'normalized', 'OuterPosition', [0.05 0.05 0.90 0.85]);

    %% Command-window output

    fprintf("\n--- FRICTION RESULTS ---\n");
    fprintf("Restart mode: %s\n", restartMode);
    fprintf("Interpolation method: %s\n", interpMethod);
    fprintf("Time interval after restart reconstruction: [%.6g, %.6g]\n", t_clean(1), t_clean(end));
    fprintf("Estimated period: T = %.6g\n", T_est);
    fprintf("Uniform dt: %.6g\n", dt_u);
    fprintf("Global Ub min percentile: %.8g\n", Ub_min_global);
    fprintf("Global Ub max percentile: %.8g\n", Ub_max_global);
    fprintf("Ub low threshold:  %.8g\n", Ub_low_global);
    fprintf("Ub high threshold: %.8g\n", Ub_high_global);
    fprintf("Detected cycle minima: %d\n", numel(idx_mins));
    fprintf("Detected %.4g%%-%.4g%% correction windows: %d\n", 100*highFraction, 100*lowFraction, size(detected_windows,1));
    fprintf("Theoretical cut: t_cut = %.6g\n", t_cut);
    fprintf("Steady-state start: t_start_ss = %.6g\n", t_start_ss);
    fprintf("Steady-state end:   t_end_ss   = %.6g\n", t_end_ss);
    fprintf("Number of full cycles in steady-state window: %g\n", N_full_cycles);
    fprintf("Scale factor for Cf in comparison figure: %.8g\n", scaleCf);
    fprintf("Mean abs Cf, uniform steady-state:    %.8g\n", meanCf_abs_uniform_ss);
    fprintf("Mean signed Cf, uniform steady-state: %.8g\n", meanCf_signed_uniform_ss);
    fprintf("Mean signed Cf, trapz steady-state:   %.8g\n\n", meanCf_signed_uniform_ss_trapz);
end

end

%% ========================================================================
% LOCAL FUNCTION: restart reconstruction for multicolumn time series
%% ========================================================================

function [t_clean, X_clean] = clean_restarts_matrix(t, X, restartMode)

restartMode = lower(string(restartMode));

switch restartMode
    case "last"
        startIdx = 1;

        for i = 2:numel(t)
            if t(i) <= t(i-1)
                startIdx = i;
            end
        end

        t_clean = t(startIdx:end);
        X_clean = X(startIdx:end,:);

        keep = [true; diff(t_clean) > 0];

        t_clean = t_clean(keep);
        X_clean = X_clean(keep,:);

    case {"merge","overwrite"}
        restartIdx = find(diff(t) <= 0) + 1;

        segStarts = [1; restartIdx(:)];
        segEnds = [restartIdx(:)-1; numel(t)];

        t_clean = [];
        X_clean = [];

        timeTolFactor = 1e-12;

        for k = 1:numel(segStarts)
            tk = t(segStarts(k):segEnds(k));
            Xk = X(segStarts(k):segEnds(k),:);

            if numel(tk) >= 2
                keepk = [true; diff(tk) > 0];
                tk = tk(keepk);
                Xk = Xk(keepk,:);
            end

            if isempty(tk)
                continue;
            end

            if isempty(t_clean)
                t_clean = tk;
                X_clean = Xk;
                continue;
            end

            switch restartMode
                case "merge"
                    tol = timeTolFactor * max(1, abs(t_clean(end)));
                    mask_new = tk > (t_clean(end) + tol);

                    t_clean = [t_clean; tk(mask_new)];
                    X_clean = [X_clean; Xk(mask_new,:)];

                case "overwrite"
                    firstNewTime = tk(1);
                    tol = timeTolFactor * max(1, abs(firstNewTime));
                    mask_keep_old = t_clean < (firstNewTime - tol);

                    t_clean = [t_clean(mask_keep_old); tk];
                    X_clean = [X_clean(mask_keep_old,:); Xk];
            end
        end

        if numel(t_clean) >= 2
            keep = [true; diff(t_clean) > 0];
            t_clean = t_clean(keep);
            X_clean = X_clean(keep,:);
        end
end

end

%% ========================================================================
% LOCAL FUNCTION: dominant-period estimation
%% ========================================================================

function [T_est, f_est] = estimate_dominant_period(t, y, dt)

t_tmp = (t(1):dt:t(end)).';
y_tmp = interp1(t, y, t_tmp, "pchip", "extrap");

y_detr = y_tmp - mean(y_tmp, "omitnan");
N = numel(y_detr);

if N < 3
    error("Not enough points to estimate the dominant period.");
end

w = 0.5*(1 - cos(2*pi*(0:N-1)'/(N-1)));
Y = fft(y_detr .* w);

P2 = abs(Y/N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

f = (0:floor(N/2))'/(N*dt);
P1(1) = 0;

[~, kmax] = max(P1);
f_est = f(kmax);

if f_est <= 0
    warning("Could not estimate a dominant period. T_est is set to NaN.");
    T_est = NaN;
else
    T_est = 1/f_est;
end

end

%% ========================================================================
% LOCAL FUNCTION: uniform time-step selection
%% ========================================================================

function dt_u = choose_uniform_dt(dt_uniform_user, dt_med, T_est, minPtsCycle)

if isempty(dt_uniform_user)
    dt_u = dt_med;
else
    dt_u = dt_uniform_user;
end

if isfinite(T_est) && T_est > 0
    dt_u = min(dt_u, T_est/minPtsCycle);
end

end

%% ========================================================================
% LOCAL FUNCTION: deceleration/backflow window detection
%% ========================================================================

function [detected_windows, negative_mask] = detect_deceleration_windows(t, Ub_smooth, Ub_low, Ub_high)

negative_mask = false(size(t));
detected_windows = [];

in_high_state = Ub_smooth(1) >= Ub_high;
in_deceleration = false;
i_start = NaN;

for i = 2:numel(t)

    if ~in_deceleration && Ub_smooth(i) >= Ub_high
        in_high_state = true;
    end

    if ~in_deceleration && in_high_state
        crossed_down_high = (Ub_smooth(i-1) >= Ub_high) && (Ub_smooth(i) < Ub_high);
        is_decreasing = Ub_smooth(i) < Ub_smooth(i-1);

        if crossed_down_high && is_decreasing
            in_deceleration = true;
            i_start = i;
        end
    end

    if in_deceleration && Ub_smooth(i) <= Ub_low
        i_end = i;

        detected_windows = [detected_windows; t(i_start), t(i_end), Ub_low, Ub_high];

        in_deceleration = false;
        in_high_state = false;
        i_start = NaN;
    end
end

if in_deceleration && ~isnan(i_start)
    detected_windows = [detected_windows; t(i_start), t(end), Ub_low, Ub_high];
end

for k = 1:size(detected_windows,1)
    negative_mask = negative_mask | (t >= detected_windows(k,1) & t <= detected_windows(k,2));
end

end

%% ========================================================================
% LOCAL FUNCTION: steady-state window from cycle minima
%% ========================================================================

function [t_cut, t_start_ss, t_end_ss, idx_ss, N_full_cycles] = build_steady_state_window(t, idx_mins, T_est, Ncycles_remove)

if isfinite(T_est) && T_est > 0
    t_cut = t(1) + Ncycles_remove*T_est;
    t_mins = t(idx_mins);

    idx_min_after_cut = find(t_mins >= t_cut, 1, "first");

    if isempty(idx_min_after_cut)
        t_start_ss = t_cut;
        t_end_ss = t(end);
        idx_ss = t >= t_start_ss;
        N_full_cycles = NaN;
    else
        t_start_ss = t_mins(idx_min_after_cut);
        idx_last_min = find(t_mins > t_start_ss, 1, "last");

        if isempty(idx_last_min)
            t_end_ss = t(end);
            idx_ss = t >= t_start_ss;
            N_full_cycles = 0;
        else
            t_end_ss = t_mins(idx_last_min);
            idx_ss = (t >= t_start_ss) & (t < t_end_ss);
            N_full_cycles = idx_last_min - idx_min_after_cut;
        end
    end
else
    t_cut = t(1);
    t_start_ss = t(1);
    t_end_ss = t(end);
    idx_ss = true(size(t));
    N_full_cycles = NaN;
end

end

%% ========================================================================
% LOCAL FUNCTION: fallback cycle-minima detection based on period estimate
%% ========================================================================

function idx_mins = detect_minima_by_period(t, y, T_est)

t = t(:);
y = y(:);

centres = t(1):T_est:t(end);
idx_mins = [];

for k = 1:numel(centres)-1
    idx_window = find(t >= centres(k) & t < centres(k+1));

    if numel(idx_window) < 5
        continue;
    end

    [~, rel_min] = min(y(idx_window));
    idx_mins = [idx_mins; idx_window(rel_min)];
end

idx_mins = unique(idx_mins);

end

%% ========================================================================
% LOCAL FUNCTION: percentile without requiring toolboxes
%% ========================================================================

function q = percentile_local(x, p)

x = sort(x(isfinite(x(:))));

if isempty(x)
    q = NaN;
    return;
end

if p <= 0
    q = x(1);
    return;
end

if p >= 100
    q = x(end);
    return;
end

pos = 1 + (numel(x)-1)*p/100;
i_low = floor(pos);
i_high = ceil(pos);

if i_low == i_high
    q = x(i_low);
else
    w = pos - i_low;
    q = (1-w)*x(i_low) + w*x(i_high);
end

end

%% ========================================================================
% LOCAL FUNCTION: common LaTeX axis formatting
%% ========================================================================

function format_axes_latex(FS)

grid on;

ax = gca;
ax.FontSize = FS;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1.2;

end
