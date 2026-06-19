function out = process_nusselt_timeseries(filename, varargin)
% PROCESS_NUSSELT_TIMESERIES
%
% Reads a temporal Nusselt-number signal, reconstructs possible restarts,
% interpolates the signal onto a uniform time grid, removes the prescribed
% number of initial cycles and computes cycle-averaged Nusselt numbers.
%
% The steady-state averaging window is built from minimum to minimum in
% order to average complete physical cycles.
%
% Example:
%   out = process_nusselt_timeseries("Nusselt", ...
%       "Ncycles_remove", 1, ...
%       "restart_mode", "overwrite", ...
%       "dt_uniform", [], ...
%       "plot", true);
%
% Name-value inputs:
%   "Ncycles_remove"   Number of initial cycles to remove. Default: 1
%   "dt_uniform"       Uniform time step. If empty, median(diff(t_clean)) is used.
%   "plot"             If true, diagnostic figures are generated. Default: true
%   "min_points_cycle" Minimum points per cycle. Default: 50
%   "spline_method"    Interpolation method: "spline" or "pchip". Default: "spline"
%   "restart_mode"     Restart treatment: "merge", "overwrite" or "last". Default: "merge"
%
% Main output:
%   out.meanNu_uniform_ss

close all

%% -------------------- Parse inputs --------------------

p = inputParser;
p.addRequired("filename", @(s)ischar(s) || isstring(s));
p.addParameter("Ncycles_remove", 1, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter("dt_uniform", [], @(x)isnumeric(x) && (isempty(x) || (isscalar(x) && x>0)));
p.addParameter("plot", true, @(x)islogical(x) && isscalar(x));
p.addParameter("min_points_cycle", 50, @(x)isnumeric(x) && isscalar(x) && x>=10);
p.addParameter("spline_method", "spline", @(s)any(strcmpi(string(s), ["spline","pchip"])));
p.addParameter("restart_mode", "merge", @(s)any(strcmpi(string(s), ["merge","overwrite","last"])));
p.parse(filename, varargin{:});

Ncycles_remove  = p.Results.Ncycles_remove;
dt_uniform_user = p.Results.dt_uniform;
doPlot          = p.Results.plot;
minPtsCycle     = p.Results.min_points_cycle;
interpMethod    = lower(string(p.Results.spline_method));
restartMode     = lower(string(p.Results.restart_mode));

%% -------------------- Load data --------------------

M = readmatrix(filename);

if size(M,2) < 2
    error("The input file must have at least two columns: time and Nusselt number.");
end

t = M(:,1);
Nu = M(:,2);

mask = isfinite(t) & isfinite(Nu);
t = t(mask);
Nu = Nu(mask);

if numel(t) < 5
    error("The Nusselt signal is too short after removing NaNs/Infs.");
end

%% -------------------- Step 1: Handle restarts --------------------

[t_clean, Nu_clean] = clean_restarts_signal(t, Nu, restartMode);

if numel(t_clean) < 5
    error("The Nusselt signal is too short after restart reconstruction.");
end

%% -------------------- Step 2: Estimate dominant period --------------------

dt_med = median(diff(t_clean));
[T_est, f_est] = estimate_dominant_period(t_clean, Nu_clean, dt_med);

%% -------------------- Step 3: Uniform time grid --------------------

dt_u = choose_uniform_dt(dt_uniform_user, dt_med, T_est, minPtsCycle);

t_u = (t_clean(1):dt_u:t_clean(end)).';
Nu_u = interp1(t_clean, Nu_clean, t_u, char(interpMethod));

%% -------------------- Step 4: Remove transient cycles and keep full cycles --------------------

[t_cut, t_start_ss, t_end_ss, idx_ss, N_full_cycles] = ...
    build_steady_state_window(t_u, Nu_u, T_est, dt_u, Ncycles_remove);

t_ss = t_u(idx_ss);
Nu_ss = Nu_u(idx_ss);

%% -------------------- Step 5: Averages --------------------

meanNu_clean_weighted = trapz(t_clean, Nu_clean) / (t_clean(end) - t_clean(1));
meanNu_uniform = mean(Nu_u, "omitnan");

if numel(t_ss) >= 2
    meanNu_uniform_ss = mean(Nu_ss, "omitnan");
else
    meanNu_uniform_ss = NaN;
end

%% -------------------- Output struct --------------------

out = struct();

out.filename = string(filename);
out.restart_mode = restartMode;

out.t_clean = t_clean;
out.Nu_clean = Nu_clean;

out.dt_uniform = dt_u;
out.t_u = t_u;
out.Nu_u = Nu_u;

out.f_est = f_est;
out.T_est = T_est;
out.Ncycles_remove = Ncycles_remove;

out.t_cut = t_cut;
out.t_start_ss = t_start_ss;
out.t_end_ss = t_end_ss;
out.N_full_cycles = N_full_cycles;

out.t_ss = t_ss;
out.Nu_ss = Nu_ss;

out.meanNu_clean_weighted = meanNu_clean_weighted;
out.meanNu_uniform = meanNu_uniform;
out.meanNu_uniform_ss = meanNu_uniform_ss;

%% -------------------- Plots --------------------

if doPlot
    figure;
    plot(t, Nu, '-', 'DisplayName', 'Raw'); hold on;
    plot(t_clean, Nu_clean, '-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Clean (%s)', restartMode));

    xlabel('t');
    ylabel('Nu');
    grid on;
    legend('Location','best');
    title('Raw and reconstructed Nusselt signal');

    figure;
    plot(t_clean, Nu_clean, '.', 'DisplayName', 'Clean (non-uniform)'); hold on;
    plot(t_u, Nu_u, '-', 'DisplayName', ...
        sprintf('Uniform interp (%s), dt=%.3g', interpMethod, dt_u));

    if isfinite(T_est)
        xline(t_cut, '--', sprintf('cut @ %.3g', t_cut), 'HandleVisibility', 'off');
        xline(t_start_ss, '--', sprintf('start ss @ %.3g', t_start_ss), 'HandleVisibility', 'off');
        xline(t_end_ss, '--', sprintf('end ss @ %.3g', t_end_ss), 'HandleVisibility', 'off');
    end

    xlabel('t');
    ylabel('Nu');
    grid on;
    legend('Location','best');
    title('Uniform interpolation and steady-state averaging window');

    fprintf("\n--- NUSSELT RESULTS ---\n");
    fprintf("Restart mode: %s\n", restartMode);
    fprintf("Time interval after restart reconstruction: [%.6g, %.6g] (N=%d)\n", ...
        t_clean(1), t_clean(end), numel(t_clean));
    fprintf("Estimated period: T = %.6g (f = %.6g)\n", T_est, f_est);
    fprintf("Uniform dt: %.6g (N=%d)\n", dt_u, numel(t_u));
    fprintf("Theoretical cut: t_cut = %.6g\n", t_cut);
    fprintf("Steady-state start: t_start_ss = %.6g\n", t_start_ss);
    fprintf("Steady-state end:   t_end_ss   = %.6g\n", t_end_ss);
    fprintf("Number of full cycles in steady-state window: %g\n", N_full_cycles);
    fprintf("Number of points in steady-state window: %d\n", numel(t_ss));
    fprintf("Mean Nu, clean weighted trapz: %.8g\n", meanNu_clean_weighted);
    fprintf("Mean Nu, uniform full signal:  %.8g\n", meanNu_uniform);
    fprintf("Mean Nu, uniform steady-state: %.8g\n\n", meanNu_uniform_ss);
end

end

%% ========================================================================
% LOCAL FUNCTION: restart reconstruction for one scalar time series
%% ========================================================================

function [t_clean, y_clean] = clean_restarts_signal(t, y, restartMode)

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
        y_clean = y(startIdx:end);

        keep = [true; diff(t_clean) > 0];

        t_clean = t_clean(keep);
        y_clean = y_clean(keep);

    case {"merge","overwrite"}
        restartIdx = find(diff(t) <= 0) + 1;

        segStarts = [1; restartIdx(:)];
        segEnds = [restartIdx(:)-1; numel(t)];

        t_clean = [];
        y_clean = [];

        timeTolFactor = 1e-12;

        for k = 1:numel(segStarts)
            tk = t(segStarts(k):segEnds(k));
            yk = y(segStarts(k):segEnds(k));

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

            switch restartMode
                case "merge"
                    tol = timeTolFactor * max(1, abs(t_clean(end)));
                    mask_new = tk > (t_clean(end) + tol);

                    t_clean = [t_clean; tk(mask_new)];
                    y_clean = [y_clean; yk(mask_new)];

                case "overwrite"
                    firstNewTime = tk(1);
                    tol = timeTolFactor * max(1, abs(firstNewTime));
                    mask_keep_old = t_clean < (firstNewTime - tol);

                    t_clean = [t_clean(mask_keep_old); tk];
                    y_clean = [y_clean(mask_keep_old); yk];
            end
        end

        if numel(t_clean) >= 2
            keep = [true; diff(t_clean) > 0];
            t_clean = t_clean(keep);
            y_clean = y_clean(keep);
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

w = hann(N);
Y = fft(y_detr .* w);

P2 = abs(Y/N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

f = (0:floor(N/2))'/(N*dt);

if numel(f) < 3
    error("Not enough frequency resolution to estimate the dominant period.");
end

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
% LOCAL FUNCTION: steady-state window from cycle minima
%% ========================================================================

function [t_cut, t_start_ss, t_end_ss, idx_ss, N_full_cycles] = ...
    build_steady_state_window(t, y, T_est, dt, Ncycles_remove)

if isfinite(T_est) && T_est > 0
    t_cut = t(1) + Ncycles_remove*T_est;

    y_det = smoothdata(y, "movmean", max(5, round(0.05*T_est/dt)));
    minSep_samples = max(5, round(0.6*T_est/dt));

    t_mins = t(islocalmin(y_det, "MinSeparation", minSep_samples));

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
