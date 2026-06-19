function [x, Re, u_bulk_average, info] = generate_Re_profile_ndim_multicycle_asym_raisedcos_double_hold( ...
    Re_i, Re_o, alpha, beta, t1, t3, n, Nsamp_per_T)
%--------------------------------------------------------------------------
% Generates a smooth periodic raised-cosine profile with double hold:
%
%   lower HOLD (t1) + acceleration (t2) + upper HOLD (t3) + deceleration (t4)
%
% repeated over n cycles.
%
% Inputs:
%   Re_i          minimum Reynolds number
%   Re_o          maximum Reynolds number
%   alpha         parameter that defines t2
%   beta          parameter that defines t4
%   t1            duration of the lower hold
%   t3            duration of the upper hold
%   n             number of cycles
%   Nsamp_per_T   number of samples per period
%
% Outputs:
%   x               total advective time
%   Re              Reynolds-number profile
%   u_bulk_average  normalized profile Re/<Re>
%   info            structure containing additional information
%
% Ramp definition:
%   Re_m = (Re_o + Re_i)/2
%   t2   = (Re_o - Re_i)/(alpha*Re_m)
%   t4   = (Re_o - Re_i)/(beta *Re_m)
%--------------------------------------------------------------------------

%---------------- Checks ----------------
assert(n >= 1 && floor(n) == n, ...
    'n must be an integer >= 1.');
assert(Nsamp_per_T >= 50 && floor(Nsamp_per_T) == Nsamp_per_T, ...
    'Nsamp_per_T must be an integer >= 50.');
assert(isfinite(Re_i) && isfinite(Re_o) && Re_o > Re_i, ...
    'The condition Re_o > Re_i must be satisfied.');
assert(isfinite(alpha) && alpha > 0 && isfinite(beta) && beta > 0, ...
    'alpha and beta must be > 0.');
assert(isfinite(t1) && t1 >= 0, ...
    't1 must be >= 0.');
assert(isfinite(t3) && t3 >= 0, ...
    't3 must be >= 0.');

%---------------- Durations ----------------
Re_m = 0.5 * (Re_o + Re_i);
t2   = (Re_o - Re_i) / (alpha * Re_m);
t4   = (Re_o - Re_i) / (beta  * Re_m);
T    = t1 + t2 + t3 + t4;

assert(isfinite(t2) && t2 > 0, 'Invalid t2.');
assert(isfinite(t4) && t4 > 0, 'Invalid t4.');
assert(isfinite(T)  && T  > 0, 'Invalid period T.');

%---------------- Time grid for one cycle ----------------
dt   = T / Nsamp_per_T;
tau1 = linspace(0, T, Nsamp_per_T + 1);   % includes tau = T

%---------------- Construction of one cycle ----------------
Re1 = zeros(size(tau1));

% 1) Lower HOLD: 0 <= tau <= t1
m1 = (tau1 <= t1);
Re1(m1) = Re_i;

% 2) Acceleration: t1 < tau <= t1+t2
m2 = (tau1 > t1) & (tau1 <= t1 + t2);
s  = (tau1(m2) - t1) / t2;   % s in (0,1]
Re1(m2) = Re_i + 0.5 * (Re_o - Re_i) .* (1 - cos(pi * s));

% 3) Upper HOLD: t1+t2 < tau <= t1+t2+t3
m3 = (tau1 > t1 + t2) & (tau1 <= t1 + t2 + t3);
Re1(m3) = Re_o;

% 4) Deceleration: t1+t2+t3 < tau <= T
m4 = (tau1 > t1 + t2 + t3);
sd = (tau1(m4) - (t1 + t2 + t3)) / t4;   % sd in (0,1]
Re1(m4) = Re_o - 0.5 * (Re_o - Re_i) .* (1 - cos(pi * sd));

% Exact closure for numerical robustness
Re1(1)   = Re_i;
Re1(end) = Re_i;

%---------------- Concatenate n cycles without duplicating periodic final point ----------------
Ntotal = n * Nsamp_per_T;

x  = zeros(1, Ntotal);
Re = zeros(1, Ntotal);

% First cycle: take all points except the final tau = T point
x(1:Nsamp_per_T)  = tau1(1:end-1);
Re(1:Nsamp_per_T) = Re1(1:end-1);

% Remaining cycles: also avoid duplicating the final point of each cycle
for k = 2:n
    idx0 = (k-1) * Nsamp_per_T;
    xk   = (k-1) * T + tau1;

    x(idx0+1 : idx0+Nsamp_per_T)  = xk(1:end-1);
    Re(idx0+1 : idx0+Nsamp_per_T) = Re1(1:end-1);
end

%---------------- Derived quantities ----------------
Re_average     = mean(Re);
u_bulk_average = Re / Re_average;

%---------------- Info ----------------
info = struct();
info.Re_m          = Re_m;
info.t1            = t1;
info.t2            = t2;
info.t3            = t3;
info.t4            = t4;
info.T             = T;
info.alpha         = alpha;
info.beta          = beta;
info.dt_table      = dt;
info.Nsamp_per_T   = Nsamp_per_T;
info.n_cycles      = n;
info.Re_i          = Re_i;
info.Re_o          = Re_o;
info.Re_avg        = Re_average;
info.slope_max_up   = 0.5 * (Re_o - Re_i) * pi / t2;
info.slope_max_down = 0.5 * (Re_o - Re_i) * pi / t4;
info.slope_ratio    = info.slope_max_up / info.slope_max_down;  % = t4/t2

%---------------- Prints ----------------
fprintf('\n=== Asymmetric raised-cosine profile with double hold ===\n');
fprintf('Re_i = %.6f, Re_o = %.6f, Re_m = %.6f\n', Re_i, Re_o, Re_m);
fprintf('alpha = %.6g, beta = %.6g, alpha/beta = %.6f\n', alpha, beta, alpha/beta);
fprintf('t1 = %.6f, t2 = %.6f, t3 = %.6f, t4 = %.6f, T = %.6f\n', t1, t2, t3, t4, T);
fprintf('n = %d, Nsamp_per_T = %d, dt = %.6e\n', n, Nsamp_per_T, dt);
fprintf('min(Re) = %.6f, max(Re) = %.6f, Re_avg = %.9f\n', min(Re), max(Re), Re_average);
fprintf('max slope up   = %.6f\n', info.slope_max_up);
fprintf('max slope down = %.6f\n', info.slope_max_down);
fprintf('ratio (up/down)= %.6f\n', info.slope_ratio);
fprintf('======================================================\n');

%---------------- File saving ----------------
filename = 'Profile_Data.txt';
data_to_save = [x', u_bulk_average'];
writematrix(data_to_save, filename, 'Delimiter', '\t');
fprintf('Data saved to: %s (t_adv, u_b/<u_b>)\n', filename);

filename_avg = 'Reynolds_Average.txt';
writematrix(Re_average, filename_avg, 'Delimiter', ' ');
fprintf('Re_avg saved to: %s\n', filename_avg);

end
