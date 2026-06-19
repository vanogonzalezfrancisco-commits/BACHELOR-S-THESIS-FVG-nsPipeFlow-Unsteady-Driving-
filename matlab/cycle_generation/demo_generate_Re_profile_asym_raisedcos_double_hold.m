function demo_generate_Re_profile_asym_raisedcos_double_hold()
% DEMO_GENERATE_RE_PROFILE_ASYM_RAISEDCOS_DOBLE_HOLD
%
% Demonstration script for generating an asymmetric Reynolds-number profile
% using raised-cosine transitions and optional lower and upper holding stages.
%
% This script defines the cycle parameters, calls
% generate_Re_profile_ndim_multicycle_asym_raisedcos_doble_hold, and produces
% diagnostic plots of the Reynolds-number profile and the normalized bulk
% velocity forcing signal.
%
% The generated forcing signal is written by the profile-generation function
% as Profile_Data.txt, which contains the non-dimensional advective time and
% the normalized bulk velocity u_b/<u_b> required by nsPipeFlow.
%
% The script also writes Resumen_Caso.txt with the main numerical parameters
% of the selected cycle.

clear; clc; close all;

%% ---------- DEFAULT PARAMETERS ----------
Re_i  = 3200;       
Re_o  = 9600;
alpha = 0.01;        % t2 = 1/alpha
beta  = 0.15;        % t4 = 1/beta
t1    = 20;          % Lower hold duration
t3    = 5;           % Upper hold duration
n     = 15;          % Number of cycles

Nsamp_per_T = 8000;  % Points per cycle for long T, increase to 12000
save_figs = false;

%% ---------- FUNCTION CALL ----------
[x, Re, u_bulk_avg, info] = generate_Re_profile_ndim_multicycle_asym_raisedcos_double_hold( ...
    Re_i, Re_o, alpha, beta, t1, t3, n, Nsamp_per_T);

T   = info.T;
t2  = info.t2;
t3  = info.t3;
t4  = info.t4;
Re_avg = mean(Re);

N = numel(x);

%% ---------- SUMMARY PRINTING ----------
fprintf('\n=== Case summary ===\n');
fprintf('Re_i = %.1f, Re_o = %.1f, Re_m = %.1f\n', Re_i, Re_o, info.Re_m);
fprintf('alpha = %.4f (d(Re/Re_m)/dt_adv), beta = %.4f (d(Re/Re_m)/dt_adv)\n', alpha, beta);
fprintf('t1 = %.4f, t2 = %.4f, t3 = %.4f, t4 = %.4f   -->  T = %.4f (t_adv)\n', ...
    info.t1, info.t2, info.t3, info.t4, info.T);
fprintf('n = %d cycles, N = %d samples (approximately %d per period)\n', n, N, round(N/n));

filename_val = 'Resumen_Caso.txt';
fileID = fopen(filename_val, 'w');
fprintf(fileID, '=== Numerical values of the case summary ===\n');

fprintf(fileID, 'Re_i : %.1f\n', Re_i);
fprintf(fileID, 'Re_o : %.1f\n', Re_o);
fprintf(fileID, 'Re_m : %.1f\n', info.Re_m);
fprintf(fileID, 'Re_avg : %.11f\n', Re_avg);
fprintf(fileID, '\n');
fprintf(fileID, 'alpha : %.4f\n', alpha);
fprintf(fileID, 'beta : %.4f\n', beta);
fprintf(fileID, '\n');
fprintf(fileID, 't1 : %.4f\n', info.t1);
fprintf(fileID, 't2 : %.4f\n', info.t2);
fprintf(fileID, 't3 : %.4f\n', info.t3);
fprintf(fileID, 't4 : %.4f\n', info.t4);
fprintf(fileID, 'T : %.4f\n', info.T);
fprintf(fileID, '\n');
fprintf(fileID, 'n_cycles : %d\n', n);
fprintf(fileID, 'N_samples : %d\n', N);
fprintf(fileID, 'Samples_per_period : %d\n', Nsamp_per_T);

fclose(fileID);
fprintf('Case summary saved to: %s\n', filename_val);

%% ---------- PLOTS ----------
idx1 = (x >= 0) & (x <= T + 1e-12);
x1   = x(idx1);
Re1  = Re(idx1);
u1   = u_bulk_avg(idx1);

% Correct phase markers for double hold
vmarks = [t1, t1+t2, t1+t2+t3, T];

%% ---------- Figure 1: Re(t) + ideal profiles ----------
figure('Name','Re(t_adv) - One period','Color','w');
plot(x1, Re1, 'LineWidth', 0.9, 'DisplayName', 'Re(t)'); hold on;

Re_m = info.Re_m;

t_acc = linspace(t1, t1+t2, 200);
Re_acc_ideal = Re_i + alpha*Re_m*(t_acc - t1);

t_dec = linspace(t1+t2+t3, T, 200);
Re_dec_ideal = Re_o - beta*Re_m*(t_dec - (t1+t2+t3));

plot(t_acc, Re_acc_ideal, '--', 'Color', [0 0.6 0], 'LineWidth', 1.5, 'DisplayName', 'Ideal acceleration');
plot(t_dec, Re_dec_ideal, '-.', 'Color', [1 0 0], 'LineWidth', 1.5, 'DisplayName', 'Ideal deceleration');

yl = ylim;
for vm = vmarks
    plot([vm vm], yl, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
ylim(yl); grid on; box on;
xlabel('t_{adv} (non-dimensional)');
ylabel('Re');
title(sprintf('Raised-cosine double hold: t1=%.3g, t2=%.3g, t3=%.3g, t4=%.3g, T=%.3g', ...
    t1, t2, t3, t4, T));
legend('Location','best');

%% ---------- Figure 2: u_b/<u_b> ----------
figure('Name','u_b/<u_b> - One period','Color','w');
plot(x1, u1, 'LineWidth', 1.2, 'DisplayName', 'u_b/<u_b>'); hold on;

yl = ylim;
for vm = vmarks
    plot([vm vm], yl, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
ylim(yl); grid on; box on;
xlabel('t_{adv} (non-dimensional)');
ylabel('u_b/\langle u_b\rangle');
title('Dimensionless velocity (mean = 1)');
legend('Location','best');

%% ---------- Figure 3: full signal ----------
figure('Name','Full signal (n periods)','Color','w');
plot(x, u_bulk_avg, '-o', 'MarkerEdgeColor', 'r', 'LineWidth', 1.0, 'MarkerSize', 1);
grid on; box on;
xlabel('t_{adv} (non-dimensional)');
ylabel('u_b/\langle u_b\rangle');
title(sprintf('Full signal (%d cycles, T=%.3g, N=%d)', n, T, N));

%% ---------- OPTIONAL SAVING ----------
if save_figs
    f1 = findobj('Type','figure','Name','Re(t_adv) - One period');
    f2 = findobj('Type','figure','Name','u_b/<u_b> - One period');
    f3 = findobj('Type','figure','Name','Full signal (n periods)');
    if ~isempty(f1), exportgraphics(f1, 'Re_one_period.png', 'Resolution',200); end
    if ~isempty(f2), exportgraphics(f2, 'u_over_umean_one_period.png', 'Resolution',200); end
    if ~isempty(f3), exportgraphics(f3, 'u_over_umean_full_signal.png', 'Resolution',200); end
    fprintf('Figures saved: Re_one_period.png, u_over_umean_one_period.png, u_over_umean_full_signal.png\n');
end

end
