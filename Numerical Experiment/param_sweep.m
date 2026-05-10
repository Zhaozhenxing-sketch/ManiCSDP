% param_sweep.m — Full-factorial joint sweep for G1 (real), Option A
% Variables: sigma0, sigma_min, tau1, tau2, TR_maxinner (4 levels)
% Fixed:     p0=1, tol=1e-8, theta=1e-3, alpha=0.1, TR_maxiter=6,
%            line_search=1, gama=2, sigma_max=1e7, AL_maxiter=100, delta=15
%            + delta supplement {5,8,12} at best sigma config
%
% Constraints: sigma_min < sigma0, tau2 > tau1
% Est. runtime: ~195 configs x 3 repeats x ~35s ≈ 5.7 hours

clear; clc;
cd(fileparts(mfilename('fullpath')));
addpath(genpath('d:/Mine/Mani/ManiSDP-matlab-main/ManiSDP-matlab-main/src'));
addpath(genpath('d:/Mine/Mani/ManiSDP-matlab-main/ManiSDP-matlab-main/manopt7.0'));
addpath(genpath('d:/Mine/Mani/ManiCSDP/src/basicfunction'));

%% Build G1 real problem
fprintf('=== Building G1 real problem ===\n');
L = Laplacian(append('../Gset/', "G1", '.txt'));
C = -1/4*sparse(L);
[At, b, c, K] = unitdiag_constraints(C);
[At_r, b_r, c_r, K_r] = convertCtoR(At', b, c, K);
fprintf('n_r=%d, m_r=%d\n\n', K_r.s, length(b_r));

%% Generate configs programmatically
sigma0_vals = [0.01, 0.1, 0.5, 2.0, 5.0];
sigmin_vals = [1e-3, 1e-2, 1e-1];
tau1_vals   = [1e-3, 1e-2];
tau2_vals   = [5e-3, 1e-2, 5e-2];
trin_vals   = [10, 20, 40, 80];   % expanded: added 10 and 20
delta_main  = 15;
delta_extra = [5, 8, 12];

% Pre-count valid combinations to preallocate
n_main = 0;
for s0 = sigma0_vals
    for sm = sigmin_vals
        if sm >= s0, continue; end
        for t1 = tau1_vals
            for t2 = tau2_vals
                if t2 <= t1, continue; end
                n_main = n_main + numel(trin_vals);
            end
        end
    end
end
N_total = n_main + numel(delta_extra);
configs = cell(N_total, 7);
k = 0;

% Prioritise known-good sigma region first so early results are meaningful
priority_s0 = [2.0, 0.5, 5.0, 0.1, 0.01];
for s0 = priority_s0
    for sm = sigmin_vals
        if sm >= s0, continue; end
        for t1 = tau1_vals
            for t2 = tau2_vals
                if t2 <= t1, continue; end
                for trin = trin_vals
                    k = k + 1;
                    label = sprintf('s0=%.2g sm=%.0e t1=%.0e t2=%.0e in=%d', ...
                                    s0, sm, t1, t2, trin);
                    configs(k,:) = {label, s0, sm, t1, t2, trin, delta_main};
                end
            end
        end
    end
end

% Delta supplement at known-best sigma config
for dlt = delta_extra
    k = k + 1;
    label = sprintf('s0=2.0 sm=1e-1 t1=1e-2 t2=5e-2 in=40 dlt=%d', dlt);
    configs(k,:) = {label, 2.0, 1e-1, 1e-2, 5e-2, 40, dlt};
end

N = size(configs, 1);
fprintf('Total configs: %d  (x3 repeats ≈ %.1f h)\n\n', N, N*3*35/3600);

n_repeats = 3;
times = nan(N, n_repeats);
objs  = nan(N, n_repeats);
convs = false(N, n_repeats);

% Resume from checkpoint if it exists and matches current sweep
i_start = 1;
chk_file = 'param_sweep_checkpoint.mat';
if exist(chk_file, 'file')
    chk = load(chk_file);
    if isfield(chk,'N') && chk.N == N && isfield(chk,'i')
        times   = chk.times;
        objs    = chk.objs;
        convs   = chk.convs;
        i_start = chk.i + 1;
        fprintf('Checkpoint found — resuming from config %d/%d\n\n', i_start, N);
    else
        fprintf('Checkpoint found but layout mismatch — starting fresh\n\n');
    end
end

fprintf('=== Full-factorial sweep: %d configs x %d repeats ===\n\n', N, n_repeats);

for i = i_start:N
    label = configs{i,1};
    s0    = configs{i,2};
    sm    = configs{i,3};
    t1    = configs{i,4};
    t2    = configs{i,5};
    trin  = configs{i,6};
    dlt   = configs{i,7};

    fprintf('--- [%3d/%d] %s ---\n', i, N, label);

    for r = 1:n_repeats
        options_r = struct();
        options_r.p0          = 1;
        options_r.tol         = 1e-8;
        options_r.theta       = 1e-3;
        options_r.delta       = dlt;
        options_r.alpha       = 0.1;
        options_r.TR_maxiter  = 6;
        options_r.TR_maxinner = trin;
        options_r.line_search = 1;
        options_r.sigma0      = s0;
        options_r.sigma_min   = sm;
        options_r.sigma_max   = 1e7;
        options_r.gama        = 2;
        options_r.AL_maxiter  = 100;
        options_r.tau1        = t1;
        options_r.tau2        = t2;

        try
            [~, obj, data] = ManiSDP(At_r, b_r, c_r, K_r, options_r);
            times(i,r) = data.time;
            objs(i,r)  = obj;
            convs(i,r) = (data.status == 0);
            fprintf('  rep%d: %5.1fs  obj=%.4f  conv=%d\n', ...
                    r, data.time, obj, data.status==0);
        catch ME
            fprintf('  rep%d ERROR: %s\n', r, ME.message);
            if isnan(times(i,r))   % only mark failed if ManiSDP itself crashed
                times(i,r) = inf;
            end
        end
    end
    t_mean = mean(times(i,:), 'omitnan');
    t_std  = std(times(i,:), 0, 'omitnan');
    fprintf('  >>> mean=%5.1fs  std=%.1fs\n\n', t_mean, t_std);

    % Checkpoint: save after every config so results survive early termination
    save('param_sweep_checkpoint.mat', 'configs', 'times', 'objs', 'convs', 'i', 'N', 'n_repeats');
end

%% Summary sorted by mean time
mean_t = mean(times, 2, 'omitnan');
std_t  = std(times, 0, 2, 'omitnan');
[~, idx] = sort(mean_t);

conv_mask = any(convs, 2);
if any(conv_mask)
    best_t = min(mean_t(conv_mask));
else
    best_t = inf;
end

fprintf('\n==================== FULL-FACTORIAL SUMMARY ====================\n');
fprintf('%-52s  %7s  %5s  %5s\n', 'Label', 'Mean(s)', 'Std', 'Conv%');
fprintf('%s\n', repmat('-', 1, 76));
for ii = 1:N
    i = idx(ii);
    conv_pct = 100 * sum(convs(i,:)) / n_repeats;
    if conv_pct == 0
        marker = ' [FAIL]';
    elseif abs(mean_t(i) - best_t) < 0.1
        marker = ' <-- BEST';
    else
        marker = '';
    end
    fprintf('%-52s  %7.1f  %5.1f  %4.0f%%%s\n', ...
            configs{i,1}, mean_t(i), std_t(i), conv_pct, marker);
end
fprintf('%s\n', repmat('=', 1, 76));
baseline_label = 's0=2 sm=1e-1 t1=1e-2 t2=5e-2 in=40';
baseline_idx = strcmp(configs(:,1), baseline_label);
if any(baseline_idx)
    fprintf('Baseline (%s): %.1fs\n', baseline_label, mean_t(baseline_idx));
end
fprintf('Best converged: %.1fs\n', best_t);

% Save final results to a timestamped file (never overwritten)
ts = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
final_file = sprintf('param_sweep_final_%s.mat', ts);
save(final_file, 'configs', 'times', 'objs', 'convs', 'mean_t', 'std_t', 'idx');
fprintf('Final results saved to %s\n', final_file);
