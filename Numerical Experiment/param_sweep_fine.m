% param_sweep_fine.m — Fine-grained sweep around Phase-3 winner
% Winner: sigma0=2, sigma_min=1e-1, tau1=1e-3, tau2=1e-2, in=40, delta=15, TR_maxiter=6
%
% Main grid (108 configs): sigma0 x sigma_min x tau1 x tau2
%   sigma0:    {1.5, 2.0, 2.5, 3.0}
%   sigma_min: {5e-2, 1e-1, 2e-1}  (all < min sigma0=1.5)
%   tau1:      {5e-4, 1e-3, 2e-3}
%   tau2:      {7e-3, 1e-2, 1.5e-2}
%   Fixed: TR_maxinner=40, delta=15, TR_maxiter=6
%
% Supplements (15 configs) at winner (s0=2, sm=0.1, t1=1e-3, t2=1e-2):
%   TR_maxinner: {30, 35, 45, 50}
%   delta:       {12, 13, 14, 16, 17, 18}
%   TR_maxiter:  {3, 4, 5, 8, 10}
%
% Est: 123 configs x 3 repeats x ~31s = 3.2h

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

%% Build config list
sigma0_vals = [1.5, 2.0, 2.5, 3.0];
sigmin_vals = [5e-2, 1e-1, 2e-1];
tau1_vals   = [5e-4, 1e-3, 2e-3];
tau2_vals   = [7e-3, 1e-2, 1.5e-2];

% Pre-count main configs
n_main = 0;
for s0 = sigma0_vals
    for sm = sigmin_vals
        if sm >= s0, continue; end
        for t1 = tau1_vals
            for t2 = tau2_vals
                if t2 <= t1, continue; end
                n_main = n_main + 1;
            end
        end
    end
end

% Supplement configs (TR_maxinner / delta / TR_maxiter at winner)
trin_supp    = [30, 35, 45, 50];
delta_supp   = [12, 13, 14, 16, 17, 18];
triter_supp  = [3, 4, 5, 8, 10];
n_supp = numel(trin_supp) + numel(delta_supp) + numel(triter_supp);

N_total = n_main + n_supp;
configs = cell(N_total, 8);   % label,s0,sm,t1,t2,trin,delta,triter
k = 0;

% Priority: winner sigma first
priority_s0 = [2.0, 2.5, 1.5, 3.0];
for s0 = priority_s0
    for sm = sigmin_vals
        if sm >= s0, continue; end
        for t1 = tau1_vals
            for t2 = tau2_vals
                if t2 <= t1, continue; end
                k = k + 1;
                label = sprintf('s0=%.1f sm=%.0e t1=%.0e t2=%.0e', s0, sm, t1, t2);
                configs(k,:) = {label, s0, sm, t1, t2, 40, 15, 6};
            end
        end
    end
end

% TR_maxinner supplement
for trin = trin_supp
    k = k + 1;
    label = sprintf('SUPP_trin=%d', trin);
    configs(k,:) = {label, 2.0, 1e-1, 1e-3, 1e-2, trin, 15, 6};
end

% delta supplement
for dlt = delta_supp
    k = k + 1;
    label = sprintf('SUPP_delta=%d', dlt);
    configs(k,:) = {label, 2.0, 1e-1, 1e-3, 1e-2, 40, dlt, 6};
end

% TR_maxiter supplement
for triter = triter_supp
    k = k + 1;
    label = sprintf('SUPP_triter=%d', triter);
    configs(k,:) = {label, 2.0, 1e-1, 1e-3, 1e-2, 40, 15, triter};
end

N = k;
fprintf('Main configs: %d  Supplements: %d  Total: %d\n', n_main, n_supp, N);
fprintf('Est runtime: %dx3x31s = %.1fh\n\n', N, N*3*31/3600);

%% Init result arrays
n_repeats = 3;
times = nan(N, n_repeats);
objs  = nan(N, n_repeats);
convs = false(N, n_repeats);

% Resume from checkpoint
i_start = 1;
chk_file = 'param_sweep_fine_checkpoint.mat';
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

fprintf('=== Fine sweep: %d configs x %d repeats ===\n\n', N, n_repeats);

%% Main loop
for i = i_start:N
    label  = configs{i,1};
    s0     = configs{i,2};
    sm     = configs{i,3};
    t1     = configs{i,4};
    t2     = configs{i,5};
    trin   = configs{i,6};
    dlt    = configs{i,7};
    triter = configs{i,8};

    fprintf('--- [%3d/%d] %s ---\n', i, N, label);

    for r = 1:n_repeats
        opts = struct();
        opts.p0          = 1;
        opts.tol         = 1e-8;
        opts.theta       = 1e-3;
        opts.delta       = dlt;
        opts.alpha       = 0.1;
        opts.TR_maxiter  = triter;
        opts.TR_maxinner = trin;
        opts.line_search = 1;
        opts.sigma0      = s0;
        opts.sigma_min   = sm;
        opts.sigma_max   = 1e7;
        opts.gama        = 2;
        opts.AL_maxiter  = 100;
        opts.tau1        = t1;
        opts.tau2        = t2;

        try
            [~, obj, data] = ManiSDP(At_r, b_r, c_r, K_r, opts);
            times(i,r) = data.time;
            objs(i,r)  = obj;
            convs(i,r) = (data.status == 0);
            fprintf('  rep%d: %5.1fs  obj=%.4f  conv=%d\n', ...
                    r, data.time, obj, data.status==0);
        catch ME
            fprintf('  rep%d ERROR: %s\n', r, ME.message);
            times(i,r) = inf;
        end
    end
    t_mean = mean(times(i,:), 'omitnan');
    t_std  = std(times(i,:), 0, 'omitnan');
    fprintf('  >>> mean=%5.1fs  std=%.1fs\n\n', t_mean, t_std);

    save(chk_file, 'configs', 'times', 'objs', 'convs', 'i', 'N', 'n_repeats');
end

%% Summary
mean_t = mean(times, 2, 'omitnan');
std_t  = std(times, 0, 2, 'omitnan');
[~, idx] = sort(mean_t);

conv_mask = any(convs, 2);
best_t = inf;
if any(conv_mask)
    best_t = min(mean_t(conv_mask));
end

fprintf('\n==================== FINE SWEEP SUMMARY ====================\n');
fprintf('%-48s  %7s  %5s  %5s\n', 'Label', 'Mean(s)', 'Std', 'Conv%');
fprintf('%s\n', repmat('-', 1, 70));
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
    fprintf('%-48s  %7.1f  %5.1f  %4.0f%%%s\n', ...
            configs{i,1}, mean_t(i), std_t(i), conv_pct, marker);
end
fprintf('%s\n', repmat('=', 1, 70));
fprintf('Phase-3 winner (s0=2 sm=1e-1 t1=1e-3 t2=1e-2 in=40): reference 31.2s\n');
fprintf('Best converged: %.1fs\n', best_t);

%% Save final results
ts = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
final_file = sprintf('param_sweep_fine_final_%s.mat', ts);
save(final_file, 'configs', 'times', 'objs', 'convs', 'mean_t', 'std_t', 'idx');
fprintf('Final results saved to %s\n', final_file);
