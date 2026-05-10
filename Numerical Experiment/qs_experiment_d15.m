% qs_experiment_d15.m
% Experiment: QS problem, d=15, seeds 1/2/3
% Each seed: ManiCSDP x3, ManiSDP x3 (data generated once per seed)

addpath(genpath('D:\Mine\Mani\ManiCSDP\src'));
addpath(genpath('D:\Mine\Mani\ManiCSDP\manopt'));
addpath(genpath('D:\Mine\Mani\ManiSDP-matlab-main\ManiSDP-matlab-main\src'));

out_file = 'D:\Mine\Mani\ManiCSDP\Numerical Experiment\qs_d15_results.txt';
if exist(out_file, 'file'), delete(out_file); end
diary(out_file);
diary on;

d     = 15;
N     = (d+1)*(2*d+1);
n_runs = 3;
seeds  = [1, 2, 3];

fprintf('=======================================================\n');
fprintf('QS Experiment: d=%d, N=%d\n', d, N);
fprintf('Started: %s\n', datestr(now));
fprintf('Runs per solver per seed: %d\n', n_runs);
fprintf('=======================================================\n\n');

for si = 1:length(seeds)
    seed = seeds(si);
    fprintf('\n###################################################\n');
    fprintf('SEED = rng(%d)\n', seed);
    fprintf('###################################################\n\n');

    rng(seed);
    Q_quartic = randn(N) + 1i*randn(N);
    Q_quartic = (Q_quartic + Q_quartic') / 2;

    fprintf('--- Generating problem data ---\n');
    t0 = tic;
    [A_sdp, b_sdp, c_sdp, K] = qsmom_complex(Q_quartic);
    [At_r, b_r, c_r, K_r]    = convertCtoR(A_sdp, b_sdp, c_sdp, K);
    fprintf('Data ready: complex n=%d m=%d | real n=%d m=%d (%.2fs)\n\n', ...
        K.s, length(b_sdp), K_r.s, length(b_r), toc(t0));

    % ---- ManiCSDP runs ----
    fprintf('=== ManiCSDP (seed=%d) ===\n', seed);
    for run = 1:n_runs
        fprintf('\n--- ManiCSDP Run %d/%d (seed=%d) ---\n', run, n_runs, seed);
        options             = struct();
        options.sigma0      = 1;
        options.AL_maxiter  = 500;
        options.p0          = 4;
        options.tol         = 1e-8;
        options.tau1        = 1e-2;
        options.tau2        = 1e-2;
        options.delta       = 15;
        options.TR_maxinner = 20;
        options.TR_maxiter  = 4;
        try
            tic;
            [~, obj_mani, data_mani] = ManiCSDP(A_sdp', b_sdp, c_sdp, K, options);
            t_mani = toc;
            emani = max([data_mani.gap, data_mani.pinf, data_mani.dinf]);
            fprintf('>>> ManiCSDP Run %d: optimum=%.8f, eta=%.1e, time=%.2fs\n', ...
                    run, obj_mani, emani, t_mani);
        catch ME
            fprintf('>>> ManiCSDP Run %d ERROR: %s\n', run, ME.message);
        end
    end

    % ---- ManiSDP runs ----
    fprintf('\n=== ManiSDP (seed=%d) ===\n', seed);
    for run = 1:n_runs
        fprintf('\n--- ManiSDP Run %d/%d (seed=%d) ---\n', run, n_runs, seed);
        options_r          = struct();
        options_r.p0       = 2;
        options_r.sigma0   = 1;
        try
            tic;
            [~, obj_r, data_r] = ManiSDP(At_r, b_r, c_r, K_r, options_r);
            t_r = toc;
            emani_r = max([data_r.gap, data_r.pinf, data_r.dinf]);
            fprintf('>>> ManiSDP Run %d: optimum=%.8f, eta=%.1e, time=%.2fs\n', ...
                    run, obj_r, emani_r, t_r);
        catch ME
            fprintf('>>> ManiSDP Run %d ERROR: %s\n', run, ME.message);
        end
    end
end

fprintf('\n=======================================================\n');
fprintf('Experiment completed: %s\n', datestr(now));
fprintf('=======================================================\n');
diary off;
