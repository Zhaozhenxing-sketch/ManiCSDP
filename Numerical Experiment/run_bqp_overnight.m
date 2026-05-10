% Overnight BQP experiment runner
%
% Phase 1: d=35, seeds 1-5  -> ManiCSDP_unitdiag + ManiSDP_complexunitdiag
% Phase 2: d=25, seeds 1,2,4 -> SeDuMi only
% Phase 3: d=30, seeds 1,2,3 -> SeDuMi only
%
% Launch from PowerShell via run_bqp_overnight_launcher.ps1

%% ---- Paths -------------------------------------------------------
addpath(genpath('D:\Mine\Mani\ManiCSDP'));
addpath(genpath('D:\Mine\Mani\ManiSDP-matlab-main\ManiSDP-matlab-main'));
addpath(genpath('D:\matlab\cvx\sedumi'));

%% ---- Output directory --------------------------------------------
base    = fileparts(mfilename('fullpath'));
out_dir = fullfile(base, 'results_bqp_overnight');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

warning('off', 'all');
fprintf('========== OVERNIGHT RUN START: %s ==========\n', datestr(now, 31));
fprintf('Output directory: %s\n\n', out_dir);

validate_required_functions({'bqpmom_complex_nounitdiag','bqpmom_complex', ...
    'convertCtoR','ManiCSDP_unitdiag','ManiSDP_complexunitdiag','sedumi'});

%% ---- Phase 1: d=35, Mani solvers only ----------------------------
fprintf('\n========== Phase 1: d=35, ManiCSDP_unitdiag + ManiSDP_complexunitdiag ==========\n');
for seed = [1, 2, 3, 4, 5]
    run_case(@() run_mani_case(35, seed), out_dir, sprintf('bqp_mani_d35_rng%d.txt', seed));
end

%% ---- Phase 2: d=25, SeDuMi only ----------------------------------
fprintf('\n========== Phase 2: d=25, SeDuMi ==========\n');
for seed = [1, 2, 4]
    run_case(@() run_sedumi_case(25, seed), out_dir, sprintf('bqp_sedumi_d25_rng%d.txt', seed));
end

%% ---- Phase 3: d=30, SeDuMi only ----------------------------------
fprintf('\n========== Phase 3: d=30, SeDuMi ==========\n');
for seed = [1, 2, 3]
    run_case(@() run_sedumi_case(30, seed), out_dir, sprintf('bqp_sedumi_d30_rng%d.txt', seed));
end

fprintf('\n========== OVERNIGHT RUN END: %s ==========\n', datestr(now, 31));

%% ===================================================================
%  Local functions
%% ===================================================================

function run_case(case_fun, out_dir, log_name)
    log_file = fullfile(out_dir, log_name);
    if exist(log_file, 'file'), delete(log_file); end
    diary(log_file); diary on;
    cleanup = onCleanup(@() diary('off'));
    fprintf('Log: %s\n', log_file);
    fprintf('Started: %s\n', datestr(now, 31));
    try
        case_fun();
        fprintf('Status: SUCCESS\n');
    catch ME
        fprintf(2, 'Status: ERROR\n');
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
    fprintf('Finished: %s\n', datestr(now, 31));
    fprintf('========== CASE END ==========\n');
    drawnow;
end

function run_mani_case(d, seed)
    fprintf('========== BQP Mani  d=%d  rng(%d) ==========\n', d, seed);
    rng(seed);
    Q = (randn(d) + 1i*randn(d)) / sqrt(2);  Q = (Q + Q') / 2;
    e = (randn(d, 1) + 1i*randn(d, 1)) / sqrt(2);

    [At_cnd, b_cnd, c_cnd, K_cnd] = bqpmom_complex_nounitdiag(Q, e);

    fprintf('\n--- ManiCSDP_unitdiag ---\n');
    opts = mani_options();
    [~, fval_c, data_c] = ManiCSDP_unitdiag(At_cnd', b_cnd, c_cnd, K_cnd, opts);
    fprintf('ManiCSDP_unitdiag: f = %.8f, eta = %.1e, t = %.2fs\n', ...
        fval_c, max([data_c.gap, data_c.pinf, data_c.dinf]), data_c.time);

    [At_r_nd, b_r_nd, c_r_nd, K_r_nd] = convertCtoR(At_cnd, b_cnd, c_cnd, K_cnd);

    fprintf('\n--- ManiSDP_complexunitdiag ---\n');
    [~, fval_r, data_r] = ManiSDP_complexunitdiag(At_r_nd, b_r_nd, c_r_nd, K_r_nd, opts);
    fprintf('ManiSDP_complexunitdiag: f = %.8f, eta = %.1e, t = %.2fs\n', ...
        fval_r, max([data_r.gap, data_r.pinf, data_r.dinf]), data_r.time);

    fprintf('\n--- Summary  d=%d rng(%d) ---\n', d, seed);
    fprintf('ManiCSDP_unitdiag:       f = %.8f, t = %.4fs\n', fval_c, data_c.time);
    fprintf('ManiSDP_complexunitdiag: f = %.8f, t = %.4fs\n', fval_r, data_r.time);
end

function run_sedumi_case(d, seed)
    fprintf('========== BQP SeDuMi  d=%d  rng(%d) ==========\n', d, seed);
    rng(seed);
    Q = (randn(d) + 1i*randn(d)) / sqrt(2);  Q = (Q + Q') / 2;
    e = (randn(d, 1) + 1i*randn(d, 1)) / sqrt(2);

    [At_c, b_c, c_c, K_c] = bqpmom_complex(Q, e);
    N_c = K_c.s;

    K_sedumi          = struct();
    K_sedumi.s        = N_c;
    K_sedumi.scomplex = 1;
    pars.maxiter      = 300;
    pars.fid          = 1;  % show SeDuMi output to diary

    maxc_c = max(abs(c_c));
    fprintf('N (matrix size) = %d,  maxc = %.4e\n', N_c, maxc_c);

    t0 = tic;
    [x_sed, y_sed, info_sed] = sedumi(At_c', b_c, c_c/maxc_c, K_sedumi, pars);
    t_sed = toc(t0);

    if info_sed.numerr == 0
        fprintf('SeDuMi OK, iter = %d, time = %.2fs\n', info_sed.iter, t_sed);
    else
        fprintf('SeDuMi numerr = %d, time = %.2fs\n', info_sed.numerr, t_sed);
    end

    obj_sed  = real(c_c' * x_sed);
    pinf_sed = norm(At_c * x_sed - b_c) / max(1, norm(b_c));
    S_sed    = reshape(c_c - maxc_c * At_c' * y_sed, N_c, N_c);
    S_sed    = (S_sed + S_sed') / 2;
    dS_sed   = real(eig(full(S_sed)));
    dinf_sed = max(0, -min(dS_sed)) / (1 + max(dS_sed));
    gap_sed  = abs(c_c'*x_sed - b_c'*y_sed*maxc_c) / ...
               (1 + abs(c_c'*x_sed) + abs(b_c'*y_sed*maxc_c));

    fprintf('SeDuMi: f = %.8f, pinf = %.2e, dinf = %.2e, gap = %.2e, t = %.2fs\n', ...
        obj_sed, pinf_sed, dinf_sed, gap_sed, t_sed);
end

function opts = mani_options()
    opts.tol         = 1e-8;
    opts.AL_maxiter  = 1000;
    opts.p0          = 2;
    opts.sigma0      = 1e-3;
    opts.sigma_min   = 1e-2;
    opts.sigma_max   = 1e7;
    opts.gama        = 2;
    opts.tau1        = 1;
    opts.tau2        = 1;
    opts.theta       = 1e-3;
    opts.delta       = 8;
    opts.alpha       = 0.1;
    opts.TR_maxinner = 25;
    opts.TR_maxiter  = 4;
    opts.line_search = 0;
end

function validate_required_functions(required)
    fprintf('=== Function path preflight ===\n');
    missing = {};
    for k = 1:numel(required)
        fn = required{k};
        resolved = which(fn);
        if isempty(resolved)
            missing{end+1} = fn; %#ok<AGROW>
            fprintf(2, 'MISSING: %s\n', fn);
        else
            fprintf('%s -> %s\n', fn, resolved);
        end
    end
    if ~isempty(missing)
        error('run_bqp_overnight:MissingFunction', ...
            'Missing: %s', strjoin(missing, ', '));
    end
    fprintf('=== Preflight OK ===\n\n');
end
