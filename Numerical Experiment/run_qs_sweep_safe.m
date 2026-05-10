% Robust QS sweep runner.
% Phase 1: d = 45, 50, 55, 60 with rng seeds 1, 2, 3.
% Phase 2: d = 65, 70 with rng seeds 1, 2, 3, after phase 1 is attempted.
%
% Each d/seed case writes one diary log. Each solver is protected by its own
% try/catch, so a solver failure does not stop later solvers or cases.
% MOSEK is intentionally skipped because mosekopt caused a MATLAB-level
% access violation on d = 55, rng(1); such MEX crashes cannot be caught.

base = fileparts(mfilename('fullpath'));
repo_root = fullfile(base, '..');
manisdp_root = 'D:\Mine\Mani\ManiSDP-matlab-main\ManiSDP-matlab-main';
sedumi_root = 'D:\matlab\cvx\sedumi';

addpath(genpath(repo_root));
addpath(genpath(manisdp_root));
addpath(genpath(sedumi_root));

out_dir = fullfile(base, 'results_qs_sweep_safe');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

warning('off', 'all');

fprintf('========== QS SWEEP SAFE RUN START ==========\n');
fprintf('Started at: %s\n', datestr(now, 31));
fprintf('Output directory: %s\n', out_dir);
print_function_preflight();

phase1_d = [45, 50, 55, 60];
phase2_d = [65, 70];
seeds = [1, 2, 3];

run_phase('phase1', phase1_d, seeds, out_dir);
run_phase('phase2', phase2_d, seeds, out_dir);

fprintf('Finished at: %s\n', datestr(now, 31));
fprintf('========== QS SWEEP SAFE RUN END ==========\n');

function run_phase(phase_name, d_values, seeds, out_dir)
    fprintf('\n========== %s START ==========\n', upper(phase_name));
    for d = d_values
        for seed = seeds
            log_name = sprintf('qs_d%d_rng%d.txt', d, seed);
            run_case(@() run_qs_case(d, seed), out_dir, log_name);
        end
    end
    fprintf('========== %s END ==========\n\n', upper(phase_name));
end

function run_case(case_fun, out_dir, log_name)
    log_file = fullfile(out_dir, log_name);
    if exist(log_file, 'file')
        fprintf('Skipping existing log file: %s\n', log_file);
        return;
    end

    diary(log_file);
    diary on;
    cleanup = onCleanup(@() diary('off'));

    fprintf('Log file: %s\n', log_file);
    fprintf('Case started at: %s\n', datestr(now, 31));
    try
        case_fun();
        fprintf('Case status: ATTEMPTED\n');
    catch ME
        fprintf(2, '\nCase status: FATAL ERROR\n');
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
    fprintf('Case finished at: %s\n', datestr(now, 31));
    fprintf('========== CASE END ==========\n');
    drawnow;
end

function run_qs_case(d, seed)
    fprintf('========== QS complex d=%d rng(%d) ==========\n', d, seed);

    A_sdp = [];
    b_sdp = [];
    c_sdp = [];
    K = [];
    At_r = [];
    b_r = [];
    c_r = [];
    K_r = [];

    try
        rng(seed);
        N = (d + 1) * (d + 2) / 2;
        fprintf('Generated quartic matrix size N = %d\n', N);
        Q_quartic = randn(N) + 1i*randn(N);
        Q_quartic = (Q_quartic + Q_quartic') / 2;

        tic;
        [A_sdp, b_sdp, c_sdp, K] = qsmom_complex(Q_quartic);
        fprintf('qsmom_complex: SUCCESS, K.s = %d, m = %d, time = %.2fs\n', ...
            K.s, numel(b_sdp), toc);
    catch ME
        fprintf(2, '\nqsmom_complex: ERROR\n');
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        fprintf('Skipping all solvers for d=%d rng(%d) because model generation failed.\n', d, seed);
        return;
    end

    run_mani_csdp(A_sdp, b_sdp, c_sdp, K);
    run_sedumi(A_sdp, b_sdp, c_sdp, K);

    try
        fprintf('\n--- convertCtoR ---\n');
        tic;
        [At_r, b_r, c_r, K_r] = convertCtoR(A_sdp, b_sdp, c_sdp, K);
        fprintf('convertCtoR: SUCCESS, K_r.s = %d, m = %d, time = %.2fs\n', ...
            K_r.s, numel(b_r), toc);
    catch ME
        fprintf(2, '\nconvertCtoR: ERROR\n');
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        fprintf('Skipping real solvers for d=%d rng(%d) because conversion failed.\n', d, seed);
        return;
    end

    run_mani_sdp(At_r, b_r, c_r, K_r);

    fprintf('\n=== Summary marker: QS d=%d rng(%d) attempted first three solvers; MOSEK skipped ===\n', d, seed);
end

function run_mani_csdp(A_sdp, b_sdp, c_sdp, K)
    fprintf('\n--- ManiCSDP ---\n');
    try
        options = struct();
        options.sigma0 = 1;
        options.AL_maxiter = 500;
        options.p0 = 1;
        options.tol = 1e-8;
        options.tau1 = 1e-2;
        options.tau2 = 5e-2;
        options.delta = 6;
        options.TR_maxinner = 20;
        options.TR_maxiter = 8;
        options.line_search = 1;

        tic;
        [~, obj_mani, data_mani] = ManiCSDP(A_sdp', b_sdp, c_sdp, K, options);
        tmani = toc;
        emani = max_eta(data_mani);
        fprintf('ManiCSDP: SUCCESS, optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', ...
            obj_mani, emani, tmani);
    catch ME
        fprintf(2, '\nManiCSDP: ERROR\n');
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

function run_sedumi(A_sdp, b_sdp, c_sdp, K)
    fprintf('\n--- SeDuMi ---\n');
    try
        K_sedumi = struct();
        K_sedumi.s = K.s;
        K_sedumi.scomplex = 1;
        pars = struct();
        pars.maxiter = 200;
        pars.fid = 0;

        tic;
        [x_sed, y_sed, info_sed] = sedumi(A_sdp', b_sdp, c_sdp, K_sedumi, pars);
        t_sedumi = toc;
        obj_sed = real(c_sdp' * x_sed);
        pinf = norm(A_sdp*x_sed - b_sdp) / max(1, norm(b_sdp));
        S_mat = reshape(c_sdp - A_sdp' * y_sed, K.s, K.s);
        dS = eig((S_mat + S_mat') / 2);
        dinf = max(0, -min(real(dS))) / (1 + max(real(dS)));
        gap = abs(c_sdp'*x_sed - b_sdp'*y_sed) / ...
            (1 + abs(c_sdp'*x_sed) + abs(b_sdp'*y_sed));
        fprintf('SeDuMi: DONE, numerr = %d, iter = %d, time = %.2fs\n', ...
            info_sed.numerr, info_sed.iter, t_sedumi);
        fprintf('SeDuMi: optimum = %.8f, pinf = %.2e, dinf = %.2e, gap = %.2e\n', ...
            obj_sed, pinf, dinf, gap);
    catch ME
        fprintf(2, '\nSeDuMi: ERROR\n');
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

function run_mani_sdp(At_r, b_r, c_r, K_r)
    fprintf('\n--- ManiSDP ---\n');
    try
        options_r = struct();
        options_r.p0 = 1;
        options_r.sigma0 = 1;
        options_r.sigma_min = 1e-1;
        options_r.theta = 1e-2;
        options_r.delta = 6;
        options_r.tau1 = 1e-2;
        options_r.tau2 = 5e-2;
        options_r.TR_maxinner = 20;
        options_r.TR_maxiter = 8;
        options_r.line_search = 1;

        tic;
        [~, obj_r, data_r] = ManiSDP(At_r, b_r, c_r, K_r, options_r);
        t_manisdp = toc;
        emani_r = max_eta(data_r);
        fprintf('ManiSDP: SUCCESS, optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', ...
            obj_r, emani_r, t_manisdp);
    catch ME
        fprintf(2, '\nManiSDP: ERROR\n');
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

function eta = max_eta(data)
    vals = [];
    if isfield(data, 'gap'), vals(end+1) = data.gap; end %#ok<AGROW>
    if isfield(data, 'pinf'), vals(end+1) = data.pinf; end %#ok<AGROW>
    if isfield(data, 'dinf'), vals(end+1) = data.dinf; end %#ok<AGROW>
    if isempty(vals)
        eta = NaN;
    else
        eta = max(vals);
    end
end

function print_function_preflight()
    required = { ...
        'qsmom_complex', ...
        'ManiCSDP', ...
        'sedumi', ...
        'convertCtoR', ...
        'ManiSDP' ...
    };

    fprintf('\n=== Function path preflight ===\n');
    missing = {};
    for k = 1:numel(required)
        fn = required{k};
        resolved = which(fn);
        if isempty(resolved)
            missing{end + 1} = fn; %#ok<AGROW>
            fprintf(2, 'MISSING: %s\n', fn);
        else
            fprintf('%s -> %s\n', fn, resolved);
        end
    end
    if ~isempty(missing)
        error('run_qs_sweep_safe:MissingFunction', ...
            'Missing required functions: %s', strjoin(missing, ', '));
    end
    fprintf('=== Function path preflight END ===\n\n');
end
