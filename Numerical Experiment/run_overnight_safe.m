% Robust overnight test runner.
% BQP: d = 15 and rng seeds 1, 2, 3.
% QS: skipped for this short rerun.
%
% Start from PowerShell with a MATLAB logfile to preserve the full outer
% command-window transcript as well:
% matlab -logfile "D:\Mine\Mani\ManiCSDP\Numerical Experiment\results_bqp_d15_AL1000\matlab_full.log" -batch "addpath(genpath('D:\Mine\Mani\ManiCSDP')); addpath(genpath('D:\Mine\Mani\ManiSDP-matlab-main\ManiSDP-matlab-main')); cd('D:\Mine\Mani\ManiCSDP\Numerical Experiment'); run_overnight_safe"

base = fileparts(mfilename('fullpath'));
repo_root = fullfile(base, '..');
manisdp_root = 'D:\Mine\Mani\ManiSDP-matlab-main\ManiSDP-matlab-main';

addpath(genpath(repo_root));
addpath(genpath(manisdp_root));

out_dir = fullfile(base, 'results_bqp_d15_AL1000');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

warning('off', 'all');

fprintf('========== SAFE OVERNIGHT RUN START ==========\n');
fprintf('Started at: %s\n', datestr(now, 31));
fprintf('Output directory: %s\n', out_dir);
validate_required_functions();

for seed = [1, 2, 3]
    run_case(@() run_bqp_case(15, seed), out_dir, sprintf('bqp_d15_rng%d.txt', seed));
end

fprintf('Finished at: %s\n', datestr(now, 31));
fprintf('========== SAFE OVERNIGHT RUN END ==========\n');

function run_case(case_fun, out_dir, log_name)
    log_file = fullfile(out_dir, log_name);
    if exist(log_file, 'file')
        delete(log_file);
    end

    diary(log_file);
    diary on;
    cleanup = onCleanup(@() diary('off'));

    fprintf('Log file: %s\n', log_file);
    fprintf('Case started at: %s\n', datestr(now, 31));
    try
        case_fun();
        fprintf('Case status: SUCCESS\n');
    catch ME
        fprintf(2, '\nCase status: ERROR\n');
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
    fprintf('Case finished at: %s\n', datestr(now, 31));
    fprintf('========== CASE END ==========\n');
    drawnow;
end

function run_bqp_case(d, seed)
    fprintf('========== BQP complex d=%d rng(%d) ==========\n', d, seed);

    rng(seed);
    Q = (randn(d) + 1i*randn(d)) / sqrt(2);
    Q = (Q + Q') / 2;
    e = (randn(d, 1) + 1i*randn(d, 1)) / sqrt(2);

    [At_cnd, b_cnd, c_cnd, K_cnd] = bqpmom_complex_nounitdiag(Q, e);

    fprintf('\n=== ManiCSDP_unitdiag (complex SDP) ===\n');
    opts_c = bqp_options();
    [~, fval_c, data_c] = ManiCSDP_unitdiag(At_cnd', b_cnd, c_cnd, K_cnd, opts_c);
    fprintf('ManiCSDP_unitdiag: f = %.8f, eta = %.1e, t = %.2fs\n', ...
        fval_c, max([data_c.gap, data_c.pinf, data_c.dinf]), data_c.time);

    fprintf('\n=== convertCtoR ===\n');
    [At_r_nd, b_r_nd, c_r_nd, K_r_nd] = convertCtoR(At_cnd, b_cnd, c_cnd, K_cnd);

    fprintf('\n=== ManiSDP_complexunitdiag (real RSDP) ===\n');
    opts_r = bqp_options();
    [~, fval_r, data_r] = ManiSDP_complexunitdiag(At_r_nd, b_r_nd, c_r_nd, K_r_nd, opts_r);
    fprintf('ManiSDP_complexunitdiag: f = %.8f, eta = %.1e, t = %.2fs\n', ...
        fval_r, max([data_r.gap, data_r.pinf, data_r.dinf]), data_r.time);

    fprintf('\n=== Summary: BQP d=%d rng(%d) ===\n', d, seed);
    fprintf('ManiCSDP_unitdiag:       f = %.8f, t = %.4fs\n', fval_c, data_c.time);
    fprintf('ManiSDP_complexunitdiag: f = %.8f, t = %.4fs\n', fval_r, data_r.time);
end

function opts = bqp_options()
    opts.tol = 1e-8;
    opts.AL_maxiter = 1000;
    opts.p0 = 2;
    opts.sigma0 = 1e-3;
    opts.sigma_min = 1e-2;
    opts.sigma_max = 1e7;
    opts.gama = 2;
    opts.tau1 = 1;
    opts.tau2 = 1;
    opts.theta = 1e-3;
    opts.delta = 8;
    opts.alpha = 0.1;
    opts.TR_maxinner = 25;
    opts.TR_maxiter = 4;
    opts.line_search = 0;
end

function validate_required_functions()
    required = { ...
        'bqpmom_complex_nounitdiag', ...
        'convertCtoR', ...
        'ManiCSDP_unitdiag', ...
        'ManiSDP_complexunitdiag' ...
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
        error('run_overnight_safe:MissingFunction', ...
            'Missing required functions: %s', strjoin(missing, ', '));
    end
    fprintf('=== Function path preflight OK ===\n\n');
end
