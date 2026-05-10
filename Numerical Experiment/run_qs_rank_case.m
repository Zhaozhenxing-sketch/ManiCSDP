function run_qs_rank_case(d, seed, solver_name)
%RUN_QS_RANK_CASE Run one QS case for SeDuMi or MOSEK and report solution rank.

rank_tol = 1e-5;
solver_name = lower(string(solver_name));

fprintf('========== QS rank case START ==========\n');
fprintf('d = %d, rng(%d), solver = %s\n', d, seed, solver_name);
fprintf('rank tolerance: eigenvalue > %.1e\n', rank_tol);
fprintf('Started at: %s\n', datestr(now, 31));

add_required_paths();
print_function_preflight(solver_name);

rng(seed);
N = (d + 1) * (d + 2) / 2;
fprintf('Generated quartic matrix size N = %d\n', N);
Q_quartic = randn(N) + 1i * randn(N);
Q_quartic = (Q_quartic + Q_quartic') / 2;

tic;
[A_sdp, b_sdp, c_sdp, K] = qsmom_complex(Q_quartic);
fprintf('qsmom_complex: SUCCESS, K.s = %d, m = %d, time = %.2fs\n', ...
    K.s, numel(b_sdp), toc);

switch solver_name
    case "sedumi"
        run_sedumi_rank(A_sdp, b_sdp, c_sdp, K, rank_tol);
    case "mosek"
        run_mosek_rank(A_sdp, b_sdp, c_sdp, K, rank_tol);
    otherwise
        error('run_qs_rank_case:UnknownSolver', ...
            'Unknown solver_name: %s. Expected sedumi or mosek.', solver_name);
end

fprintf('Finished at: %s\n', datestr(now, 31));
fprintf('========== QS rank case END ==========\n');
end

function run_sedumi_rank(A_sdp, b_sdp, c_sdp, K, rank_tol)
fprintf('\n--- SeDuMi rank run ---\n');
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
pinf = norm(A_sdp * x_sed - b_sdp) / max(1, norm(b_sdp));
S_mat = reshape(c_sdp - A_sdp' * y_sed, K.s, K.s);
S_mat = (S_mat + S_mat') / 2;
dS = real(eig(full(S_mat)));
dinf = max(0, -min(dS)) / (1 + max(dS));
gap = abs(c_sdp' * x_sed - b_sdp' * y_sed) / ...
    (1 + abs(c_sdp' * x_sed) + abs(b_sdp' * y_sed));

X_sed = reshape(x_sed, K.s, K.s);
X_sed = (X_sed + X_sed') / 2;
ev = sort(real(eig(full(X_sed))), 'descend');
rank_x = sum(ev > rank_tol);

fprintf('SeDuMi: DONE, numerr = %d, iter = %d, time = %.2fs\n', ...
    info_sed.numerr, info_sed.iter, t_sedumi);
fprintf('SeDuMi: optimum = %.12g, pinf = %.3e, dinf = %.3e, gap = %.3e\n', ...
    obj_sed, pinf, dinf, gap);
print_rank_summary('SeDuMi complex X', ev, rank_x, rank_tol);
end

function run_mosek_rank(A_sdp, b_sdp, c_sdp, K, rank_tol)
fprintf('\n--- convertCtoR for MOSEK ---\n');
tic;
[At_r, b_r, c_r, K_r] = convertCtoR(A_sdp, b_sdp, c_sdp, K);
fprintf('convertCtoR: SUCCESS, K_r.s = %d, m = %d, time = %.2fs\n', ...
    K_r.s, numel(b_r), toc);

fprintf('\n--- MOSEK rank run ---\n');
max_time = 10000;
prob = convert_sedumi2mosek(At_r, b_r, c_r, K_r);
param = struct();
param.MSK_DPAR_OPTIMIZER_MAX_TIME = max_time - 10;
param.MSK_IPAR_LOG = 0;

tic;
[rcode, res] = mosekopt('minimize echo(0)', prob, param);
t_mosek = toc;

if rcode ~= 0
    fprintf('MOSEK: FAILED, response code = %d, time = %.2fs\n', rcode, t_mosek);
    if isfield(res, 'rmsg')
        fprintf('MOSEK message: %s\n', res.rmsg);
    end
    return;
end

K_mosek = struct();
K_mosek.s = K_r.s;
[X_mosek, y_mosek, S_mosek, mobj] = recover_mosek_sol_blk(res, K_mosek);
if isempty(mobj) || isempty(X_mosek)
    fprintf('MOSEK: returned empty solution, time = %.2fs\n', t_mosek);
    return;
end

X_real = full(X_mosek{1});
X_real = (X_real + X_real') / 2;
ev = sort(real(eig(X_real)), 'descend');
rank_real = sum(ev > rank_tol);

x_mos = X_real(:);
pinf = norm(At_r' * x_mos - b_r) / max(1, norm(b_r));
by = b_r' * y_mosek;
gap = abs(mobj(1) - by) / (1 + abs(mobj(1)) + abs(by));
S_mat = full(S_mosek{1});
S_mat = (S_mat + S_mat') / 2;
dS = real(eig(S_mat));
dinf = max(0, -min(dS)) / (1 + max(dS));

fprintf('MOSEK: DONE, time = %.2fs\n', t_mosek);
fprintf('MOSEK: optimum = %.12g, pinf = %.3e, dinf = %.3e, gap = %.3e\n', ...
    mobj(1), pinf, dinf, gap);
print_rank_summary('MOSEK real lifted X', ev, rank_real, rank_tol);
fprintf('MOSEK complex-equivalent rank estimate = %.1f (real lifted rank / 2)\n', ...
    rank_real / 2);
end

function print_rank_summary(label, ev, rank_x, rank_tol)
fprintf('%s rank(> %.1e) = %d\n', label, rank_tol, rank_x);
fprintf('%s eigenvalue max = %.12e, min = %.12e\n', label, max(ev), min(ev));
head_count = min(12, numel(ev));
tail_count = min(6, numel(ev));
fprintf('%s top %d eigenvalues:\n', label, head_count);
fprintf(' %.12e', ev(1:head_count));
fprintf('\n');
fprintf('%s bottom %d eigenvalues:\n', label, tail_count);
fprintf(' %.12e', ev(end - tail_count + 1:end));
fprintf('\n');
end

function add_required_paths()
repo_root = 'D:\Mine\Mani\ManiCSDP';
manisdp_root = 'D:\Mine\Mani\ManiSDP-matlab-main\ManiSDP-matlab-main';
sedumi_root = 'D:\matlab\cvx\sedumi';
mosek_root = 'D:\Mosek\11.1\toolbox\r2019b';

addpath(genpath(repo_root));
addpath(genpath(manisdp_root));
addpath(genpath(sedumi_root));
addpath(genpath(mosek_root));
end

function print_function_preflight(solver_name)
required = {'qsmom_complex'};
if solver_name == "sedumi"
    required = [required, {'sedumi'}];
elseif solver_name == "mosek"
    required = [required, {'convertCtoR', 'convert_sedumi2mosek', ...
        'mosekopt', 'recover_mosek_sol_blk'}];
end

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
    error('run_qs_rank_case:MissingFunction', ...
        'Missing required functions: %s', strjoin(missing, ', '));
end
fprintf('=== Function path preflight END ===\n\n');
end
