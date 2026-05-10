% run_qs_all_solvers.m
% QS complex 实验：d=5,10,15 (四求解器) + d=20,25 (仅 ManiCSDP + ManiSDP)
% Seeds: rng(1), rng(2), rng(3)
% 每个 (d, seed) 组合保存为独立 txt 文件

base = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(base, '..', 'src')));
addpath(genpath(fullfile(base, '..', 'manopt')));
addpath(genpath('D:\Mine\Mani\ManiSDP-matlab-main\ManiSDP-matlab-main\src'));
addpath('D:\matlab\cvx\sedumi');
addpath('D:\Mosek\11.1\toolbox\r2019b');

out_dir = fullfile(base, 'results_qs_all_solvers');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

warning('off', 'all');

seeds = [1, 2, 3];

% ============================================================
%  Phase 1: d=5,10,15 — 四个求解器全测
% ============================================================
for d = [5, 10, 15]
    N = (2*d+1)*(d+1);

    for seed = seeds
        fname = fullfile(out_dir, sprintf('qs_d%d_rng%d.txt', d, seed));
        if exist(fname, 'file'), delete(fname); end
        diary(fname); diary on;

        fprintf('========================================\n');
        fprintf('QS complex  d=%d  N=%d  rng(%d)\n', d, N, seed);
        fprintf('Started: %s\n', datetime("now"));
        fprintf('========================================\n\n');

        rng(seed);
        Q_quartic = randn(N) + 1i*randn(N);
        Q_quartic = (Q_quartic + Q_quartic') / 2;

        % --- 生成约束矩阵 ---
        fprintf('--- 生成约束矩阵 ---\n');
        t0 = tic;
        [A_sdp, b_sdp, c_sdp, K] = qsmom_complex(Q_quartic);
        t_gen = toc(t0);
        fprintf('qsmom_complex: n=%d, m=%d (%.2fs)\n', K.s, length(b_sdp), t_gen);

        [At_r, b_r, c_r, K_r] = convertCtoR(A_sdp, b_sdp, c_sdp, K);
        fprintf('convertCtoR:   n_r=%d, m_r=%d\n\n', K_r.s, length(b_r));

        % ---- ManiCSDP ----
        fprintf('--- ManiCSDP (complex SDP) ---\n');
        opts_c = struct();
        opts_c.sigma0      = 1;
        opts_c.AL_maxiter  = 500;
        opts_c.p0          = 1;
        opts_c.tol         = 1e-8;
        opts_c.tau1        = 1e-2;
        opts_c.tau2        = 5e-2;
        opts_c.delta       = 6;
        opts_c.TR_maxinner = 20;
        opts_c.TR_maxiter  = 4;
        opts_c.line_search = 1;
        obj_mani = NaN; emani = NaN; t_mani = NaN;
        try
            tic;
            [~, obj_mani, data_mani] = ManiCSDP(A_sdp', b_sdp, c_sdp, K, opts_c);
            t_mani = toc;
            emani = max([data_mani.gap, data_mani.pinf, data_mani.dinf]);
            fprintf('ManiCSDP: optimum=%.8f, eta=%.1e, time=%.2fs\n', obj_mani, emani, t_mani);
        catch ME
            fprintf('ManiCSDP ERROR: %s\n', ME.message);
        end

        % ---- ManiSDP ----
        fprintf('\n--- ManiSDP (real SDP) ---\n');
        opts_r = struct();
        opts_r.p0          = 1;
        opts_r.sigma0      = 1;
        opts_r.sigma_min   = 1e-1;
        opts_r.theta       = 1e-2;
        opts_r.delta       = 6;
        opts_r.tau1        = 1e-2;
        opts_r.tau2        = 1e-2;
        opts_r.TR_maxinner = 20;
        opts_r.TR_maxiter  = 4;
        opts_r.line_search = 1;
        obj_r = NaN; emani_r = NaN; t_msdp = NaN;
        try
            tic;
            [~, obj_r, data_r] = ManiSDP(At_r, b_r, c_r, K_r, opts_r);
            t_msdp = toc;
            emani_r = max([data_r.gap, data_r.pinf, data_r.dinf]);
            fprintf('ManiSDP: optimum=%.8f, eta=%.1e, time=%.2fs\n', obj_r, emani_r, t_msdp);
        catch ME
            fprintf('ManiSDP ERROR: %s\n', ME.message);
        end

        % ---- SeDuMi ----
        fprintf('\n--- SeDuMi (complex SDP) ---\n');
        obj_sed = NaN; t_sed = NaN;
        try
            K_sedumi        = struct();
            K_sedumi.s      = K.s;
            K_sedumi.scomplex = 1;
            pars_sed.maxiter = 200;
            pars_sed.fid    = 0;
            tic;
            [x_sed, y_sed, info_sed] = sedumi(A_sdp', b_sdp, c_sdp, K_sedumi, pars_sed);
            t_sed = toc;
            obj_sed = real(c_sdp' * x_sed);
            pinf_s = norm(A_sdp*x_sed - b_sdp) / max(1, norm(b_sdp));
            S_mat  = reshape(c_sdp - A_sdp' * y_sed, K.s, K.s);
            dS     = eig((S_mat + S_mat') / 2);
            dinf_s = max(0, -min(real(dS))) / (1 + max(real(dS)));
            gap_s  = abs(c_sdp'*x_sed - b_sdp'*y_sed) / ...
                     (1 + abs(c_sdp'*x_sed) + abs(b_sdp'*y_sed));
            fprintf('SeDuMi: optimum=%.8f, pinf=%.2e, dinf=%.2e, gap=%.2e, time=%.2fs', ...
                    obj_sed, pinf_s, dinf_s, gap_s, t_sed);
            if info_sed.numerr ~= 0
                fprintf('  [numerr=%d]', info_sed.numerr);
            end
            fprintf('\n');
        catch ME
            fprintf('SeDuMi ERROR: %s\n', ME.message);
        end

        % ---- MOSEK ----
        fprintf('\n--- MOSEK (real SDP) ---\n');
        obj_mos = NaN; t_mos = NaN;
        try
            prob = convert_sedumi2mosek(At_r, b_r, c_r, K_r);
            param_mos.MSK_DPAR_OPTIMIZER_MAX_TIME = 10000;
            param_mos.MSK_IPAR_LOG = 0;
            tic;
            [rcode, res] = mosekopt('minimize echo(0)', prob, param_mos);
            t_mos = toc;
            if rcode == 0
                K_mosek.s = K_r.s;
                [X_mosek, y_mosek, S_mosek, mobj] = recover_mosek_sol_blk(res, K_mosek);
                if ~isempty(mobj)
                    x_mos  = X_mosek{1}(:);
                    pinf_m = norm(At_r'*x_mos - b_r) / max(1, norm(b_r));
                    by     = b_r' * y_mosek;
                    gap_m  = abs(mobj(1)-by) / (1+abs(mobj(1))+abs(by));
                    dS     = eig(S_mosek{1});
                    dinf_m = max(0, -min(dS)) / (1 + max(dS));
                    obj_mos = mobj(1);
                    fprintf('MOSEK: optimum=%.8f, pinf=%.2e, dinf=%.2e, gap=%.2e, time=%.2fs\n', ...
                            obj_mos, pinf_m, dinf_m, gap_m, t_mos);
                else
                    fprintf('MOSEK 返回空解。\n');
                end
            else
                fprintf('MOSEK 求解失败 (rcode=%d)\n', rcode);
                if isfield(res, 'rmsg'), fprintf('  %s\n', res.rmsg); end
            end
        catch ME
            fprintf('MOSEK ERROR: %s\n', ME.message);
        end

        % ---- 汇总 ----
        fprintf('\n=== SUMMARY  d=%d  rng(%d) ===\n', d, seed);
        fprintf('Solver      optimum           eta       time\n');
        fprintf('ManiCSDP    %.8f    %.1e  %.2fs\n', obj_mani, emani,   t_mani);
        fprintf('ManiSDP     %.8f    %.1e  %.2fs\n', obj_r,    emani_r, t_msdp);
        fprintf('SeDuMi      %.8f              %.2fs\n', obj_sed, t_sed);
        fprintf('MOSEK       %.8f              %.2fs\n', obj_mos, t_mos);
        fprintf('========================================\n\n');

        diary off;
        fprintf('Saved: %s\n', fname);
    end
end

% ============================================================
%  Phase 2: d=20,25 — 仅 ManiCSDP + ManiSDP
% ============================================================
for d = [20, 25]
    N = (2*d+1)*(d+1);

    for seed = seeds
        fname = fullfile(out_dir, sprintf('qs_d%d_rng%d.txt', d, seed));
        if exist(fname, 'file'), delete(fname); end
        diary(fname); diary on;

        fprintf('========================================\n');
        fprintf('QS complex  d=%d  N=%d  rng(%d)\n', d, N, seed);
        fprintf('Started: %s\n', datetime("now"));
        fprintf('========================================\n\n');

        rng(seed);
        Q_quartic = randn(N) + 1i*randn(N);
        Q_quartic = (Q_quartic + Q_quartic') / 2;

        % --- 生成约束矩阵 ---
        fprintf('--- 生成约束矩阵 ---\n');
        t0 = tic;
        [A_sdp, b_sdp, c_sdp, K] = qsmom_complex(Q_quartic);
        t_gen = toc(t0);
        fprintf('qsmom_complex: n=%d, m=%d (%.2fs)\n', K.s, length(b_sdp), t_gen);

        [At_r, b_r, c_r, K_r] = convertCtoR(A_sdp, b_sdp, c_sdp, K);
        fprintf('convertCtoR:   n_r=%d, m_r=%d\n\n', K_r.s, length(b_r));

        % ---- ManiCSDP ----
        fprintf('--- ManiCSDP (complex SDP) ---\n');
        opts_c = struct();
        opts_c.sigma0      = 1;
        opts_c.AL_maxiter  = 500;
        opts_c.p0          = 1;
        opts_c.tol         = 1e-8;
        opts_c.tau1        = 1e-2;
        opts_c.tau2        = 5e-2;
        opts_c.delta       = 6;
        opts_c.TR_maxinner = 20;
        opts_c.TR_maxiter  = 4;
        opts_c.line_search = 1;
        obj_mani = NaN; emani = NaN; t_mani = NaN;
        try
            tic;
            [~, obj_mani, data_mani] = ManiCSDP(A_sdp', b_sdp, c_sdp, K, opts_c);
            t_mani = toc;
            emani = max([data_mani.gap, data_mani.pinf, data_mani.dinf]);
            fprintf('ManiCSDP: optimum=%.8f, eta=%.1e, time=%.2fs\n', obj_mani, emani, t_mani);
        catch ME
            fprintf('ManiCSDP ERROR: %s\n', ME.message);
        end

        % ---- ManiSDP ----
        fprintf('\n--- ManiSDP (real SDP) ---\n');
        opts_r = struct();
        opts_r.p0          = 1;
        opts_r.sigma0      = 1;
        opts_r.sigma_min   = 1e-1;
        opts_r.theta       = 1e-2;
        opts_r.delta       = 6;
        opts_r.tau1        = 1e-2;
        opts_r.tau2        = 1e-2;
        opts_r.TR_maxinner = 20;
        opts_r.TR_maxiter  = 4;
        opts_r.line_search = 1;
        obj_r = NaN; emani_r = NaN; t_msdp = NaN;
        try
            tic;
            [~, obj_r, data_r] = ManiSDP(At_r, b_r, c_r, K_r, opts_r);
            t_msdp = toc;
            emani_r = max([data_r.gap, data_r.pinf, data_r.dinf]);
            fprintf('ManiSDP: optimum=%.8f, eta=%.1e, time=%.2fs\n', obj_r, emani_r, t_msdp);
        catch ME
            fprintf('ManiSDP ERROR: %s\n', ME.message);
        end

        % ---- 汇总 ----
        fprintf('\n=== SUMMARY  d=%d  rng(%d) ===\n', d, seed);
        fprintf('Solver      optimum           eta       time\n');
        fprintf('ManiCSDP    %.8f    %.1e  %.2fs\n', obj_mani, emani,   t_mani);
        fprintf('ManiSDP     %.8f    %.1e  %.2fs\n', obj_r,    emani_r, t_msdp);
        fprintf('========================================\n\n');

        diary off;
        fprintf('Saved: %s\n', fname);
    end
end

fprintf('\n所有实验完成。结果保存在: %s\n', out_dir);
