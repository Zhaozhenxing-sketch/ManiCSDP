% run_G55_G57.m — 批量对 G55, G56, G57 运行三个求解器并保存输出
clear; clc;
script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);
addpath(genpath(fullfile(script_dir, '..')));   % ManiCSDP 全部子目录（绝对路径）

graphs = {'G55', 'G56', 'G57'};
result_dir = fullfile(pwd, 'results_G55_G57');
if ~exist(result_dir, 'dir'), mkdir(result_dir); end

for gi = 1:length(graphs)
    gname = graphs{gi};

    diary_file = fullfile(result_dir, sprintf('output_%s.txt', gname));
    diary(diary_file);
    diary on;

    fprintf('==================== %s ====================\n', gname);

    %% 加载图
    L = Laplacian(append('../Gset/', gname, '.txt'));
    C = -1/4*sparse(L);
    [At, b, c, K] = unitdiag_constraints(C);
    n = K.s;
    m = length(b);
    fprintf('图 %s: n=%d, m=%d\n', gname, n, m);

    %% SeDuMi
    fprintf('\n调用 SeDuMi (复数 SDP, 最大迭代 300 次)\n');
    K_sedumi        = struct();
    K_sedumi.s      = n;
    K_sedumi.scomplex = 1;
    pars_sedumi.maxiter = 300;
    pars_sedumi.fid     = 0;
    tic;
    try
        [x_sedumi, y_sedumi, info_sedumi] = sedumi(At, b, c, K_sedumi, pars_sedumi);
        t_sedumi = toc;
        if info_sedumi.numerr == 0
            fprintf('SeDuMi 成功完成，用时 %.2f 秒 (迭代 %d 次)\n', t_sedumi, info_sedumi.iter);
        else
            fprintf('SeDuMi 数值问题 (numerr = %d)，用时 %.2f 秒\n', info_sedumi.numerr, t_sedumi);
        end
        obj_sedumi = real(c' * x_sedumi);
        pinf_sed   = norm(At' * x_sedumi - b) / max(1, norm(b));
        S_mat      = reshape(c - At * y_sedumi, n, n);
        dS_sed     = eig((S_mat + S_mat') / 2);
        dinf_sed   = max(0, -min(real(dS_sed))) / (1 + max(real(dS_sed)));
        gap_sed    = abs(c'*x_sedumi - b'*y_sedumi) / ...
                     (1 + abs(c'*x_sedumi) + abs(b'*y_sedumi));
        fprintf('  最优值 = %.8f, pinf = %.2e, dinf = %.2e, gap = %.2e\n', ...
                obj_sedumi, pinf_sed, dinf_sed, gap_sed);
    catch ME
        fprintf('SeDuMi 出错: %s\n', ME.message);
    end

    %% ManiCSDP_onlyunitdiag
    fprintf('\n调用 ManiCSDP_onlyunitdiag\n');
    options_cu              = struct();
    options_cu.p0           = 2;
    options_cu.tol          = 1e-8;
    options_cu.theta        = 1e-1;
    options_cu.delta        = 15;
    options_cu.alpha        = 0.5;
    options_cu.TR_maxiter   = 40;
    options_cu.TR_maxinner  = 100;
    options_cu.line_search  = 0;
    tic;
    [X_cu, obj_cu, data_cu] = ManiCSDP_onlyunitdiag(C, options_cu);
    t_cu = toc;
    fprintf('ManiCSDP_onlyunitdiag: optimum = %0.8f, dinf = %0.1e, time = %0.2fs\n', ...
            obj_cu, data_cu.dinf, t_cu);

    %% convertCtoR + ManiSDP_onlycomplexunitdiag
    [At_r, b_r, c_r, K_r] = convertCtoR(At', b, c, K);
    n_r = K_r.s;
    C_r = reshape(c_r, n_r, n_r);
    fprintf('\n约束转换完成，开始求解\n');
    fprintf('调用 ManiSDP_onlycomplexunitdiag\n');
    options_pob             = struct();
    options_pob.p0          = 2;
    options_pob.tol         = 1e-8;
    options_pob.theta       = 1e-1;
    options_pob.delta       = 15;
    options_pob.alpha       = 0.5;
    options_pob.TR_maxiter  = 40;
    options_pob.TR_maxinner = 100;
    options_pob.line_search = 0;
    tic;
    [X_pob, obj_pob, data_pob] = ManiSDP_onlycomplexunitdiag(C_r, options_pob);
    t_pob = toc;
    fprintf('ManiSDP_onlycomplexunitdiag: optimum = %0.8f, dinf = %0.1e, time = %0.2fs\n', ...
            obj_pob, data_pob.dinf, t_pob);

    fprintf('\n==================== %s 完成 ====================\n\n', gname);
    diary off;
    fprintf('输出已保存至: %s\n', diary_file);
end

fprintf('\n所有图运行完毕，结果保存在 %s\n', result_dir);
