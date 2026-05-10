% run_G55_G57_sedumi.m — 仅对 G55, G56, G57 运行 SeDuMi 并保存输出
clear; clc;
script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);
addpath(genpath(fullfile(script_dir, '..')));
addpath(genpath('d:/matlab/cvx/sedumi'));

graphs = {'G55', 'G56', 'G57'};
result_dir = fullfile(script_dir, 'results_G55_G57');
if ~exist(result_dir, 'dir'), mkdir(result_dir); end

for gi = 1:length(graphs)
    gname = graphs{gi};

    diary_file = fullfile(result_dir, sprintf('output_%s_sedumi.txt', gname));
    diary(diary_file);
    diary on;

    fprintf('==================== %s ====================\n', gname);

    L = Laplacian(append('../Gset/', gname, '.txt'));
    C = -1/4*sparse(L);
    [At, b, c, K] = unitdiag_constraints(C);
    n = K.s;
    fprintf('图 %s: n=%d\n', gname, n);

    fprintf('调用 SeDuMi (复数 SDP, 最大迭代 300 次)\n');
    K_sedumi          = struct();
    K_sedumi.s        = n;
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

    fprintf('==================== %s 完成 ====================\n\n', gname);
    diary off;
    fprintf('输出已保存至: %s\n', diary_file);
end

fprintf('所有图运行完毕，结果保存在 %s\n', result_dir);
