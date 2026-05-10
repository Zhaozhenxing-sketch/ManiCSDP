% Overnight test runner
% BQP: d=15 and d=20, seeds 1-3 (ManiCSDP_unitdiag + ManiSDP_complexunitdiag)
% QS:  d=20,          seeds 1-3 (ManiCSDP        + ManiSDP)

base = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(base, '..', 'src')));
addpath(genpath(fullfile(base, '..', 'manopt')));
addpath(genpath('D:\Mine\Mani\ManiSDP-matlab-main\ManiSDP-matlab-main\src\primal'));

out_dir = fullfile(base, 'results_overnight');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

warning('off', 'all');

%% ============================================================
%  BQP tests
%% ============================================================
for d = [15, 20]
    for seed = [1, 2, 3]
        fname = fullfile(out_dir, sprintf('bqp_d%d_rng%d.txt', d, seed));
        if exist(fname, 'file'), delete(fname); end
        diary(fname); diary on;

        fprintf('========== BQP complex  d=%d  rng(%d) ==========\n', d, seed);

        rng(seed);
        Q = (randn(d) + 1i*randn(d)) / sqrt(2);
        Q = (Q + Q') / 2;
        e = (randn(d, 1) + 1i*randn(d, 1)) / sqrt(2);

        [At_cnd, b_cnd, c_cnd, K_cnd] = bqpmom_complex_nounitdiag(Q, e);

        %--- ManiCSDP_unitdiag ---
        fprintf('\n=== ManiCSDP_unitdiag (complex SDP) ===\n');
        oc.tol         = 1e-8;
        oc.AL_maxiter  = 50;
        oc.p0          = 2;
        oc.sigma0      = 1e-3;
        oc.sigma_min   = 1e-2;
        oc.sigma_max   = 1e7;
        oc.gama        = 2;
        oc.tau1        = 1;
        oc.tau2        = 1;
        oc.theta       = 1e-3;
        oc.delta       = 8;
        oc.alpha       = 0.1;
        oc.TR_maxinner = 25;
        oc.TR_maxiter  = 4;
        oc.line_search = 0;

        [~, fval_c, data_c] = ManiCSDP_unitdiag(At_cnd', b_cnd, c_cnd, K_cnd, oc);
        fprintf('ManiCSDP_unitdiag: f = %.8f, eta = %.1e, t = %.2fs\n', ...
                fval_c, max([data_c.gap, data_c.pinf, data_c.dinf]), data_c.time);

        %--- convertCtoR ---
        [At_r_nd, b_r_nd, c_r_nd, K_r_nd] = convertCtoR(At_cnd, b_cnd, c_cnd, K_cnd);

        %--- ManiSDP_complexunitdiag ---
        fprintf('\n=== ManiSDP_complexunitdiag (real RSDP) ===\n');
        or.tol         = 1e-8;
        or.AL_maxiter  = 50;
        or.p0          = 2;
        or.sigma0      = 1e-3;
        or.sigma_min   = 1e-2;
        or.sigma_max   = 1e7;
        or.gama        = 2;
        or.tau1        = 1;
        or.tau2        = 1;
        or.theta       = 1e-3;
        or.delta       = 8;
        or.alpha       = 0.1;
        or.TR_maxinner = 25;
        or.TR_maxiter  = 4;
        or.line_search = 0;

        [~, fval_r, data_r] = ManiSDP_complexunitdiag(At_r_nd, b_r_nd, c_r_nd, K_r_nd, or);
        fprintf('ManiSDP_complexunitdiag: f = %.8f, eta = %.1e, t = %.2fs\n', ...
                fval_r, max([data_r.gap, data_r.pinf, data_r.dinf]), data_r.time);

        fprintf('\n=== Summary: d=%d rng(%d) ===\n', d, seed);
        fprintf('ManiCSDP_unitdiag:        f = %.8f, t = %.4fs\n', fval_c, data_c.time);
        fprintf('ManiSDP_complexunitdiag:  f = %.8f, t = %.4fs\n', fval_r, data_r.time);
        fprintf('========== DONE ==========\n\n');

        diary off;
        fprintf('Saved: %s\n', fname);
    end
end

%% ============================================================
%  QS tests  (d=20, seeds 1-3)
%% ============================================================
d_qs = 20;
N_qs = (d_qs + 1) * (2*d_qs + 1);

for seed = [1, 2, 3]
    fname = fullfile(out_dir, sprintf('qs_d%d_rng%d.txt', d_qs, seed));
    if exist(fname, 'file'), delete(fname); end
    diary(fname); diary on;

    fprintf('========== QS complex  d=%d  rng(%d) ==========\n', d_qs, seed);

    rng(seed);
    Q_quartic = randn(N_qs) + 1i*randn(N_qs);
    Q_quartic = (Q_quartic + Q_quartic') / 2;

    [A_sdp, b_sdp, c_sdp, K] = qsmom_complex(Q_quartic);

    %--- ManiCSDP ---
    fprintf('\n=== ManiCSDP (complex SDP) ===\n');
    opts_c = struct();
    opts_c.sigma0      = 1;
    opts_c.AL_maxiter  = 500;
    opts_c.p0          = 4;
    opts_c.tol         = 1e-8;
    opts_c.tau1        = 1e-2;
    opts_c.tau2        = 1e-2;
    opts_c.delta       = 15;
    opts_c.TR_maxinner = 20;
    opts_c.TR_maxiter  = 6;

    [~, obj_mani, data_mani] = ManiCSDP(A_sdp', b_sdp, c_sdp, K, opts_c);
    emani = max([data_mani.gap, data_mani.pinf, data_mani.dinf]);
    fprintf('ManiCSDP: f = %.8f, eta = %.1e, t = %.2fs\n', obj_mani, emani, data_mani.time);

    %--- convertCtoR + ManiSDP ---
    [At_r, b_r, c_r, K_r] = convertCtoR(A_sdp, b_sdp, c_sdp, K);

    fprintf('\n=== ManiSDP (real RSDP) ===\n');
    opts_r = struct();
    opts_r.p0     = 1;
    opts_r.sigma0 = 1;
    opts_r.tau2   = 1e-2;
    opts_r.delta  = 6;

    [~, obj_r, data_r] = ManiSDP(At_r, b_r, c_r, K_r, opts_r);
    emani_r = max([data_r.gap, data_r.pinf, data_r.dinf]);
    fprintf('ManiSDP: f = %.8f, eta = %.1e, t = %.2fs\n', obj_r, emani_r, data_r.time);

    fprintf('\n=== Summary: QS d=%d rng(%d) ===\n', d_qs, seed);
    fprintf('ManiCSDP:  f = %.8f, t = %.4fs\n', obj_mani, data_mani.time);
    fprintf('ManiSDP:   f = %.8f, t = %.4fs\n', obj_r, data_r.time);
    fprintf('========== DONE ==========\n\n');

    diary off;
    fprintf('Saved: %s\n', fname);
end

fprintf('\nAll overnight tests completed.\n');
