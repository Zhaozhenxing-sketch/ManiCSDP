function [At, b, c, K] = mcmom_complex(C_small)
    % 1. 修正变量数：C_small 是 n x n，变量数就是 n
    n = size(C_small, 1); 
    N = 2*n + 1;           % 基底: [1, z_1...z_n, conj(z_1)...conj(z_n)]
    K.s = N;
    
    % 2. 修正目标矩阵 C 的位置：跳过第一行第一列的常数项
    % C_small 对应 z_i * conj(z_j)，放入块 (2:n+1, 2:n+1)
    C_large = sparse(N, N);
    C_large(2:n+1, 2:n+1) = C_small;
    c = C_large(:);
    
    % 3. 构造约束 (保留完整对角线约束)
    % 估计非零元数量：对角线 N 个，等价类 n*(n-1)/2 对
    estimated_nnz = N + 4 * n * (n-1);
    
    % 这里我们直接构造 At，维度为 N^2 x m
    % 这样 ManiCSDP 内部执行 A = At' 后，A 就是 m x N^2，匹配 x (N^2 x 1)
    rows = zeros(estimated_nnz, 1);
    cols = zeros(estimated_nnz, 1);
    vals = zeros(estimated_nnz, 1);
    
    m_cnt = 0; % 约束计数器
    nnz_ptr = 0;

    % --- A. 对角线约束: X(i,i) = 1 ---
    fprintf('正在生成对角线约束...\n');
    for i = 1:N
        m_cnt = m_cnt + 1;
        nnz_ptr = nnz_ptr + 1;
        rows(nnz_ptr) = (i-1)*N + i; % 矩阵线性索引 X(i,i)
        cols(nnz_ptr) = m_cnt;
        vals(nnz_ptr) = 1;
    end
    b_diag = ones(N, 1);

    % --- B. 等价类约束: X(i+1, j+1) = X(n+j+1, n+i+1) ---
    fprintf('正在生成等价类约束 (n=%d)...\n', n);
    for i = 1:n
        for j = (i+1):n
            % 物理位置索引
            idx1 = (j)*N + (i+1);      % (i+1, j+1)
            idx1_s = (i)*N + (j+1);    % (j+1, i+1) 对称位
            idx2 = (n+i)*N + (n+j+1);  % (n+j+1, n+i+1)
            idx2_s = (n+j)*N + (n+i+1);% (n+i+1, n+j+1) 对称位
            
            % 实部相等约束: 0.5(X_idx1 + X_idx1_s) - 0.5(X_idx2 + X_idx2_s) = 0
            m_cnt = m_cnt + 1;
            curr_idxs = [idx1; idx1_s; idx2; idx2_s];
            rows(nnz_ptr+1:nnz_ptr+4) = curr_idxs;
            cols(nnz_ptr+1:nnz_ptr+4) = m_cnt;
            vals(nnz_ptr+1:nnz_ptr+4) = [0.5; 0.5; -0.5; -0.5];
            nnz_ptr = nnz_ptr + 4;
            
            % 虚部相等约束: -0.5i(X_idx1 - X_idx1_s) - (-0.5i)(X_idx2 - X_idx2_s) = 0
            m_cnt = m_cnt + 1;
            rows(nnz_ptr+1:nnz_ptr+4) = curr_idxs;
            cols(nnz_ptr+1:nnz_ptr+4) = m_cnt;
            vals(nnz_ptr+1:nnz_ptr+4) = [-0.5i; 0.5i; 0.5i; -0.5i];
            nnz_ptr = nnz_ptr + 4;
        end
    end
    
    % 构造最终的 At (N^2 x m)
    At = sparse(rows(1:nnz_ptr), cols(1:nnz_ptr), vals(1:nnz_ptr), N^2, m_cnt);
    b = [b_diag; zeros(m_cnt - N, 1)];
    
    fprintf('约束生成完毕。At 维度: %d x %d (N^2 x m)\n', size(At,1), size(At,2));
end