% Generate the second-order SOS-SDP relaxation for the binary quadratic program:
% Min  x'*Q*x + x'*e
% s.t. x_i^2 = 1, i = 1,...,n.
% Output sedumi format data.
% dAAt: the diagnal of A*A'

function [A, b, dAAt, mb] = bqpsos(Q, e, n)
sp = get_basis(n, 4);
ind = true(size(sp, 2), 1);
ind(sum(sp>1)> 0) = false;
sp = sp(:,ind);
mb = nchoosek(n+2, 2) - n;
lsp = size(sp, 2);
row = ones(mb^2, 1);
col = zeros(mb^2, 1);
val = ones(mb^2, 1);
dAAt = zeros(lsp, 1);
dAAt(1) = mb;
col(1:mb) = (0:(mb-1))*mb+(1:mb);
ind = mb + 1;
for i = 1:mb
    for j = i+1:mb
        bi = mod(sp(:,i) + sp(:,j), 2);
        locb = bfind(sp, lsp, bi, n);
        row(ind) = locb;
        row(ind+1) = locb;
        col(ind) = (i-1)*mb+j;
        col(ind+1) = (j-1)*mb+i;
        dAAt(locb) = dAAt(locb) + 2;
        ind = ind + 2;
    end
end
A = sparse(row,col,val,lsp,mb^2);

b = zeros(lsp, 1);
b(1) = trace(Q);
b(2:n+1) = e;
b(n+2:(n+1)*(n+2)/2-n) = 2*Q(triu(ones(n),1)~=0);
end
