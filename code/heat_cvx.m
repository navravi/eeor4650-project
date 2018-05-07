function U = heat_cvx(x, y, t, bc, ic)
    m = length(x);
    n = length(y);
    k = length(t);
    dx = (x(end)-x(1))/(m-1);
    dy = (y(end)-y(1))/(n-1);
    dt = (t(end)-t(1))/(k-1);

    e = ones(m*n, 1);
    A = spdiags([e/dy^2 e/dx^2 -2*e*(dx^-2+dy^-2) e/dx^2 e/dy^2], ...
        [-m -1 0 1 m], m*n, m*n);           % boundaries don't matter

    cvx_begin
        variable U(m,n,k);
        variable F1(m,n,k-1);
        variable F2(m,n,k-1);
        variable R(m,n,k-1);
        minimize norm(reshape(R(2:end-1,2:end-1,:), [], 1));
        subject to
            U(1,:,2:end) == repmat(bc{1}, 1, 1, k-1);
            U(end,:,2:end) == repmat(bc{2}, 1, 1, k-1);
            U(2:end-1,1,2:end) == repmat(bc{3}(2:end-1)', 1, 1, k-1);
            U(2:end-1,end,2:end) == repmat(bc{4}(2:end-1)', 1, 1, k-1);
            U(:,:,1) == ic;
            F1 == reshape(A*reshape(U(:,:,1:end-1), m*n, k-1), m, n, k-1);
            F2 == reshape(A*reshape(U(:,:,2:end), m*n, k-1), m, n, k-1);
            R == 0.5*dt*(F1+F2) - (U(:,:,2:end)-U(:,:,1:end-1));
    cvx_end
end
