function U = poisson_cvx(x, y, bc, f)
    m = length(x);
    n = length(y);
    dx = (x(end)-x(1)) / (m-1);
    dy = (y(end)-y(1)) / (n-1);

    e = ones(m*n, 1);
    A = spdiags([e/dy^2 e/dx^2 -2*e*(dx^-2+dy^-2) e/dx^2 e/dy^2], ...
        [-m -1 0 1 m], m*n, m*n);           % boundaries don't matter
    F = f(x', y);

    cvx_begin
        variable U(m, n);
        variable R(m, n);
        minimize norm(reshape(R(2:end-1, 2:end-1), [], 1));
        subject to
            U(1, :) == bc{1};
            U(end, :) == bc{2};
            U(2:end-1, 1) == bc{3}(2:end-1)';
            U(2:end-1, end) == bc{4}(2:end-1)';
            R == reshape(A*U(:), m, n) - F;
    cvx_end
end
