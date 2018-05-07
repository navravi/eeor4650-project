function U = poisson_linsys(x, y, bc, f)
    % Determine mesh parameters
    m = length(x)-2;
    n = length(y)-2;
    dx = (x(end)-x(1)) / (m+1);
    dy = (y(end)-y(1)) / (n+1);

    % Construct sparse Laplacian operator matrix
    Im = speye(m);
    In = speye(n);
    Em = sparse(2:m, 1:m-1, 1, m, m);
    En = sparse(2:n, 1:n-1, 1, n, n);
    T = (Em+Em')/dx^2 - 2*Im*(dx^-2+dy^-2);
    S = En+En';
    D = Im/dy^2;
    A = kron(In, T) + kron(S, D);

    % Initialize solution with boundary conditions
    U = zeros(m+2, n+2);
    U(1, :) = bc{1};
    U(end, :) = bc{2};
    U(:, 1) = bc{3}';
    U(:, end) = bc{4}';

    % Construct source vector with boundary correction terms
    F = f(x(2:end-1)', y(2:end-1));
    F(1, :) = F(1, :) - bc{1}(2:end-1)/dx^2;
    F(end, :) = F(end, :) - bc{2}(2:end-1)/dx^2;
    F(:, 1) = F(:, 1) - bc{3}(2:end-1)'/dy^2;
    F(:, end) = F(:, end) - bc{4}(2:end-1)'/dy^2;

    % Solve system
    U(2:end-1, 2:end-1) = reshape(A \ reshape(F, m*n, 1), m, n);
end
