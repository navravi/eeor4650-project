x = linspace(0, 3*pi, 40);
y = linspace(0, 2*pi, 25);
bc = {zeros(1, length(y)), zeros(1, length(y)), sin(x), sin(x)};
f = @(x, y) -2*sin(x)*cos(y);
U_true = sin(x)'*cos(y);

U_linsys = poisson_linsys(x, y, bc, f);
err_linsys = max(abs(U_linsys(:)-U_true(:)));

U_cvx = poisson_cvx(x, y, bc, f);
err_cvx = max(abs(U_cvx(:)-U_true(:)));

figure
[xx, yy] = meshgrid(x, y);
s = surf(xx, yy, U_cvx');

light               % add a light
lighting gouraud    % preferred lighting for a curved surface
axis equal off      % set axis equal and remove axis
view(40,30)         % set viewpoint
camzoom(1.5)        % zoom into scene
