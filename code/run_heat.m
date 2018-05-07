x = linspace(0, 3*pi, 40);
y = linspace(0, 2*pi, 25);
t = linspace(0, 1.6, 17);
bc = {zeros(1, length(y)), zeros(1, length(y)), sin(x), sin(x)};
ic = sin(x)'*cos(y);

U = heat_cvx(x, y, t, bc, ic);

figure
[xx, yy] = meshgrid(x, y);
s = surf(xx, yy, U(:,:,1)');

light               % add a light
lighting gouraud    % preferred lighting for a curved surface
axis equal off      % set axis equal and remove axis
view(40,30)         % set viewpoint
camzoom(1.5)        % zoom into scene

M(length(t)) = struct('cdata', [], 'colormap', []);
for i = 1:length(t)
    s.ZData = U(:,:,i)';
    M(i) = getframe;
end

movie(M, 1, 8)
