clear all; close all; clc

% Solves u'' = 1 + e^(2x), u(0) = 0 = u(1)

epsilon = 50;
f   = @(x,center) ( exp(-epsilon.*((x-center).^2)) );
D1f = @(x,center) ( -2.*epsilon.*(x-center).*f(x,center) );
D2f = @(x,center) ( 2.*epsilon.*(2.*epsilon.*((x-center).^2)-1).* ...
                    f(x,center) );

rhs = @(x) ( 1 + exp(2.*x) );

hold on;
pts = linspace(0,1);
colors = 'rgkyc';
for N=[4 8 16 23];
    colloc_pts = linspace(0,1,N);
    tmp = repmat(colloc_pts, N, 1);
    colloc_mat = D2f(tmp',tmp);
    coef = colloc_mat\[0;rhs(colloc_pts(2:end-1))';0];

    u_numeric = @(x) ( f(x,colloc_pts)*coef );
    plot(pts, arrayfun(u_numeric, pts), colors(1));
    colors = colors(2:end);
end
u_analytic = @(x) ( 0.25.*((2.*pts.^2)-exp(2).*pts-pts+exp(2.*pts)-1) );
plot(pts, u_analytic(pts));
title('Solution to u'' = 1+exp(2x), u(0) = 0 = u(1)');
legend('4 elems', ...
       '8 elems', ...
       '16 elems', ...
       '23 elems', ...
       'analytic')
