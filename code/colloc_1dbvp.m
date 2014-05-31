clear all; close all; clc

% Solves u'' = 1 + e^(2x), u(0) = 0 = u(1)

epsilon = 350;
K   = @(x,center) ( exp(-epsilon.*((x-center).^2)) );
D1K = @(x,center) ( -2.*epsilon.*(x-center).*K(x,center) );
D2K = @(x,center) ( 2.*epsilon.*(2.*epsilon.*((x-center).^2)-1).* ...
                    K(x,center) );

rhs = @(x) ( 1 + exp(2.*x) );
u_analytic = @(x) ( 0.25.*((2.*x.^2)-exp(2).*x-x+exp(2.*x)-1) );

pts = linspace(0,1);
subplot(1,2,1);
hold on;
colors = 'grcmyk';

%% Calculate numerical solution
for N=[4 8 12 16 20 24];
    colloc_pts = linspace(0,1,N);
    
    % translated basis
    tmp = repmat(colloc_pts, N, 1);
    colloc_mat = [D2K(tmp',tmp)   zeros(N,2); 
                  K(0, colloc_pts) 1 0;  % we need 1 + x term
                  K(1, colloc_pts) 1 1]; % for BCs
    coef = colloc_mat\[rhs(colloc_pts)';0;0];
    
    % newton basis
    lambda = zeros(N,1);
    for i=2:N;
        lambda(i) = rhs(colloc_pts(i)) - lambda(i-1);
    end
    

    u_numeric = @(x) ( [K(x,colloc_pts) 1 x]*coef );
    plot(pts, arrayfun(u_numeric, pts), colors(1));
    colors = colors(2:end);
end

plot(pts, u_analytic(pts));
title('Solution to u\prime\prime = 1+exp(2x), u(0) = 0 = u(1)');
legend('4 elems', ...
       '8 elems', ...
       '12 elems', ...
       '16 elems', ...
       '20 elems', ...
       '24 elems', ...
       'analytic');

%% Calculate condition numbers of collocation matrices
trans_cond = zeros(100);
newt_cond = zeros(100);
subplot(1,2,2);
for N=1:100;
    colloc_pts = linspace(0,1,N);
    tmp = repmat(colloc_pts, N, 1);
    trans_cond(N) = cond([D2K(tmp',tmp)  zeros(N,2); 
                        K(0, colloc_pts) 1 0;
                        K(1, colloc_pts) 1 1]);
end
semilogy(trans_cond);
title('Condition Number of collocation matrices for N points');
legend('translated basis');
xlabel('N');
