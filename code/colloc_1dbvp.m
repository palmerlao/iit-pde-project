function [] = colloc_1dbvp()
    clear all; close all; clc

    % Solves u'' = 1 + e^(2x), u(0) = 0 = u(1)

    epsilon = 400;
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

    %% Calculate and plot numerical solution
    for N=[12 20 28];
        colloc_pts = linspace(0,1,N);

        tmp = repmat(colloc_pts, N, 1);
        KM = K(tmp',tmp);
        D2KM = D2K(tmp',tmp);
        
        % translated basis
        colloc_mat = [D2KM             zeros(N,2); 
                      K(0, colloc_pts) 1 0;  % we need 1 + x term
                      K(1, colloc_pts) 1 1]; % for BCs
        coef = colloc_mat\[rhs(colloc_pts)';0;0];
        u_numeric = @(x) ( [K(x,colloc_pts) 1 x]*coef );
        plot(pts, arrayfun(u_numeric, pts), colors(1));
        
        [B, V] = calculate_beta_v(KM, N, colloc_pts, K);
        D2V_ = B\D2KM; % maybe bad if D2KM is ill-cond.
        colloc_mat = [D2V_'   zeros(N,2);
                      V(:,1)' 1 0;
                      V(:,2)' 1 1];
        coef = colloc_mat\[rhs(colloc_pts)';0;0];

        plot(colloc_pts, coef'*[V;ones(1,N);colloc_pts], [colors(2) '*:']);

        colors = colors(3:end);
    end

    plot(pts, u_analytic(pts));
    title('Solution to u\prime\prime = 1+exp(2x), u(0) = 0 = u(1)');
    legend('12, translate', ...
           '12, Newton', ...
           '20, translate', ...
           '20, Newton', ...
           '28, translate', ...
           '28, Newton', ...
           'analytic');

    %% Calculate condition numbers of collocation matrices
    trans_cond = zeros(1,100);
    newt_cond = zeros(1,100);
    for N=2:100;
        colloc_pts = linspace(0,1,N);
        tmp = repmat(colloc_pts, N, 1);
        KM = K(tmp',tmp);
        D2KM = D2K(tmp',tmp);
        
        trans_cond(N) = cond([D2KM  zeros(N,2); 
                              K(0, colloc_pts) 1 0;
                              K(1, colloc_pts) 1 1]);
        [B, V] = calculate_beta_v(KM, N, colloc_pts, K);
        D2V_ = B\D2KM; % maybe bad if D2KM is ill-cond.
        newt_cond(N) = cond([D2V_'   zeros(N,2);
                             V(:,1)' 1 0;
                             V(:,2)' 1 1]);
    end

    subplot(1,2,2);
    hold on;
    semilogy(log(trans_cond), 'b*-');
    semilogy(log(newt_cond), 'go-');
    title('log log of condition number of collocation matrices for N points');
    legend('translated basis', 'Newton basis');
    ylabel('log of log of condition number');
    xlabel('N');
end

function [B, V] = calculate_beta_v(KM, N, xs, K)
% here we alternately calculate betas and vs since they depend on each
% other (see paper for relationship). B should be lower triangular,
% and V should be unit upper triangular s.t. KM = B*V
    B = zeros(N,N);
    V = eye(N);
    for c=1:(N-1)
        for i=c:N
            B(i,c) = calculate_single_beta(B,V,i,c,K,xs);
        end
        % not sure if this is bad when KM is ill-cond.
        V(1:c,c+1) = B(1:c,1:c)\KM(1:c,c+1); 
    end
    B(N,N) = calculate_single_beta(B,V,N,N,K,xs);
end

function [res] = calculate_single_beta(B, V, i, j, K, xs)
    res = K(xs(j), xs(i));
    for k=1:j-1
        res = res - B(i,k).*V(k,j);
    end
end
