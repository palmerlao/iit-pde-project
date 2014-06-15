function [] = colloc_1dbvp()
    clear all; close all; clc
    format compact %remove blank lines from output

    % Solves u'' = 1 + e^(2x), u(0) = 0 = u(1)

    rhs = @(x) ( 1 + exp(2.*x) );
    u_analytic = @(x) ( 0.25.*((2.*x.^2)-exp(2).*x-x+exp(2.*x)-1) );

    pts = linspace(0,1);
    
    Ns = ceil(1.4.^(1:17));
    num_Ns=numel(Ns);
    
    %% Calculate condition numbers of collocation matrices and
    %  calculate error of numerical solutions

    trans_err_cond = zeros(2,num_Ns);
    newt_err_cond  = zeros(2,num_Ns);
    newt2_err_cond = zeros(2,num_Ns);
    
    for i=1:num_Ns;
       N=Ns(i);
        
        epsilon = (N/8).^2;
        K   = @(x,center) ( exp(-epsilon.*((x-center).^2)) );
        D1K = @(x,center) ( -2.*epsilon.*(x-center).*K(x,center) );
        D2K = @(x,center) ( 2.*epsilon.*(2.*epsilon.*((x-center).^2)-1).* ...
                            K(x,center) );

        colloc_pts = linspace(0,1,N);
        tmp = repmat(colloc_pts,N,1);
        KM = K(tmp',tmp);
        D2KM = D2K(tmp',tmp);
        KM_evals = K( repmat(pts',1,size(colloc_pts,2)), repmat(colloc_pts,size(pts,2),1));

        % usual basis
        colloc_mat = [D2KM             zeros(N,2); 
                      K(0, colloc_pts) 1 0;
                      K(1, colloc_pts) 1 1];
        coef = colloc_mat\[rhs(colloc_pts)';0;0;];

        trans_err_cond(1,i) = norm(([KM_evals ones(size(pts,2),1) pts']*coef) ...
                                          -u_analytic(pts)',Inf);
        trans_err_cond(2,i) = cond(colloc_mat);

        % Creating a Newton basis for span{ K(\cdot, x_1), ... }
        [B, V] = calculate_beta_v(KM, N, colloc_pts, K);
        D2V = B\D2KM; % maybe bad if D2KM is ill-cond.
        colloc_mat = [D2V'    zeros(N,2);
                      V(:,1)' 1 0;
                      V(:,N)' 1 1];
        coef = colloc_mat\[rhs(colloc_pts)';0;0;];

        newt_err_cond(1,i) = norm(([(B\KM_evals')' ones(size(pts,2),1) pts']*coef) ...
                                  -u_analytic(pts)',Inf);
        newt_err_cond(2,i) = cond(colloc_mat);

        % Creating a Newton basis for span{ LK(\cdot, x_1), ... }        
        [B, D2V] = calculate_beta_v(D2KM, N, colloc_pts, D2K);
        V = B\KM;
        colloc_mat = [D2V'    zeros(N,2);
                      V(:,1)' 1 0;
                      V(:,N)' 1 1];
        coef = colloc_mat\[rhs(colloc_pts)';0;0;];

        newt2_err_cond(1,i) = norm(([(B\KM_evals')' ones(size(pts,2),1) pts']*coef) ...
                                   -u_analytic(pts)',Inf);
        newt2_err_cond(2,i) = cond(colloc_mat);
        
    end    

    subplot(1,2,1);  
    loglog(Ns, trans_err_cond(1,:), 'b*-');
    hold on;  
    loglog(Ns, newt_err_cond(1,:), 'go-');
    loglog(Ns, newt2_err_cond(1,:), 'r+-');

    title('Maximum error on 100 evenly spaced pts, when \epsilon_n=n^2/16');
    legend('Usual basis',  ...
           'Newton basis', ...
           'differentiated kernel Newton basis');
    ylabel('Error');
    xlabel('# of collocation points');


    subplot(1,2,2);
    semilogy(Ns, trans_err_cond(2,:), 'b*-');
    hold on;
    semilogy(Ns, newt_err_cond(2,:), 'go-');
    semilogy(Ns, newt2_err_cond(2,:), 'r+-');
    title('condition number of collocation matrices for N points');
    ylabel('condition number');
    xlabel('N');
end

function [B, V] = chol_calc_beta_v(KM, N, xs, K)
    B = chol(KM);
    V = B';
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
