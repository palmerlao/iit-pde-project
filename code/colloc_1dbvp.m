clear all; close all; clc
format compact %remove blank lines from output

% Solves u'' = 1 + e^(2x), u(0) = 0 = u(1)

rhs = @(x) ( 1 + exp(2.*x) );
u_analytic = @(x) ( 0.25.*((2.*x.^2)-exp(2).*x-x+exp(2.*x)-1) );

pts = linspace(0,1);

Ns = ceil(1.4.^(1:9));
num_Ns=numel(Ns);

%% Calculate condition numbers of collocation matrices and
%  calculate error of numerical solutions

trans_err_cond = zeros(2,num_Ns);
newt_err_cond  = zeros(2,num_Ns);
newt2_err_cond = zeros(2,num_Ns);
newt3_err_cond = zeros(2,num_Ns);

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
    colloc_mat = [D2KM(2:end-1,:);
                  K(0,colloc_pts);
                  K(1, colloc_pts)];
    coef = colloc_mat\[rhs(colloc_pts(2:end-1))';0;0;];

    trans_err_cond(1,i) = norm(([KM_evals]*coef)-u_analytic(pts)',Inf);
    trans_err_cond(2,i) = cond(colloc_mat);

    % Creating a Newton basis for span{ K(\cdot, x_1), ... }
    %    [B, V] = calculate_beta_v(KM);
    [B,V] = calculate_beta_v(KM);
    D2V = B\D2KM; % maybe bad if D2KM is ill-cond.
    colloc_mat = [D2V(:,2:end-1)';
                  V(:,1)';
                  V(:,N)'];
    coef = colloc_mat\[rhs(colloc_pts(2:end-1))';0;0;];

    newt_err_cond(1,i) = norm(([(B\KM_evals')']*coef) ...
                              -u_analytic(pts)',Inf);
    if any(isnan(colloc_mat))
        dbstop
    end
    
    newt_err_cond(2,i) = cond(colloc_mat);

    % Creating a Newton basis for span{ LK(\cdot, x_1), ... }        
    [B, D2V] = calculate_beta_v(D2KM);
    V = B\KM;
    colloc_mat = [D2V(:,2:end-1)';
                  V(:,1)';
                  V(:,N)'];
    coef = colloc_mat\[rhs(colloc_pts(2:end-1))';0;0;];

    newt2_err_cond(1,i) = norm(([(B\KM_evals')']*coef) ...
                               -u_analytic(pts)',Inf);
    newt2_err_cond(2,i) = cond(colloc_mat);

    V = calculate_newton_basis(KM)';
    B = V';
    D2V = B\D2KM; % maybe bad if D2KM is ill-cond.
    colloc_mat = [D2V(:,2:end-1)';
                  V(:,1)';
                  V(:,N)'];
    coef = colloc_mat\[rhs(colloc_pts(2:end-1))';0;0;];

    newt3_err_cond(1,i) = norm(([(B\KM_evals')']*coef) ...
                              -u_analytic(pts)',Inf);
    newt3_err_cond(2,i) = cond(colloc_mat);
    
end    

subplot(1,2,1);  
loglog(Ns, trans_err_cond(1,:), 'b*-');
hold on;
loglog(Ns, newt_err_cond(1,:), 'go-');
loglog(Ns, newt2_err_cond(1,:), 'r+-');
loglog(Ns, newt3_err_cond(1,:), 'yd-');

title('Maximum error on 100 evenly spaced pts, when \epsilon_n=n^2/16');
legend('Usual basis',  ...
       'Newton basis', ...
       'differentiated kernel Newton basis', ...
       '2011 Newton basis');
ylabel('Error');
xlabel('# of collocation points');


subplot(1,2,2);
semilogy(Ns, trans_err_cond(2,:), 'b*-');
hold on;
semilogy(Ns, newt_err_cond(2,:), 'go-');
semilogy(Ns, newt2_err_cond(2,:), 'r+-');
semilogy(Ns, newt3_err_cond(1,:), 'yd-');
title('condition number of collocation matrices for N points');
ylabel('condition number');
xlabel('N');
