clear all; close all; clc
format compact %remove blank lines from output
warning('off','MATLAB:nearlySingularMatrix'); % suppress cond. warnings

% Solves u'' = 1 + e^(2x), u(0) = 0 = u(1)

rhs = @(x) ( 1 + exp(2.*x) );
u_analytic = @(x) ( 0.25.*((2.*x.^2)-exp(2).*x-x+exp(2.*x)-1) );

pts = linspace(0,1);

default_tol_mult = 10;
tol_mult = default_tol_mult;
maxitr = 25;

Ns = ceil(1.4.^[1:17 20]);
Ns(Ns==218) = 217; %218 is magic in a bad way
zs_used = zeros(1,numel(Ns));
num_Ns=numel(Ns);

%% Calculate condition numbers of collocation matrices and
%  calculate error of numerical solutions

trans_err_cond = nan(2,num_Ns);
newt_err_cond  = nan(2,num_Ns);
newt2_err_cond = nan(2,num_Ns);
newt3_err_cond = nan(2,num_Ns);
newt4_err_cond = nan(2,num_Ns);
figure(1);

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

    %% usual basis
    colloc_mat = [D2KM(2:end-1,:);
                  K(0,colloc_pts);
                  K(1, colloc_pts)];
    coef = colloc_mat\[rhs(colloc_pts(2:end-1))';0;0;];

    trans_err_cond(1,i) = norm(([KM_evals]*coef)-u_analytic(pts)',Inf);
    trans_err_cond(2,i) = cond(colloc_mat);

    %% Creating a Newton basis for span{ K(\cdot, x_1), ... }
    [B,V] = calculate_beta_v(KM);
    D2V = B\D2KM;
    colloc_mat = [D2V(:,2:end-1)';
                  V(:,1)';
                  V(:,N)'];
    coef = colloc_mat\[rhs(colloc_pts(2:end-1))';0;0;];

    newt_err_cond(1,i) = norm(([(B\KM_evals')']*coef) ...
                              -u_analytic(pts)',Inf);
    newt_err_cond(2,i) = cond(colloc_mat);

    %% Creating a Newton basis for span{ LK(\cdot, x_1), ... }
    [B, D2V] = calculate_beta_v(D2KM);
    V = B\KM;
    colloc_mat = [D2V(:,2:end-1)';
                  V(:,1)';
                  V(:,N)'];
    coef = colloc_mat\[rhs(colloc_pts(2:end-1))';0;0;];

    newt2_err_cond(1,i) = norm(([(B\KM_evals')']*coef) ...
                               -u_analytic(pts)',Inf);
    newt2_err_cond(2,i) = cond(colloc_mat);

    %% Adaptive strategy
    tol_too_strict = true;
    tol_mult = default_tol_mult;
    itr = 1;
    while tol_too_strict && itr<maxitr
        tol_too_strict = false;
        [B, zminds] = calculate_newton_basis(KM,tol_mult);
        V = B';
        D2V = B\D2KM;
        colloc_mat = [D2V(:,2:end-1)';
                      V(:,1)';
                      V(:,N)'];
        coef = colloc_mat\[rhs(colloc_pts(2:end-1))';0;0;];

        newt3_err_cond(1,i) = norm(([(B\KM_evals')']*coef) ...
                                   -u_analytic(pts)',Inf);
        try
            newt3_err_cond(2,i) = cond(colloc_mat);
        catch exception
            disp('could not calculate cond');
            tol_too_strict = true;
            tol_mult = tol_mult + 1;
        end
        
        zs_used(i) = numel(zminds);
        itr = itr + 1;
    end
end    

subplot(1,3,1);  
loglog(Ns, trans_err_cond(1,:), 'b*-');
hold on;
loglog(Ns, newt_err_cond(1,:), 'go-');
loglog(Ns, newt2_err_cond(1,:), 'ro-');
loglog(Ns, newt3_err_cond(1,:), 'm+-');

title('Maximum error on 100 evenly spaced pts, when \epsilon_n=n/8');
legend('Usual basis',  ...
       'overconstrained Newton basis', ...
       'differentiated Newton basis', ...
       '2011 Newton basis')
ylabel('Error');
xlabel('# of collocation points');

subplot(1,3,2);
semilogy(Ns, trans_err_cond(2,:), 'b*-');
hold on;
loglog(Ns, newt_err_cond(2,:), 'go-');
loglog(Ns, newt2_err_cond(2,:), 'ro-');
loglog(Ns, newt3_err_cond(2,:), 'm+-');
title('condition number of collocation matrices for N points');
ylabel('condition number');
xlabel('N');

subplot(1,3,3)
plot(Ns, zs_used./Ns,'m+-');
title(['percentage of points used in adaptive strategy, tol. mult. = ' num2str(tol_mult)]);
