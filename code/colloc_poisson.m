clear all; close all; clc;
format compact %remove blank lines from output
warning('off','MATLAB:nearlySingularMatrix'); % suppress cond. warnings

% Lu = f on \Omega = (0,1)x(0,1), where L is the Laplacian
% u = g on bdry
% g also happens to be the analytic soln.
f = @(x) (-5.*pi.^2.*sin(pi.*x(:,1)).*cos(2.*pi.*x(:,2)));
g = @(x) (sin(pi.*x(:,1)).*cos(2.*pi.*x(:,2)));

Ns = [3:2:25];
eN = 100;
num_Ns=numel(Ns);

tol_mult = 30;

trans_err_cond = zeros(2,num_Ns);
newt_err_cond  = zeros(2,num_Ns);
newt2_err_cond = zeros(2,num_Ns);
newt3_err_cond = zeros(2,num_Ns);

for i=1:num_Ns

    N = Ns(i);
    epsilon = (N./8);
    G = @(r) ( exp(-(epsilon.*r).^2) );
    LG = @(r) ( (2.*(epsilon.^2).*G(r)).*((2.*(epsilon.*r).^2)-1) + (-2.*(epsilon.^2).*G(r)) );

    [cx,cy] = meshgrid(linspace(0,1,N));
    cx=cx(:);
    cy=cy(:);

    [ex,ey] = meshgrid(linspace(0,1,eN));
    ex=ex(:);
    ey=ey(:);

    tmp = repmat(cx,1,N.^2);
    cdiff_1 = tmp'-tmp;
    tmp = repmat(cy,1,N.^2);
    cdiff_2 = tmp'-tmp;
    crad_mat = sqrt(cdiff_1.^2+cdiff_2.^2);

    tmp  = repmat(ex,1,N.^2);
    tmp_ = repmat(cx',eN.^2,1);
    ediff_1 = tmp-tmp_;
    tmp  = repmat(ey,1,N.^2);
    tmp_ = repmat(cy',eN.^2,1);
    ediff_2 = tmp-tmp_;
    erad_mat = sqrt(ediff_1.^2+ediff_2.^2);
    
    KM  = G(crad_mat);
    EKM = G(erad_mat);
    LKM = LG(crad_mat);

    dirichlet_pts_ind = [1:N N.*(2:N-1) N.*(1:N-2)+1 (N.*(N-1)+1):N.^2];
    interior_pts_ind = reshape(1:N.^2,N,N);
    interior_pts_ind = interior_pts_ind(2:end-1,2:end-1);
    interior_pts_ind = interior_pts_ind(:);
    
    colloc_mat = [ LKM(interior_pts_ind,:);
                   KM(dirichlet_pts_ind,:)];
    rhs = [ f([cx(interior_pts_ind) cy(interior_pts_ind)]);
            g([cx(dirichlet_pts_ind) cy(dirichlet_pts_ind)]) ];
    coef = colloc_mat\rhs;
    
    trans_err_cond(1,i) = norm((EKM*coef)-g([ex ey]),Inf);
    trans_err_cond(2,i) = cond(colloc_mat);
    
    [B,VM] = calculate_beta_v(KM);
    LVM = B\LKM;
    colloc_mat = [LVM(:,interior_pts_ind)';
                  VM(:,dirichlet_pts_ind)'];
    coef = colloc_mat\rhs;

    newt_err_cond(1,i) = norm(([(B\EKM')']*coef)-g([ex ey]),Inf);
    newt_err_cond(2,i) = cond(colloc_mat);
    
    [B, LVM] = calculate_beta_v(LKM);
    VM = B\KM;
    colloc_mat = [LVM(:,interior_pts_ind)';
                  VM(:,dirichlet_pts_ind)'];
    coef = colloc_mat\rhs;

    newt2_err_cond(1,i) = norm(([(B\EKM')']*coef)-g([ex ey]),Inf);
    newt2_err_cond(2,i) = cond(colloc_mat);

    [B, zminds] = calculate_newton_basis(KM,tol_mult);
    VM = B';
    LVM = B\LKM;
    colloc_mat = [LVM(:,interior_pts_ind)';
                  VM(:,dirichlet_pts_ind)'];
    coef = colloc_mat\rhs;

    newt3_err_cond(1,i) = norm(([(B\EKM')']*coef)-g([ex ey]),Inf);
    if all(all(isnan(colloc_mat))) || all(all(isinf(colloc_mat)))
        newt3_err_cond(2,i) = NaN;
    else
        newt3_err_cond(2,i) = cond(colloc_mat);
    end
    zs_used(i) = numel(zminds);
    
    %    subplot(1,num_Ns+1,i);
    %    surfc(reshape(ex,eN,eN),reshape(ey,eN,eN),reshape(EKM*coef,eN,eN));
end
%subplot(1,num_Ns+1,num_Ns+1)
%surfc(reshape(ex,eN,eN),reshape(ey,eN,eN),reshape(g([ex ey]),eN,eN));
subplot(1,3,1);
loglog(Ns.^2,trans_err_cond(1,:),'b*-');
hold on;
loglog(Ns.^2, newt_err_cond(1,:), 'go-');
loglog(Ns.^2, newt2_err_cond(1,:), 'r+-');
loglog(Ns.^2, newt3_err_cond(1,:), 'md-');
title('Maximum error on 100 evenly spaced pts, when \epsilon_n=n/8');
legend('Usual basis',  ...
       'Newton basis', ...
       'differentiated kernel Newton basis', ...
       '2011 Newton basis');
ylabel('Error');
xlabel('# of collocation points');


subplot(1,3,2);
semilogy(Ns.^2,trans_err_cond(2,:),'b*-');
hold on;
semilogy(Ns.^2, newt_err_cond(2,:), 'go-');
semilogy(Ns.^2, newt2_err_cond(2,:), 'r+-');
semilogy(Ns.^2, newt3_err_cond(2,:), 'md-');
title('condition number of collocation matrices for N points');
ylabel('condition number');
xlabel('# of collocation points');

subplot(1,3,3)
plot(Ns.^2, zs_used./(Ns.^2),'b*-');
title(['percentage of points used in adaptive strategy, tol. mult. = ' num2str(tol_mult)]);
