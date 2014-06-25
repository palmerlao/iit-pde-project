clear all; close all; clc
format compact %remove blank lines from output

% Lu = f on \Omega = (0,1)x(0,1), where L is the Laplacian
% u = g on bdry
% g also happens to be the analytic soln.
f = @(x) (-5.*pi.^2.*sin(pi.*x(:,1)).*cos(2.*pi.*x(:,2)));
g = @(x) (sin(pi.*x(:,1)).*cos(2.*pi.*x(:,2)));

Ns = [3 4 5];
num_Ns=numel(Ns);

for i=1:num_Ns
    N=Ns(i);

    epsilon = 15;
    K   = @(x,center) ( exp(-epsilon.*(norm(x-center).^2)) );
    D1K = @(x,center,ind) ( -2.*epsilon.*(x(ind)-center(ind)).*K(x(ind),center(ind)) );
    D2K = @(x,center,ind) ( 2.*epsilon.*(2.*epsilon.*((x(ind)-center(ind)).^2)-1).* ...
                        K(x(ind),center(ind)) );

    [ex, ey] = meshgrid(linspace(0,1,10));
    pts = [ex(:) ey(:)];
    mesh = zeros(N,N,2);
    [mesh(:,:,1), mesh(:,:,2)] = meshgrid(linspace(0,1,N));
    interior_pts  = reshape(mesh(2:end-1,2:end-1,:),(N-2).^2,2,1);
    dirichlet_pts = [ reshape(mesh(:,1,:),N,2,1);
                      reshape(mesh(N,2:end-1,:),N-2,2,1);
                      reshape(mesh(1,2:end-1,:),N-2,2,1);
                      reshape(mesh(:,N,:),N,2,1); ];
    dirichlet_pts_ind = [1:N N.*(2:N-1) N.*(1:N-2)+1 (N.*(N-1)+1):N.^2];
    interior_pts_ind = reshape(1:N.^2,N,N);
    interior_pts_ind = interior_pts_ind(2:end-1,2:end-1);
    colloc_pts = reshape(mesh,N.^2,2,1);

    KM  = zeros(N.^2,N.^2);
    LKM = zeros(N.^2,N.^2);
    EKM = zeros(size(pts,1),N.^2);

    % there's probably a better way to do this
    for j=1:N.^2
        for k=1:N.^2
            KM(j,k)  = K(colloc_pts(j,:), colloc_pts(k,:));
            LKM(j,k) = D2K(colloc_pts(j,:), colloc_pts(k,:),1) + D2K(colloc_pts(j,:), colloc_pts(k,:),2);
        end
    end
    for j=1:size(pts,1)
        for k=1:N.^2
            EKM(j,k) = K(pts(j,:), colloc_pts(k,:));
        end
    end

    rhs = [f(colloc_pts(interior_pts_ind,:));
           g(colloc_pts(dirichlet_pts_ind,:))];
    % usual basis
    colloc_mat = [LKM(interior_pts_ind,:);
                  KM(dirichlet_pts_ind,:)];
    coef = colloc_mat\rhs;
    
    subplot(2,2,i);
    meshc(ex,ey,reshape(EKM*coef,10,10));

    trans_err_cond(1,i) = norm((EKM*coef)-g(pts),Inf);
    trans_err_cond(2,i) = cond(colloc_mat);

    % Creating a Newton basis for span{ K(\cdot, x_1), ... }
    [B, VM] = calculate_beta_v(KM);
    LVM = B\LKM; % maybe bad if D2KM is ill-cond.
    colloc_mat = [LVM(:,interior_pts_ind)'
                  VM(:,dirichlet_pts_ind)'];
    coef = colloc_mat\rhs;

    newt_err_cond(1,i) = norm(([(B\EKM')']*coef)-g(pts),Inf);
    newt_err_cond(2,i) = cond(colloc_mat);

    % Creating a Newton basis for span{ LK(\cdot, x_1), ... }        
    [B, LVM] = calculate_beta_v(LKM);
    VM = B\KM;
    colloc_mat = [LVM(:,interior_pts_ind)';
                  VM(:,dirichlet_pts_ind)'];
    coef = colloc_mat\rhs;

    newt2_err_cond(1,i) = norm(([(B\EKM')']*coef)-g(pts),Inf);
    newt2_err_cond(2,i) = cond(colloc_mat);
    
end

subplot(2,2,4);
meshc(ex,ey,reshape(g(pts),10,10));

figure
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
