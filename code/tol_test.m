clear all; close all; clc
format compact %remove blank lines from output
warning('off','MATLAB:nearlySingularMatrix'); % suppress cond. warnings

% Solves u'' = 1 + e^(2x), u(0) = 0 = u(1)

rhs = @(x) ( 1 + exp(2.*x) );
u_analytic = @(x) ( 0.25.*((2.*x.^2)-exp(2).*x-x+exp(2.*x)-1) );

pts = linspace(0,1);

Ntols = 10;
tol_mult = linspace(20,100,Ntols);

Ns = ceil(1.4.^[1:17 20]);
Ns(Ns==218) = 217; %218 is magic in a bad way
zs_used = zeros(Ntols,numel(Ns));
num_Ns=numel(Ns);

close all;
cs = colormap;
csn = round(linspace(1,64,Ntols));

%% Calculate condition numbers of collocation matrices and
%  calculate error of numerical solutions

trans_err_cond = nan(2,num_Ns);
newt_err  = nan(Ntols,num_Ns);
newt_cond = nan(Ntols,num_Ns);
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

    %% Adaptive strategy
    for j=1:Ntols
        [B, zminds] = calculate_newton_basis(KM,tol_mult(j));
        V = B';
        D2V = B\D2KM;
        colloc_mat = [D2V(:,2:end-1)';
                      V(:,1)';
                      V(:,N)'];
        coef = colloc_mat\[rhs(colloc_pts(2:end-1))';0;0;];

        newt_err(j,i) = norm(([(B\KM_evals')']*coef) ...
                                  -u_analytic(pts)',Inf);
        try
            newt_cond(j,i) = cond(colloc_mat);
        catch
            %qdsfaffasd
        end
        zs_used(j,i) = numel(zminds);
    end
end    


loglog(Ns, trans_err_cond(1,:), 'b*-');
hold on;
for i=1:Ntols
    loglog(Ns,newt_err(i,:),'Color',cs(csn(i),:),'Marker','d');
end
title('Maximum error on 100 evenly spaced pts, when \epsilon_n=n^2/64');
ylabel('Error');
xlabel('# of collocation points');
colorbar


figure;
semilogy(Ns, trans_err_cond(2,:), 'b*-');
hold on;
for i=1:Ntols
    semilogy(Ns,newt_cond(i,:),'Color',cs(csn(i),:),'Marker','d');
end
title('condition number of collocation matrices for N points');
ylabel('condition number');
xlabel('N');
colorbar

figure;
hold on;
for i=1:Ntols
    plot(Ns,zs_used(i,:)./Ns,'Color',cs(csn(i),:),'Marker','d');
end
title(['percentage of points used in adaptive strategy, tol. mult.='  num2str(tol_mult)]);
colorbar
