clear all, close all, clc;

global N D2KM KM KM_evals

init_cond = @(x) (2.*sin(pi.*x./2) - sin(pi.*x) + 4.*sin(2.*pi.*x));
u_analytic = @(x,t) ( 2.*sin(pi.*x./2).*exp(-(pi/4).^2.*t) ...
                      - sin(pi.*x).*exp(-(pi./2).^2.*t) ...
                      - sin(2.*pi.*x).*exp(-pi.^2.*t));

N = 5;
epsilon = 0.5;
K   = @(x,center) ( exp(-epsilon.*((x-center).^2)) );
D1K = @(x,center) ( -2.*epsilon.*(x-center).*K(x,center) );
D2K = @(x,center) ( 2.*epsilon.*(2.*epsilon.*((x-center).^2)-1).* ...
                    K(x,center) );

colloc_pts = linspace(0,2,N);
pts = linspace(0,2);
tmp = repmat(colloc_pts,N,1);
KM = K(tmp',tmp);
D2KM = D2K(tmp',tmp);
KM_evals = K( repmat(pts',1,size(colloc_pts,2)), repmat(colloc_pts,size(pts,2),1));

ic = init_cond(colloc_pts);
disc_time = [0:0.007:5];

[t,u_numeric] = ode23(@mol_1dheat_helper,disc_time,ic);
size(u_numeric)
size(colloc_pts)

%% plot movies of analytic soln
figure;
pts = linspace(0,2);
title('solution to 4u_t = u_xx');
for i=1:length(disc_time)
    plot(pts,u_analytic(pts,disc_time(i)),colloc_pts,u_numeric(i,:),'gd:');
    axis([0 2 0 5]);
    M(i) = getframe;
end

movie(M,1,30);
