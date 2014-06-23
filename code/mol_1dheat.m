clear all, close all, clc;

global N KM KM_evals D2KM D2KM_evals

u_analytic = @(x,t) ( 2.*sin(pi.*x./2).*exp(-(pi/4).^2.*t) ...
                      - sin(pi.*x).*exp(-(pi./2).^2.*t) ...
                      - sin(2.*pi.*x).*exp(-(pi).^2.*t));

N = 15;
epsilon = 0.25;
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
D2KM_evals = D2K( repmat(pts',1,size(colloc_pts,2)),repmat(colloc_pts,size(pts,2),1));
D = D2KM/KM;

ic = u_analytic(colloc_pts,0);
dt = 0.01;
disc_time = 0:dt:3;
u_numeric = zeros(length(disc_time),length(ic));
u_numeric(1,:) = ic;

% $$$ for i=2:length(disc_time)
% $$$     u_numeric(i,:) = u_numeric(i-1,:) + (0.25.*dt.*D*u_numeric(i-1,:)')';
% $$$     u_numeric(i,[1,end]) = [0 0];
% $$$ end

[t,u_numeric] = ode45(@(t,y) (0.25.*D*y).*[0;ones(N-2,1);0],disc_time,ic);

%% plot movies of analytic soln
figure;
pts = linspace(0,2);
title('solution to 4u_t = u_xx');
for i=1:length(disc_time)
    plot(pts,u_analytic(pts,disc_time(i)),colloc_pts,u_numeric(i,:),'gd--');
    axis([0 2 -1 5]);
    M(i) = getframe;
end

movie(M,1,30);
