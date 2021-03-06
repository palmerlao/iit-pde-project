clear all, close all, clc;

u_analytic = @(x,t) ( 2.*sin(pi.*x./2).*exp(-(pi/4).^2.*t) ...
                      - sin(pi.*x).*exp(-(pi./2).^2.*t) ...
                      - sin(2.*pi.*x).*exp(-(pi).^2.*t));

N = 25;
epsilon = (N/8);
K   = @(x,center) ( exp(-epsilon.*((x-center).^2)) );
D1K = @(x,center) ( -2.*epsilon.*(x-center).*K(x,center) );
D2K = @(x,center) ( 2.*epsilon.*(2.*epsilon.*((x-center).^2)-1).* ...
                    K(x,center) );

colloc_pts = linspace(0,2,N);
pts = linspace(0,2);
tmp = repmat(colloc_pts,N,1);

KM = K(tmp',tmp);
D2KM = D2K(tmp',tmp);

[B, zminds, was_inc] = calculate_newton_basis(KM,10);
B = B(zminds,:);
VM = B';
D2VM = B\D2KM(zminds,zminds);

D_usual = D2KM/KM;
D_newton = D2VM/VM;

ic = u_analytic(colloc_pts,0);
dt = 0.01;
disc_time = 0:dt:3;
u_usual = zeros(length(disc_time),length(ic));
u_usual(1,:) = ic;
u_newton = zeros(length(disc_time),length(ic));
u_newton(1,:) = ic;
u_newton = u_newton(:,was_inc(zminds));

% explicit euler
for i=2:length(disc_time)
    u_usual(i,:) = u_usual(i-1,:) + (0.25.*dt.*D_usual*u_usual(i-1,:)')';
    u_usual(i,[1,end]) = [0 0];

    u_newton(i,:) = u_newton(i-1,:) + (0.25.*dt.*D_newton*u_newton(i-1,:)')';
    u_newton(i,[1,end]) = [0 0];

end

%[t,u_usual] = ode45(@(t,y) (0.25.*D_usual*y).*[0;ones(N-2,1);0],disc_time,ic);

%% plot movies of analytic soln
figure;

for i=1:length(disc_time)
    plot(pts,u_analytic(pts,disc_time(i)), ...
         colloc_pts,u_usual(i,:),'gd:', ...
         colloc_pts(was_inc),u_newton(i,:),'rd');
    axis([0 2 -1 5]);
    M(i) = getframe;
end

movie(M,1,30);
