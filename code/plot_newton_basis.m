clear all; close all; clc;

colloc_pts = [0.225995, 0.272585, 0.368313, 0.779455];
epsilon = 1;
N = 4;
K = @(x,center) ( exp(-epsilon.*((x-center).^2)) );
pts = linspace(-4,4);
cs = 'bgmy';

Kf1 = @(x) (exp(-(-colloc_pts(1)+x).^2));
Kf2 = @(x) (exp(-(-colloc_pts(2)+x).^2));
Kf3 = @(x) (exp(-(-colloc_pts(3)+x).^2));
Kf4 = @(x) (exp(-(-colloc_pts(4)+x).^2));

Nf1 = @(x) (Kf1(x));
Nf2 = @(x) (1.47751.*Kf4(x) - 1.08767.*Kf1(x));
Nf3 = @(x) (-3.30077.*Kf4(x) + 12.2844.*Kf3(x) - 9.60817.*Kf1(x));
Nf4 = @(x) (8.20995.*Kf4(x) - 154.831.*Kf3(x) +393.471.*Kf2(x) - 246.935.*Kf1(x));

hold on;
subplot(1,2,1);

plot(pts, Nf1(pts),cs(1), ...
     pts, Nf2(pts),cs(2), ...
     pts, Nf3(pts),cs(3), ...
     pts, Nf4(pts),cs(4) );
title('Dr. Erickson"s analytical solns');

tmp = repmat(colloc_pts,N,1);
KM = K(tmp',tmp);
KM_evals = K( repmat(pts',1,size(colloc_pts,2)), repmat(colloc_pts, ...
                                                  size(pts,2),1));
B = calculate_newton_basis(KM);
V = B';
V_evals = B\KM_evals';
subplot(1,2,2);
plot(pts, V_evals(1,:),cs(1), ...
     pts, V_evals(2,:),cs(2), ...
     pts, V_evals(3,:),cs(3), ...
     pts, V_evals(4,:),cs(4));
title('mine')

