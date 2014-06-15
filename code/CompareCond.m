clear all; close all; clc;
set(0,'defaultLineLineWidth',3) %thick lines



subplot(1, 2, 1);

Nstart = 1;
Nend = 100;

condPlotK = zeros(1, Nend - Nstart);
condPlotKTilda = zeros(1, Nend - Nstart);
condPlotL = zeros(1, Nend - Nstart);
condPlotLTilda = zeros(1, Nend - Nstart);
condPlotV = zeros(1, Nend - Nstart);
condPlotP = zeros(1, Nend - Nstart);
f=@(x) 1+exp(2*x);
u_analytic = @(x) ( 0.25.*((2.*x.^2)-exp(2).*x-x+exp(2.*x)-1) );
errL=zeros(1,Nend-Nstart);
errP=zeros(1,Nend-Nstart);
for N = Nstart : 1 : Nend
   %fix K and take derivatives (taken from Palmer's code)
   %epsilon = 200;
   epsilon = (N/4).^2;
   K = @(x, point) exp(-epsilon.*(x-point).^2);
   %D1K = @(x,center) ( -2.*epsilon.*(x-center).*K(x,center) );
   D2K = @(x,center) ( 2.*epsilon.*(2.*epsilon.*((x-center).^2)-1).*K(x,center));

    samplePoints = linspace(0, 1, N);
    fsample=f(samplePoints)';
    usample=u_analytic(samplePoints)';
    temp = repmat(samplePoints, N, 1);
    %Kernel matrix
    KMatrix = K(temp', temp);
    
    KMatrixTilda = [K(temp', temp) ones(N, 1) samplePoints';
                    ones(1, N) 0 0;
                    samplePoints 0 0];

    %For this example, L = 2nd derivative
    LMatrix = D2K(temp', temp);
    LMatrixTilda = [D2K(temp', temp) zeros(N, 2);
               K(0, samplePoints) 1 0;
               K(1, samplePoints) 1 1];
    cond(KMatrix);
    condPlotK(N-Nstart+1) = cond(KMatrix);
    condPlotKTilda(N-Nstart+1) = cond(KMatrixTilda);
    condPlotL(N-Nstart+1) = cond(LMatrix);
    condPlotLTilda(N-Nstart+1) = cond(LMatrixTilda);
    
    %keyboard
    uL=[KMatrix ones(N,1) samplePoints']*(LMatrixTilda\[fsample;0;0]);
    errL(N-Nstart+1)=max(abs(usample-uL));
    
    
    %%Newton Basis
    [B, VMatrix] = calculate_beta_v(KMatrix, N, samplePoints, K);
    condPlotV(N-Nstart+1) = cond(VMatrix);
    
    %should work since L is linear operator (double check)
    [B, PMatrix] = calculate_beta_v(LMatrix, N, samplePoints, D2K);
    condPlotP(N-Nstart+1) = cond(PMatrix);
     PMatrixTilda = [PMatrix zeros(N, 2);
               VMatrix(:,1)' 1 0;
               VMatrix(:,end)' 1 1];

    uP=[VMatrix ones(N,1) samplePoints']*(PMatrixTilda\[fsample;0;0]);
    errP(N-Nstart+1)=max(abs(usample-uP));


end

    semilogy(Nstart:1:Nend, condPlotK, ':b')
    hold on;
    semilogy(Nstart:1:Nend, condPlotKTilda, 'b')
    semilogy(Nstart:1:Nend, condPlotL, ':g')
    semilogy(Nstart:1:Nend, condPlotLTilda, 'g')
    semilogy(Nstart:1:Nend, condPlotV, ':r')
    semilogy(Nstart:1:Nend, condPlotP, ':c')
    
    title('Condition Number vs. Number of points sampled');
    legend('K Matrix', 'K~ Matrix', 'L Matrix', 'L~ Matrix', 'V Matrix', 'P Matrix', 'Location', 'NorthWest');
    ylabel('Condition Number');
    xlabel('Number of points sampled');


   subplot(1, 2, 2);
    hold on;
    plot(Nstart:1:Nend, (condPlotK./condPlotL), ':g')
    plot(Nstart:1:Nend, (condPlotKTilda./condPlotLTilda), 'g')
    plot(Nstart:1:Nend, (condPlotV./condPlotP), ':r')
    
    title('Ratio of Condition Numbers');
    legend('K/L', 'K~/L~', 'V/P', 'Location', 'NorthWest');
    ylabel('Ratios');
    xlabel('Number of points sampled');
    
    figure
    semilogy(Nstart:1:Nend,errL,'g',Nstart:1:Nend,errP,'r')

    %see what ratios are 
    %KLRatio = condPlotK (100)/condPlotL(100)
    %KLTildaRatio = condPlotKTilda (100)/condPlotLTilda(100)
    %VPRatio = condPlotV (100)/condPlotP(100)
    condPlotV(end)
