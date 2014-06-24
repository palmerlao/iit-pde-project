
clear all; close all; clc;
set(0,'defaultLineLineWidth',3) %thick lines



subplot(1, 2, 1);

count = 1;
NStart = 0;
NEnd = 7;

NVector = floor(2.^(NStart:1:NEnd));

condPlotK = zeros(1, NEnd - NStart);
condPlotKTilda = zeros(1, NEnd - NStart);
condPlotL = zeros(1, NEnd - NStart);
condPlotLTilda = zeros(1, NEnd - NStart);
condPlotV = zeros(1, NEnd - NStart);
condPlotVTilda = zeros(1, NEnd - NStart);
condPlotP = zeros(1, NEnd - NStart);
condPlotPTilda = zeros(1, NEnd - NStart);
f=@(x) (x-1).^2;
u_analytic = @(x) (x.^4)/12 - (x.^3)/3 + (x.^2)/2;
errL=zeros(1,NEnd-NStart);
errP=zeros(1,NEnd-NStart);
for N = NVector
   %keyboard
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
               zeros(2, N) [1 0; 0 1]];
    cond(KMatrix);
    condPlotK(count) = cond(KMatrix);
    condPlotKTilda(count) = cond(KMatrixTilda);
    condPlotL(count) = cond(LMatrix);
    condPlotLTilda(count) = cond(LMatrixTilda);
    
    %keyboard
    uL=[KMatrix ones(N,1) samplePoints']*(LMatrixTilda\[fsample;0;0]);
    errL(count)=max(abs(usample-uL));
    
    
    %%Newton Basis
    [B, VMatrix] = calculate_beta_v(KMatrix, N, samplePoints, K);
    condPlotV(count) = cond(VMatrix);
    
    VMatrixTilda = [VMatrix ones(N, 1) samplePoints';
                    ones(1, N) 0 0;
                    samplePoints 0 0];
    
    condPlotVTilda(count) = cond(VMatrixTilda);
    
    %should work since L is linear operator 
    [B, PMatrix] = calculate_beta_v(LMatrix, N, samplePoints, D2K);
    condPlotP(count) = cond(PMatrix);
     PMatrixTilda = [PMatrix zeros(N, 2);
               zeros(2, N) [1 0; 0 1]];

    condPlotPTilda(count) = cond(PMatrixTilda);
    uP=[VMatrix ones(N,1) samplePoints']*(PMatrixTilda\[fsample;0;0]);
    errP(count)=max(abs(usample-uP));
    
    

count = count+1;
end

    semilogy(NVector, (condPlotK), ':b')
    hold on;
    semilogy(NVector, (condPlotKTilda), 'b')
    semilogy(NVector, (condPlotL), ':g')
    semilogy(NVector, (condPlotLTilda), 'g')
    semilogy(NVector, (condPlotV), ':r')
    semilogy(NVector, (condPlotVTilda), 'r')
    semilogy(NVector, (condPlotP), ':c')
    semilogy(NVector, (condPlotPTilda), 'c')
    
    title('Condition Number vs. Number of points sampled');
    legend('K Matrix', 'K~ Matrix', 'L Matrix', 'L~ Matrix', 'V Matrix','V~ Matrix', 'P Matrix', 'P~ Matrix', 'Location', 'NorthWest');
    ylabel('Condition Number');
    xlabel('Number of points sampled');


   subplot(1, 2, 2);
    hold on;
    plot(NVector, (condPlotK./condPlotL), ':g')
    plot(NVector, (condPlotKTilda./condPlotLTilda), 'g')
    plot(NVector, (condPlotV./condPlotP), ':r')
    
    title('Ratio of Condition Numbers');
    legend('K/L', 'K~/L~', 'V/P', 'Location', 'NorthWest');
    ylabel('Ratios');
    xlabel('Number of points sampled');
    
    %keyboard
    figure
    subplot(1, 2, 1);
    plot(samplePoints, usample, 'b', samplePoints, uL, 'g', samplePoints, uP, 'c')
    
    title('Analytic and Numerical functions');
    legend('Analytic Solution', 'Kernel Basis Solution', 'Newton Basis Solution');
    
    subplot(1, 2, 2);
    semilogy(NVector,errL,'g',NVector,errP,'c')
    
    title('Absolute Error');
    legend('Error with Kernel Basis', 'Error with Newton Basis', 'Location', 'Northwest');
    ylabel('Absolute Error');
    xlabel('Number of points sampled');





