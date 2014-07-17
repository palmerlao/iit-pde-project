clear all; close all; clc;
set(0,'defaultLineLineWidth',3) % thick lines


count = 1;
NStart = 1;
NEnd = 200;

NVector = NStart:1:NEnd;

f=@(x) 3 + 2*x;


errS = zeros(1, NEnd - NStart);
errN = zeros(1, NEnd - NStart);
count = 1;
for N = NVector
    epsilon = (N/3).^2;
   K = @(x, point) exp(-epsilon.*(x-point).^2);
   % D1K = @(x,point) ( -2.*epsilon.*(x-point).*K(x,point) );
   D2K = @(x,point) ( 2.*epsilon.*(2.*epsilon.*((x-point).^2)-1).*K(x,point));

    samplePoints = linspace(0, 1, N);
    fsample=f (samplePoints)';
    temp = repmat(samplePoints, N, 1);
    % Kernel matrix
    KMatrix = K(temp', temp);
    
    KMatrixTilda = [K(temp', temp) ones(N, 1) samplePoints';
                    ones(1, N) 0 0;
                    samplePoints 0 0];

    uK=[KMatrix ones(N,1) samplePoints']*(KMatrixTilda\[fsample;0;0]);
    
    fAppx = ([KMatrix ones(N,1) samplePoints']*(KMatrixTilda\[fsample;0;0]));
    %% Newton Basis
    
     [B, VMatrix] = calculate_beta_v(KMatrix);
    
     VMatrixTilda = [VMatrix ones(N, 1) samplePoints';
                    ones(1, N) 0 0;
                    samplePoints 0 0];
                
     fAppxN = [KMatrix ones(N,1) samplePoints']*(VMatrixTilda\[fsample; 0; 0]);
     %% Errors
     
     errS(count) = max(abs(fAppx-fsample));
     errN(count) = max(abs(fAppxN-fsample));
     
     count = count + 1;
end


%keyboard;
plot(samplePoints, fsample, 'b', samplePoints, fAppx, 'g', samplePoints, fAppxN, 'c');

figure
plot(NVector, errS, 'g', NVector, errN, 'c')
xlabel('number of points')
ylabel('absolute error')

