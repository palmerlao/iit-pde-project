function [] = main()
    clear all; close all; clc

    % Solves u'' = 1 + e^(2x), u(0) = 0 = u(1)

    epsilon = 1;
    K   = @(x,center) ( exp(-epsilon.*((x-center).^2)) );
    D1K = @(x,center) ( -2.*epsilon.*(x-center).*K(x,center) );
    D2K = @(x,center) ( 2.*epsilon.*(2.*epsilon.*((x-center).^2)-1).* ...
                        K(x,center) );

    rhs = @(x) ( 1 + exp(2.*x) );
    u_analytic = @(x) ( 0.25.*((2.*x.^2)-exp(2).*x-x+exp(2.*x)-1) );

    N = 6;
    colloc_pts = linspace(0,1,N);
    pts = linspace(0,1);
    
    tmp = repmat(colloc_pts,N,1);
    KM = K(tmp',tmp);
    
    [B, V] = calculate_beta_v(KM, N, colloc_pts, K);

end

function [B, V] = calculate_beta_v(KM, N, xs, K)
    
    B = tril(-999.*ones(N,N));
    V = B' - diag(diag(B)) + diag(ones(1,N));
    for c=1:(N-1)
        for i=c:N
            B(i,c) = calculate_single_beta(B,V,i,c,K,xs);
        end
        V(1:c,c+1) = B(1:c,1:c)\KM(1:c,c+1);
    end
    B(N,N) = calculate_single_beta(B,V,N,N,K,xs);
end

function [res] = calculate_single_beta(B, V, i, j, K, xs)
    res = K(xs(j), xs(i));
    for k=1:j-1
        if ( V(k,j) == -999 )
            fprintf(['V(%d,%d) is uninitialized, but accessed ' ...
                     'anyway!\n'], k, i);
        end
        if ( B(i,k) == -999 )
            fprintf(['B(%d,%d) is uninitialized, but accessed ' ...
                     'anyway!\n'], j, k);
        end
        res = res - B(i,k).*V(k,j);
    end
end
