function [B, V] = calculate_beta_v(KM, N, xs, K)
% here we alternately calculate betas and vs since they depend on each
% other (see paper for relationship). B should be lower triangular,
% and V should be unit upper triangular s.t. KM = B*V
    B = zeros(N,N);
    V = eye(N);
    for c=1:(N-1)
        for i=c:N
            B(i,c) = calculate_single_beta(B,V,i,c,K,xs);
        end
        % not sure if this is bad when KM is ill-cond.
        V(1:c,c+1) = B(1:c,1:c)\KM(1:c,c+1); 
    end
    B(N,N) = calculate_single_beta(B,V,N,N,K,xs);
end

function [res] = calculate_single_beta(B, V, i, j, K, xs)
    res = K(xs(j), xs(i));
    for k=1:j-1
        res = res - B(i,k).*V(k,j);
    end
end
