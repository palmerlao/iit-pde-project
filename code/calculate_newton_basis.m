function [V] = calculate_newton_basis(KM)
    N = size(KM,1);
    V = zeros(N,N);
    w = zeros(N,1);
    u = zeros(N,1);
    zminds = zeros(N,1);
    z = diag(KM);

    % get first column of value matrix
    zminds(1) = 1;
    V(:,1) = KM(:,1)./ sqrt(max(abs(z)));
    w = w + V(:,1).^2;
    for i=2:N
        if any(isnan(V))
            disp(['found NaN, N=' num2str(size(KM,1)) ',iteration ' ...
                 num2str(i)]);
        end
        
        [zm, zmind] = max(abs(z-w));
        if zm < eps*10
            disp(['power fcn below eps*10, only using ' num2str(i-1) ...
                 'pts']);
            V = V(:,1:(i-1));
            break
        end
        
        V(:,i) = (KM(:,zmind) - V(:,1:i)*V(zmind,1:i)') ...
                  ./ sqrt(abs(z(zmind)-w(zmind)));
        w = w + V(:,i).^2;
        zminds(i) = zmind;
    end
    
end
