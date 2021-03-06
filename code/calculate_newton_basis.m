function [V, varargout] = calculate_newton_basis(KM,tol_mult)
    N = size(KM,1);
    V = zeros(N,N);
    w = zeros(N,1);
    u = zeros(N,1);
    zminds = zeros(N,1);
    was_included = false(N,1);
    z = diag(KM);

    % get first column of value matrix
    zminds(1) = 1;
    V(:,1) = KM(:,1)./ sqrt(max(abs(z)));
    w = w + V(:,1).^2;
    was_included(1) = true;
    for i=2:N
        
        [zm, zmind] = max(abs(z-w));
        if zm < eps*tol_mult
            disp(['power fcn below eps*' num2str(tol_mult)  ', only using ' num2str(i-1) ...
                 '/' num2str(size(KM,1)) 'pts']);
            zminds = zminds(1:(i-1));
            V = V(:,1:(i-1));
            break
        end
        
        V(:,i) = (KM(:,zmind) - V(:,1:i)*V(zmind,1:i)') ...
                  ./ sqrt(abs(z(zmind)-w(zmind)));
        w = w + V(:,i).^2;
        zminds(i) = zmind;
        was_included(zmind) = true;
        
    end
    if nargout==2
        varargout = {zminds};
    end
    if nargout==3
        varargout = {zminds, was_included};
    end
    
end
