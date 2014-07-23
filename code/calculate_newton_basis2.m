function [V, varargout] = calculate_newton_basis2(KM,LKM,tol_mult)
    N = size(KM,1);
    V = zeros(N,N);
    w = zeros(N,1);
    u = zeros(N,1);
    zminds = zeros(N,1);
    was_included = false(N,1);
    z = diag(KM);
    
    pf = zeros(N,1);

    % get first column of value matrix
    zminds(1) = 1;
    V(:,1) = KM(:,1)./ sqrt(max(abs(z)));
    w = w + V(:,1).^2;
    was_included(1) = true;
    for i=2:N
        
        for j=1:N
            u = V(j,1:i-1);
            pf(j) = -2.*u*LKM(j,1:i-1)' + u*KM(1:i-1,1:i-1)*u';
        end
        
        pf(was_included) = 0;
        
        [zm, zmind] = max(pf);
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
