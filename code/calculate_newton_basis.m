function [V] = calculate_newton_basis(KM)
    N = size(KM,1);
    V = zeros(N,N);
    w = zeros(N,1);
    u = zeros(N,1);

    % get first column of value matrix
    V(:,1) = KM(:,1)./ sqrt(abs(KM(1,1)));
    z = diag(KM);
    w = w + V(:,1).^2;
    for i=2:N
        V(:,i) = (KM(i,:)' - V(:,1:(i-1))*V(i,1:(i-1))') ...
                  ./ sqrt(z(i)-w(i));
        w = w + V(:,i).^2;    
    end
end
