function [V] = calculate_newton_basis(KM)
    N = size(KM,1);
    V = zeros(N,N);
    w = zeros(N,1);
    u = zeros(N,1);

    z = diag(KM);
    % get first column of value matrix
    V(:,1) = KM(:,1)./ sqrt(max(abs(z)));
    w = w + V(:,1).^2;
    for i=2:N
        [zm, zmind] = max(z-w);
        V(:,i) = (KM(:,zmind) - V(:,1:i)*V(zmind,1:i)') ...
                  ./ sqrt(abs(z(zmind)-w(zmind)));
        w = w + V(:,i).^2;
    end
end
