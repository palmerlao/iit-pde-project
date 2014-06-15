function [B, V] = chol_calc_beta_v(KM, N, xs, K)
    B = chol(KM);
    V = B';
end
