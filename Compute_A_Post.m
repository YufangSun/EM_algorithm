function [p_k_n, cost_n] = Compute_A_Post(y,miu,R, prior, n, k)
% p = Compute_A_Post(y,miu,R, pi, n, k)
% y is the measurement data
% miu and R is parameter of Gaussian distribution given x_k
% The function computes P(x_n = k|Y = y, theta)

[M,N]= size(y);
[temp ,K] = size(miu);
p = zeros(1,K);

for i=1:K
    miu_k = miu(:,i);
    R_k = R(:,(i-1)*M+1:i*M);
    pi_k = prior(i);
    scaler = pi_k/(sqrt(det(R_k))*(2*pi)^(M/2));

    p(i) = p(i) + scaler*exp(-1/2*(y(:,n)-miu_k)'*(inv(R_k))*(y(:,n)-miu_k));
end % for k

p_k_n = p(k)/sum(p);
cost_n = sum(p);