function [miu, R, prior, cost] = EM_iter(y, miu_init, R_init, prior_init, K)
% [miu, R, prior] = EM_iter(y, mui, R, prior, num_iter) performs EM
% Perform one iteration of EM algorithm on data y, given initial conditions
[M,N]=size(y);

miu = miu_init;
R= R_init;
prior= prior_init;

%perform one iteration
%E step
for n = 1:N
    for k=1:K
    [p(n,k),cost_local(n)] = Compute_A_Post(y,miu,R, prior, n, k);
    end
end

log(cost_local);

cost = -sum(log(cost_local))
for k = 1:K
    p_k = p(:,k);
    N_k = sum(p_k);
    t1_k = y*p_k;
    t2_k = y*diag(p_k)*y';
    
    %M Step
    miu(:,k) = t1_k/N_k;
    R(:,(k-1)*M+1:k*M)=t2_k/N_k-t1_k*t1_k'/(N_k)^2;
    prior(k)= N_k/N;
end %for k