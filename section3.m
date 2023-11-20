% EM lab section 3
clear all
% Generate Data

mk_data;
y=x';
[M,N]=size(y);

%Initalize theta
K = 3;
prior = 1/K*ones(1,K);
miu = y(:,1:K);
R = zeros(2,2*K);
for k=1:K
    R(:,(k-1)*M+1:k*M)=eye(2)
end
%R = [eye(2), eye(2), eye(2)];

% for iter = 1:20
% %E step
% for k = 1:K
% for n=1:N
% p(n,k) = Compute_A_Post(y,miu,R, prior, n, k);
% end
% end
%
% for k = 1:K
% p_k = p(:,k);
% N_k = sum(p_k)
% t1_k = y*p_k
% t2_k = y*diag(p_k)*y'
%
% %M Step
% miu(:,k) = t1_k/N_k
% R(:,(k-1)*M+1:k*M)=t2_k/N_k-t1_k*t1_k'/(N_k)^2
% prior(k)= N_k/N
%
% end %for k
% end %for iter
for iter = 1:20
    [miu, R, prior, cost] = EM_iter(y, miu, R, prior, K)
end