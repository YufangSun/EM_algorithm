% EM lab section 4

clear all
% Generate Data
mk_data;
y=x';
[M,N]=size(y);

%Initalize theta
K = 9;
prior = 1/K*ones(1,K);
miu = y(:,1:K);
R = zeros(2,K);
for k=1:K
    R(:,(k-1)*M+1:k*M)=eye(2);
end

for k = K:-1:1
    L =k*(1+M+(M+1)*M/2)-1;
    for iter = 1:20
    % run EM algorithm
    [miu, R, prior, cost_k] = EM_iter(y, miu, R, prior, k);
    cost(k, iter) = cost_k + 0.5*L*log(M*N);
    end %for iter

    % merge cluster
    [l,m, miu, R, prior] = merge_cluster(y, miu, R, prior)
end %for k