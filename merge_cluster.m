function [l,m, miu_new, R_new, prior_new] = merge_cluster(y, miu, R, prior)
% [l,m, miu_new, R_new, prior_new] = merge_cluster(y, miu, R, prior)
% return the same size parameters, set the tails to zero

[M,N]=size(y);
K = length(prior);

d = zeros(K,K); % distance matrix
for l=1:K
    for m = l+1:K
        miu_l = miu(:,l);
        miu_m = miu(:,m);
        R_l = R(:,(l-1)*M+1:l*M);
        R_m = R(:,(m-1)*M+1:m*M);
        miu_lm = (prior(l)*miu_l+prior(m)*miu(m))/(prior(l)+prior(m));
        R_lm = (prior(l)*(R_l+(miu_l-miu_lm)*(miu_l-miu_lm)') + prior(m)*(R_m+(miu_m-miu_lm)*(miu_m-miu_lm)'))/(prior(l)+prior(m));
        d(l,m)= N*prior(l)/2*log(det(R_lm)/det(R_l)) + N*prior(m)/2*log(det(R_lm)/det(R_m));
    end
end

% find the min distance by two 1D search
d1 = d+d'+10000000*diag(ones(1,K));
[min_row, I]=min(d1);
[min_all, l] = min(min_row);
m = I(l);

prior_new = zeros(size(prior,1), K-1);
miu_new = zeros(size(miu,1),K-1);
R_new = zeros(size(R,1),(K-1)*M);

% update the parameter
for i=1:K-1
    if (i<max(l,m))
        prior_new(i)=prior(i);
        miu_new(:,i) = miu(:,i);
        R_new(:,(i-1)*M+1:i*M)=R(:,(i-1)*M+1:i*M);
    end %if

    if (i==min(l,m))
        prior_new(i)= prior(l)+prior(m);
        miu_l = miu(:,l);
        miu_m = miu(:,m);
        R_l = R(:,(l-1)*M+1:l*M);
        R_m = R(:,(m-1)*M+1:m*M);
        miu_lm = (prior(l)*miu_l+prior(m)*miu(m))/(prior(l)+prior(m));
        R_lm = (prior(l)*(R_l+(miu_l-miu_lm)*(miu_l-miu_lm)') + prior(m)*(R_m+(miu_m-miu_lm)*(miu_m-miu_lm)'))/(prior(l)+prior(m));
    end %if i==min(l,m)

    if (i>=max(l,m))
        prior_new(i)=prior(i+1);
        miu_new(:,i) = miu(:,i+1);
        R_new(:,(i-1)*M+1:i*M)=R(:,i*M+1:(i+1)*M);
    end %if i>max(l,m)
end %for i