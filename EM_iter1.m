function [miu_0, R_0, prior_0, miu_1, R_1, prior_1, miu_2, R_2, prior_2] = ...
EM_iter1(data, miu_0_init, R_0_init, prior_0_init, miu_1_init, R_1_init, prior_1_init, miu_2_init, R_2_init, prior_2_init)

%load data;
y=data';
[data_num,data_dim]=size(data);
cluster_num = 3;

miu_0 = miu_0_init;
miu_1 = miu_1_init;
miu_2 = miu_2_init;
R_0 = R_0_init;
R_1 = R_1_init;
R_2 = R_2_init;
prior_0 = prior_0_init;
prior_1 = prior_1_init;
prior_2 = prior_2_init;

P_0=zeros(1,data_num);
P_1=zeros(1,data_num);
P_2=zeros(1,data_num);
f_0=zeros(1,data_num);
f_1=zeros(1,data_num);
f_2=zeros(1,data_num);
P=zeros(1,data_num);

t_1_0 = zeros(2,1);
t_1_1 = zeros(2,1);
t_1_2 = zeros(2,1);

t_2_0 = zeros(2,2);
t_2_1 = zeros(2,2);
t_2_2 = zeros(2,2);

for n = 1: data_num
    P_0(n) = ((2*pi)^(-data_dim/2)) * ((det(R_0))^(-0.5)) * exp(-0.5*(y(:,n)-miu_0)'*(det(R_0)^(-1))*(y(:,n)-miu_0)) * prior_0;
    P_1(n) = ((2*pi)^(-data_dim/2)) * ((det(R_1))^(-0.5)) * exp(-0.5*(y(:,n)-miu_1)'*(det(R_1)^(-1))*(y(:,n)-miu_1)) * prior_1;
    P_2(n) = ((2*pi)^(-data_dim/2)) * ((det(R_2))^(-0.5)) * exp(-0.5*(y(:,n)-miu_2)'*(det(R_2)^(-1))*(y(:,n)-miu_2)) * prior_2;

    P(n)=P_0(n)+P_1(n)+P_2(n);

    f_0(n) = P_0(n)/P(n);
    f_1(n) = P_1(n)/P(n);
    f_2(n) = P_2(n)/P(n);
end

%E-step
N_0 = sum(f_0);
N_1 = sum(f_1);
N_2 = sum(f_2);

for n = 1:data_num
    t_1_0 = t_1_0 + y(:,n)*f_0(n);
    t_1_1 = t_1_1 + y(:,n)*f_1(n);
    t_1_2 = t_1_2 + y(:,n)*f_2(n);
end

for n = 1:data_num
    t_2_0 = t_2_0 + y(:,n)*(y(:,n))'*f_0(n);
    t_2_1 = t_2_1 + y(:,n)*(y(:,n))'*f_1(n);
    t_2_2 = t_2_2 + y(:,n)*(y(:,n))'*f_2(n);
end

%M-step
miu_0 = t_1_0 / N_0;
miu_1 = t_1_1 / N_1;
miu_2 = t_1_2 / N_2;

prior_0 = N_0/data_num;
prior_1 = N_1/data_num;
prior_2 = N_2/data_num;

R_0 = t_2_0/N_0 - t_1_0*t_1_0'/(N_0^2);
R_1 = t_2_1/N_1 - t_1_1*t_1_1'/(N_1^2);
R_2 = t_2_2/N_2 - t_1_2*t_1_2'/(N_2^2);