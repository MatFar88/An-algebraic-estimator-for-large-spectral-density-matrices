%% clearvars -except f p T c r excl N_tot rho delta deltabis tau alpha imm same_rank sparse_s sparse_diff coef_lag new_method past_method k_1 var_1 var_2 pert_soft prec_pert no_pert sign_diff

clear;
p=100 %% dimension
T=1000 %% sample size
c=2 %% condition number of the low rank component
N_tot=100 %% number of replicates
r=2 %% latent rank
delta=1 %% magnitude of residual covariances VS their theoretical maximum
deltabis=0.7 %% threshold to set surviving residual elements
tau=3 %% scale parameter VS p
alpha=0.9 %% latent variance parameter
imm=1i; %% imaginary unit
same_rank=1 %% 0 is wrong. Set to 1.
sparse_s=2 %% 0 for fixed residual spectral matrix across frequencies
sparse_diff=1 %% 1 to change the residual pattern across coefficient matrices
coef_lag=[0.5 0.5 0] %% squared lag coefficients
new_method=1; %% spectral computation method
past_method=0;
k_1=-0.2 %% var parameter
var_1=0 %% if var1 desired
var_2=0 %% if var2 desired
no_pert=1; %% 0 to perturb latent eigenvalues
pert_soft=1; %% 1 for soft perturbation, 0 for extreme perturbation
prec_pert=0.1; %% parameter controlling the amount of perturbation
f=0:1/12:0.5 %% sampled frequencies (include 0.5 for tecnical reasons)
%sparse_diff=1
%clear Gamma_pre Eig_low
%%

% clearvars -except p T c r excl N_tot rho delta deltabis tau alpha imm same_rank sparse_s sparse_diff coef_lag new_method past_method k_1 var_1 var_2 var_1_com var_2_com pert_soft prec_pert no_pert
% ne, find!

var_1_com=0;
var_2_com=0;

%%

R=zeros(p,p);
V=zeros(p,p);
E=zeros(p,p);
v=zeros(p,1);
K=rand(p)*eye(p);
rank(K);
E(:,1)=K(:,1);
%for i=1:p
 %   for j=(i+1):p
%R(i,j)=dot(K(j,:),E(i,:))/norm(E(i,:));
%    end;
%end;
for j=2:p
for i=1:(j-1)
    R(i,j)=dot(K(:,j),E(:,i))/(norm(E(:,i))^2);
end;
for i=1:(j-1)
V(:,i)=R(i,j)*E(:,i);
end;
    for h=1:p
    v(h)=sum(V(h,1:(j-1)));
    end;
   E(:,j)=K(:,j)-v;
end;

rank(E)
E'*E

for i=1:p
E(:,i)=E(:,i)/norm(E(:,i));
end;
rank(E)
E
E'*E

v=randperm(p,r);
v
E_r=E(:,v)
rank(E_r)

%T=[256];
%r=[round(rho*p)];%%function of T!!!p=[(T/64) T/32 T/16 T/8 T/4 (T/2)]
lambda=zeros(r, length(c));
Lambda=zeros(r,r,length(c));
condn=zeros(1,length(c));
normFrL=zeros(1,length(c));
    for i=1:length(c);
    a=[1-c(i) r-1];
    f_past=[(r-1)+c(i) sum(1:(r-2))];
    A=[a;f_past];
    B=[0 tau*alpha*p];
    sol(:,i)=linsolve(A, B');
    lambda(r,i)=sol(1,i);
    for q=2:r 
    lambda(r-q+1,i)=lambda(r,i)+((q-1))*sol(2,i);
    end;
    Lambda(:,:,i)=diag(lambda(:,i));
    condn(i)=cond(Lambda(:,:,i),2);
    end;
Lambda
condn

%Lambda(2,2)=Lambda(2,2)+1/30;
%Lambda(3,3)=Lambda(3,3)-1/30;

Lambda

for i=1:length(c);
    normFrL(i)=norm(Lambda(:,:,i),'fro');
end;
normFrL

i=1;
j=1;
B_new=zeros(p);
B=E_r*Lambda*E_r';

rank(B)
trace(B)
diag_B=svds(B,rank(B))
diag_B(rank(B))/diag_B(1)
%[B_U B_D]=svd(B)

%%

if same_rank==0
    n_lag=length(coef_lag)-1;
end
if same_rank==1
    n_lag=0;
end

T_A=(n_lag+1)*p;

R=zeros(T_A,T_A);
V=zeros(T_A,T_A);
E=zeros(T_A,T_A);
v=zeros(T_A,1);
K=rand(T_A)*eye(T_A);
rank(K);
E(:,1)=K(:,1);
%for i=1:T_A
 %   for j=(i+1):T_A
%R(i,j)=dot(K(j,:),E(i,:))/norm(E(i,:));
%    end;
%end;
for j=2:T_A
for i=1:(j-1)
    R(i,j)=dot(K(:,j),E(:,i))/(norm(E(:,i))^2);
end;
for i=1:(j-1)
V(:,i)=R(i,j)*E(:,i);
end;
    for h=1:T_A
    v(h)=sum(V(h,1:(j-1)));
    end;
   E(:,j)=K(:,j)-v;
end;

rank(E)
E'*E

for i=1:T_A
E(:,i)=E(:,i)/norm(E(:,i));
end;
rank(E)
E
E'*E

%v=randperm(p,r);
%v
E_r=E(:,1:(n_lag+1))
rank(E_r)

E_r_past=E_r
v=randperm(p,r);
v
E_r=E(:,v)

%E_r=E(:,1:(r))
%rank(E_r)

if same_rank==0
clear A_delay
for lag=0:n_lag
A_delay(:,:,lag+1)=E_r(lag*p+1:(lag+1)*p,:);
end;
end

if var_1_com==0 && var_2_com==0
%%
n_boh=1;
if same_rank==1
n_lag=find(coef_lag==0)-2%length(coef_lag)-1;
lag_0=0;
clear A_delay Coef_lag
for lag=0:n_lag
    
    
    if no_pert==1
    A_delay(:,:,lag+1)=sqrt(coef_lag(lag+1))*E_r(lag_0*p+1:(lag_0+1)*p,:);
    end
    
    if no_pert==0
        
    if pert_soft==0
        
    vec_coef=repmat(r*(coef_lag(lag+1)),r,1);
    vec_lag = gamrnd(repmat(vec_coef, 1, n_boh),1);
    vec_sum =sum(vec_lag);
    %vec_lag = vec_lag./repmat(vec_sum, size(vec_lag, 1), 1);
    vec_lag=vec_lag./(vec_sum)*(r*(coef_lag(lag+1)));
    if lag<n_lag
    Coef_lag(:,:,lag+1)=diag(vec_lag);
    end
    if lag==n_lag
        for i=1:r
    Coef_lag(i,i,lag+1)=sum((coef_lag))-sum(Coef_lag(i,i,:));
        end
    end
    
    end
    
    if pert_soft==1
    if lag<n_lag
        pert_coef=(-prec_pert*coef_lag(lag+1)+2*prec_pert*coef_lag(lag+1)*rand(1,r)); 
    for i=1:r
    Coef_lag(i,i,lag+1)=coef_lag(lag+1)+pert_coef(i);
    end
    end
    if lag==n_lag
    for i=1:r
    Coef_lag(i,i,lag+1)=sum((coef_lag))-sum(Coef_lag(i,i,:));
    end
    end
    
    end
    
    A_delay(:,:,lag+1)=E_r(lag_0*p+1:(lag_0+1)*p,:)*sqrt(Coef_lag(:,:,lag+1));
    end
end;
end
%%

size(E_r)
size(A_delay)

for lag=0:n_lag
for ncomp=1:r
norm_comp(lag+1,ncomp)=norm(A_delay(:,ncomp,lag+1));
end;
end

for n_comp=1:r
norm_r(n_comp)=norm(E_r(:,n_comp));
end

%%
for lag_n=n_lag:-1:0
for lag=0:lag_n
    u_lag=n_lag-lag_n
    AA_delay_all(:,:,lag+1,u_lag+1)=A_delay(:,:,lag+u_lag+1)*Lambda(1:r,1:r)*A_delay(:,:,lag+1)';
end
end

%past_method=0;
if past_method==1
for lag=0:2
    AA_delay(:,:,lag+1)=A_delay(:,:,lag+1)*Lambda(1:r,1:r)*A_delay(:,:,lag+1)';
end

for lag=0:1
    AA_delay_1(:,:,lag+1)=A_delay(:,:,lag+2)*Lambda(1:r,1:r)*A_delay(:,:,lag+1)';
end

for lag=0:0
    AA_delay_2(:,:,lag+1)=A_delay(:,:,lag+3)*Lambda(1:r,1:r)*A_delay(:,:,lag+1)';
end

size(AA_delay)
for lag=0:n_lag
r_lag(lag+1)=rank(AA_delay(:,:,lag+1));
end

size(AA_delay_1)
for lag=0:n_lag-1
r_lag_1(lag+1)=rank(AA_delay_1(:,:,lag+1));
end



for i=1:p
    for j=1:p
        Sigma_0(i,j)=sum(AA_delay(i,j,:));
        Sigma_1(i,j)=sum(AA_delay_1(i,j,:));
        Sigma_2(i,j)=sum(AA_delay_2(i,j,:));
    end
end
rank(Sigma_0)
rank(Sigma_1)
rank(Sigma_2)

size(Sigma_0)
rank(Sigma_0)
eigs(Sigma_0,p)
svds(Sigma_0,p)
svds(Sigma_1,p)
svds(Sigma_2,p)
trace(Sigma_0)
trace(Sigma_1)
trace(Sigma_2)
trace(Sigma_0)+trace(Sigma_1)+trace(Sigma_2)
cond(Sigma_0)
cond(Sigma_1)
cond(Sigma_2)

Sigma_0(1,2)-Sigma_0(2,1)


end

for i=1:p
    for j=1:p
        for lag=0:n_lag
Sigma(i,j,lag+1)=sum(AA_delay_all(i,j,:,lag+1));
        end
    end
end

Sigma_0=Sigma(:,:,1)
if n_lag>0
Sigma_1=Sigma(:,:,2)
end
if n_lag>1
Sigma_2=Sigma(:,:,3)
end

end

if var_1_com==1 || var_2_com==1
   Sigma_0=B; 
end

%B_past=B;
%B=E_r*diag(lambda)*E_r';
%%

for i=1:1
    if sparse_diff==0
%for g=1:length(f)

for i=1:1%length(c)
bb=diag(Sigma_0(:,:,i));
end

l=zeros(p,1);
for i=1:p
l(i,1)=sum(bb)/alpha*(1-alpha);
end;

n_boh=1;

clear v

y = gamrnd(repmat(l, 1, n_boh),1);
v =sum(y);
v = y./repmat(v, size(y, 1), 1);

v=sum(bb)/(alpha)*(1-alpha)*v;
%v=sort(v,'descend');
sum(v);

clear A

%bb=diag(B);
bbord=sort(bb,'descend');
aaord=sort(v,'descend');
A=zeros(p,p);
for i=1:p
for j=1:p    
    if Sigma_0(i,i)==bbord(j)
       A(i,i)=aaord(j);
    end;
end;
end;
    
%A=zeros(p);
%for i=1:p
%    A(i,i)=v(i);
%end;
%A
    
f_past2(1,1)=0;
for i=2:p
f_past2(i,i)=sum(p-(i-1):p-1);
end;

for i=1:p
    for j=i+1:p
        f_past2(i,j)=f_past2(i,j-1)+1;
    end;
end;

%%THERE!!!

u=zeros(p);
h=zeros(p);
for i=1:p
    for j=i+1:p
        u(i,j)=random('unif',0,delta*sqrt(A(i,i))*sqrt(A(j,j)));
        %h(i,j)=random('unif',-u(i,j),u(i,j));
        unif(i,j)=random('unif',0,u(i,j));
        %h(i,j)=sign(MA_low(i,j,g))*unif(i,j)+1i*sign(MA_low(i,j,g))*unif(i,j);
        h(i,j)=sign(Sigma_0(i,j))*random('unif',0,u(i,j));
        %h(i,j)=sign(B(i,j))*random('unif',-sqrt(v(i))*sqrt(v(j)),sqrt(v(i))*sqrt(v(j)));
    end;
end;

iter=1;
for i=1:p
    for j=i+1:p
        resid(iter)=h(i,j);
        iter=iter+1;
    end;
end;

%gamma=0.2;
numvar=p*(p-1)/2;
%s=round(gamma*numvar);

order=sort(abs(resid),'descend');
%yes=order(1:s);
%sum(yes)/sum(order)
maxcov=max(order);

%sp_thr=1.5*mean(order);

%max(yes);
%min(yes);

P_B=B*inv(B'*B)*B';
N=eye(p);
for i=1:p
beta(i)=norm(P_B*N(:,i));
end;
inc=max(beta)/(p)
%inc=max(beta)

deg_max=p/(108*inc)%p
eigB=sort(svds(B,p),'descend');

min_allowed=sqrt(delta)*tau*p*inc^2*eigB(r)/(deg_max);
%sp_thr=0.3*delta*eigB(r)/(p*deg_max)

%eigB(r)/((2*inc)^2)%max
%eigB(r)/((inc)^2)%min
%eigB(r)/deg_max

sp_thr=(deltabis)*maxcov;
maxcov;

%sp_thr=tau*p*(1-alpha)/(s*inc)
%sp_thr=sp_thr*delta

%sp_thr=delta*eigB(r)/(s*inc);
%%OK!
%trace(M)*(1-alpha)/(2*s);

for i=1:length(order)
if order(i)<min_allowed
    break;
end;
end;
if i>1
s_max=i-1;
else s_max=0;
%sp_max=order(s_max);
perc_max=s_max/numvar;
end

for i=1:length(order)
if order(i)<sp_thr
    break;
end;
end;
i-1;
s_pre=i-1;

perc_proposed_pre=s_pre/numvar;

s=min(s_pre,s_max);
perc_proposed=s/numvar;

if s>0
%for i=1:p
thr=order(s);
k=randsample(numvar,s);
%gamma=s/numvar;
%end;
%sp_thr=thr;
%else sp_thr=order(1)+1;
end

%A=zeros(p);
for i=1:p
    for j=i+1:p
    %for q=1:s
    if abs(h(i,j))>=sp_thr%yes(s);
       %ge(abs(h(i,j),order(s))==1; 
       A(i,j)=h(i,j);
    else A(i,j)=0;
    end;
       zs(i,j)=sign(A(i,j));
    %end;
    end;
end;
A 

%A=zeros(p);
%A(1,:)=[v(1) a(1:9)];
%for i=2:p
%    A(i,:)=[zeros(1,i-1) v(i) a(sum(p-i+1:p-1)+1:sum(p-i:p-1))];
%end;

for i=1:p
    for j=i+1:p
        A(j,i)=A(i,j);
        zs(j,i)=zs(i,j);
    end;
end;

A
trace(A)
rank(A)
cond(A)
svds(A)

aa=diag(A);
bb=diag(B);
max(aa);
max(bb);
plot(sort(aa,'descend'),'Color','m')
line(1:p,sort(bb,'descend'),'Color','r')
title('Diagonal elements of S* (magenta) and L* (red)');
xlabel('order');
ylabel('diagonal elements');
diffab=bb-aa;
plot(diffab);
cntab=0;
for i=1:length(diffab)
    if diffab(i)<0
       cntab=cntab+1;
end;
end;
cntab

%MA_sparse(:,:,g)=A;

P_B=B*inv(B'*B)*B';
N=eye(p);
for i=1:p
beta(i)=norm(P_B*N(:,i));
end;
inc=max(beta);
for i=1:p
v=0;
    for j=1:p
        if A(i,j)~=0
            v=v+1;
        end;
    end;
deg(i)=v;
v=0;
end;
deg_max_real=max(deg);
sort(deg,'descend')
deg_max_real*inc;
deg_max
s_max
s
s_pre
ident=inc*deg_max_real/p^2
1/108

    end
    
    if sparse_diff==1
    
        for lag=0:max(find(ne(coef_lag,0)))-1%n_lag
        
for i=1:1%length(c)
bb=diag(Sigma_0(:,:,i));
end

l=zeros(p,1);
for i=1:p
l(i,1)=sum(bb)/alpha*(1-alpha)*coef_lag(lag+1);
end;

n_boh=1;

y = gamrnd(repmat(l, 1, n_boh),1);
v =sum(y);
v = y./repmat(v, size(y, 1), 1);

v=sum(bb)/(alpha)*(1-alpha)*v*coef_lag(lag+1);
%v=sort(v,'descend');
sum(v);

clear A

%bb=diag(B);
bbord=sort(bb,'descend');
aaord=sort(v,'descend');
A=zeros(p,p);
for i=1:p
for j=1:p    
    if Sigma_0(i,i)==bbord(j)
       A(i,i)=aaord(j);
    end;
end;
end;
    
%A=zeros(p);
%for i=1:p
%    A(i,i)=v(i);
%end;
%A
    
f_past2(1,1)=0;
for i=2:p
f_past2(i,i)=sum(p-(i-1):p-1);
end;

for i=1:p
    for j=i+1:p
        f_past2(i,j)=f_past2(i,j-1)+1;
    end;
end;

%%THERE!!!

u=zeros(p);
h=zeros(p);
for i=1:p
    for j=i+1:p
        u(i,j)=random('unif',0,delta*sqrt(A(i,i))*sqrt(A(j,j)));
        %h(i,j)=random('unif',-u(i,j),u(i,j));
        unif(i,j)=random('unif',0,u(i,j));
        %h(i,j)=sign(MA_low(i,j,g))*unif(i,j)+1i*sign(MA_low(i,j,g))*unif(i,j);
        h(i,j)=sign(Sigma_0(i,j))*random('unif',0,u(i,j));
        %h(i,j)=sign(B(i,j))*random('unif',-sqrt(v(i))*sqrt(v(j)),sqrt(v(i))*sqrt(v(j)));
    end;
end;

iter=1;
for i=1:p
    for j=i+1:p
        resid(iter)=h(i,j);
        iter=iter+1;
    end;
end;

%gamma=0.2;
numvar=p*(p-1)/2;
%s=round(gamma*numvar);

order=sort(abs(resid),'descend');
%yes=order(1:s);
%sum(yes)/sum(order)
maxcov=max(order);

%sp_thr=1.5*mean(order);

%max(yes);
%min(yes);

P_B=B*inv(B'*B)*B';
N=eye(p);
for i=1:p
beta(i)=norm(P_B*N(:,i));
end;
inc=max(beta)/(p)
%inc=max(beta)

deg_max=p/(108*inc)%p
eigB=sort(svds(B,p),'descend');

min_allowed=sqrt(delta)*tau*p*inc^2*eigB(r)/(deg_max);
%sp_thr=0.3*delta*eigB(r)/(p*deg_max)

%eigB(r)/((2*inc)^2)%max
%eigB(r)/((inc)^2)%min
%eigB(r)/deg_max

sp_thr=(deltabis)*maxcov;
maxcov;

%sp_thr=tau*p*(1-alpha)/(s*inc)
%sp_thr=sp_thr*delta

%sp_thr=delta*eigB(r)/(s*inc);
%%OK!
%trace(M)*(1-alpha)/(2*s);

for i=1:length(order)
if order(i)<min_allowed
    break;
end;
end;
if i>1
s_max=i-1;
else s_max=0;
%sp_max=order(s_max);
perc_max=s_max/numvar;
end

for i=1:length(order)
if order(i)<sp_thr
    break;
end;
end;
i-1;
s_pre=i-1;

perc_proposed_pre=s_pre/numvar;

s=min(s_pre,s_max);
perc_proposed=s/numvar;

if s>0
%for i=1:p
thr=order(s);
k=randsample(numvar,s);
%gamma=s/numvar;
%end;
%sp_thr=thr;
%else sp_thr=order(1)+1;
end

%A=zeros(p);
for i=1:p
    for j=i+1:p
    %for q=1:s
    if abs(h(i,j))>=sp_thr%yes(s);
        %ge(abs(h(i,j),order(s))==1; 
       A(i,j)=h(i,j);
    else A(i,j)=0;
    end;
           zs(i,j)=sign(A(i,j));
    %end;
    end;
end;
A 

%A=zeros(p);
%A(1,:)=[v(1) a(1:9)];
%for i=2:p
%    A(i,:)=[zeros(1,i-1) v(i) a(sum(p-i+1:p-1)+1:sum(p-i:p-1))];
%end;

for i=1:p
    for j=i+1:p
        A(j,i)=A(i,j);
        zs(j,i)=zs(i,j);
    end;
end;

A
trace(A)
rank(A)
cond(A)
svds(A)

aa=diag(A);
bb=diag(B);
max(aa);
max(bb);
plot(sort(aa,'descend'),'Color','m')
line(1:p,sort(bb,'descend'),'Color','r')
title('Diagonal elements of S* (magenta) and L* (red)');
xlabel('order');
ylabel('diagonal elements');
diffab=bb-aa;
plot(diffab);
cntab=0;
for i=1:length(diffab)
    if diffab(i)<0
       cntab=cntab+1;
end;
end;
cntab

%MA_sparse(:,:,g)=A;

P_B=B*inv(B'*B)*B';
N=eye(p);
for i=1:p
beta(i)=norm(P_B*N(:,i));
end;
inc=max(beta);
for i=1:p
v=0;
    for j=1:p
        if A(i,j)~=0
            v=v+1;
        end;
    end;
deg(i)=v;
v=0;
end;
deg_max_real=max(deg);
sort(deg,'descend')
deg_max_real*inc;
deg_max
s_max
s
s_pre
ident=inc*deg_max_real/p^2
1/108
ident_lag(:,:,lag+1)=ident;

        A_lag(:,:,lag+1)=A;
        zs_all(:,:,lag+1)=zs;
        s_top_lag(:,:,lag+1)=s_pre;
        end

    end
end
%end

%%

if sparse_diff==0
ident_lag=ident
end

if sparse_diff==1
for i=1:p
    for j=1:p
        A(i,j)=sum(A_lag(i,j,:));
    end
end

for lag=0:n_lag
   trace_A_lag(lag+1)=trace(A_lag(:,:,lag+1));
end

trace(A)
end

P_B=B*inv(B'*B)*B';
N=eye(p);
for i=1:p
beta(i)=norm(P_B*N(:,i));
end;
inc=max(beta);
for i=1:p
v=0;
    for j=1:p
        if A(i,j)~=0
            v=v+1;
        end;
    end;
deg(i)=v;
v=0;
end;
deg_max_real=max(deg);
sort(deg,'descend')
deg_max_real*inc;
deg_max
s_max
s
s_pre
ident=inc*deg_max_real/p^2
1/108

s_top_all=0;
for i=1:p
    for j=i+1:p
    if ne(A(i,j),0)==1 
        s_top_all=s_top_all+1;
    end
    end
end
s_top_all

if sparse_diff==1
for lag=0:n_lag
    s_true=0;
for i=1:p
    for j=i+1:p
    if ne(A_lag(i,j,lag+1),0)==1 
        s_true=s_true+1;
    end
    end
end
    s_true_lag(lag+1)=s_true
end
s_true_lag
end
%%

if sparse_diff ~=1
[U_A D_A]=svds(A,p);
end

if sparse_diff==0

for lag=0:n_lag
A_delay_S(:,:,lag+1)=sqrt(coef_lag(lag+1))*U_A;%sqrt(D_A)
end

end

if sparse_diff==1

for lag=0:n_lag
[U_A_pre D_A_pre]=svds(A_lag(:,:,lag+1),p);
% for i=1:p
%     for j=1:p
%         if U_A_pre(i,j)<1e-08
%            U_A_pre(i,j)=0;
%         end
% end
% end
U_A_lag(:,:,lag+1)=U_A_pre;
D_A_lag(:,:,lag+1)=D_A_pre;
A_delay_S(:,:,lag+1)=U_A_pre;%sqrt(D_A)
end

end

sumSigma=0;
for i=1:p
for j=(i+1):p
    sumSigma=sumSigma+abs(Sigma_0(i,j));
end;
end;
sumSigma
max(sum(Sigma_0))

sumB=0;
for i=1:p
for j=(i+1):p
    sumB=sumB+abs(B(i,j));
end;
end;
sumB
max(sum(B))

sumA=0;
for i=1:p
for j=(i+1):p
    sumA=sumA+abs(A(i,j));
end;
end;
sumA
max(sum(A))

tott=sumA+sumB
rapptrue=sumA/tott

%%

for h=1:length(f)

i=1;
    
if var_1_com==1 
   j=f(h);%j=(h-1)/T;%%
   E_top=E_r*E_r';%horzcat(E_r, zeros(p,p-r));
   %Lambda_top=E_r*vertcat([Lambda(1:r,1:r)*(1-2*(k_1))*(1-k_1^2) zeros(r,p-r)],zeros(p-r,p))*E_r';
   Lambda_top=E_r*Lambda(1:r,1:r)*(1-abs(k_1))*E_r';
   %(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))*(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))'
   %inv((eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))*(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))')
   MA_low(:,:,h,i)=inv((eye(p,p)+E_top*k_1*exp(-imm*2*pi*j)))*Lambda_top*inv((eye(p,p)+E_top*k_1*exp(-imm*2*pi*j))');
end

if var_2_com==1
   j=f(h);%j=(h-1)/T;%%
   E_top=E_r*E_r';%horzcat(E_r, zeros(p,p-r));
   %Lambda_top=E_r*vertcat([Lambda(1:r,1:r)*(1-2*(k_1))*(1-k_1^2) zeros(r,p-r)],zeros(p-r,p))*E_r';
   Lambda_top=E_r*Lambda(1:r,1:r)*(1-2*(k_1))*(1-k_1^2)*E_r';
   %(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))*(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))'
   %inv((eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))*(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))')
   MA_low(:,:,h,i)=inv((eye(p,p)+E_top*2*k_1*exp(-imm*2*pi*j)+E_top*k_1^2*exp(-imm*2*pi*2*j)))*Lambda_top*inv((eye(p,p)+E_top*2*k_1*exp(-imm*2*pi*j)+E_top*k_1^2*exp(-imm*2*pi*2*j))');
end

end

if var_1_com==0 && var_2_com==0

    %%
%imm=1i
%f=0:1/T:(0.5-1/T);
%spec_diag=zeros(p,length(f));
Gamma=zeros(p,r,length(f));
%b1=0.5;

for h=1:length(f)
    
    %j=f(h);%j=(h-1)/T;%%
   for lag=0:n_lag 
       j=f(h);
   Gamma_pre(:,:,h,lag+1)=A_delay(:,:,lag+1)*exp(-imm*2*pi*lag*j);
   for i1=1:p
       for i2=1:r
           Gamma(i1,i2,h)=sum(Gamma_pre(i1,i2,h,:));
       end
   end
   end
   
G0=A_delay(:,:,1);
if n_lag>0
G1=A_delay(:,:,2);
end
if n_lag>1
G2=A_delay(:,:,3);
end
end
   
if past_method==1
for h=1:length(f)
    j=f(h);%j=(h-1)/T;%%
    %for q=1:p
%spec_diag(q,h)=(1+b1*exp(-imm*h))*(1+b1*exp(-imm*h))';
    %end;
Gamma(:,:,h)=G0+G1*exp(-imm*2*pi*j)+G2*exp(-imm*2*pi*j*2);
end;
%spec
Gamma;

G0=A_delay(:,:,1);
if n_lag>0
G1=A_delay(:,:,2);
end
if n_lag>1
G2=A_delay(:,:,3);
end
end

%%

%MA_spec=zeros(p,p,length(f),length(c));
%traceMA_spec=zeros(length(f),length(c));
%condMA_spec=zeros(length(f),length(c));
%normMA_spec=zeros(length(f),length(c));

end

%for i=1:length(c)
i=1;
for h=1:length(f)
    if var_1_com==0 && var_2_com==0
    MA_low(:,:,h,i)=Gamma(:,:,h)*Lambda(:,:,i)*Gamma(:,:,h)';
    end
    traceMA_low(h,i)=trace(MA_low(:,:,h,i));
    eig_low=svds(MA_low(:,:,h),r);
    Eig_low(:,h)=eig_low;
    %condMA_low(h,i)=cond(MA_low(:,:,h))
    condMA_low(h,i)=eig_low(1)/eig_low(r);
    r_low(h,i)=rank(MA_low(:,:,h,i));
    normMA_low(h,i)=norm(MA_low(:,:,h,i),'fro');
end;
%end;
MA_low;
traceMA_low;
r_low;
condMA_low;
normMA_low;

%MA_low=MA_spec;
%traceMA_low=traceMA_spec;
%r_low=r_spec;
%condMA_low=condMA_spec;
%normMA_low=normMA_spec;

%%
if sparse_s==2
%%
if sparse_diff==0

for lag_n=n_lag:-1:0
for lag=0:lag_n
    u_lag=n_lag-lag_n
    AA_delay_S_all(:,:,lag+1,u_lag+1)=A_delay_S(:,:,lag+u_lag+1)*D_A*A_delay_S(:,:,lag+1)';
end
end

end

if sparse_diff==1
    
for lag_n=n_lag:-1:0
for lag=0:lag_n
    u_lag=n_lag-lag_n
    AA_delay_S_all(:,:,lag+1,u_lag+1)=A_delay_S(:,:,lag+u_lag+1)*sqrt(D_A_lag(:,:,u_lag+1))*sqrt(D_A_lag(:,:,lag+1)')*A_delay_S(:,:,lag+1)';
end
end

end

if past_method==1
%if sparse_diff==0
for lag=0:2
    AA_delay_S(:,:,lag+1)=A_delay_S(:,:,lag+1)*D_A*A_delay_S(:,:,lag+1)';
end

for lag=0:1
    AA_delay_1_S(:,:,lag+1)=A_delay_S(:,:,lag+2)*D_A*A_delay_S(:,:,lag+1)';
end

for lag=0:0
    AA_delay_2_S(:,:,lag+1)=A_delay_S(:,:,lag+3)*D_A*A_delay_S(:,:,lag+1)';
end

%end

if sparse_diff==1
for lag=0:2
    AA_delay_S(:,:,lag+1)=A_delay_S(:,:,lag+1)*D_A_lag(:,:,lag+1)*A_delay_S(:,:,lag+1)';
end

for lag=0:1
    AA_delay_1_S(:,:,lag+1)=A_delay_S(:,:,lag+2)*sqrt(D_A_lag(:,:,lag+2))*sqrt(D_A_lag(:,:,lag+1)')*A_delay_S(:,:,lag+1)';
end

for lag=0:0
    AA_delay_2_S(:,:,lag+1)=A_delay_S(:,:,lag+3)*sqrt(D_A_lag(:,:,lag+3))*sqrt(D_A_lag(:,:,lag+1)')*A_delay_S(:,:,lag+1)';
end

end

end

if past_method==1
size(AA_delay)
for lag=0:n_lag
r_lag(lag+1)=rank(AA_delay_S(:,:,lag+1));
end

size(AA_delay)
for lag=0:n_lag-1
r_lag_1(lag+1)=rank(AA_delay_1_S(:,:,lag+1));
end

for i=1:p
    for j=1:p
        Sigma_0_S(i,j)=sum(AA_delay_S(i,j,:));
        Sigma_1_S(i,j)=sum(AA_delay_1_S(i,j,:));
        Sigma_2_S(i,j)=sum(AA_delay_2_S(i,j,:));
    end
end
rank(Sigma_0_S)
rank(Sigma_1_S)
rank(Sigma_2_S)

size(Sigma_0_S)
rank(Sigma_0_S)
eigs(Sigma_0_S,p)
svds(Sigma_0_S,p)
svds(Sigma_1_S,p)
svds(Sigma_2_S,p)
trace(Sigma_0_S)
trace(Sigma_1_S)
trace(Sigma_2_S)
trace(Sigma_0_S)+trace(Sigma_1_S)+trace(Sigma_2_S)
cond(Sigma_0_S)
cond(Sigma_1_S)
cond(Sigma_2_S)

Sigma_0_S(1,2)-Sigma_0_S(2,1)

end


%Sigma=Sigma_eps+Sigma_chi
%Sigma_all=Sigma+Sigma_chi


for i=1:p
    for j=1:p
        for lag=0:n_lag
        Sigma_eps(i,j,lag+1)=sum(AA_delay_S_all(i,j,:,lag+1));
        end
    end
end

end


%Sigma_tot=Sigma+Sigma_eps;


%if sparse_diff==1
%    A=Sigma_0_S;
%   [U_A,D_A] =svds(Sigma_0_S)
%end

%B_past=B;
%B=E_r*diag(lambda)*E_r';


if var_1_com==1
Sigma_0=B;
Sigma_1=B*k_1;
end

if var_1==1
   Sigma_0_S=A;
   Sigma_1_S=A*k_1;%%OK!
end

E_top=E_r*E_r';
if var_2_com==1
   Sigma_0=B; 
   Sigma_1=Sigma_0*(E_top'*2*k_1)*inv(eye(p)-E_top'*k_1^2);
   Sigma_2=Sigma_1*(E_top'*2*k_1)+Sigma_0*(E_top'*k_1^2);
end


if var_2==1
   U_top=U_A*U_A';
   Sigma_0_S=A; 
   Sigma_1_S=A*(U_top'*2*k_1)*inv(eye(p)-U_top'*k_1^2);
   Sigma_2_S=Sigma_1_S*(U_top'*2*k_1)+Sigma_0_S*(U_top'*k_1^2);
end


%%

%[U_A D_A]=svds(A,p)

for h=1:length(f)

    i=1;
    
if var_1==1
   j=f(h);%j=(h-1)/T;%%
   U_top=U_A*U_A';%horzcat(E_r, zeros(p,p-r));
   D_top=U_A*D_A*(1-abs(k_1))*U_A';
   %(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))*(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))'
   %inv((eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))*(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))')
   MA_res(:,:,h,i)=inv((eye(p,p)+U_top*k_1*exp(-imm*2*pi*j)))*D_top*inv((eye(p,p)+U_top*k_1*exp(-imm*2*pi*j))');
end

if var_2==1
   j=f(h);%j=(h-1)/T;%%
   U_top=U_A*U_A';%horzcat(E_r, zeros(p,p-r));
   D_top=U_A*D_A*(1-2*(k_1))*(1-k_1^2)*U_A';
   %(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))*(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))'
   %inv((eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))*(eye(p,p)-E_top*k_1*exp(-imm*2*pi*j))')
   MA_res(:,:,h,i)=inv((eye(p,p)+U_top*2*k_1*exp(-imm*2*pi*j)+U_top*k_1^2*exp(-imm*2*pi*2*j)))*D_top*inv((eye(p,p)+U_top*2*k_1*exp(-imm*2*pi*j)+U_top*k_1^2*exp(-imm*2*pi*2*j))');
end

end

if var_1==0 && var_2==0

if sparse_s==0
MA_sparse=A;
A_delay_S(:,:,1)=U_A;
Sigma_eps=A;
n_lag=0;
end

if sparse_s==2
    
G0_res=A_delay_S(:,:,1);

if n_lag>0
G1_res=A_delay_S(:,:,2);
end

if n_lag>1
G2_res=A_delay_S(:,:,3);
end


%%
Gamma_res=zeros(p,p,length(f));
for h=1:length(f)
    %j=f(h);%j=(h-1)/T;%%
    j=f(h);
   for lag=0:n_lag 
       if sparse_diff==0
       Gamma_res_pre(:,:,h,lag+1)=A_delay_S(:,:,lag+1)*exp(-imm*2*pi*lag*j);
       end
       if sparse_diff==1
       Gamma_res_pre(:,:,h,lag+1)=A_delay_S(:,:,lag+1)*sqrt(D_A_lag(:,:,lag+1))*exp(-imm*2*pi*lag*j);
       end
   
   for i1=1:p
       for i2=1:p
           Gamma_res(i1,i2,h)=sum(Gamma_res_pre(i1,i2,h,:));
       end
   end
   end
end

if past_method==1
%imm=1i
%f=0:1/T:(0.5-1/T);
%spec_diag=zeros(p,length(f));
Gamma_res=zeros(p,p,length(f));
%b1=0.5;
for h=1:length(f)
    j=f(h);%%j=(h-1)/T;%%
    %for q=1:p
%spec_diag(q,h)=(1+b1*exp(-imm*h))*(1+b1*exp(-imm*h))';
    %end;
Gamma_res(:,:,h)=G0_res+G1_res*exp(-imm*2*pi*j)+G2_res*exp(-imm*2*pi*j*2);
end;
%spec
Gamma_res;
end

%MA_res=zeros(p,p,length(f),length(c));
%traceMA_res=zeros(length(f),length(c));
%condMA_res=zeros(length(f),length(c));
%normMA_res=zeros(length(f),length(c));

end
end

%for i=1:length(c)
i=1;
for h=1:length(f)
    if var_1==0 && var_2==0
        if sparse_diff==0
    MA_res(:,:,h,i)=Gamma_res(:,:,h)*D_A*Gamma_res(:,:,h)';
        end
        if sparse_diff==1
    MA_res(:,:,h,i)=Gamma_res(:,:,h)*Gamma_res(:,:,h)';
        end
    end
end

for h=1:length(f)
for i=1:p
    for j=1:p
        if sparse_diff==1
        zs_tot(i,j)=sum(abs(zs_all(i,j,:)));
        else
        zs_tot(i,j)=zs(i,j);    
        end
        %if i~=j && zs_tot(i,j)==0
        %MA_res(i,j,h)=0;
        %end
    end
end
end

sum(sum(abs(zs)/2))
if sparse_diff==1
sum(sum(abs(zs_all)/2))
s_top_lag
end

        
i=1;
for h=1:length(f)
    traceMA_res(h,i)=trace(MA_res(:,:,h,i));
    eig_sparse=svds(MA_res(:,:,h,i),p);
    Eig_res(:,h)=eig_sparse;
    %Eig_low(:,h)=eig_low;
    %condMA_res(h,i)=cond(MA_res(:,:,h))
    condMA_res(h,i)=eig_sparse(1)/eig_sparse(p);
    r_res(h,i)=rank(MA_res(:,:,h,i));
    normMA_res(h,i)=norm(MA_res(:,:,h,i),'fro');
end;
%end;
traceMA_res;
condMA_res;
normMA_res;    

%%

%f=0:1/T:(0.5-1/T)

if var_1==0 && var_2==0
if sparse_s==0
   for h=1:length(f)
       MA_res(:,:,h)=MA_sparse; 
   end
end
end

i=1;
for h=1:length(f)
    MA_spec_tot(:,:,h,i)=MA_low(:,:,h,i)+MA_res(:,:,h,i);
    alpha_th(h)=trace(MA_low(:,:,h,i))/trace(MA_spec_tot(:,:,h,i));
    traceMA_spec(h,i)=trace(MA_spec_tot(:,:,h,i));
    eig_spec=svds(MA_spec_tot(:,:,h),p);
    Eig_spec(:,h)=eig_spec;
    alpha_spec(h)=sum(eig_spec(1:r))/sum(eig_spec);
    %condMA_spec(h,i)=cond(MA_spec(:,:,h))
    condMA_spec(h,i)=eig_spec(1)/eig_spec(r);
    r_spec(h,i)=rank(MA_spec_tot(:,:,h,i));
    normMA_spec(h,i)=norm(MA_spec_tot(:,:,h,i),'fro');
end;
    
%MA_low=MA_spec;
%MA_spec=MA_spec_tot;

%for i=1:length(c)
%for h=1:length(f)
%    MA_spec(:,:,h,i)=MA_low(:,:,h,i)+MA_sparse;
%end;
%end; 

i=1;
for h=1:length(f)
r_tot(h,i)=rank(MA_spec_tot(:,:,h,i));
cond_tot(h,i)=cond(MA_spec_tot(:,:,h,i));
trace_tot(h,i)=trace(MA_spec_tot(:,:,h,i));
%r_low(h,i)=rank(MA_low(:,:,h,i));
end;
c

i=1;
for h=1:length(f)
S_11(h,i)=MA_spec_tot(1,1,h,i);
S_22(h,i)=MA_spec_tot(p,p,h,i);
S_12(h,i)=MA_spec_tot(1,p,h,i);
%end;
end;


subplot(1,3,1);
plot(f,S_11(:,1));
title('First component')
xlabel('frequency')
subplot(1,3,2);
plot(f,S_22(:,1));
title('Second component')
xlabel('frequency')
subplot(1,3,3);
plot(f,(S_12(:,1)).^2);
title('Squared cross spectrum')
xlabel('frequency')
%subplot(1,3,3)
%plot(f,S_11(:,1))

MA_spec=MA_spec_tot;
traceMA_spec
condMA_spec
normMA_spec
alpha_spec
Eig_spec
alpha_low=alpha_th
%%
% 
% subplot(1,3,1)
% plot(0:0.01:0.49,S_11(:,1));
% subplot(1,3,2)
% plot(0:0.01:0.49,S_22(:,1));
% subplot(1,3,3)
% plot(0:0.01:0.49,S_12(:,1));

for h=1:length(f)
    alpha_true(h)=trace(MA_low(:,:,h))/trace(MA_spec_tot(:,:,h));
    
    fac=0;
    for i=1:p
    for j=(i+1):p
    fac=fac+abs(MA_low(i,j,h));
    end;
    end;

    plus=0;
    for i=1:p
    for j=(i+1):p
    plus=plus+abs(MA_res(i,j,h));
    end;
    end;

    
    count_zeros=0;
    for i=1:p
    for j=(i+1):p
        if abs(MA_res(i,j,h))>1e-08
        count_zeros=count_zeros+1; 
        pos_res(i,j,h)=1;
        Res_out(count_zeros)=MA_res(i,j,h);
        end
    end;
    end;
    
    min_res_out(h)=min(abs(Res_out));
    
    nz_res(h)=count_zeros;
    
    lsum_true(h)=sum(diag(MA_low(:,:,h)));
    ssum_true(h)=sum(diag(MA_res(:,:,h)));
    tot_true(h)=lsum_true(h)+ssum_true(h);
    diagtot_true(h)=lsum_true(h)+tot_true(h);
    rappvar_true(h)=lsum_true(h)/(tot_true(h));
    corrtot_true(h)=fac+plus;
    rapptrue_tot(h)=plus/(fac+plus);
    
    
end

%%

alpha_true
tot_true
rappvar_true
corrtot_true
rappcorr_true=rapptrue_tot
nz_res
s_pre
s_top_all 

sum(abs(MA_res)>0.0001)

% for i=1:length(c)
% for h=1:length(f)
%     traceMA_low(h,i)=trace(MA_low(:,:,h,i));
%     eig_low=svds(MA_low(:,:,h),r);
%     %condMA_spec(h,i)=cond(MA_spec(:,:,h))
%     condMA_low(h,i)=eig_low(1)/eig_low(r);
%     r_low(h,i)=rank(MA_low(:,:,h,i));
%     normMA_low(h,i)=norm(MA_low(:,:,h,i),'fro');
% end;
% end;

MA_low;
traceMA_low
condMA_low
normMA_low
%rank_ave

subplot(1,1,1)
plot(f,Eig_low(1,:),'Color', 'red')
xlabel('frequency')
line(f,Eig_low(2,:),'Color', 'magenta')
if r>2
line(f,Eig_low(3,:),'Color', 'yellow')
end
line(f,Eig_res(1,:),'Color', 'blue')
line(f,min_res_out,'Color', 'green')
%if r>2
legend('\lambda_1(L)','\lambda_2(L)','||S||','min_{mag}(S)')
%end
if r>2
legend('\lambda_1(L)','\lambda_2(L)','\lambda_3(L)','||S||','min_{mag}(S)')
end
title('Setting features')

plot(f,rappvar_true,'Color', 'red')
line(f,rapptrue_tot,'Color', 'magenta')
legend('\beta','\zeta')
title('Structural features')
xlabel('frequency')
%%

%%
%sum(eig(A)<0)
%s
%s/numvar
numvar=p*(p-1)/2
s_lag=repmat(s_pre,1,length(f))
s_lag/numvar

s_top_lag/numvar

s_top_all/numvar

%A=eye(p)

%N_tot=100;

%sparse_s=1;
%sparse_s=2;
%T=100;

% A_delay_S(:,:,3)=zeros(p,p)
% A_delay(:,:,3)=zeros(p,r)
% D_A_lag(:,:,3)=zeros(p,p)

% %%%X=zeros(T,p);
% Xt=zeros(T,p);
% %mu=zeros(1,p);
% %sigma=zeros(1,p);
% F=zeros(T,p,length(c));
% %Store=zeros(T,p,length(c),N);
% %Kolm=zeros(T,p,length(c));
% %%!!!
% %Xf=zeros((length(f)),1);%%-1!!!
% epsilon=zeros(T,r);
% for n=1:N_tot
% for i=1:length(c)%% -->scegli il c!! sì!!!
%     %for q=1:p
%     %mu(1,q)=mean(X(:,q));
%     %sigma(1,q)=sqrt(B_yes(q,q,i));%%!!sì!!!
%     %%!!!sarà per radice di T!! sì!!!
%     %end;
%     %%bene perchè periodicob,èò.!! sì!!!
%     %epsilon(:,1)=sqrtM(:,:,i)*normrnd(zeros(p,1),diag(eye(p)));
% if new_method==0
% if sparse_s==1    
%     epsilon(1,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r));%%!! %%!!!X(k,:) via!!!
%     for h=1:length(f)
%         xi(:,h)=mvnrnd(zeros(p,1),T*MA_res(:,:,h));%by T?
%     end
%     for h=length(f)+1:2*length(f)
%         xi(:,h)=xi(:,length(f)+1+(length(f)-h));
%     end
%     xi_new(:,1:length(f))=xi(:,length(f)+1:2*length(f));
%     xi_new(:,length(f)+1:2*length(f))=xi(:,1:length(f));
%     Xt_res=ifft(xi_new','symmetric');
% end
% 
% if sparse_s==2 && var_1==0 && var_2==0
%     if sparse_diff==0
%     epsilon(1,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r));%%!! %%!!!X(k,:) via!!!
%     res_eps(1,:)=mvnrnd(zeros(1,p),D_A);    
%    %[U_A D_A]=svds(A);
%    Xt_res(1,:)=res_eps(1,:)*A_delay_S(:,:,1)';
%    res_eps(2,:)=mvnrnd(zeros(1,p),D_A);
%    Xt_res(2,:)=res_eps(2,:)*A_delay_S(:,:,1)'+res_eps(1,:)*A_delay_S(:,:,2)';
%    for k=3:T
%        res_eps(k,:)=mvnrnd(zeros(1,p),D_A);
%        Xt_res(k,:)=res_eps(k,:)*A_delay_S(:,:,1)'+res_eps(k-1,:)*A_delay_S(:,:,2)'+res_eps(k-2,:)*A_delay_S(:,:,3)';
%    end
%    end
%     if sparse_diff==1
%     epsilon(1,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r));%%!! %%!!!X(k,:) via!!!
%     res_eps(1,:)=mvnrnd(zeros(1,p),eye(p));    
%    %[U_A D_A]=svds(A);
%    Xt_res(1,:)=res_eps(1,:)*A_delay_S(:,:,1)'*sqrt(D_A_lag(:,:,1));
%    res_eps(2,:)=mvnrnd(zeros(1,p),eye(p));
%    Xt_res(2,:)=res_eps(2,:)*A_delay_S(:,:,1)*sqrt(D_A_lag(:,:,1))+res_eps(1,:)*A_delay_S(:,:,2)'*sqrt(D_A_lag(:,:,2));
%    for k=3:T
%        res_eps(k,:)=mvnrnd(zeros(1,p),eye(p));
%        Xt_res(k,:)=res_eps(k,:)*A_delay_S(:,:,1)'*sqrt(D_A_lag(:,:,1))+res_eps(k-1,:)*A_delay_S(:,:,2)'*sqrt(D_A_lag(:,:,2))+res_eps(k-2,:)*A_delay_S(:,:,3)'*sqrt(D_A_lag(:,:,3));
%    end
%    end   
% end
% 
%     if var_1_com==0 && var_2_com==0
%     epsilon(1,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r));%%!! %%!!!X(k,:) via!!!
%     %%Xt_res=mvnrnd(zeros(1,p),A);
%     Xt(1,:)=epsilon(1,:)*A_delay(:,:,1)'+Xt_res(1,:);
%     epsilon(2,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r));
%     Xt(2,:)=epsilon(2,:)*A_delay(:,:,1)'+epsilon(1,:)*A_delay(:,:,2)'+Xt_res(2,:);
%     for k=3:T
%     epsilon(k,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r));
%     Xt(k,:)=epsilon(k-2,:)*A_delay(:,:,3)'+epsilon(k-1,:)*A_delay(:,:,2)'+epsilon(k,:)*A_delay(:,:,1)'+Xt_res(k,:);%%!! %%!!!X(k,:) via!!!
%     %%bene perchè periodico!! sì!!!    
%     end
%     end
%     
%     if sparse_s==0
%        %Xt_res=mvnrnd(zeros(1,p),A);
%     epsilon(1,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r));%%!! %%!!!X(k,:) via!!!
%     res_eps(1,:)=mvnrnd(zeros(1,p),A);       
%     Xt(1,:)=epsilon(1,:)*A_delay(:,:,1)'+mvnrnd(zeros(1,p),A);
%     epsilon(2,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r));
%     Xt(2,:)=epsilon(2,:)*A_delay(:,:,1)'+epsilon(1,:)*A_delay(:,:,2)'+mvnrnd(zeros(1,p),A);
%     for k=3:T
%     %    for q=1:p
%     %       Xf(1,1)=B_yes(q,q,1,i);%%h=1-->f=0!!!sì!!!
%     %        for h=2:length(k_n)%%-1!!
%     %            j=(h-1)/T;
%     %            Xf(h,1)=exp(imm*k*j)*B_yes(q,q,h,i)+exp(imm*k*(-j))*B_yes(q,q,h,i);
%     %            %%è giusto??!!
%                 %%perchè primi alti??? dipende dai coseni!! sì!!! MF!!
%      %       end;
%      %   X(k,q)=sum(Xf);   
%       %  end;
%     %mu=zeros(1,p);
%     %for q=1:p
%     %mu(1,q)=mean(X(:,q));
%     %sigma(1,q)=sqrt(B_yes(q,q,i));%%!!sì!!!
%     %%!!!sarà per radice di T!! sì!!!
%     %end;
%     epsilon(k,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r));
%     Xt(k,:)=epsilon(k-2,:)*A_delay(:,:,3)'+epsilon(k-1,:)*A_delay(:,:,2)'+epsilon(k,:)*A_delay(:,:,1)'+mvnrnd(zeros(1,p),A);%%!! %%!!!X(k,:) via!!!
%     %%bene perchè periodico!! sì!!!
%     end;
%     end
% end
% 
% %%
%     
%     if new_method==1 && var_1==0 && var_2==0
%             for k=3:T
%                 
%                 if k<n_lag+1
%                 u_yes=n_lag+1-k;
%                 else u_yes=0;
%                 end
%                 
%                 if sparse_diff==0
%                 res_eps(k,:)=mvnrnd(zeros(1,p),D_A);   
%                 for lag=0:n_lag-u_yes
%                 res_eps_yes(k,lag+1,:)=res_eps(k-lag,:)*A_delay_S(:,:,lag+1)';
%                 end
%                 end;
%                 
%                                 
%                 if sparse_diff==1
%                 res_eps(k,:)=mvnrnd(zeros(1,p),eye(p));   
%                 for lag=0:n_lag-u_yes
%                 res_eps_yes(k,lag+1,:)=res_eps(k-lag,:)*A_delay_S(:,:,lag+1)'*sqrt(D_A_lag(:,:,lag+1));
%                 end
%                 end;%%ok, go!
%                 
%                 for j=1:p
%                 Xt_res(k,j)=sum(res_eps_yes(k,:,j));   
%                 end
%             end
%     end
% 
% if var_1==1
%     
%    res_eps(1,:)=mvnrnd(zeros(1,p),D_A*(1-abs(k_1))); 
%    Xt_res(1,:)=res_eps(1,:)*U_A;
%    
%    for k=2:T
%    res_eps(k,:)=mvnrnd(zeros(1,p),D_A*(1-abs(k_1)));
%    Xt_res(k,:)=-k_1*Xt_res(k-1,:)*U_A*U_A'+res_eps(k,:)*U_A;
%    end
%    
% end
% 
% 
% 
% if var_2==1
%    %h_1=1/k_1;
%    res_eps(1,:)=mvnrnd(zeros(1,p),D_A*(1-2*(k_1))*(1-k_1^2)); 
%    Xt_res(1,:)=res_eps(1,:)*U_A;
%    res_eps(2,:)=mvnrnd(zeros(1,p),D_A*(1-2*(k_1))*(1-k_1^2)); 
%    Xt_res(2,:)=-Xt_res(1,:)*U_A*U_A'*2*k_1+res_eps(2,:)*U_A;
%       for k=3:T
%    res_eps(k,:)=mvnrnd(zeros(1,p),D_A*(1-2*(k_1))*(1-k_1^2));
%    Xt_res(k,:)=-Xt_res(k-2,:)*U_A*U_A'*k_1^2-Xt_res(k-1,:)*U_A*U_A'*2*k_1+res_eps(k,:)*U_A;
%       end
% end
% 
%     if new_method==1 && var_1_com==0 && var_2_com==0
%                 for k=1:T
%                 %%
%                 epsilon(k,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r));%%!! %%!!!X(k,:) via!!!
% 
%                 for lag=0:n_lag-u_yes
%                 epsilon_yes(k,lag+1,:)=epsilon(k-lag,:)*A_delay(:,:,lag+1)';
%                 end%%OK. step by step!
%                 
%                 for j=1:p
%                 Xt_comp(k,j)=sum(epsilon_yes(k,:,j));
%                 Xt(k,j)=Xt_comp(k,j)+Xt_res(k,j);   
%                 end
%                 end
%     end
% 
% if var_1_com==1
%     
%    epsilon(1,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r)*(1-abs(k_1)));
%    Xt_comp(1,:)=epsilon(1,:)*E_r';
%    
%    for k=2:T
%        
%    epsilon(k,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r)*(1-abs(k_1)));
%    Xt_comp(k,:)=-k_1*Xt_comp(k-1,:)*E_r*E_r'+epsilon(k,:)*E_r';
%    
%    Xt(k,:)=Xt_res(k,:)+Xt_comp(k,:);
%    
%    end
%    
% end
% 
% 
% if var_2_com==1
%    epsilon(1,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r)*(1-2*(k_1))*(1-k_1^2));
%    Xt_comp(1,:)=epsilon(1,:)*E_r';
%    Xt(1,:)=Xt_comp(1,:)+Xt_res(1,:);
%    epsilon(2,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r)*(1-2*(k_1))*(1-k_1^2));
%    Xt_comp(2,:)=-Xt_comp(1,:)*E_r*E_r'*2*k_1+epsilon(2,:)*E_r';
%    Xt(2,:)=Xt_comp(2,:)+Xt_res(2,:);
%    for k=3:T
%    epsilon(k,:)=mvnrnd(zeros(1,r),Lambda(1:r,1:r)*(1-2*(k_1))*(1-k_1^2));
%    Xt_comp(k,:)=-Xt_comp(k-2,:)*E_r*E_r'*k_1^2'-Xt_comp(k-1,:)*E_r*E_r'*2*k_1+epsilon(k,:)*E_r';
%    Xt(k,:)=Xt_comp(k,:)+Xt_res(k,:);
%    end
% end
% 
% %Kolm(:,:,i)=X(:,:);    
% F(:,:,i)=Xt(:,:); %%plot!! +raw periodogram!! sì!!!
% F_res(:,:,i)=Xt_res(:,:);
% %%T=256!! then 10000 N!
% end;
% %Kolm;
% F;
% Store(:,:,:,n)=F(:,:,:);
% Store_res(:,:,:,n)=F_res(:,:,:);
% for qq=1:p
% var_ratio(n,qq)=var(Store_res(:,qq,1,n))/var(Store(:,qq,1,n));
% end
% end;
% Store;
% %plot(Kolm(:,:,3));
% subplot(1,1,1);
% plot(Store(:,5,1,5))
% %plot(Store(:,:,3,2)-%Kolm(:,:,3));
% size(Store)
% 
% sum(abs(Xt_res))/sum(abs(Xt))
% plot(Store(:,1),'red')
% xlabel('time')
% title('Common and idiosyncratic components')
% line(1:T,Xt_res(:,1),'Color','blue')
% legend('common','idio')
% var(Xt_res(:,1))/var(Store(:,1))
% var_ratio
% mean(mean(var_ratio))
% %%
% 
% i=1;
% f_orig=f%(T/2+1:T);
% 
% lippi=1;%1
% dan=0;
% bart=1;%1
% ave=0;
% if bart==1
%    dan=0;
% end
% if bart==0
%    dan=1;
% end
% if bart==1
% w=sqrt(T);
% end
% if bart==0
%    w=1; 
% end
% %f=f_orig(1:10:length(f))
% %f=f_orig(m+1:10:2*m+1)
% %f=f_orig(1+1:2:length(f_orig))
% %f_past=f;
% %f=1/T:(1/T):(0.5);
% %f=0:1/T:(0.5-1/T);
% %m=[(T/2)/50 (T/2)/25 (T/2)/12.5 (T/2)/6.25 (T/2)/3.125];%%%
% %m=[(T/2)/25 (T/2)/12.5 (T/2)/6.25 (T/2)/3.125 (T/2)/3.125*2];%%%
% %% [(T/2)/32] (T/2)/16 (T/2)/8 (T/2)/4 (T/2)/2 %%;%%T=64 (T/2)/32=1!!
% %S_Slep=zeros(p,p,length(f), length(m),length(c));
% %S_Sine=zeros(p,p,length(f), length(m),length(c));
% %S_ave=zeros(p,p,length(f), length(m),length(c));
% Store_per=zeros(p,p,length(f),N_tot);
% m=w;
% Store_ave=zeros(p,p,length(f),length(m),N_tot);
% %W=zeros(1,length(m));
% %TimeBW=zeros(1);
% %x=zeros(T,1);
% x1=zeros(T,1);
% x2=zeros(T,1);
% u_Slep=zeros(T,1);
% u_Sine=zeros(T,1);
% %u_ave=zeros(T,1);
% add1=zeros(T,1);
% add2=zeros(T,1);
% %Fs=T;
% %T=2*(length(f)-1);
% wv1=zeros(T,length(f));
% wv2=zeros(T,length(f));
% for h=1:length(f)
%     j=f(h);
%     for t=1:T
%         wv1(t,h)=exp(-imm*(j)*2*pi*t);
%         wv2(t,h)=exp(+imm*(j)*2*pi*t);
%     end;
% end;
% clear Store_per
% %i=5;%length(c)!!!
% for n=1:N_tot
%     %%FIX i HERE!!!!!
%     for h=1:length(f)
%         for q1=1:p
%             for q2=q1:p
%                 %j=(h-1)/T;
%                 %%ora prende le serie!! sì!!!
%                 x1(:,1)=Store(:,q1,i,n); x2(:,1)=Store(:,q2,i,n);
%                 %%un ciclo in T e uno in k!! sì!!
%                 %for k=1:K
%                 %for t=1:T
%                 add1=(1/(sqrt(T)))*x1.*wv1(:,h);%%2pi=1!!sì!!!
%                 add2=(1/(sqrt(T)))*x2.*wv2(:,h);
%                 %end;
%                 %addK(k,1)=sum(add1)*sum(add2);
%                 %end;
%                 %S_ave(q1,q2,h,b,i)
%                 Store_per(q1,q2,h,n)=sum(add1)*sum(add2);%(1/K)*sum(addK);
%                 if q1~= q2
%                     Store_per(q2,q1,h,n)=Store_per(q1,q2,h,n)';
%                 end;
%             end;
%         end;
%     end;
% end;
% % Store_per_new(:,:,1:(m),:)=Store_per(:,:,m+2:2*m+1,:);
% % Store_per_new(:,:,m+1:2*m+1+m,:)=Store_per(:,:,1:2*m+1,:);
% % Store_per_new(:,:,3*m+2:4*m+1,:)=Store_per(:,:,1:m,:);
% Store_per_new(:,:,1:length(f)-1,:)=Store_per(:,:,length(f):-1:2,:);
% Store_per_new(:,:,length(f):2*length(f)-1,:)=Store_per(:,:,:,:);
% Store_per_new(:,:,2*length(f):3*length(f)-2,:)=Store_per_new(:,:,1:length(f)-1,:);
% size(Store_per_new)%3m-2i
% 
% cov_plus=0;
% %round(nthroot(T,3));%floor((2*m+1)/4)
% i=1;
% clear Store_per_new_weight Store_ave cov_Store_k Add_cov_Store
% for b=1:length(m) 
%     %K(b)=2*w(b)+1;%%<--questa è l'equivalenza!! sì!!!
%     %%punto chiave!! sì!!!%%
%     for n=1:N_tot
%                        %mu_Store=mean(Store(:,:,i,n));
%  
%                        for l=1:p
%                            mu_Store(l)=mean(Store(:,l,i,n));
%                        end
%                        
%                        for k=1:T
%                        Store_cent(k,:,i,n)=Store(k,:,i,n)-mu_Store;
%                        end
%                        
%                        if bart==1
%                        if lippi==1 
%                           %clear spec_lippi
%                        spec_lippi(:,:,:,b,n)=spectral_x(Store_cent(:,:,i,n),p,length(f)-1,length(f)-1);
%                        end
%                        end
%                        
%         for h=1:1:length(f)
%                     j=f(h);
%                     if bart==1
%                        if lippi==0
%                        for k=length(f):-1:-length(f)+2
%                            if k>=1
%                            cov_Store_k(:,:,k+length(f)-1)=(Store_cent(k:T,:,i,n))'*(Store_cent(1:T+1-k,:,i,n))/(T-k);
%                            end
%                            if k<1
%                            k_neg=abs(2-k);
%                            if cov_plus==1
%                            cov_Store_k(:,:,k+length(f)-1)=(Store_cent(1:T+1-k_neg,:,i,n))'*(Store_cent(k_neg:T,:,i,n))/(T-k_neg);
%                            end
%                            if cov_plus==0
%                            cov_Store_k(:,:,k+length(f)-1)=cov_Store_k(:,:,length(f)+1-k)';
%                            end                           
%                            end
%                        end
%                        for k=-length(f)+1:1:length(f)-1
%                        Add_cov_Store(:,:,k+length(f))=cov_Store_k(:,:,k+length(f))*(1-(abs(k)/length(f)))*exp(-2*pi*imm*j*k);
%                        %exp(-sqrt(-1)*(-m:m)'*(-2*pi*h/h:2*pi/h:2*pi*h/h))
%                        end
%                        for q1=1:p
%                        for q2=q1:p
%                            Store_ave(q1,q2,h,b,n)=sum(Add_cov_Store(q1,q2,:));
%                            if q2>q1
%                            Store_ave(q2,q1,h,b,n)=Store_ave(q1,q2,h,b,n)';
%                            end
%                        end
%                        end
%                        end
%                     end              
%   
%             %Store_per(:,:,h,b,n)=S_ave(:,:,h,b,i);%%!!
%             %i1=(m+h)-(w);%max(1,h-(m(b)-1));
%             %i2=(m+h)+(w);
%             
%             for q1=1:p
%                 for q2=q1:p
%                     x1=Store(:,q1,i,n); x2=Store(:,q2,i,n);
%                     %%Slepian Tapers
%                     %%un ciclo in T e uno in k!! sì!!
%                     %%Averaged periodogram
%                     
%                     
%                     %i1=max(1,h-(m(b)-1));
%                     %i2=min(length(f),h+(m(b)-1));
%                     
%                     if dan==1                                            
%                         i1=(length(f)+h)-(w);%max(1,h-(m(b)-1));
%                         i2=(length(f)+h)+(w);
%                         
%                         %i1=max((length(f)+h)-(w),length(f));%max(1,h-(m(b)-1));
%                         %i2=min((length(f)+h)+(w),2*length(f)+1);
%                             
%                         for ii=i1:i2
%                             Store_per_new_weight(q1,q2,ii,n)=(1-(abs(ii-(length(f)+h))/(w+1)))*Store_per_new(q1,q2,ii,n);
%                         end
%                         
%                         Store_ave(q1,q2,h,b,n)=sum(Store_per_new_weight(q1,q2,i1:i2,n));
%                         Store_ave(q2,q1,h,b,n)=Store_ave(q1,q2,h,b,n)';
%                     end
%                     
%                     if ave==1
% %                         base_f=(length(f)-1);
% %                         i1=max(base_f+1,base_f+h-(w-1));
% %                         i2=min(base_f+length(f),base_f+h+(w-1));                            
%                         
%                         Store_per_new_weight=Store_per_new;
%                         
%                         Store_ave(q1,q2,h,b,n)=(1/(i2-i1+1))*sum(Store_per_new_weight(q1,q2,i1:i2,n));
%                         Store_ave(q2,q1,h,b,n)=Store_ave(q1,q2,h,b,n)';
%                     end
%                     
%                     %if trng==1
% 
%                     %%basta questo!!sì!!!
%                     %end
%                     
%                 end;
%             end;
%             if lippi==1;
%             Store_ave(:,:,h,b,n)=spec_lippi(:,:,h+length(f)-1,b,n);
%             end
%             eigs_now=svds(Store_ave(:,:,h,b,n),p);
%             eig_ave(:,h,n)=eigs_now;
%             cond_ave(h,n)=eigs_now(1)/eigs_now(r);
%             alpha_ave(h,n)=sum(eigs_now(1:r))/sum(eigs_now);
%             rank_ave(h,n)=rank(Store_ave(:,:,h,b,n));
%         end
%     end;
% end;
% Lambda(:,:,i)
% %MA_spec(:,:,15,i)
% %Store_ave(:,:,15,1,2)
% 
% eig_ave%5,1,:,1:r)'
% cond_ave(:,1,:)
% alpha_ave(:,1,:)
% alpha_low
% %subplot(1,2,1)
% %plot(alpha_ave(50,1,:))
% %subplot(1,2,2);
% %plot(f,S_11(:,1));
% % subplot(1,2,2);
% % plot(f,S_22(:,1));
% plot(alpha_ave(:,100,1))
% mean(alpha_ave(:,100,1))
% %lines(1:length(f),alpha_true)
% 
% 
% Store_ave_m=zeros(p,p,length(f),length(m));
% %for i=5:5%%length(c)
%     for b=1:length(m)
%         for h=1:length(f)
%             for q1=1:p
%             for q2=1:p
%                 Store_ave_m(q1,q2,h,b)=sum(Store_ave(q1,q2,h,b,:))/N_tot;
%             end
%             end;
%         end;
%     end;
% %end;
% Lambda(:,:,i)
% %MA_spec(:,:,17,i)%%
% %Store_ave_m(:,:,17,1)
% 
% for h=1:length(f)
%     err_per(h)=norm(Store_per(:,:,h)-MA_spec(:,:,h),'fro');
% end
% 
% for n=1:N_tot
% for h=1:length(f)
%     err_ave(h,n)=norm(Store_ave(:,:,h,b,n)-MA_spec(:,:,h),'fro');
% end
% end
% 
% subplot(1,3,1)
% plot(f,mean(alpha_ave'),'Color', 'red')
% %legend('\beta')
% title('Input \beta')
% xlabel('frequency')
% 
% subplot(1,3,2)
% plot(f,mean(cond_ave'),'Color', 'red')
% title('Input c')
% %axis('frequency')
% xlabel('frequency')
% 
% subplot(1,3,3)
% plot(f,mean(rank_ave'),'Color', 'red')
% title('Input r')
% xlabel('frequency')
% 
% %C
% %Store_ave(:,:,h_yes,3)
% 
% %err_ave(h_yes,3)-TL(fin1,fin2)
% %%OK!
% 
% %%
% jump=1;
% if jump==0
% for ncomp=1:min(p,T)
% for n=1:N_tot
% for b=1:length(m)
% for h=1:length(f)
%     [U_spec_ave(:,:,h,b,n) D_spec_ave(:,:,h,b,n)]=svds(Store_ave(:,:,h,b,n),min(T,p));
% end
% end
% end
% 
% for n=1:N_tot
% for b=1:length(m)
% for lag=0:2
%     for h2=1:length(f)
%         j2=(h2-1)/T;
%     for h=1:length(f)
%         j1=(h-1)/T;
%         ADD2_ave(:,:,h,b,n)=U_spec_ave(:,ncomp,h,b,n)*exp(1i*lag*j1*2*pi);
%     end
%     for i=1:p
%         for j=1:1
%         add2_ave(i,j)=sum(ADD2_ave(i,j,:,b,n));
%         end
%     end
%     ADD_k2_ave(:,:,lag+1,h2,b,n)=add2_ave*exp(1i*lag*j2*2*pi);
%     end
% end
% 
% for h=1:length(f)
% for i=1:p
%     for j=1:1
%         MA_spec_tot_dyn2_ave(i,j,h,b,n)=sum(ADD_k2_ave(i,j,:,h,b,n));
%     end
% end
% end
% 
% for h=1:length(f)
% %Dynamic_eig(ncomp,h,b,n)=MA_spec_tot_dyn(:,1,h,b,n)'*MA_spec_tot(:,:,h,b,n)*MA_spec_tot_dyn2(:,1,h,b,n);
% Dynamic_eig_def_ave(ncomp,h,b,n)=U_spec_ave(:,ncomp,h,b,n)'*D_spec_ave(ncomp,ncomp,h,b,n)*U_spec_ave(:,ncomp,h,b,n);
% end
% end
% end
% end
% 
% for n=1:N_tot
% for b=1:length(m)
% for h=1:length(f)
% %for n=1:N_tot
% %for b=1:length(m)
% Add_Def_ave(:,:,h,b,n)=U_spec_ave(:,1:min(T,p),h,b,n)*D_spec_ave(1:min(T,p),1:min(T,p),h,b,n)*U_spec_ave(:,1:min(T,p),h,b,n)';
% Add_TOT_ave(:,:,h,b,n)=U_spec_ave(:,:,h,b,n)*D_spec_ave(:,:,h,b,n)*U_spec_ave(:,:,h,b,n)';
% %end
% %end
% end
% 
% 
% %for n=1:N_tot
% %for b=1:length(m)
% for i=1:p
%     for j=1:p
%         ADD_yes_ave(i,j,b,n)=sum(Add_Def_ave(i,j,:,b,n));
%         ADD_TOT_ave(i,j,b,n)=sum(Add_TOT_ave(i,j,:,b,n));
%     end
% end
% 
% [U_yes_ave(:,:,b,n), D_yes_ave(:,:,b,n)]=svds(ADD_yes_ave(:,:,b,n),min(T,p));
% 
% end
% end
% 
% 
% for ncomp=1:10
% for n=1:N_tot
%     dyn1(ncomp,n)=D_yes_ave(ncomp,ncomp,1,n);
% end
% end
% 
% ksdensity(dyn1(1,:))
% col_letter=['y', 'm', 'c', 'r', 'g', 'b', 'w' ,'k'];
% for ncomp=2:9
% line(1:N_tot,ksdensity(dyn1(ncomp,:)),'Color',(col_letter(ncomp-1)))
% end
% kstest(dyn1(ncomp,:))
% 
% for ncomp=1:10
% subplot(2,5,ncomp)
% ksdensity(dyn1(ncomp,:))
% end
% 
% rank_ave
% end
% %%clear
% 
% %%
% 
% Sigma_0
% Sigma_1
% if var_1~=1
% Sigma_2
% end
% 
% if past_method==1
% Sigma_0_S
% Sigma_1_S
% Sigma_2_S
% else
% Sigma_eps
% end
% %%
% 
% zs_tot
% sum(abs(MA_res)>0.001)
% sum(abs(MA_res)>0.001)-1
% sum(sum(abs(MA_res)>0.001)-1)/2
% if sparse_diff==1
% s_top_lag
% end
% s_top_all
% ident_lag
% 
% p
% T
% rappvar_true
% rapptrue_tot
% 
% delta
% deltabis
% 
% MA_low
% MA_res
% %%

normMA_low
normMA_res

r
s_top_all

alpha_true
rappcorr_true

subplot(1,1,1)
plot(f,Eig_low(1,:),'Color', 'red')
xlabel('frequency')
line(f,Eig_low(2,:),'Color', 'magenta')
if r>2
line(f,Eig_low(3,:),'Color', 'yellow')
end
line(f,Eig_res(1,:),'Color', 'blue')
line(f,min_res_out,'Color', 'green')
%if r>2
legend('\lambda_1(L)','\lambda_2(L)','||S||','min_{mag}(S)')
%end
if r>2
legend('\lambda_1(L)','\lambda_2(L)','\lambda_3(L)','||S||','min_{mag}(S)')
end
title('Setting features')