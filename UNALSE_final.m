% UNALSE_final.m carries out UNALSE computational routine, introduced in
%   
% Barigozzi, M. and Farnè, M. (2021), 'An algebraic estimator 
% for large spectral density matrices.
%
% The INPUT arguments are: 
% 'Data': T (time points) times p (variables) data matrix
% 'm': smoothing parameter (to be passed to spectral_x: defaults to [sqrt(T)])
% 'r_thr_input': input latent eigenvalue threshold grid parameter across frequencies
% 's_thr_input': input sparsity threshold grid parameter across frequencies
%
% The OUTPUT argument is a list containing the following fields:
% 'freq': the used frequencies in radians
% 'L': latent spectral density matrix estimates across frequencies
% 'S': residual spectral density matrix estimates across frequencies
% 'Sigma': overall spectral density matrix estimates across frequencies
% 'psi_opt': optimal eigenvalue threshold at each frequency according to MC criterion
% 'rho_opt': optimal sparsity threshold at each frequency according to MC criterion
% 'psi_ind': grid index of the selected eigenvalue threshold by MC criterion
% 'rho_ind': grid index of the selected sparsity threshold by MC criterion
% 'rank': estimated rank across frequencies
% 'nonzeros_prop': estimated proporton of residual non-zeros across frequencies
% 'latent_var_prop': estimated proportion of latent variance across frequencies
% 'res_cov_prop': estimated proportion of residual
% covariance across frequencies
% 'rank_all': estimated rank across all threshold pairs and frequencies
% 'nonzeros_all': estimated number of nonzeros across all threshold pairs and frequencies
% 'pos_def_ind': 0 if the overall or the residual estimate across frequencies
% is not positive definite
% 'neg_res_eig': number of negative eigenvalues in the residual estimate
% across frequencies
% 'neg_all_eig': number of negative eigenvalues in the overall estimate
% across frequencies

function[Out_UNALSE]=UNALSE_final(Data,m,r_thr_input,s_thr_input,num_thr)

% sample size and dimension
[T,p] = size(Data);

% r_thr_input=1
% s_thr_input=0.5
% num_thr=10

% total number of off-diagonal elements
numvar=p*(p-1)/2;

% smoothing parameter
m=min(m,floor(sqrt(T)));

% Bartlett's input estimates
[Sigma_X,f_orig_rad] = spectral_x(Data, p, m, m);

% stored on positive frequencies (starting from 0)
Store_ave=Sigma_X(:,:,m+1:2*m+1);

% extract relevant frequencies
f_top=f_orig_rad(m+1:2*m+1);

% maximum number of iterations
N_max=1000;

% loop in frequencies
for hh=1:length(f_top)
    
        h_yes=hh;
        
        % input
        M=(Store_ave(:,:,h_yes,1));
        M_star=(Store_ave(:,:,h_yes,1));
        
        % grid for psi
        
        alpha_M=sqrt(sqrt(r_thr_input(hh)/p));
        psi_orig=(sqrt(p/T)/(2*alpha_M)):1/num_thr*(sqrt(p/T)/(2*alpha_M)):(sqrt(p/T)/(alpha_M));
        gamma=s_thr_input(hh)*(p^(4/8)/p:1/num_thr*(p^(3/4)/p-p^(4/8)/p):p^(3/4)/p);

        % loop in sparsity thresholds
        
        for t1=1:num_thr
            
            % loop in eigenvalue thresholds
            
            for t2=1:num_thr
                
                % algorithm convergence statistics
                
                clear arr al criterion
                arr=zeros(1,N_max);
                al=zeros(1,N_max);
                criterion=zeros(1,N_max);
                k=1;
                
                arr(1)=1;
                criterion(1)=1;
                al(1)=1;
                
                % initializers
                
                L_Thr=diag(diag(M))/2;
                S_Thr=diag(diag(M))/2;
                E=zeros(p,p);
                Y=L_Thr;
                Z=L_Thr;
                
                % singular value thresholding plus soft thresholding
                while k<N_max && abs(criterion(k))>1.0e-02
                    
                    % differential of the smooth part
                    L=Y-1/2*(Y+Z-M);
                    
                    % SVT of the differential
                    [U,D]=svds(L,rank(L));
                    d_coef= gini(repmat(1,rank(L),1),diag(D));
                    psi=(1/d_coef)*psi_orig*trace(D)/p;
                    sparse_real=gamma*max(psi_orig);
                    
                    % SVT checks
                    
                    if t1*t2==1
                        clear D_Thr
                    end
                    
                    if rank(D)>0
                       d=1;
                    end
                    
                    if t1*t2>1 && d==1
                        clear D_Thr
                    end
                    
                    % eigenvalue thresholding
                    for i=1:rank(D)
                        D_Thr(i,i)=max(D(i,i)-psi(t2),0);
                    end;
                    
                    r_Thr=rank(D_Thr);
                    rank_Thr(t1,t2)=r_Thr;
                    L_pre=L_Thr;
                    % low rank current estimate
                    L_Thr=U(:,1:r_Thr)*D_Thr(1:r_Thr,1:r_Thr)*U(:,1:r_Thr)';
                    % convergence criterion for L
                    add1=norm(L_Thr-L_pre,'fro')/(1+norm(L_pre,'fro'));
                    
                    % 
                    S=Z-1/2*(Y+Z-M);
                    S_pre=S_Thr;
                    M_pre=L_pre+S_pre;
                    
                    % soft thresholding
                    for i=1:(p)
                        for j=i+1:p
                            S_Thr(i,j)=(S(i,j))/abs(S(i,j))*max(abs(S(i,j))-sparse_real(t1),0);
                        end;
                     end
                    for i=1:p
                            S_Thr(i,i)=S(i,i);
                    end;
                    for i=2:p
                            for j=1:(i-1)
                                S_Thr(i,j)=S_Thr(j,i);
                            end;
                    end;
                       
                    % convergence criterion for S
                    add2=norm(S_Thr-S_pre,'fro')/(1+norm(S_pre,'fro'));
                    M_star=L_Thr+S_Thr;
                    E=M_star-M;
                    k=k+1;
                    arr(k)=norm(E)/norm(M);
                    
                    % overall convergence criterion
                    criterion(k)=add1+add2;
                    
                    % solution updates
                    al(k)=(1+sqrt(1+4*al(k-1)^2))/2;
                    Y=L_Thr+((al(k-1)-1)/al(k))*(L_Thr-L_pre);
                    Z=S_Thr+((al(k-1)-1)/al(k))*(S_Thr-S_pre);
                    
                end
                
                % unshrinkage of thresholded eigenvalues
                for i=1:r_Thr
                    D_Thr(i,i)=D_Thr(i,i)+psi(t2);
                end;
                L_Thr=U(:,1:r_Thr)*D_Thr(1:r_Thr,1:r_Thr)*U(:,1:r_Thr)';
                for i=1:p
                    S_Thr(i,i)=M_star(i,i)-L_Thr(i,i);
                end;
                    
                % non zeros count
                v=0;
                for i=1:(p-1)
                    for j=(i+1):p
                        if S_Thr(i,j)~=0
                           v=v+1;
                        end;
                    end;
                end;
                    
                % statistics for each threshold pair
                
                Ar(t1,t2)=arr(k);
                Al(t1,t2)=al(k);
                Crit(t1,t2)=criterion(k);
                K(t1,t2)=k;
                nz(t1,t2)=v;
                Low(:,:,t1,t2)=L_Thr;
                Sparse(:,:,t1,t2)=S_Thr;
                UU(:,1:r_Thr,t1,t2)=U(:,1:r_Thr);
                rank_Thr(t1,t2)=r_Thr;
                Sigma_hat(:,:,t1,t2)=L_Thr+S_Thr;
                linf_s(t1,t2)=norm(S_Thr,Inf);
                l2_s(t1,t2)=norm(L_Thr);
                scale(t1,t2)=sparse(t1)/psi(t2);
                t_L(t1,t2)=trace(Low(:,:,t1,t2));
                t_S(t1,t2)=trace(Sparse(:,:,t1,t2));
                t_TOT(t1,t2)=trace(Low(:,:,t1,t2))+trace(Sparse(:,:,t1,t2));
                rappvar(t1,t2)=t_L(t1,t2)/t_TOT(t1,t2);
                fac=0;
                for i=1:p
                    for j=(i+1):p
                        fac=fac+abs(L_Thr(i,j));
                    end
                end;
                plus=0;
                for i=1:p
                    for j=(i+1):p
                        plus=plus+abs(S_Thr(i,j));
                    end;
                end;
                lsum=sum(diag(L_Thr));
                ssum=sum(diag(S_Thr));
                diagtot(t1,t2)=lsum+ssum;
                tott(t1,t2)=plus+fac;
                rappcorr(t1,t2)=plus/tott(t1,t2);
                
                % count negative eigenvalues
                
                negeigSp(t1,t2)=sum(eig((Sparse(:,:,t1,t2)))<0);
                negeigSigma(t1,t2)=sum(eig((Sigma_hat(:,:,t1,t2)))<0);
                               
                % MC: criterion for threshold selection
                
                if linf_s(t1,t2)<l2_s(t1,t2) && negeigSp(t1,t2)==0 && negeigSigma(t1,t2)==0
                   pos_def(t1,t2)=1;
                   MC(t1,t2)=max(linf_s(t1,t2)/(scale(t1,t2)*(1-rappvar(t1,t2))*trace(M)),(rank_Thr(t1,t2)*l2_s(t1,t2))/(trace(M)*rappvar(t1,t2)));
                   difflinfyes(t1,t2)=MC(t1,t2)-linf_s(t1,t2)/(scale(t1,t2)*(1-rappvar(t1,t2))*trace(M));
                   diffl2yes(t1,t2)=MC(t1,t2)-(rank_Thr(t1,t2)*l2_s(t1,t2))/(trace(M)*rappvar(t1,t2));
                else
                   pos_def(t1,t2)=0;
                   MC(t1,t2)=norm(M,'Inf');    
                end
            end
        end
        
       % threshold pair selection
        [minmin1 ind_1]=min(MC);
        [minmin2 ind_2]=min(min(MC));
       
       % selected thresholds indicators
       fin1=ind_1(ind_2);
       fin2=ind_2;
       
       % selected thresholds at each frequency
       Psi(hh)=psi(fin2);
       Rho(hh)=sparse_real(fin1);
          
       rho_ind(hh)=fin1;
       psi_ind(hh)=fin2;
       
       % solutions at each frequency
       Low_all(:,:,hh)=Low(:,:,fin1,fin2);
       Sparse_all(:,:,hh)=Sparse(:,:,fin1,fin2);
       Sigma_hat_all(:,:,hh)=Sigma_hat(:,:,fin1,fin2);
        
       % algebraic statistics
       r_all(hh)=rank_Thr(fin1,fin2);
       nz_all(hh)=nz(fin1,fin2);
       
       rappvar_all(hh)=rappvar(fin1,fin2);
       rappcorr_all(hh)=rappcorr(fin1,fin2);
              
       rank_Thr_all(:,:,hh)=rank_Thr;
       nz_Thr_all(:,:,hh)=nz;
       
       % negative eigenvalues
       negeigSp_all(hh)=negeigSp(fin1,fin2);
       negeigSigma_all(hh)=negeigSigma(fin1,fin2);
       pos_def_ind(hh)=pos_def(fin1,fin2);
       %
      
       % frequency index display
       
       string = ['frequency ',num2str(hh), ' out of ',num2str(length(f_top))];
       disp(string)
end

% solution output list
Out_UNALSE={'freq',f_top,'L',Low_all,'S',Sparse_all,'Sigma',Sigma_hat_all,'psi_opt',Psi,'rho_opt',Rho,'psi_ind',psi_ind,'rho_ind',rho_ind,'rank',r_all,'nonzeros_prop',nz_all/numvar,'latent_var_prop',rappvar_all,'res_cov_prop',rappcorr_all,'rank_all',rank_Thr_all,'nonzeros_all',nz_Thr_all,'pos_def_ind',pos_def_ind,'neg_res_eig',negeigSp_all,'neg_all_eig',negeigSigma_all};

end
