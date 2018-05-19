function VBS_syn_main()
close all;clc;clear all;
nmiset=zeros(11,1);
times=1;
for p=0:0 %10
    for i=1:times
%% Data input
    [A,Q] = GenerateSignAYang(4,32,32,0.5,0.5,p*0.05);
%      A=load('war.mat');
%      A=A.g_A;
    imagesc(A)
    k_min=4; k_max=4;
    evidenceset=[];s_tao=[];s_rho=[];s_mu=[];s_evidence=-10^10;
    k_t=0;
    tic
for k=k_min:k_max
    k_t=k_t+1; timess=1;
    [k_tao,k_rho,k_mu,k_evidence]=ssbm_times(A,k,timess);
    if k==k_min
        s_tao=k_tao;
        s_evidence=k_evidence;
        s_rho=k_rho;
        s_mu=k_mu;
    end
    evidenceset(k_t)=k_evidence;
    if k_evidence>s_evidence
        s_tao=k_tao; s_rho=k_rho;s_mu=k_mu;s_evidence=k_evidence;
    end
end
toc
% figure;plot(likeset,'-o');xlabel('iterations');ylabel('Cost function');title('Change of Cost funtion')
Z = comp_Z(s_tao);
[TA,C_number,S_C,index] =Draw(A,Z,1);
%% Compute NMI
    % generate the label of nodes
    Z=full(Z); z=zeros(length(Z),1);
    for i=1:size(Z,2)
        for j=1:length(Z)
            if Z(j,i)==1 z(j)=i; end
        end
    end
    % compute nmi
    label=[];
    n_num=32;%set the number of each block
    for i=1:4
        label((i-1)*n_num+1:i*n_num)=i;        
    end
%     label=load('SPP_c.txt');
    acc=compute_NMI(z,label);    
    disp(['NMI = ',num2str(acc)]);
%---
    nmiset(p+1)=nmiset(p+1)+acc;
    nmiset'
    end
end
nmiset=nmiset/times;

function [t_tao,t_rho,t_mu,t_evidence]=ssbm_times(A,k,times)
    t_evidence=-10^10;
    t_tao=[];t_rho=[];t_mu=[];
    for it=1:times
        it;
        [t_taot,t_rhot,t_mut,t_evidencet]=ssbm(A,k);
        if it==1
           t_tao=t_taot;
           t_evidence=t_evidencet;
           t_rho=t_rhot;
           t_mu=t_mut;
        end
        if t_evidence<t_evidencet
            t_evidence=t_evidencet;
            t_tao=t_taot;
            t_rho=t_rhot;
            t_mu=t_mut;
        end   
    end

function [tao,rho,mu,evidence]=ssbm(A,k)
%% initial
   n=length(A);
   rho_0=ones(1,k)*1/2;% non-informative Jeffreys prior distribution
   mu_0=ones(k,k,3)*1/2; % non-informative Jeffreys prior distribution
   evidence=0;
   evidencet=0;
   torf=1;   
   rho=ones(1,k)*0.5;
   mu=ones(k,k,3)*0.5; %1-positive link 2-no link 3-nagetive link
   e_p=sum(sum(sum(A+abs(A))))/(2*(length(A)*(length(A)-1)))
   e_n=sum(sum(sum(-A+abs(A))))/(2*(length(A)*(length(A)-1)))
   e_z=1-e_p-e_n
   m_e=max([e_p,e_n,e_z]);

   for i=1:k
       for j=1:k
           if i==j               
               mu(i,j,1)=0.4;
               mu(i,j,2)=0.5;
               mu(i,j,3)=0.1;
           else
               mu(i,j,1)=0.1;
               mu(i,j,2)=0.5;
               mu(i,j,3)=0.2;
           end
       end
   end
%     tao=zeros(n,k);
%     tao=kmeans_com(A,k);    
    tao=rand(n,k);
    tao=tao./(sum(tao,2)*ones(1,k));
    taot=tao;

%%-initializing parameters end

while torf
   

%% compute the approximate distribution of Z 
for i=1:n
    for r=1:k
        % tao(i,r)=exp(psi(rho(r))-psi(sum(rho)));
        tao(i,r)=psi(rho(r))-psi(sum(rho));
        %  disp(['tao:',num2str(tao(i,r))])
        for j=1:n
            if j~=i
               for rr=1:k
                   if A(i,j)==1
                       % tao(i,r)=tao(i,r)*exp(tao(j,rr)*(psi(mu(r,rr,1))-psi(sum(mu(r,rr,:)))));
                       tao(i,r)=tao(i,r)+tao(j,rr)*(psi(mu(r,rr,1))-psi(sum(mu(r,rr,:))));
                       % disp(['psi1:',num2str(psi(mu(r,rr,1))-psi(sum(mu(r,rr,:))))])
                   elseif A(i,j)==0
                       % tao(i,r)=tao(i,r)*exp(tao(j,rr)*(psi(mu(r,rr,2))-psi(sum(mu(r,rr,:)))));
                       tao(i,r)=tao(i,r)+tao(j,rr)*(psi(mu(r,rr,2))-psi(sum(mu(r,rr,:))));
                   else
                       % tao(i,r)=tao(i,r)*exp(tao(j,rr)*(psi(mu(r,rr,3))-psi(sum(mu(r,rr,:)))));
                       tao(i,r)=tao(i,r)+tao(j,rr)*(psi(mu(r,rr,3))-psi(sum(mu(r,rr,:))));
                       % disp(['tao-1:',num2str(tao(i,r))])
                   end
                end
            end
        end
    end
    % tao(i,:)
    if max(tao(i,:))-min(tao(i,:))>700        
        [mv,mi]=max(tao(i,:));
        tao(i,:)=0;
        tao(i,mi)=1;
%     else
%         tao(i,:)=tao(i,:)-min(tao(i,:))+1;
    end
    tao(i,:)=exp(tao(i,:));
    tao(i,:)=tao(i,:)/sum(tao(i,:));
end
%% compute the approximate distribution of omega  
   rho=rho_0+sum(tao);

%% compute the approximate distribution of pi
for r=1:k
    for rr=1:k
        if r~=rr
          mu(r,rr,1)=mu_0(r,rr,1);
          mu(r,rr,2)=mu_0(r,rr,2);
          mu(r,rr,3)=mu_0(r,rr,3);
          for i=1:n
             for j=1:n
                 if i~=j
                    mu(r,rr,1)=mu(r,rr,1)+tao(i,r)*tao(j,rr)*(A(i,j)==1);
                    mu(r,rr,2)=mu(r,rr,2)+tao(i,r)*tao(j,rr)*(A(i,j)==0);
                    mu(r,rr,3)=mu(r,rr,3)+tao(i,r)*tao(j,rr)*(A(i,j)==-1);
                 end
             end
          end
        else
            mu(r,rr,1)=mu_0(r,rr,1);
            mu(r,rr,2)=mu_0(r,rr,2);
            mu(r,rr,3)=mu_0(r,rr,3);
            for i=1:n
               for j=i+1:n
                    mu(r,rr,1)=mu(r,rr,1)+tao(i,r)*tao(j,rr)*(A(i,j)==1);
                    mu(r,rr,2)=mu(r,rr,2)+tao(i,r)*tao(j,rr)*(A(i,j)==0);
                    mu(r,rr,3)=mu(r,rr,3)+tao(i,r)*tao(j,rr)*(A(i,j)==-1);
               end
            end
        end
    end
end 

%% compute evidence
   evidence1=gammaln(sum(rho_0))-gammaln(sum(rho));
   for r=1:k
       evidence1=evidence1+gammaln(rho(r))-gammaln(rho_0(r));
   end
   evidence2=0;
   for r=1:k
       for rr=r:k
           tmp=0;
           tmp=gammaln(sum(mu_0(r,rr,:)))-gammaln(sum(mu(r,rr,:)));
           tmp=tmp+gammaln(mu(r,rr,1))+gammaln(mu(r,rr,2))+gammaln(mu(r,rr,3))-(gammaln(mu_0(r,rr,1))+gammaln(mu_0(r,rr,2))+gammaln(mu_0(r,rr,3)));
           evidence2=evidence2+tmp;
       end
   end
   evidence3=0;
   for i=1:n
       for r=1:k
           evidence3=evidence3+tao(i,r)*log(tao(i,r));
       end
   end
   evidence=evidence1+evidence2-evidence3;
   %% if convergence based on evidence   
%    li_vary=abs(abs(evidence)-abs(evidencet))
%    if li_vary<10^-1
%        torf=0;
%    end
%    evidencet=evidence;
   %% if convergence based on tao
   tao_vary=sum(sum(abs(tao-taot)))
   if tao_vary<10^1
       torf=0;
   end
   taot=tao;
end


function Z = comp_Z(Gamma)
Z = Gamma;
[L,K]= size(Gamma);
for i=1:L
    [m,l] = max(Z(i,:));
    Z(i,:)=zeros(1,K);
    Z(i,l)= 1;   
end
   
function [B,C_number,S_C,index] = Draw(A,Z,flag)
[n,K] = size(Z);
index = [];
S_C = [];
% for clustering a bipartite network
index_col = ones(1,n);
index_col = find(index_col);
for k=1:K
    cluster = find(Z(:,k)==1);
    s_c = length(cluster);
    index = [index; cluster];
    S_C = [S_C s_c];
end
B = A(index,index);
C_number = length(find(S_C>0));
% if flag==1
%  disp(['The number of cluster is:' num2str(C_number)]);
% end
 disp(['The number of cluster is:' num2str(C_number)] );
% S_C
% figure;
if flag==1
    imagesc(B);
else
    spy(B);
end
d=0;
for k = 1:K-1
    d = d+S_C(k);
    if S_C(k)==0
        continue;
    end
    line([0.5,n+0.5],[d+0.5,d+0.5],'Color','r');
    line([d+0.5,d+0.5],[0.5,n+0.5],'Color','r');
end

function draw_random_net(A,S_C)
    K = length(S_C);
    n = length(A);
    d=0;
    spy(A)
    for k = 1:K-1
      d = d+S_C(k);
      if S_C(k)==0
        continue;
      end
      line([0.5,n+0.5],[d+0.5,d+0.5],'Color','r');
      line([d+0.5,d+0.5],[0.5,n+0.5],'Color','r');
    end
       