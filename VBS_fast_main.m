function VBS_fast_main()
close all;clc;clear all;
nmiset=zeros(11,1);
times=1;
for p=0:0 %10
    for i=1:times
%% Data input
% [A,Q] = GenerateSignAYang(4,32,32,0.5,p*0.05,0.5);
A=load('4_10000.dat');
A=sparse(A(:,1),A(:,2),A(:,3),10000,10000);
%     [A]=generate_RN_signed_large()

%      A=load('war_max.mat');
%      A=A.A;
%      imagesc(A)
% A=load('symwiki_nodegree_0.mat');
% A=A.c;
% A=load('Slashdot081106_no0.mat');
% A=A.A;
% spy(A)
% A=load('Gahuku-Gama Subtribes of Highland New Guinea_1.txt');
    k_min=2; k_max=10;
    evidenceset=[];s_tao=[];s_rho=[];s_mu=[];s_evidence=[];
    k_t=0;
tic
for k=k_min:k_max
    k_t=k_t+1
    timess=1;
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
Z = comp_Z(s_tao);
size(Z,2)
Z(:,find(sum(Z)==0))=[];
size(Z,2)

%---
    nmiset(p+1)=nmiset(p+1)+acc;
    nmiset'
    [TA,C_number,S_C,index] =Draw(A,Z,0);
    end    
end
evidenceset
nmiset=nmiset/times

function [t_tao,t_rho,t_mu,t_evidence]=ssbm_times(A,k,times)
    t_evidence=[];
    t_tao=[];t_rho=[];t_mu=[];
    for it=1:times
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
   rho_0=ones(1,k)*0.1;% non-informative Jeffreys prior distribution
   mu_0=ones(k,k,3)*1/2; % non-informative Jeffreys prior distribution
   evidence=0;
   torf=1;   
   rho=ones(1,k)*0.5;
   mu=ones(k,k,3)*0.5; %%1-positive link 2-no link 3-nagetive link

   % generate edge matix
   [re,ce]=find(A==1); % positive edge
   edgep=[re,ce];
   [re,ce]=find(A==-1); %negative edge
   edgen=[re,ce];
   for i=1:k
       for j=1:k
           if i==j               
               mu(i,j,1)=0.4;
               mu(i,j,2)=0.5;%0.5
               mu(i,j,3)=0.1;%0.1
           else
               mu(i,j,1)=0.1;%0.1;
               mu(i,j,2)=0.5;%0.5
               mu(i,j,3)=0.2;%
           end
       end
   end 
    tao=rand(n,k);
    tao=tao./(sum(tao,2)*ones(1,k));
    taot=tao;
    tao_sum=sum(tao);
    nadp=cell(1,n);
    nadn=cell(1,n);
    for i=1:n
        nadp{i}=find(A(i,:)==1);
        nadn{i}=find(A(i,:)==-1);
    end
%%-initializing parameters end
tt=0; % set end condition
evidencet=0;
while torf
   tt=tt+1
   %% compute the approximate distribution of Z
   psimu=psi(mu);
   psirho=0;   
   for i=1:k
       for j=1:k
          psimu(i,j,4)=psi(sum(mu(i,j,:)));
       end
       psirho(i)=psi(rho(i))-psi(sum(rho));
   end
   
   for i=1:n
     tao_t=tao(i,:);
     tp=sum(tao(nadp{i},:),1);
     tn=sum(tao(nadn{i},:),1);
     tz=tao_sum-tp-tn-tao(i,:);

     for r=1:k
        tao(i,r)=psirho(r);
        for rr=1:k
            tao(i,r)=tao(i,r)+tp(rr)*(psimu(r,rr,1)-psimu(r,rr,4))...
                             +tn(rr)*(psimu(r,rr,3)-psimu(r,rr,4))...
                             +tz(rr)*(psimu(r,rr,2)-psimu(r,rr,4));
        end
     end  % disp(['tao(i,:),i=',num2str(i),':',num2str(tao(i,:))])
     mi= min(tao(i,:));
     if max(tao(i,:))-mi>700
        [mv,mind]=max(tao(i,:));
        tao(i,:)=0;
        tao(i,mind)=1;
     else if mi<-700
          tao(i,:)=tao(i,:)-(mi+700);
         end         
     end
     tao(i,:)=exp(tao(i,:))/sum(exp(tao(i,:)));  % disp(['exp tao(i,:),i=',num2str(i),':',num2str(tao(i,:))])
     %% updating the related parameters
     tao_dif=tao(i,:)-tao_t;
     tao_sum=tao_sum+tao_dif;
   end

%% compute the approximate distribution of omega  
   rho=rho_0+sum(tao);
   rho(find(rho<10^-300))=10^-300; 
 
%% compute the approximate distribution of pi
for r=1:k
    for rr=r:k   
        mu(r,rr,1)=0;
        mu(r,rr,2)=0;
        mu(r,rr,3)=0;        
        for i=1:size(edgep,1)
            mu(r,rr,1)=mu(r,rr,1)+tao(edgep(i,1),r)*tao(edgep(i,2),rr);
        end        
        for i=1:size(edgen,1)
             mu(r,rr,3)=mu(r,rr,3)+tao(edgen(i,1),r)*tao(edgen(i,2),rr);
        end
        if r~=rr
            mu(r,rr,1)=mu(r,rr,1)+mu_0(r,rr,1);
            mu(r,rr,3)=mu(r,rr,3)+mu_0(r,rr,3);
            mu(r,rr,2)=mu_0(r,rr,2)+sum(tao(:,r))*sum(tao(:,rr))-sum(tao(:,r).*tao(:,rr))-mu(r,rr,1)-mu(r,rr,3)+mu_0(r,rr,1)+mu_0(r,rr,3);
            mu(r,rr,find(mu(r,rr,:)<10^-300))=10^-300; 
            mu(rr,r,1)=mu(r,rr,1);
            mu(rr,r,2)=mu(r,rr,2);
            mu(rr,r,3)=mu(r,rr,3);
        else
            mu(r,rr,1)=mu(r,rr,1)/2+mu_0(r,rr,1);
            mu(r,rr,3)=mu(r,rr,3)/2+mu_0(r,rr,3);
            mu(r,rr,2)=mu_0(r,rr,2)+(sum(tao(:,r))*sum(tao(:,rr))-sum(tao(:,r).*tao(:,rr)))/2-mu(r,rr,1)-mu(r,rr,3)+mu_0(r,rr,1)+mu_0(r,rr,3);
            mu(r,rr,find(mu(r,rr,:)<10^-300))=10^-300;
        end        
    end    
end

%% compute evidence
   evidence1=gammaln(sum(rho_0))-gammaln(sum(rho))+sum(gammaln(rho)-gammaln(rho_0));
%    for r=1:k
%        evidence1=evidence1+gammaln(rho(r))-gammaln(rho_0(r));
%    end  
   evidence2=0;
   for r=1:k
       for rr=r:k
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
%    evidence3=sum(sum(tao.*log(tao)));
   evidence=evidence1+evidence2-evidence3;
   %% if convergence based on evidence   
   li_vary=abs(abs(evidence)-abs(evidencet));
   if li_vary<10^-5
       torf=0;
   end
   evidencet=evidence;

   %% if convergence based on the times of iteration
   if tt>8
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
       
