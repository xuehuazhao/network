function [Z]=SISN_main(A,k_min,k_max)
% input: A is the Adjacency matrix of signed network;
%        k_min is the minimum of the number of communities;
%        k_max is the maximum of the number of communities;
% output:Z is the indicating vector of nodes belonging to communities;

% Written by Xuehua Zhao in 2017. 
% Please cite Xuehua Zhao, Bo Yang, Xueyan Liu, and Huiling Chen. Statistical 
% inference for community detection in signed networks. Physical Review E, 2016, 95(4):042313. 

close all;clc;clear all;
%% initialization
likelihoodf=[];
costf=10^10;
likelihoodset=[];
costset=[];
pzf=[];
pif=[];
thetaf=[];
k_t=0;
for k=k_min:k_max
    k_t=k_t+1;    
    times=60;
    [pz,pi,theta,likelihood,cost]=msm_times(A,k,times);
    cost=cost_one(A,k,pi,theta);
    cost=-likelihood+cost;
    likelihoodset(k_t)=likelihood
    costset(k_t)=cost
    if costf>cost
        pzf=pz; pif=pi;thetaf=theta;likelihoodf=likelihood;
        costf=cost;
    end
end
Z = comp_Z(pzf);
[TA,C_number,S_C,index] =Draw(A,Z,1);

%% calculating NMI
    % generating the class label of nodes
    Z=full(Z); z=zeros(length(Z),1);
    for i=1:size(Z,2)
        for j=1:length(Z)
            if Z(j,i)==1 z(j)=i; end
        end
    end
    % calculating NMI
    label=load('D:\mprogram\mix signed model\sam_c.txt');
    acc=compute_NMI(z,label);
    disp(['NMI = ',num2str(acc)]);
%---
function [pz,pi,theta,likelihood,cost]=msm_times(A,k,times)
    likelihood=-10^10;
    pz=[];pi=[];theta=[];cost=[];
    for it=1:times
        [pzt,pit,thetat,likelihoodt,costt]=msm(A,k);
        if likelihood<likelihoodt
            likelihood=likelihoodt;
            pz=pzt;pi=pit;theta=thetat;likelihood=likelihoodt;cost=costt;
        end            
    end

function [pz,pi,theta,likelihood,cost]=msm(A,k)
%% initialization
   n=length(A);
   pz=[];pi=[];theta=[];likehihood=[];cost=[];
   pi=rand(1,k); pi=pi/sum(pi);
   likelihood=0;
   torf=1;   
 
   theta = rand(k,n,3); 
   for i=1:k
       for j=1:n
           for h=1:3
               theta(i,j,h)=rand();
           end
           theta(i,j,:)=theta(i,j,:)/(sum(theta(i,j,:)));
       end
   end
%%-init. end
while torf
%% E step£ºcalculating the posterior of Z
for i=1:n
    for r=1:k
        pz1=1;
        for j=1:n
            if A(i,j)==1 pz1=pz1*theta(r,j,1); end
            if A(i,j)==-1 pz1=pz1*theta(r,j,2); end
            if A(i,j)==0  pz1=pz1*theta(r,j,3); end
            
        end
        pz(i,r)=pz1*pi(k);
    end
    pz(i,:)=pz(i,:)/sum(pz(i,:));
end
%% M step: calculating parameters
pi=sum(pz)/n;
for r=1:k
    for j=1:n
        theta(r,j,:)=0;
        for i=1:n
            if A(i,j)==1 theta(r,j,1)=theta(r,j,1)+pz(i,r); end
            if A(i,j)==-1 theta(r,j,2)=theta(r,j,2)+pz(i,r); end
            if A(i,j)==0  theta(r,j,3)=theta(r,j,3)+pz(i,r); end        
        end
        theta(r,j,:)=theta(r,j,:)/sum(pz(:,r));
    end
end

%% calculating likelihood
   likelihoodt=0;
   for i=1:n
       tem1=0;
       for r=1:k
           tem2=1;
           for j=1:n
               if A(i,j)==1 tem2=tem2*theta(r,j,1); end
               if A(i,j)==-1 tem2=tem2*theta(r,j,2); end
               if A(i,j)==0 tem2=tem2*theta(r,j,3); end
           end
           tem1=tem1+pi(r)*tem2;
       end
       likelihoodt=likelihoodt+log(tem1);
   end
   li_vary=abs(abs(likelihood)-abs(likelihoodt));
   likelihood=likelihoodt;
   if li_vary<10^-2
       torf=0;
   end
end
%% calculating cost function
function cost=cost_one(A,k,pi,theta)
    cost=0;
    alpha=1/(size(A,1)*0.8);
    for i=1:k
        if pi(i)>alpha && i<k
            cost=cost+log(pi(i)/alpha);
        end
        for j=1:size(A,1)
            if theta(i,j,1)>alpha
                cost=cost+log(theta(i,j,1)/alpha);
            end
            if theta(i,j,2)>alpha
                cost=cost+log(theta(i,j,2)/alpha);
             end
        end        
    end
function cost=cost_two(A,k,pi,theta)
    n=size(A,1);
    cost=0;   
%      for i=1:k
%            cost=cost+0.5*log(pi(i));
%         for j=1:n
%             if theta(i,j,1)>0.001
%             cost=cost+0.5*log(theta(i,j,1));
%             end
%             if theta(i,j,2)>0.001
%             cost=cost+0.5*log(theta(i,j,2));
%             end
%             if theta(i,j,3)>0.001
%             cost=cost+0.5*log(theta(i,j,3));
%             end
%          end
%      end
     cost=cost+0.5*k*(2*k+1)*log(n)+0.5*k*(k+1)*(1-log(12));

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
figure;
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
       