function [A,Q] = GenerateSignAYang(c,n,d,pin,pn,pp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%signed network generate model based the paper of Bo Yang S(c,n,d,pin,p-,p+)
%input£ºc the number of communities£»
%      n the number of nodes of each community£»
%      d the average degree of each node£»
%      pin the probability of positive links to communities£»
%      p- the probability of negative links in communities£»
%      p+ the probability of positive links between communities£»
%output£ºsigned matix A,label Q
% c=4;n=32;d=32;pin=0.8;pn=0;pp=0;
num=c*n;  % the number of nodes
incn=ceil(d*pin);
outcn=d-incn;

A=zeros(num,num);  
for i=1:num
    ran=randperm(num);
    for j=1:num
        idegreein=sum(A(i,:)+abs(A(i,:)))/2; % the positive degree of node i
        if idegreein<incn
            nodej=ran(j);  % the label of node j
            if i~=nodej
               if nodej<=(ceil(i/n)*n)&&nodej>(ceil(i/n)-1)*n
                  jdegreein=sum(A(nodej,:)+abs(A(nodej,:)))/2; %the positive degree of node j
                  if jdegreein<incn
                      A(i,nodej)=1; A(nodej,i)=1;
                  end
               end
            end
        end        
     end
     for j=1:num
        idegreeout=-sum(A(i,:)-abs(A(i,:)))/2; %the negative degree of node i
        if idegreeout<outcn
            nodej=ran(j);  % the label of node j
            if i~=nodej
            if nodej>(ceil(i/n)*n)||nodej<=(ceil(i/n)-1)*n
                jdegreeout=-sum(A(nodej,:)-abs(A(nodej,:)))/2; % the negative degree of node j
                if jdegreeout<outcn
                    A(i,nodej)=-1;
                    A(nodej,i)=-1;
                end
            end
            end
        end
    end
end
for i=1:num
    for j=i+1:num
        if A(i,j)~=0
           if j>ceil(i/n)*n           
              if rand()<=pp
                A(i,j)=-A(i,j);
                A(j,i)=-A(j,i);
              end
           else
              if rand()<=pn
                A(i,j)=-A(i,j);
                A(j,i)=-A(j,i);
              end
          end
        end
    end
end
Q=[];
n_num=n; % set the number of each block
for i=1:c
        Q((i-1)*n_num+1:i*n_num)=i;        
end  
imagesc(A);
