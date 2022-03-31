% function [index,index1,index2]=Dheat(varargin)
function [index]=Dheat(varargin)%compute density deviation and find density peaks
D=varargin{1};

[Sup,NN,RNN,NNN,nb,A]=NaNSearching(D);
% index1=find(nb==0);
[dev]=findCenter2(D,NN,NNN,A,nb);
% index2=find(dev<(mean(dev)-2.5*std(dev)));%噪声点
index=find(dev>(mean(dev)-(0.5)*std(dev)));%核心点
end



function[Sup,NN,RNN,NNN,nb,A]=NaNSearching(varargin)
D=varargin{1};
r=1;
nb=zeros(size(D,1),1);
C=cell(size(D,1),1);
NN=cell(size(D,1),1);%初始化每个点的KNN邻居
RNN=cell(size(D,1),1);%初始化每个点的RKNN邻居
NNN=cell(size(D,1),1);%是NN和RNN的交集，也就每个点的
A=pdist2(D,D);
Numb1=0;
Numb2=0;
for ii=1:size(D,1)
   [sa,index]=sort(A(:,ii));
   C{ii}=[sa,index];
end
while(r<size(D,1))
 for kk=1:size(D,1)
     x=kk;
     y=C{x}(r+1,2);
     nb(y)= nb(y)+1;
     NN{x}=[NN{x},y];
     RNN{y}=[RNN{y},x];
 end
    Numb1=sum(nb==0);
    if Numb2~=Numb1
        Numb2=Numb1;
    else
       break; 
    end
    r=r+1;
end
for jj=1:size(D,1)
NNN{jj}=intersect(NN{jj},RNN{jj});
end
Sup=r;
end




function[dev]=findCenter2(varargin)
D=varargin{1};
NN=varargin{2};
NNN=varargin{3};
A=varargin{4};
nb=varargin{5};


r1=[];
Pr1=[];
dev=[];


CP=[];


Nei1=cell(size(D,1),1);                                              

 
for kk=1:size(D,1)
 CL(kk)=0;
 end 
    for ii=1:size(D,1)
       if ~isempty(NN{ii})
         r1(ii)=1*max(pdist2(D(ii,:),D(NN{ii},:)));
		 Nei1{ii}=find(A(:,ii)<r1(ii));         
         Pr1(ii)=size(Nei1{ii},1);
         dev(ii)=log(Pr1(ii))-(size(D,2))*log(r1(ii)+1);
		else
           r1(ii)=0;
           dev(ii)=0;
       end
    end

end
