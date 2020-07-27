function [CL,CL1,d,index,mst]=main(varargin)% IMST-LDP
D=varargin{1};
CL=zeros(size(D,1),1);
% [index,noiseindex1,noiseindex2]=Dheat(D);
% [index,CL]=Dheat(D);
[index]=Dheat(D);
D2 = D(index,:);
C=cell(size(D2,1),1);
A=pdist2(D2,D2);
mindis=zeros(size(D2,1),1);
for ii=1:size(D2,1)
   [sa,index1]=sort(A(:,ii));
   C{ii}=[sa,index1];
   mindis(ii,1) = sa(2,:);
end
maxcut=max(mindis(:,1));
[G,mst] = Kruskal(D2);
maxedge=find(mst(3,:)>maxcut);
cutnumber=size(maxedge,2);
edge=mst(3,:);
[Y,I]=sort(edge,2,'descend');
sortmst=mst(1:2,I);
sortmst(3,:)=Y;
cutmst=sortmst;
G=graph(cutmst(1,:),cutmst(2,:));
for ii=1:cutnumber
    a=cutmst(1,ii);
    b=cutmst(2,ii);
    G=rmedge(G,a,b);
end
bins=conncomp(G);
cl=bins';
CL(index,:)=cl;
a=max(CL);
CS=zeros(a,1);%the number of points in each cluster
for i=1:a
    b1=find(CL==i);
    CS(i,1)=size(b1,1);
end
[CS_1,I1]=sort(CS,'descend');
score=zeros(a,1);
for i=2:a
    score(i)=1-exp(-CS_1(i-1)/CS_1(i)/CS_1(i));
end
count=0;
test=[];
for j=1:a
    if score(j)>0.1
         b2=find(CL==I1(j));
         CL(b2)=0;
         test=[test,I1(j)];
         count=count+1;
    end
end
for ii=1:size(test,2)
    b3=find(CL>test(1,ii));
    CL(b3)=CL(b3)-1;
end
testCL=CL;
d=a-count;
CL1=CL;
D4=find(CL==0);
D5=find(CL~=0);
%allocate remaining points;
distance=pdist2(D(D4,:),D(D5,:));
for i=1:size(D4,1)
    tempdistance=distance(i,:);
    [di,la]=min(tempdistance);
    templa=D5(la);
    CL(D4(i))=CL(templa);
end
d=max(CL);
end