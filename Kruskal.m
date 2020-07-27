function [T,result]=Kruskal(varargin)
%构造最小生成树
%算法：Kruskal
D=varargin{1};
%%做连接图（邻接矩阵）
a = zeros(size(D,1));
% for i = 1:size(D,1)
%     for j = 1:size(D,1)
%         if i<j
%             a(i,j)=sqrt((D(i,1)-D(j,1))*(D(i,1)-D(j,1))+(D(i,2)-D(j,2))*(D(i,2)-D(j,2)));
%         end
%     end
% end
a=pdist2(D,D);
G = graph(a,'upper');  
[T,pred] = minspantree(G,'Method','sparse');
result = [ T.Edges.EndNodes(:,1), T.Edges.EndNodes(:,2), T.Edges.Weight]';
end%构造最小生成树