function [D,D_E] = Geodisic(Y,K) 
% function [D] = Geodisic(Y, K, dim) 
% 计算测地距离

%  Input:
%       Y: 每一行代表一个点，如Y=[x1;x2;...;x_N]
%       K: 选取的K个近邻点
%  Output:
%       D: 形成的测地距离矩阵
%       D_E: 形成的K个欧氏距离
%       yy： 这K个最小欧氏距离的标签

%%%%% Step 0: 初始化 %%%%%
D=squareform(pdist(Y));  % Y每一行表示一个点，生成每两个点之间的欧氏距离矩阵
%% 计算K个欧氏距离 
% S=Y';
% xx=S(:,1);
% x_train=S(:,2:end);
% y_train=2:size(S,2);
% yy=knn(xx,x_train,y_train,K);
D_E=D;

%% 计算K个测地距离
N = size(D,1); 
INF =  1000*max(max(D))*N;  % 使无穷大有意义
% dims = dim;
% comp = 1;  % Assume one component.
%Y.coords = cell(length(dims),1); 
% R = zeros(1,length(dims)); 
%%%%% Step 1: 构建邻近图 %%%%%
[tmp, ind] = sort(D);  % tmp是排好序的矩阵，ind是tmp在原矩阵D下的坐标

% 将超过K个较大的距离值设为无穷大
for i=1:N
    D(i,ind((2+K):end,i)) = inf; 
end; 

D = min(D,D');    % 确保距离矩阵对称
%%%%% Step 2: 计算最短距离，即测地距离 %%%%%
%% 此步相当于Dijkstra算法
for k=1:N
     D = min(D,repmat(D(:,k),[1 N])+repmat(D(k,:),[N 1])); 
end



%%%%% Remove outliers from graph %%%%%
% n_connect = sum(~(D==INF));        %% number of points each point connects to
% [tmp, firsts] = min(D==INF);       %% first point each point connects to
% [comps, I, J] = unique(firsts);    %% represent each connected component once
% size_comps = n_connect(comps);     %% size of each connected component
% [tmp, comp_order] = sort(size_comps);  %% sort connected components by size
% comps = comps(comp_order(end:-1:1));    
% size_comps = size_comps(comp_order(end:-1:1)); 
% n_comps = length(comps);               %% number of connected components
% if (comp>n_comps)                
%      comp=1;                              %% default: use largest component
% end
% Y.index = find(firsts==comps(comp)); 
% D = D(Y.index, Y.index); 