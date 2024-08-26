function [D,D_E] = Geodisic(Y,K) 
% function [D] = Geodisic(Y, K, dim) 
% �����ؾ���

%  Input:
%       Y: ÿһ�д���һ���㣬��Y=[x1;x2;...;x_N]
%       K: ѡȡ��K�����ڵ�
%  Output:
%       D: �γɵĲ�ؾ������
%       D_E: �γɵ�K��ŷ�Ͼ���
%       yy�� ��K����Сŷ�Ͼ���ı�ǩ

%%%%% Step 0: ��ʼ�� %%%%%
D=squareform(pdist(Y));  % Yÿһ�б�ʾһ���㣬����ÿ������֮���ŷ�Ͼ������
%% ����K��ŷ�Ͼ��� 
% S=Y';
% xx=S(:,1);
% x_train=S(:,2:end);
% y_train=2:size(S,2);
% yy=knn(xx,x_train,y_train,K);
D_E=D;

%% ����K����ؾ���
N = size(D,1); 
INF =  1000*max(max(D))*N;  % ʹ�����������
% dims = dim;
% comp = 1;  % Assume one component.
%Y.coords = cell(length(dims),1); 
% R = zeros(1,length(dims)); 
%%%%% Step 1: �����ڽ�ͼ %%%%%
[tmp, ind] = sort(D);  % tmp���ź���ľ���ind��tmp��ԭ����D�µ�����

% ������K���ϴ�ľ���ֵ��Ϊ�����
for i=1:N
    D(i,ind((2+K):end,i)) = inf; 
end; 

D = min(D,D');    % ȷ���������Գ�
%%%%% Step 2: ������̾��룬����ؾ��� %%%%%
%% �˲��൱��Dijkstra�㷨
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