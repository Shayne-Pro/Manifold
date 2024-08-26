function [D_E,D_G,H,S]=nonlinearity_score(X,K)
% 计算非线性得分函数S
%
%  Input:
%       X: 数据集X=[x1,x2,...,xN],每个xi都是一个样本，尺寸是D*N
%       K: K个近邻点
%  Output:
%       D_E: 欧式距离矩阵
%       H: K邻近的标签
%       S: 非线性得分函数

%%%%%%%%%%%%%%% 对输入的X去掉相同的列 %%%%%%%%%%%%%%%
% X_tp=unique(X','rows');
% X=X_tp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,X_size]=size(X);
small_distance=2*10^(-6);

%%% 计算欧式距离矩阵D_E，维数是N*N %%%
for i=1:X_size
    for j=1:X_size
        D_E(i,j)=norm(X(:,i)-X(:,j));
    end
end
for i=1:X_size
    for j=1:X_size
        if D_E(i,j)==0
           D_E(i,j)=small_distance;
        end
    end
end

%%% 计算测地距离矩阵D_G，维数是N*N %%%
n_size_step=2;
% INF=1000*max(max(D_E))*size(D_E,1); 
%  增加K值以保证单连通图
while (1)    
   [D_G,~]= Geodisic(X',K);
    if (1)
        aa=(D_G==inf);
        [tmp, firsts] = min(D_G==inf);     %% first point each point connects to
        [comps, I, J] = unique(firsts);    %% first point in each connected component
        n_comps = length(comps);           %% number of connected components
        % size_comps = sum((repmat(firsts,n_comps,1)==(comps'*ones(1,X_size)))');                                        
    end
%     disp([' --> fun_GetHDCMLPs Step 1: n_size: ' num2str(n_size) ' and Number of comps: ' num2str(n_comps)]);
    if ( n_comps == 1 )
        break;
    else
        K = K + n_size_step;
    end
end
for i=1:X_size
    for j=1:X_size
        if D_G(i,j)==0
           D_G(i,j)=small_distance;
        end
    end
end
%%% 计算比率矩阵R %%%
R=D_G./D_E;

%%% 计算K邻近矩阵H %%%
X_label=1:X_size;
[dist,neighbor]=top_K_neighbors(X,X,X_label,K+1);
dist(1,:)=[];
neighbor(1,:)=[];
H=neighbor;
 
%%% 计算非线性得分函数S %%%
sum_X=sum(sum(R));
S=sum_X/(X_size*X_size);
%S=mean(R(:));