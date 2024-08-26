function L = Compute_L(X,sigma1)
% 计算测试集的Laplace图 （一个测试集被划分成m个测试集根据稀疏子空间聚类）
% Input:
%     X: 大小d*p*m。测试集被划分成m个后求得的m个子空间 X=[X1,X2,...,Xm],其中Xi的大小为d*p,p为Grassmann流形G(p,d)
%     sigma1: Gassian核函数的系数
% Output:
%     L: 测试集内部的Laplace图

m = size(X,3);

K_X = grassmann_proj(X);

%% 先计算相似性矩阵W
W = zeros(size(K_X,1),size(K_X,2));
for i = 1:size(W,1)
    for j = 1:size(W,2)
        dis_i_j = K_X(i,i)+K_X(j,j)-2*K_X(i,j);
        W(i,j) = exp(-dis_i_j/sigma1);
    end
end
%% 计算D
D = zeros(size(W,1),size(W,2));
for i = 1:size(D,1)
    D(i,i) = sum(W(i,:));
end
%% 计算L
L = D-W;
















































