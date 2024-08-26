function B = Compute_B(X,A,Train_label)
% 计算测试集的Laplace图 （一个测试集被划分成m个测试集根据稀疏子空间聚类）
% Input:
%     X: 大小d*p*m。测试集被划分成m个后求得的m个子空间 X=[X1,X2,...,Xm],其中Xi的大小为d*p,p为Grassmann流形G(p,d)
%     A: 大小d*p*N。字典
%     C: 一个图像集被分成C簇
% Output:
%     B: m*N。测试集与字典之间的距离，B的每一行是测试集被分成m个后其中一个与N个字典的距离

% m = size(X,3);
% N = size(A,3);
% 
% K_XA = grassmann_proj(X,A);
% K_X = grassmann_proj(X);
% K_A = grassmann_proj(A);
% 
% B = zeros(size(K_XA,1),size(K_XA,2));
% for i = 1:size(B,1)
%     for j = 1:size(B,2)
%         dis_i_j = K_X(j,j)+K_A(i,i)-2*K_XA(i,j);
%         B(i,j) = exp(dis_i_j/sigma2);
%     end
% end


B = zeros(size(A,3),size(X,3));
Number_Of_Classes = max(Train_label);
for i = 1 : size(B,2)
    X_i = X(:,:,i);
    sum_all = 0;
    for j = 1 : size(A,3)
        sum_all = sum_all + sum(sum((X_i*X_i'-A(:,:,j)*A(:,:,j)').^2));
    end
    for j = 1 : Number_Of_Classes 
        classIndex = (Train_label==j);
        class_j = find(classIndex==1);
        sum_Xi = 0;
        for k = 1 : length(class_j)
            sum_Xi = sum_Xi + sum(sum((X_i*X_i'-A(:,:,class_j(k))*A(:,:,class_j(k))').^2));
        end
        B(class_j,i) = sum_Xi/sum_all;
    end
end

end