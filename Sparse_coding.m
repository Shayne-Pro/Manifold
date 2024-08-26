function Y = Sparse_coding(X,A,Y_Init,alpha,beta,lambda,miu,v,L,B)
% solve X所有列的稀疏系数 其中Y_Init,miu,v是初始化的参数，alpha,beta,lambda是正则化参数，X,A,L,B已知
% Input:
%     X: 一个测试集，通过稀疏子空间聚类，分成m个测试集，再通过流形对称投影变为m个样本   m*m
%     A: 字典 m*N
%     Y_Init: 初始化的稀疏系数矩阵 N*m
%     alpha: LLC二范数前的正则化参数
%     beta: Laplace图前的正则化参数
%     lambda: y的一范数前的正则化参数
%     miu: Lagrange乘子前的初始化系数初始化
%     v: Lagrange乘子中的参数vector初始化
%     L: Laplace图: 测试集分成m个子空间，这m个子空间在流形上的相互距离
%     B: 测试集到训练集的距离

[d,m] = size(X);
[d,N] = size(A);
Y = Y_Init;
max_Iter = 1e6;
tol = 1e-6;
fobj3_temp = 0;
% for i = 1:m
%     S_i = diag(B(:,i)'*diag(B(:,i)));
%     Y_i = Y(:,i);
%     fobj3_temp = fobj3_temp+Y_i'*S_i*Y_i;
% end
% fobj_old = [norm(X-A*Y,'fro')]^2 + beta*trace(Y*L*Y') + alpha*fobj3_temp + lambda*sum(sum(abs(Y)));

% 将稀疏系数矩阵的每一列迭代更新（更新第i列时固定其他列）
%start loop
% iter = 0;
% while iter<max_Iter
%     iter = iter + 1;
    for i = 1:m
        x = X(:,i);
        Y(:,i) = Ialm(i,x,A,Y,alpha,beta,lambda,miu,v,L,B);
    end
%     fobj1 = [norm(X-A*Y,'fro')]^2;
%     fobj2 = beta*trace(Y*L*Y');
%     fobj3_temp = 0;
%     for i = 1:m
%         S_i = diag(B(:,i)'*diag(B(:,i)));
%         Y_i = Y(:,i);
%         fobj3_temp = fobj3_temp+Y_i'*S_i*Y_i;
%     end
%     fobj3 = alpha*fobj3_temp;
%     fobj4 = lambda*sum(sum(abs(Y)));
%     fobj_new = fobj1+fobj2+fobj3+fobj4;
%     stopC = abs(fobj_new-fobj_old);
%     plot_residual(iter) = stopC;
%     if iter==1 || mod(iter,50)==0 || stopC<tol
%         disp(['iter ' num2str(iter) ',stopALM=' num2str(stopC,'%2.3e')]);
%     end
%     if stopC<tol
%         break;
%     end
%     fobj_old = fobj_new;
%     fobj_search(iter) = fobj_old;
%     clear fobj1 fobj2 fobj3 fobj4 fobj_new fobj_temp;
% end
% 
% %%%% 目标函数值的变化是否收敛 %%%%
%  plot(fobj_search);
% 
% end



































