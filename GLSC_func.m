function [outLabel] = GLSC_func(X,A,Train_label,lambda1,lambda2)
% Input:
%     X: 大小d*p*m。测试集被划分成m个后求得的m个子空间 X=[X1,X2,...,Xm],其中Xi的大小为d*p,p为Grassmann流形G(p,d)
%     A: 大小d*p*N。字典
%     Y_Init: 初始化的稀疏系数矩阵 N*m
%     alpha: LLC二范数前的正则化参数
%     beta: Laplace图前的正则化参数
%     lambda: y的一范数前的正则化参数
%     miu: Lagrange乘子前的初始化系数初始化
%     v: Lagrange乘子中的参数vector初始化
%     sigma1: Gassian核函数的系数
%     sigma2: Gassian核函数的系数
% Output:
%     Y: 稀疏系数矩阵
%     qX: 流形上的测试样本
%     qA: 流形上的字典

[d,p,m] = size(X);
[d,p,N] = size(A);

%% 先计算两个相似性矩阵L,B
% L = Compute_L(X,sigma1);
% B = Compute_B(X,A,Train_label);

%% 计算核
K_XA = grassmann_proj(X,A);
K_A = grassmann_proj(A);
K_X = grassmann_proj(X);

%% 将欧式空间华为流形
[KX_U,KX_A,~] = svd(K_X);
qVx = diag(sqrt(diag(KX_A)))*KX_U';

[KA_U,KA_A,~] = svd(K_A);
qA = diag(sqrt(diag(KA_A)))*KA_U';

M_TEM = inv(qA'+0.001*eye(size(qA',1)))*K_XA*inv(qVx+0.001*eye(size(qVx,1)));
M = M_TEM';
QQ=qA'*(eye(size(M,2))-M'*M)*qA;
AA = [qVx -1*M*qA];
BBB =zeros(m+N);
BBD =zeros(m+N);
I11 = eye(m);
I22 = eye(N);
BBB(1:m,1:m)=lambda1*I11;
BBB(m+1:m+N,m+1:m+N)=lambda2*I22;
BBD =zeros(m+N);
BBD(m+1:m+N,m+1:m+N)=QQ;
BB=BBB+BBD;

dd= [ones(1,m) zeros(1,N)];

    cc=AA'*AA+BB;
	Z0=inv(cc+0.001*eye(size(cc,1)))*dd';
	zx=Z0 / (dd*Z0);
    ax=zx(1:m);
    bx=zx(m+1:m+N);
    sprintf('here is l2 norm to excute...');

% D_Inv = KA_U*diag(1./sqrt(diag(KA_A)));
% qX = D_Inv'*K_XA;
% 
% %% 稀疏表示
% Y = Sparse_coding(qX,qA,Y_Init,alpha,beta,lambda,miu,v,L,B);

outLabel=Classify_ISCRC(qVx,M*qA,Train_label,bx,ax);

























