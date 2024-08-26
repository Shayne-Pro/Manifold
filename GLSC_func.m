function [outLabel] = GLSC_func(X,A,Train_label,lambda1,lambda2)
% Input:
%     X: ��Сd*p*m�����Լ������ֳ�m������õ�m���ӿռ� X=[X1,X2,...,Xm],����Xi�Ĵ�СΪd*p,pΪGrassmann����G(p,d)
%     A: ��Сd*p*N���ֵ�
%     Y_Init: ��ʼ����ϡ��ϵ������ N*m
%     alpha: LLC������ǰ�����򻯲���
%     beta: Laplaceͼǰ�����򻯲���
%     lambda: y��һ����ǰ�����򻯲���
%     miu: Lagrange����ǰ�ĳ�ʼ��ϵ����ʼ��
%     v: Lagrange�����еĲ���vector��ʼ��
%     sigma1: Gassian�˺�����ϵ��
%     sigma2: Gassian�˺�����ϵ��
% Output:
%     Y: ϡ��ϵ������
%     qX: �����ϵĲ�������
%     qA: �����ϵ��ֵ�

[d,p,m] = size(X);
[d,p,N] = size(A);

%% �ȼ������������Ծ���L,B
% L = Compute_L(X,sigma1);
% B = Compute_B(X,A,Train_label);

%% �����
K_XA = grassmann_proj(X,A);
K_A = grassmann_proj(A);
K_X = grassmann_proj(X);

%% ��ŷʽ�ռ仪Ϊ����
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
% %% ϡ���ʾ
% Y = Sparse_coding(qX,qA,Y_Init,alpha,beta,lambda,miu,v,L,B);

outLabel=Classify_ISCRC(qVx,M*qA,Train_label,bx,ax);

























