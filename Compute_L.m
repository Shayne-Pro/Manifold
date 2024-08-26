function L = Compute_L(X,sigma1)
% ������Լ���Laplaceͼ ��һ�����Լ������ֳ�m�����Լ�����ϡ���ӿռ���ࣩ
% Input:
%     X: ��Сd*p*m�����Լ������ֳ�m������õ�m���ӿռ� X=[X1,X2,...,Xm],����Xi�Ĵ�СΪd*p,pΪGrassmann����G(p,d)
%     sigma1: Gassian�˺�����ϵ��
% Output:
%     L: ���Լ��ڲ���Laplaceͼ

m = size(X,3);

K_X = grassmann_proj(X);

%% �ȼ��������Ծ���W
W = zeros(size(K_X,1),size(K_X,2));
for i = 1:size(W,1)
    for j = 1:size(W,2)
        dis_i_j = K_X(i,i)+K_X(j,j)-2*K_X(i,j);
        W(i,j) = exp(-dis_i_j/sigma1);
    end
end
%% ����D
D = zeros(size(W,1),size(W,2));
for i = 1:size(D,1)
    D(i,i) = sum(W(i,:));
end
%% ����L
L = D-W;
















































