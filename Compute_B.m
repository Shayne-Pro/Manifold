function B = Compute_B(X,A,Train_label)
% ������Լ���Laplaceͼ ��һ�����Լ������ֳ�m�����Լ�����ϡ���ӿռ���ࣩ
% Input:
%     X: ��Сd*p*m�����Լ������ֳ�m������õ�m���ӿռ� X=[X1,X2,...,Xm],����Xi�Ĵ�СΪd*p,pΪGrassmann����G(p,d)
%     A: ��Сd*p*N���ֵ�
%     C: һ��ͼ�񼯱��ֳ�C��
% Output:
%     B: m*N�����Լ����ֵ�֮��ľ��룬B��ÿһ���ǲ��Լ����ֳ�m��������һ����N���ֵ�ľ���

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