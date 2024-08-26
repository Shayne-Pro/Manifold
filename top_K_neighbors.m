function [dist,neighbor]=top_K_neighbors(x,x_train,y_train,K)

% Input:
%      x: ��Ԥ��㼯(ÿһ�б�ʾһ����)
%      x_train: ѵ���㼯
%      y_train: ѵ���㼯��ǩ
%      K: ѡȡ��K�����ڽ��ĵ�
% Output:
%      dist: �����K�����루��С�������� 
%      neighbor: ��K��ѵ���㼯�ı�ǩ

[~,size_x]=size(x);
[~,size_xtrain]=size(x_train);
for i=1:size_x
    for j=1:size_xtrain
        temp(j,i)=norm(x(:,i)-x_train(:,j));
    end
end

for i=1:size_x
    tmp1=[temp(:,i) y_train'];  %��ÿһ�м����ǩ����������
    tmp2=sortrows(tmp1,1);
    if size(tmp2,1)>=K
        tmp3=tmp2(1:K,:);  %ȡ����K����С�ľ���
        dist(:,i)=tmp3(:,1);
        neighbor(:,i)=tmp3(:,2);
    else
        tmp3=tmp2;
        dist(:,i)=tmp3(:,1);
        neighbor(:,i)=tmp3(:,2);
    end
end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        