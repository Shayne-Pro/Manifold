function [dist,neighbor]=top_K_neighbors(x,x_train,y_train,K)

% Input:
%      x: 待预测点集(每一列表示一个点)
%      x_train: 训练点集
%      y_train: 训练点集标签
%      K: 选取的K个最邻近的点
% Output:
%      dist: 最近的K个距离（由小到大排序） 
%      neighbor: 这K个训练点集的标签

[~,size_x]=size(x);
[~,size_xtrain]=size(x_train);
for i=1:size_x
    for j=1:size_xtrain
        temp(j,i)=norm(x(:,i)-x_train(:,j));
    end
end

for i=1:size_x
    tmp1=[temp(:,i) y_train'];  %对每一列加入标签并进行排序
    tmp2=sortrows(tmp1,1);
    if size(tmp2,1)>=K
        tmp3=tmp2(1:K,:);  %取其中K个最小的距离
        dist(:,i)=tmp3(:,1);
        neighbor(:,i)=tmp3(:,2);
    else
        tmp3=tmp2;
        dist(:,i)=tmp3(:,1);
        neighbor(:,i)=tmp3(:,2);
    end
end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        