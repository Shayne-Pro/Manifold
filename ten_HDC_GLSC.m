function [mean_accuracy,var_accuracy]=ten_HDC_GLSC(sub_dim,C,K)
% 计算GLSC在Extended yaleB上的十次平均值
% Input:
%     sub_dim: Grassmann流形上的子空间维度，相当于G(p,d)中的p
%     alpha: LLC二范数前的正则化参数
%     beta: Laplace图前的正则化参数
%     lambda: y的一范数前的正则化参数
%     miu: Lagrange乘子前的初始化系数初始化
%     v: Lagrange乘子中的参数vector初始化
%     sigma1: Gassian核函数的系数
%     sigma2: Gassian核函数的系数
% Output:
%     accuracy: eth80十次交叉实验的准确率
%     mean_accuracy: 十次交叉实验的平均准确率

%% 初始化训练集，测试集
load random_index_extyaleB;
load Extended_yaleB;
num_class=28;
num_sets_per_class=3;
accuracy_rate=zeros(1,10);
for num=1:1
    %% 初始化训练集，测试集，及它们的标签
    X_Train=cell(84,1);
    X_Test=cell(168,1);
    Train_label=zeros(1,84);
    Test_label=zeros(1,168);
    for j=1:num_class
        for k=1:num_sets_per_class
            temp1=Extended_yaleB{j,random_index{num}(j,k)};
            X_Train{(j-1)*num_sets_per_class+k,1}=temp1;
            Train_label(1,(j-1)*num_sets_per_class+k)=j;
        end
    end
    for j=1:num_class
        for k=1:6
            temp2=Extended_yaleB{j,random_index{num}(j,k+3)};
            X_Test{(j-1)*6+k,1}=temp2;
            Test_label(1,(j-1)*6+k)=j;
        end
    end
    
   %% 将所有训练集进行子空间聚类
    for i = 1 : size(X_Train,1)
        trn_i = X_Train{i,1};
        if size(trn_i,2) >= sub_dim*C
            trn_cluster(:,i) = HDC_local_model(trn_i,K,C);
        else
            trn_cluster{1,i} = trn_i;
          %  trn_cluster{2:C,i} = [];
        end
%%% 可视化分簇
% for ii=1:C
% a_i=trn_cluster{ii,1};
% for j=1:size(a_i,2)
% a_ii=a_i(:,j);
% a=reshape(a_ii,20,20);
% %imshow(a);
% %subplot(4,20,j);
% end
% end       
    end
    
   %% 将所有测试集进行子空间聚类
    for i = 1 : size(X_Test,1)
        tst_i = X_Test{i,1};
        tst_cluster(:,i) = HDC_local_model(tst_i,K,C);
%         if size(tst_i,2) >= sub_dim*C
%            tst_cluster(:,i) = HDC_local_model(tst_i,K,C);
%         else
%             tst_cluster{1,i} = tst_i;
%         end
%%% 可视化分簇
% for ii=1:C
% a_i=trn_cluster{ii,1};
% for j=1:size(a_i,2)
% a_ii=a_i(:,j)/255;
% a=reshape(a_ii,20,20);
% %imshow(a);
% %subplot(4,20,j);
% end
% end       
    end 
    
   ss=1;
   %% 将所有训练集进行子空间聚类
   for i=1:size(trn_cluster,2)
       for j=1:size(trn_cluster,1)
           if size(trn_cluster{j,i},2)>=sub_dim  
             trn_subspace(:,:,ss)=Compute_Subspace(trn_cluster{j,i},sub_dim);
             Train_label(1,ss)=floor((i-1)/num_sets_per_class)+1;
             ss=ss+1;
           end
       end
   end
    
    %% 将聚类完的所有测试集提取子空间
    for i = 1 : size(tst_cluster,1)
        for j = 1 : size(tst_cluster,2)
            if size(tst_cluster{i,j},2)>=sub_dim                
               tst_subspace{1,j}(:,:,i) = Compute_Subspace(tst_cluster{i,j},sub_dim);
            end
        end
    end

    %% 对每个测试集子空间进行稀疏编码
    for i = 1 : size(tst_subspace,2)
        X = tst_subspace{1,i};
        A = trn_subspace;
        %%% 初始化参数 %%%
        m = size(A,3);
        N = size(X,3);


        lambda1=0.1;
        lambda2=0.1;
        y_hat = GLSC_func(X,A,Train_label,lambda1,lambda2);
       
        tst_label(1,i) = y_hat;
         
        clear X A m N Y_Init v y_hat;
    end
   %% 计算准确率
    l_size=length(tst_label);
    l_count=0;
    for sx=1:l_size
        if tst_label(sx)==Test_label(sx)
            l_count=l_count+1;
        end
    end

    accuracy_rate(1,num)=l_count/l_size
end

mean_accuracy=mean(accuracy_rate);
var_accuracy=var(accuracy_rate);


































