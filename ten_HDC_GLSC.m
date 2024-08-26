function [mean_accuracy,var_accuracy]=ten_HDC_GLSC(sub_dim,C,K)
% ����GLSC��Extended yaleB�ϵ�ʮ��ƽ��ֵ
% Input:
%     sub_dim: Grassmann�����ϵ��ӿռ�ά�ȣ��൱��G(p,d)�е�p
%     alpha: LLC������ǰ�����򻯲���
%     beta: Laplaceͼǰ�����򻯲���
%     lambda: y��һ����ǰ�����򻯲���
%     miu: Lagrange����ǰ�ĳ�ʼ��ϵ����ʼ��
%     v: Lagrange�����еĲ���vector��ʼ��
%     sigma1: Gassian�˺�����ϵ��
%     sigma2: Gassian�˺�����ϵ��
% Output:
%     accuracy: eth80ʮ�ν���ʵ���׼ȷ��
%     mean_accuracy: ʮ�ν���ʵ���ƽ��׼ȷ��

%% ��ʼ��ѵ���������Լ�
load random_index_extyaleB;
load Extended_yaleB;
num_class=28;
num_sets_per_class=3;
accuracy_rate=zeros(1,10);
for num=1:1
    %% ��ʼ��ѵ���������Լ��������ǵı�ǩ
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
    
   %% ������ѵ���������ӿռ����
    for i = 1 : size(X_Train,1)
        trn_i = X_Train{i,1};
        if size(trn_i,2) >= sub_dim*C
            trn_cluster(:,i) = HDC_local_model(trn_i,K,C);
        else
            trn_cluster{1,i} = trn_i;
          %  trn_cluster{2:C,i} = [];
        end
%%% ���ӻ��ִ�
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
    
   %% �����в��Լ������ӿռ����
    for i = 1 : size(X_Test,1)
        tst_i = X_Test{i,1};
        tst_cluster(:,i) = HDC_local_model(tst_i,K,C);
%         if size(tst_i,2) >= sub_dim*C
%            tst_cluster(:,i) = HDC_local_model(tst_i,K,C);
%         else
%             tst_cluster{1,i} = tst_i;
%         end
%%% ���ӻ��ִ�
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
   %% ������ѵ���������ӿռ����
   for i=1:size(trn_cluster,2)
       for j=1:size(trn_cluster,1)
           if size(trn_cluster{j,i},2)>=sub_dim  
             trn_subspace(:,:,ss)=Compute_Subspace(trn_cluster{j,i},sub_dim);
             Train_label(1,ss)=floor((i-1)/num_sets_per_class)+1;
             ss=ss+1;
           end
       end
   end
    
    %% ������������в��Լ���ȡ�ӿռ�
    for i = 1 : size(tst_cluster,1)
        for j = 1 : size(tst_cluster,2)
            if size(tst_cluster{i,j},2)>=sub_dim                
               tst_subspace{1,j}(:,:,i) = Compute_Subspace(tst_cluster{i,j},sub_dim);
            end
        end
    end

    %% ��ÿ�����Լ��ӿռ����ϡ�����
    for i = 1 : size(tst_subspace,2)
        X = tst_subspace{1,i};
        A = trn_subspace;
        %%% ��ʼ������ %%%
        m = size(A,3);
        N = size(X,3);


        lambda1=0.1;
        lambda2=0.1;
        y_hat = GLSC_func(X,A,Train_label,lambda1,lambda2);
       
        tst_label(1,i) = y_hat;
         
        clear X A m N Y_Init v y_hat;
    end
   %% ����׼ȷ��
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


































