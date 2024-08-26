function result_C = HDC_local_model(X,K,split_number)
% 构建局部线性模型
%
%  Input:
%       X: 数据集X=[x1,x2,...,xN],每个xi都是一个样本，尺寸是D*N
%       K: K个近邻点
%  Output:
%       result_C: 是一个cell，每个里面是一个Ci

[~,X_size]=size(X);

%%% 初始化 %%%
X_temp=X;
C={};
C_temp=X;
C_L=[];
C_R=[];
D={};
HH={};
ni=1;
sigma=1.2;
small_distance=2*10^(-6);
[~,D_G,H,S_init]=nonlinearity_score(X,K);
S_max=S_init;
SY=[];
%split_number=6;   %最大分成局部模型的个数

%%% 分割迭代 %%%
% while S_max>sigma
if S_max<=sigma
    C{1}=X;
end
while ni<split_number && S_max>sigma
      %%% 从测地矩阵D_G中选取距离最远的两个种子点x_index,y_index
      [x,y]=find(D_G==max(max(D_G)));
      x_index=x(1);
      y_index=y(1);
      C_L=X(:,x_index);
      C_R=X(:,y_index);
      
      %%% P_L和P_R分别表示各自的K个最近邻点
      P_L_index=H(:,x_index);
      for i=1:length(P_L_index)
          P_L(:,i)=C_temp(:,P_L_index(i));
      end
      P_R_index=H(:,y_index);
      for i=1:length(P_R_index)
          P_R(:,i)=C_temp(:,P_R_index(i));
      end
      
      C_temp(:,[x_index,y_index])=[];
      clear H;
      %%% 更新C_L，C_R，C_temp
      [C_L_temp,~,ib]=intersect(P_L',C_temp','rows');
      C_L_transform=union(C_L',C_L_temp,'rows');
      C_L=C_L_transform';
      C_temp(:,ib)=[];
          
      [C_R_temp,~,jb]=intersect(P_R',C_temp','rows');
      C_R_transform=union(C_R',C_R_temp,'rows');
      C_R=C_R_transform';
      C_temp(:,jb)=[];
      
      while isempty(C_temp)==0
          [~,size_C_tmp]=size(C_temp);
          label_L=1:size_C_tmp;
          if size_C_tmp>K 
              [~,neighbor_L]=top_K_neighbors(P_L,C_temp,label_L,K);
              H_L=neighbor_L;
              H_unique_L=unique(H_L);
              P_L=C_temp(:,H_unique_L);
          else
              P_L=[];
              %%% 计算P_L的平均值
              CL_mean=mean(C_L,2);
              CR_mean=mean(C_R,2);
              for tt=1:size(C_temp,2)
                  LTemp=norm(CL_mean-C_temp(:,tt));
                  RTemp=norm(CR_mean-C_temp(:,tt));
                  if LTemp<=RTemp
                      P_L=[P_L C_temp(:,tt)];
                  end
              end
          end
          
          label_R=1:size_C_tmp;
          if size_C_tmp>K 
              [~,neighbor_R]=top_K_neighbors(P_R,C_temp,label_R,K);
              H_R=neighbor_R;
              H_unique_R=unique(H_R);
              P_R=C_temp(:,H_unique_R);
          else
              if isempty(P_L)==1
                  P_R=C_temp;
              else
                  P_R_trans=setdiff(C_temp',P_L','rows');
                  P_R=P_R_trans';
              end
          end
          
          %%% 更新C_L,C_R,C_temp
          if isempty(P_L)~=1 
              [C_L_temp,~,ib_tp]=intersect(P_L',C_temp','rows');
              C_L_transform=union(C_L',C_L_temp,'rows');
              C_L=C_L_transform';
              C_temp(:,ib_tp)=[];
          end
          
          if isempty(P_R)~=1
              [C_R_temp,~,jb_tp]=intersect(P_R',C_temp','rows');
              C_R_transform=union(C_R',C_R_temp,'rows');
              C_R=C_R_transform';
              C_temp(:,jb_tp)=[];
          end
      end
      clear C_temp;
      %%% 对C_L C_R求S，选取最大的
      [~,D_G_L,L_H,S_L]=nonlinearity_score(C_L,K);
      [~,D_G_R,R_H,S_R]=nonlinearity_score(C_R,K);
      SY(ni)=S_L;
      C{ni}=C_L;
      D{ni}=D_G_L;
      HH{ni}=L_H;
      
      SY(ni+1)=S_R;
      C{ni+1}=C_R;
      D{ni+1}=D_G_R;
      HH{ni+1}=R_H;
      
      [SY,loc]=sort(SY);
      C=C(loc);  % 将C,D和SY都从小到大排序
      D=D(loc);
      HH=HH(loc);
      loc_max=length(loc);
      D_G=D{loc_max};
      C_temp=C{loc_max};
      H=HH{loc_max};
      S_max=SY(loc_max);
      
      ni=ni+1;
      clear C_L C_R;
      clear P_L P_R;
      
end

%%%%%%%%%%% 将尺寸小于6的扩成6 %%%%%%%%%%%%
s_C=size(C,2);
if s_C<split_number
   sm_C=size(C{s_C},2);
   difference=1+split_number-s_C;
   tmp=fix(sm_C/difference);
   for iii=s_C+1:split_number
       C{iii}=C{s_C}(:,1:tmp);
       C{s_C}(:,1:tmp)=[];
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result_C=C;


          
          
























































