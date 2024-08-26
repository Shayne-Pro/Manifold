function y = Ialm(num,x,A,Y,alpha,beta,lambda,miu,v,L,B)
% solve problem: min(y)||x-Ay||2 + beta(tr(YLY')) + alpha(y'Sy) + lambda||y||_1
% Input:
%     num: 更新的测试样本的列数i
%     x: 测试样本m*1
%     A: 字典m*N
%     Y__i: Y中除了要更新的y_i 其他的所有列
%     alpha: LLC二范数前的正则化参数
%     beta: Laplace图前的正则化参数
%     lambda: y的一范数前的正则化参数
%     miu: Lagrange乘子前的初始化系数
%     v: Lagrange乘子中的参数 m*1
%     L: Laplace图: 测试集分成m个子空间，这m个子空间在流形上的相互距离
%     B: 测试集到训练集的距离
% Output:
%     y: 稀疏系数

%% 初始化
maxIter = 1e6;
maxmiu=10^6;
tol = 1e-8;
rou = 1.1;
[d,N] = size(A); % d表示维数，N表示训练集或字典的个数
y = zeros(d,1);  % 初始化y
z = zeros(d,1);  % 初始化z
B_i = B(:,num);
S_i = diag(B_i)'*diag(B_i);
L_i = L(:,num);
L_ii = L_i(num);
Y__i = Y;  
Y__i(:,num) = [];
L_i_i = L_i;
L_i_i(num) = [];
% v = zeros(d,1);   % 初始化v
% miu =0;    % 初始化miu
I_N = eye(N);

%% 交替迭代更新
%start loop
iter = 0;
while iter<maxIter
    iter = iter + 1;    
    % 更新z
    Temp_l = A'*A + beta*L_ii*I_N + alpha*S_i + miu*I_N;
    if size(Y__i,2) >= 1
        Temp_r = A'*x - beta*Y__i*L_i_i + miu*y - 0.5*v;
    else
        Temp_r = A'*x + miu*y - 0.5*v;
    end
    z = inv(Temp_l)*Temp_r;
    clear Temp_1 Temp_r;
    % 更新y
    temp = z+(1/(2*miu))*v;
    y = soft(temp,lambda/(2*miu));
    clear temp;
    % 更新v
    v = v+miu*(z-y);
    % 更新miu
    miu = min(rou*miu,maxmiu);
    % 检查收敛条件
    stopC = norm(y-z);
    plot_residual(iter) = stopC;
    if iter==1 || mod(iter,50)==0 || stopC<tol
        disp(['iter ' num2str(iter) ',stopALM=' num2str(stopC,'%2.3e')]);
    end
    if stopC<tol
        break;
    end
end

%%% 检查是否收敛
 plot(plot_residual);
 
end




























































