function y = Ialm(num,x,A,Y,alpha,beta,lambda,miu,v,L,B)
% solve problem: min(y)||x-Ay||2 + beta(tr(YLY')) + alpha(y'Sy) + lambda||y||_1
% Input:
%     num: ���µĲ�������������i
%     x: ��������m*1
%     A: �ֵ�m*N
%     Y__i: Y�г���Ҫ���µ�y_i ������������
%     alpha: LLC������ǰ�����򻯲���
%     beta: Laplaceͼǰ�����򻯲���
%     lambda: y��һ����ǰ�����򻯲���
%     miu: Lagrange����ǰ�ĳ�ʼ��ϵ��
%     v: Lagrange�����еĲ��� m*1
%     L: Laplaceͼ: ���Լ��ֳ�m���ӿռ䣬��m���ӿռ��������ϵ��໥����
%     B: ���Լ���ѵ�����ľ���
% Output:
%     y: ϡ��ϵ��

%% ��ʼ��
maxIter = 1e6;
maxmiu=10^6;
tol = 1e-8;
rou = 1.1;
[d,N] = size(A); % d��ʾά����N��ʾѵ�������ֵ�ĸ���
y = zeros(d,1);  % ��ʼ��y
z = zeros(d,1);  % ��ʼ��z
B_i = B(:,num);
S_i = diag(B_i)'*diag(B_i);
L_i = L(:,num);
L_ii = L_i(num);
Y__i = Y;  
Y__i(:,num) = [];
L_i_i = L_i;
L_i_i(num) = [];
% v = zeros(d,1);   % ��ʼ��v
% miu =0;    % ��ʼ��miu
I_N = eye(N);

%% �����������
%start loop
iter = 0;
while iter<maxIter
    iter = iter + 1;    
    % ����z
    Temp_l = A'*A + beta*L_ii*I_N + alpha*S_i + miu*I_N;
    if size(Y__i,2) >= 1
        Temp_r = A'*x - beta*Y__i*L_i_i + miu*y - 0.5*v;
    else
        Temp_r = A'*x + miu*y - 0.5*v;
    end
    z = inv(Temp_l)*Temp_r;
    clear Temp_1 Temp_r;
    % ����y
    temp = z+(1/(2*miu))*v;
    y = soft(temp,lambda/(2*miu));
    clear temp;
    % ����v
    v = v+miu*(z-y);
    % ����miu
    miu = min(rou*miu,maxmiu);
    % �����������
    stopC = norm(y-z);
    plot_residual(iter) = stopC;
    if iter==1 || mod(iter,50)==0 || stopC<tol
        disp(['iter ' num2str(iter) ',stopALM=' num2str(stopC,'%2.3e')]);
    end
    if stopC<tol
        break;
    end
end

%%% ����Ƿ�����
 plot(plot_residual);
 
end




























































