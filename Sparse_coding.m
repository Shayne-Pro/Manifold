function Y = Sparse_coding(X,A,Y_Init,alpha,beta,lambda,miu,v,L,B)
% solve X�����е�ϡ��ϵ�� ����Y_Init,miu,v�ǳ�ʼ���Ĳ�����alpha,beta,lambda�����򻯲�����X,A,L,B��֪
% Input:
%     X: һ�����Լ���ͨ��ϡ���ӿռ���࣬�ֳ�m�����Լ�����ͨ�����ζԳ�ͶӰ��Ϊm������   m*m
%     A: �ֵ� m*N
%     Y_Init: ��ʼ����ϡ��ϵ������ N*m
%     alpha: LLC������ǰ�����򻯲���
%     beta: Laplaceͼǰ�����򻯲���
%     lambda: y��һ����ǰ�����򻯲���
%     miu: Lagrange����ǰ�ĳ�ʼ��ϵ����ʼ��
%     v: Lagrange�����еĲ���vector��ʼ��
%     L: Laplaceͼ: ���Լ��ֳ�m���ӿռ䣬��m���ӿռ��������ϵ��໥����
%     B: ���Լ���ѵ�����ľ���

[d,m] = size(X);
[d,N] = size(A);
Y = Y_Init;
max_Iter = 1e6;
tol = 1e-6;
fobj3_temp = 0;
% for i = 1:m
%     S_i = diag(B(:,i)'*diag(B(:,i)));
%     Y_i = Y(:,i);
%     fobj3_temp = fobj3_temp+Y_i'*S_i*Y_i;
% end
% fobj_old = [norm(X-A*Y,'fro')]^2 + beta*trace(Y*L*Y') + alpha*fobj3_temp + lambda*sum(sum(abs(Y)));

% ��ϡ��ϵ�������ÿһ�е������£����µ�i��ʱ�̶������У�
%start loop
% iter = 0;
% while iter<max_Iter
%     iter = iter + 1;
    for i = 1:m
        x = X(:,i);
        Y(:,i) = Ialm(i,x,A,Y,alpha,beta,lambda,miu,v,L,B);
    end
%     fobj1 = [norm(X-A*Y,'fro')]^2;
%     fobj2 = beta*trace(Y*L*Y');
%     fobj3_temp = 0;
%     for i = 1:m
%         S_i = diag(B(:,i)'*diag(B(:,i)));
%         Y_i = Y(:,i);
%         fobj3_temp = fobj3_temp+Y_i'*S_i*Y_i;
%     end
%     fobj3 = alpha*fobj3_temp;
%     fobj4 = lambda*sum(sum(abs(Y)));
%     fobj_new = fobj1+fobj2+fobj3+fobj4;
%     stopC = abs(fobj_new-fobj_old);
%     plot_residual(iter) = stopC;
%     if iter==1 || mod(iter,50)==0 || stopC<tol
%         disp(['iter ' num2str(iter) ',stopALM=' num2str(stopC,'%2.3e')]);
%     end
%     if stopC<tol
%         break;
%     end
%     fobj_old = fobj_new;
%     fobj_search(iter) = fobj_old;
%     clear fobj1 fobj2 fobj3 fobj4 fobj_new fobj_temp;
% end
% 
% %%%% Ŀ�꺯��ֵ�ı仯�Ƿ����� %%%%
%  plot(fobj_search);
% 
% end



































