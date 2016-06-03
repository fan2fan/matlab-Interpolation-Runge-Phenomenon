function [polyFun,yi] = NewtInt(x,y,xi)
n = length(x);
fx = zeros(n); fx(:,1) = y;
%����n�������n-1�ײ���
for j = 2:n
    for i = j:n
    fx(i,j) = (fx(i,j-1) - fx(i-1,j-1))/(x(i)-x(i-j+1));
    end
end
%�����������õ��Ĳ��̼�����Ӧ�Ĳ�ֵ
%����ͨ��������˵���ʽ�������Ӧ�Ĳ�ֵ������ѭ����
%   yi = [1   (xi1-x<1>)  (xi1-x<1>)(xi1-x<2>) ... (xi1-x<1>)(xi1-x<2>)...(xi1-x<n-1>)
%         1   (xi2-x<1>)  (xi2-x<1>)(xi2-x<2>) ... (xi2-x<1>)(xi2-x<2>)...(xi2-x<n-1>)
%         ...
%         1   (xim-x<1>)  (xim-x<1>)(xim-x<2>) ... (xim-x<1>)(xim-x<2>)...(xim-x<n-1>)]
%         *[b1;b2;...;bn]
xi = xi(:); %��xiת��Ϊ������
X = ones(length(xi),n);
for k = 2:n
    X(:,k) = X(:,k-1).*(xi-x(k-1));
end
yi = X*diag(fx);

%���ɶ���ʽ����
syms X b poly
Y=[1];
for k=2:n
    b=Y(end)*(X-x(k-1));%���һ��һ���۳�
    Y=[Y b];%���²�������׷��ԭ���ľ���
end
poly=Y*diag(fx);
polyFun = poly2str(sym2poly(poly),'x');
end

