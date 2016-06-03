function [polyFun,yi] = NewtInt(x,y,xi)
n = length(x);
fx = zeros(n); fx(:,1) = y;
%根据n个点计算n-1阶差商
for j = 2:n
    for i = j:n
    fx(i,j) = (fx(i,j-1) - fx(i-1,j-1))/(x(i)-x(i-j+1));
    end
end
%根据上面计算得到的差商计算相应的插值
%下面通过矩阵相乘的形式计算出对应的插值，减少循环量
%   yi = [1   (xi1-x<1>)  (xi1-x<1>)(xi1-x<2>) ... (xi1-x<1>)(xi1-x<2>)...(xi1-x<n-1>)
%         1   (xi2-x<1>)  (xi2-x<1>)(xi2-x<2>) ... (xi2-x<1>)(xi2-x<2>)...(xi2-x<n-1>)
%         ...
%         1   (xim-x<1>)  (xim-x<1>)(xim-x<2>) ... (xim-x<1>)(xim-x<2>)...(xim-x<n-1>)]
%         *[b1;b2;...;bn]
xi = xi(:); %将xi转化为列向量
X = ones(length(xi),n);
for k = 2:n
    X(:,k) = X(:,k-1).*(xi-x(k-1));
end
yi = X*diag(fx);

%生成多项式函数
syms X b poly
Y=[1];
for k=2:n
    b=Y(end)*(X-x(k-1));%最后一项一次累乘
    Y=[Y b];%将新产生的项追加原来的矩阵
end
poly=Y*diag(fx);
polyFun = poly2str(sym2poly(poly),'x');
end

