function [polyFun,yi] = NewtInt(x,y,xi)
n = length(x);
fx = zeros(n); fx(:,1) = y;

for j = 2:n
    for i = j:n
    fx(i,j) = (fx(i,j-1) - fx(i-1,j-1))/(x(i)-x(i-j+1));
    end
end

xi = xi(:);
X = ones(length(xi),n);
for k = 2:n
    X(:,k) = X(:,k-1).*(xi-x(k-1));
end
yi = X*diag(fx);

syms X b poly
Y=[1];
for k=2:n
    b=Y(end)*(X-x(k-1));
    Y=[Y b];
end
poly=Y*diag(fx);
polyFun = poly2str(sym2poly(poly),'x');
end

