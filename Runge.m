%% lec7 编程实现牛顿内插公式，并对如下程序进行n阶多项式拟合
%   f(x) = 1/(1+25x^2), x∈[-1,1]
% 1. 当xi在[-1,1]中均匀取值时，随着n的增加，N(x)是否更加接近f(x) -> linspace(-1,1,n)
% 2. 如果xi不是均匀分布，随着n的增加，N(x)对f(x)的收敛性是否发生改变-> 
%   观察初等函数形式的分布对收敛性的影响：
%   使用函数：f = @(t)(t.^2-1),t∈[0,1]; f = sqrt(t)-1,t∈[0,4];
%            f = ln(t),t∈[e^-1,e];
%            f = exp(t)-2,t∈[0,log(3)];
%            f = -cos(t),t∈[0,pi];f = sin(t);t∈[-pi/2,pi/2];
% 3. 能否找到一种xi分布，使随着n的增加，N(x)能收敛到f(x)?

close all;clear all;clc;

%% 题目1，x均匀取值时不同阶数插值
% x均匀取值实际上对应的是一次函数y=ax+b;
% n个数对应n-1阶
% 定义原函数：
f = @(x)1./(1+25*x.^2);

xi = linspace(-1,1,100);    %插值x点
x1 = linspace(-1,1,6);      %5阶
x2 = linspace(-1,1,11);     %10阶
x3 = linspace(-1,1,21);     %20阶

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(1)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('原函数','5阶插值','10阶插值','20阶插值');
xlabel('x');ylabel('y');title('x取均值时5,10,20阶牛顿插值');
axis([-1,1,-0.2,1.2]);
disp('x取均值时对应的插值多项式');
disp(['5阶： ', polyFun1]);disp(['10阶： ', polyFun2]);disp(['20阶： ', polyFun3]);

%% 题目2，x非均匀分布时不同阶数插值
%   要使x非均匀分布，可以对数学的初等函数进行尝试        
%   对于幂函数：  f = @(t)(t.^2-1),t∈[0,sqrt(2)]; f = sqrt(t)-1,t∈[0,4];
%   对于指数函数：f = exp(t)-2,t∈[0,log(3)];
%   对于对数函数：f = log(t),t∈[exp(1).^-1,exp(1)];
%   对于三角函数：f = -cos(t),t∈[0,pi]; f = sin(t),t∈[-pi/2,pi/2]


% 1.对于幂函数分布：（有两种形式）
% 1.1对于二次函数
f1 = @(t)(t.^2-1);           %x = t^2 - 1;t∈[0,sqrt(2)];
t = linspace(0,sqrt(2),100); xi = f1(t);
x1 = f1(linspace(0,sqrt(2),6));      %5阶
x2 = f1(linspace(0,sqrt(2),11));     %10阶
x3 = f1(linspace(0,sqrt(2),21));     %20阶

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(2)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('原函数','5阶插值','10阶插值','20阶插值');
xlabel('x');ylabel('y');
title(['$x=t^2-1, t \in [0, \sqrt{2}]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=t^2-1, t∈[0, sqrt(2)]');
disp(['5阶： ', polyFun1]);disp(['10阶： ', polyFun2]);disp(['20阶： ', polyFun3]);

% 1.2对于x的开方：
f2 = @(t)(sqrt(t)-1);           %f = sqrt(t)-1;t∈[0,4];
t = linspace(0,4,100); xi = f2(t);
x1 = f2(linspace(0,4,6));      %5阶
x2 = f2(linspace(0,4,11));     %10阶
x3 = f2(linspace(0,4,21));     %20阶

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(3)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('原函数','5阶插值','10阶插值','20阶插值');
xlabel('x');ylabel('y');
title(['$x=\sqrt{t}-1, t \in [0, 4]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=sqrt(t)-1, t∈[0, 4]');
disp(['5阶： ', polyFun1]);disp(['10阶： ', polyFun2]);disp(['20阶： ', polyFun3]);

% 2对于指数函数分布：
f3 = @(t)exp(t)-2;           % f = exp(t)-2,t∈[0,log(3)];
t = linspace(0,log(3),100); xi = f3(t);
x1 = f3(linspace(0,log(3),6));      %5阶
x2 = f3(linspace(0,log(3),11));     %10阶
x3 = f3(linspace(0,log(3),21));     %20阶

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(4)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('原函数','5阶插值','10阶插值','20阶插值');
xlabel('x');ylabel('y');
title(['$x=\mathrm{e}^{t}-2, t \in [0, ln3]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=e^t-2, t∈[0, ln3]');
disp(['5阶： ', polyFun1]);disp(['10阶： ', polyFun2]);disp(['20阶： ', polyFun3]);

% 3对于对数函数分布：
f4 = @(t)log(t);           % f = ln(t),t∈[e^-1,e];
t = linspace(exp(1)^-1,exp(1),100); xi = f4(t);
x1 = f4(linspace(exp(1)^-1,exp(1),6));      %5阶
x2 = f4(linspace(exp(1)^-1,exp(1),11));     %10阶
x3 = f4(linspace(exp(1)^-1,exp(1),21));     %20阶

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(5)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('原函数','5阶插值','10阶插值','20阶插值');
xlabel('x');ylabel('y');
title(['$x=ln(t), t \in [\mathrm{e}^{-1}, \mathrm{e}]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=ln(t), t∈[e-1, e]');
disp(['5阶： ', polyFun1]);disp(['10阶： ', polyFun2]);disp(['20阶： ', polyFun3]);


% 4对于三角函数分布：
% 4.1对于余弦
f5 = @(t)-cos(t);           %f = -cos(t);t∈[0,pi];
t = linspace(0,pi,100); xi = f5(t);
x1 = f5(linspace(0,pi,6));      %5阶
x2 = f5(linspace(0,pi,11));     %10阶
x3 = f5(linspace(0,pi,21));     %20阶

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(6)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('原函数','5阶插值','10阶插值','20阶插值');
xlabel('x');ylabel('y');
title(['$x=-cos(t), t \in [0, \pi]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=-cos(t), t∈[0, π]');
disp(['5阶： ', polyFun1]);disp(['10阶： ', polyFun2]);disp(['20阶： ', polyFun3]);

% 4.2对于正弦:
f6 = @(t)sin(t);           %f = sin(t);t∈[-pi/2,pi/2];
t = linspace(-pi/2,pi/2,100); xi = f6(t);
x1 = f6(linspace(-pi/2,pi/2,6));      %5阶
x2 = f6(linspace(-pi/2,pi/2,11));     %10阶
x3 = f6(linspace(-pi/2,pi/2,21));     %20阶

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(7)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('原函数','5阶插值','10阶插值','20阶插值');
xlabel('x');ylabel('y');
title(['$x=sin(t), t \in [-\frac{\pi}{2}, \frac{\pi}{2}]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=sin(t), t∈[-π/2, π/2]');
disp(['5阶： ', polyFun1]);disp(['10阶： ', polyFun2]);disp(['20阶： ', polyFun3]);


%% 题3 找出一种xi分布使N(x)收敛到f(x)
% 从上图不同的初等函数插值，可以知道三角函数能有效的使多项式收敛到原函数，阶数越高，收敛结果越好
% 实际上xi满足切比雪夫结点时能使N(x)有效收敛于f(x),此时有xk = cos((2k+1)*pi/(2(n+1))), k =n,n-1,,,0;
f7 = @(k)cos(((2*k+1)*pi)./(2*length(k)));
k = linspace(99,0,100); xi = f7(k);
x1 = f7(linspace(5,0,6));      %5阶
x2 = f7(linspace(10,0,11));     %10阶
x3 = f7(linspace(20,0,21));     %20阶

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(8)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('原函数','5阶插值','10阶插值','20阶插值');
xlabel('x');ylabel('y');
title(['$x_k=cos(\frac{(2k+1)\pi}{2(n+1)}), k=n,n-1,...,0$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('xk = cos((2k+1)*pi/(2(n+1))), k =n,n-1,,,0');
disp(['5阶： ', polyFun1]);disp(['10阶： ', polyFun2]);disp(['20阶： ', polyFun3]);