%% Newton Interpolation and Runge Phenomenon with high order interpolation
%  the original function is f(x) = 1/(1+25x^2), x∈[-1,1]

close all;clear all;clc;

% 1. when x obeys uniform distribution
% original function：
f = @(x)1./(1+25*x.^2);

xi = linspace(-1,1,100);    % original interpolation points
x1 = linspace(-1,1,6);      % order = 5
x2 = linspace(-1,1,11);     % order = 10
x3 = linspace(-1,1,21);     % order = 20

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(1)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('original','order=5','order=10','order=20');
xlabel('x');ylabel('y');title('x obeys uniform distribution ');
axis([-1,1,-0.2,1.2]);


%  Then we want to see whether different distributions of x will lead to different results
%  Simply，we can use  the elementary functions in math  as below  
%  quadratic function：   f = @(t)(t.^2-1),t∈[0,sqrt(2)]; 
%  power function:        f = sqrt(t)-1,t∈[0,4];
%  exponential function： f = exp(t)-2,t∈[0,log(3)];
%  logarithmic function： f = log(t),t∈[exp(1).^-1,exp(1)];
%  consine function：     f = -cos(t),t∈[0,pi]; 
%  sine function:         f = sin(t),t∈[-pi/2,pi/2]

% 1.quadratic_distribution
f1 = @(t)(t.^2-1);
t = linspace(0,sqrt(2),100); xi = f1(t);
x1 = f1(linspace(0,sqrt(2),6));      % order = 5
x2 = f1(linspace(0,sqrt(2),11));     % order = 10
x3 = f1(linspace(0,sqrt(2),21));     % order = 20

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(2)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('original','order=5','order=10','order=20');
xlabel('x');ylabel('y');
title(['$x=t^2-1, t \in [0, \sqrt{2}]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);

% 2.power_distribution
f2 = @(t)(sqrt(t)-1);
t = linspace(0,4,100); xi = f2(t);
x1 = f2(linspace(0,4,6));      % order = 5
x2 = f2(linspace(0,4,11));     % order = 10
x3 = f2(linspace(0,4,21));     % order = 20

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(3)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('original','order=5','order=10','order=20');
xlabel('x');ylabel('y');
title(['$x=\sqrt{t}-1, t \in [0, 4]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);

% 3.exponential_distribution
f3 = @(t)exp(t)-2;           % f = exp(t)-2,t∈[0,log(3)];
t = linspace(0,log(3),100); xi = f3(t);
x1 = f3(linspace(0,log(3),6));      % order = 5
x2 = f3(linspace(0,log(3),11));     % order = 10
x3 = f3(linspace(0,log(3),21));     % order = 20

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(4)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('original','order=5','order=10','order=20');
xlabel('x');ylabel('y');
title(['$x=\mathrm{e}^{t}-2, t \in [0, ln3]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);

% 4.logarithmic_distribution
f4 = @(t)log(t);           % f = ln(t),t∈[e^-1,e];
t = linspace(exp(1)^-1,exp(1),100); xi = f4(t);
x1 = f4(linspace(exp(1)^-1,exp(1),6));      % order = 5
x2 = f4(linspace(exp(1)^-1,exp(1),11));     % order = 10
x3 = f4(linspace(exp(1)^-1,exp(1),21));     % order = 20

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(5)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('original','order=5','order=10','order=20');
xlabel('x');ylabel('y');
title(['$x=ln(t), t \in [\mathrm{e}^{-1}, \mathrm{e}]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);

% 5.cosine_distribution
f5 = @(t)-cos(t);           %f = -cos(t);t∈[0,pi];
t = linspace(0,pi,100); xi = f5(t);
x1 = f5(linspace(0,pi,6));      % order = 5
x2 = f5(linspace(0,pi,11));     % order = 10
x3 = f5(linspace(0,pi,21));     %  order = 20

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(6)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('original','order=5','order=10','order=20');
xlabel('x');ylabel('y');
title(['$x=-cos(t), t \in [0, \pi]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);

% 5.sine_distribution
f6 = @(t)sin(t);           %f = sin(t);t∈[-pi/2,pi/2];
t = linspace(-pi/2,pi/2,100); xi = f6(t);
x1 = f6(linspace(-pi/2,pi/2,6));      % order = 5
x2 = f6(linspace(-pi/2,pi/2,11));     % order = 10
x3 = f6(linspace(-pi/2,pi/2,21));     % order = 20

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(7)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('original','order=5','order=10','order=20');
xlabel('x');ylabel('y');
title(['$x=sin(t), t \in [-\frac{\pi}{2}, \frac{\pi}{2}]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);

% Actually when x satisfies the chebyshev polynomial, the interpolation function convergents
% xk = cos((2k+1)*pi/(2(n+1))), k =n,n-1,,,0;
f7 = @(k)cos(((2*k+1)*pi)./(2*length(k)));
k = linspace(99,0,100); xi = f7(k);
x1 = f7(linspace(5,0,6));      % order = 5
x2 = f7(linspace(10,0,11));     %order = 10
x3 = f7(linspace(20,0,21));     %order = 20

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(8)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('original','order=5','order=10','order=20');
xlabel('x');ylabel('y');
title(['$x_k=cos(\frac{(2k+1)\pi}{2(n+1)}), k=n,n-1,...,0$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
