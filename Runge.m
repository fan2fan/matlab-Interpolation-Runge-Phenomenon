%% lec7 ���ʵ��ţ���ڲ幫ʽ���������³������n�׶���ʽ���
%   f(x) = 1/(1+25x^2), x��[-1,1]
% 1. ��xi��[-1,1]�о���ȡֵʱ������n�����ӣ�N(x)�Ƿ���ӽӽ�f(x) -> linspace(-1,1,n)
% 2. ���xi���Ǿ��ȷֲ�������n�����ӣ�N(x)��f(x)���������Ƿ����ı�-> 
%   �۲���Ⱥ�����ʽ�ķֲ��������Ե�Ӱ�죺
%   ʹ�ú�����f = @(t)(t.^2-1),t��[0,1]; f = sqrt(t)-1,t��[0,4];
%            f = ln(t),t��[e^-1,e];
%            f = exp(t)-2,t��[0,log(3)];
%            f = -cos(t),t��[0,pi];f = sin(t);t��[-pi/2,pi/2];
% 3. �ܷ��ҵ�һ��xi�ֲ���ʹ����n�����ӣ�N(x)��������f(x)?

close all;clear all;clc;

%% ��Ŀ1��x����ȡֵʱ��ͬ������ֵ
% x����ȡֵʵ���϶�Ӧ����һ�κ���y=ax+b;
% n������Ӧn-1��
% ����ԭ������
f = @(x)1./(1+25*x.^2);

xi = linspace(-1,1,100);    %��ֵx��
x1 = linspace(-1,1,6);      %5��
x2 = linspace(-1,1,11);     %10��
x3 = linspace(-1,1,21);     %20��

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(1)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('ԭ����','5�ײ�ֵ','10�ײ�ֵ','20�ײ�ֵ');
xlabel('x');ylabel('y');title('xȡ��ֵʱ5,10,20��ţ�ٲ�ֵ');
axis([-1,1,-0.2,1.2]);
disp('xȡ��ֵʱ��Ӧ�Ĳ�ֵ����ʽ');
disp(['5�ף� ', polyFun1]);disp(['10�ף� ', polyFun2]);disp(['20�ף� ', polyFun3]);

%% ��Ŀ2��x�Ǿ��ȷֲ�ʱ��ͬ������ֵ
%   Ҫʹx�Ǿ��ȷֲ������Զ���ѧ�ĳ��Ⱥ������г���        
%   �����ݺ�����  f = @(t)(t.^2-1),t��[0,sqrt(2)]; f = sqrt(t)-1,t��[0,4];
%   ����ָ��������f = exp(t)-2,t��[0,log(3)];
%   ���ڶ���������f = log(t),t��[exp(1).^-1,exp(1)];
%   �������Ǻ�����f = -cos(t),t��[0,pi]; f = sin(t),t��[-pi/2,pi/2]


% 1.�����ݺ����ֲ�������������ʽ��
% 1.1���ڶ��κ���
f1 = @(t)(t.^2-1);           %x = t^2 - 1;t��[0,sqrt(2)];
t = linspace(0,sqrt(2),100); xi = f1(t);
x1 = f1(linspace(0,sqrt(2),6));      %5��
x2 = f1(linspace(0,sqrt(2),11));     %10��
x3 = f1(linspace(0,sqrt(2),21));     %20��

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(2)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('ԭ����','5�ײ�ֵ','10�ײ�ֵ','20�ײ�ֵ');
xlabel('x');ylabel('y');
title(['$x=t^2-1, t \in [0, \sqrt{2}]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=t^2-1, t��[0, sqrt(2)]');
disp(['5�ף� ', polyFun1]);disp(['10�ף� ', polyFun2]);disp(['20�ף� ', polyFun3]);

% 1.2����x�Ŀ�����
f2 = @(t)(sqrt(t)-1);           %f = sqrt(t)-1;t��[0,4];
t = linspace(0,4,100); xi = f2(t);
x1 = f2(linspace(0,4,6));      %5��
x2 = f2(linspace(0,4,11));     %10��
x3 = f2(linspace(0,4,21));     %20��

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(3)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('ԭ����','5�ײ�ֵ','10�ײ�ֵ','20�ײ�ֵ');
xlabel('x');ylabel('y');
title(['$x=\sqrt{t}-1, t \in [0, 4]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=sqrt(t)-1, t��[0, 4]');
disp(['5�ף� ', polyFun1]);disp(['10�ף� ', polyFun2]);disp(['20�ף� ', polyFun3]);

% 2����ָ�������ֲ���
f3 = @(t)exp(t)-2;           % f = exp(t)-2,t��[0,log(3)];
t = linspace(0,log(3),100); xi = f3(t);
x1 = f3(linspace(0,log(3),6));      %5��
x2 = f3(linspace(0,log(3),11));     %10��
x3 = f3(linspace(0,log(3),21));     %20��

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(4)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('ԭ����','5�ײ�ֵ','10�ײ�ֵ','20�ײ�ֵ');
xlabel('x');ylabel('y');
title(['$x=\mathrm{e}^{t}-2, t \in [0, ln3]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=e^t-2, t��[0, ln3]');
disp(['5�ף� ', polyFun1]);disp(['10�ף� ', polyFun2]);disp(['20�ף� ', polyFun3]);

% 3���ڶ��������ֲ���
f4 = @(t)log(t);           % f = ln(t),t��[e^-1,e];
t = linspace(exp(1)^-1,exp(1),100); xi = f4(t);
x1 = f4(linspace(exp(1)^-1,exp(1),6));      %5��
x2 = f4(linspace(exp(1)^-1,exp(1),11));     %10��
x3 = f4(linspace(exp(1)^-1,exp(1),21));     %20��

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(5)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('ԭ����','5�ײ�ֵ','10�ײ�ֵ','20�ײ�ֵ');
xlabel('x');ylabel('y');
title(['$x=ln(t), t \in [\mathrm{e}^{-1}, \mathrm{e}]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=ln(t), t��[e-1, e]');
disp(['5�ף� ', polyFun1]);disp(['10�ף� ', polyFun2]);disp(['20�ף� ', polyFun3]);


% 4�������Ǻ����ֲ���
% 4.1��������
f5 = @(t)-cos(t);           %f = -cos(t);t��[0,pi];
t = linspace(0,pi,100); xi = f5(t);
x1 = f5(linspace(0,pi,6));      %5��
x2 = f5(linspace(0,pi,11));     %10��
x3 = f5(linspace(0,pi,21));     %20��

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(6)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('ԭ����','5�ײ�ֵ','10�ײ�ֵ','20�ײ�ֵ');
xlabel('x');ylabel('y');
title(['$x=-cos(t), t \in [0, \pi]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=-cos(t), t��[0, ��]');
disp(['5�ף� ', polyFun1]);disp(['10�ף� ', polyFun2]);disp(['20�ף� ', polyFun3]);

% 4.2��������:
f6 = @(t)sin(t);           %f = sin(t);t��[-pi/2,pi/2];
t = linspace(-pi/2,pi/2,100); xi = f6(t);
x1 = f6(linspace(-pi/2,pi/2,6));      %5��
x2 = f6(linspace(-pi/2,pi/2,11));     %10��
x3 = f6(linspace(-pi/2,pi/2,21));     %20��

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(7)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('ԭ����','5�ײ�ֵ','10�ײ�ֵ','20�ײ�ֵ');
xlabel('x');ylabel('y');
title(['$x=sin(t), t \in [-\frac{\pi}{2}, \frac{\pi}{2}]$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('x=sin(t), t��[-��/2, ��/2]');
disp(['5�ף� ', polyFun1]);disp(['10�ף� ', polyFun2]);disp(['20�ף� ', polyFun3]);


%% ��3 �ҳ�һ��xi�ֲ�ʹN(x)������f(x)
% ����ͼ��ͬ�ĳ��Ⱥ�����ֵ������֪�����Ǻ�������Ч��ʹ����ʽ������ԭ����������Խ�ߣ��������Խ��
% ʵ����xi�����б�ѩ����ʱ��ʹN(x)��Ч������f(x),��ʱ��xk = cos((2k+1)*pi/(2(n+1))), k =n,n-1,,,0;
f7 = @(k)cos(((2*k+1)*pi)./(2*length(k)));
k = linspace(99,0,100); xi = f7(k);
x1 = f7(linspace(5,0,6));      %5��
x2 = f7(linspace(10,0,11));     %10��
x3 = f7(linspace(20,0,21));     %20��

[polyFun1,N1] = NewtInt(x1,f(x1),xi);
[polyFun2,N2] = NewtInt(x2,f(x2),xi);
[polyFun3,N3] = NewtInt(x3,f(x3),xi);

figure(8)
hold on;
plot(xi,f(xi),'LineWidth',1.5); plot(xi,N1,'LineWidth',1.5); plot(xi,N2,'LineWidth',1.5); plot(xi,N3,'LineWidth',1.5);
grid on; legend('ԭ����','5�ײ�ֵ','10�ײ�ֵ','20�ײ�ֵ');
xlabel('x');ylabel('y');
title(['$x_k=cos(\frac{(2k+1)\pi}{2(n+1)}), k=n,n-1,...,0$'],'interpreter','latex','FontSize',13);
axis([-1,1,-0.2,1.2]);
disp('xk = cos((2k+1)*pi/(2(n+1))), k =n,n-1,,,0');
disp(['5�ף� ', polyFun1]);disp(['10�ף� ', polyFun2]);disp(['20�ף� ', polyFun3]);