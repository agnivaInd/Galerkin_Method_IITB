clear
close all
clc

res_equal = readmatrix('equally_spaced_interpolation.csv');
res_legendre = readmatrix('legendre_interpolation.csv');
res_lobatto = readmatrix('lobatto_interpolation.csv');

derivative_equal = readmatrix('derivative_equally_spaced.csv');
derivative_legendre = readmatrix('derivative_legendre.csv');
derivative_lobatto = readmatrix('derivative_lobatto.csv');

pts = numel(res_equal);
delx = 2/(pts-1);
x = -1:delx:1;
analytic_result = sin(pi.*x/2);
analytic_derivative = (pi/2)*cos(pi.*x/2);

figure(1)
plot(x,analytic_result,LineWidth=1.5,Color='b')
hold on
plot(x,res_equal,LineStyle="none",Marker="o",MarkerSize=6,LineWidth=1.5)
hold off
xlabel('x')
ylabel('f(x)')
legend('Analytical','Numerical - Equally Spaced Points',Location='northwest')
xlim([-1 1])
ylim([-1 1])
dim = [0.2 0.5 0.2 0.3];
str = {'N = 3'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

figure(2)
plot(x,analytic_result,LineWidth=1.5,Color='b')
hold on
plot(x,res_legendre,LineStyle="none",Marker="square",MarkerSize=6,LineWidth=1.5)
hold off
xlabel('x')
ylabel('f(x)')
legend('Analytical','Numerical - Legendre Points',Location='northwest')
xlim([-1 1])
ylim([-1 1])
dim = [0.2 0.5 0.2 0.3];
str = {'N = 3'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

figure(3)
plot(x,analytic_result,LineWidth=1.5,Color='b')
hold on
plot(x,res_lobatto,LineStyle="none",Marker="hexagram",MarkerSize=6,LineWidth=1.5)
hold off
xlabel('x')
ylabel('f(x)')
legend('Analytical','Numerical - Lobatto Points',Location='northwest')
xlim([-1 1])
ylim([-1 1])
dim = [0.2 0.5 0.2 0.3];
str = {'N = 3'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

figure(4)
plot(x,analytic_derivative,LineWidth=1.5,Color='b')
hold on
plot(x,derivative_lobatto,LineStyle="none",Marker="hexagram",MarkerSize=6,LineWidth=1.5)
hold off
xlabel('x')
ylabel('df/dx')
xlim([-1 1])
legend('Analytical Derivative','Numerical Derivative - Lobatto Points',Location='south')
dim = [0.2 0.5 0.2 0.3];
str = {'N = 3'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

N = [1,3,5];
l2_lobatto = [0.209146,0.00661555,9.16535e-05];
l2_legendre = [0.182898,0.00645124,0.00010623];
l2_equal = [0.209146,0.00826785,0.000154893];

figure(5)
plot(N,l2_equal,LineWidth=1.5)
hold on
plot(N,l2_legendre,linewidth=1.5)
hold on
plot(N,l2_lobatto,LineWidth=1.5)
legend('Equally Spaced Points','Legendre Points','Lobatto Points')
xlabel('N')
ylabel('L2 Norm')
hold off
