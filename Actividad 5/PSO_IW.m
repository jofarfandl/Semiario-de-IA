clear all
close all
clc

f = @(x,y) (x-2).^2 + (y-2).^2; % Sphere
% f = @(x,y) -20*exp(-0.2*sqrt(0.5*(x.^2 + y.^2))) - exp(0.5*(cos(2*pi*x)+cos(2*pi*y)))+20+exp(1); % Ackley
% f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2)); % Drop-Wave Function
% f = @(x,y) 10*2 + x.^2 + y.^2 -10*cos(2*pi*x) - 10*cos(2*pi*y); %Rastrigin
% f = @(x,y) ((x.^2/4000)+(y.^2/4000))-(cos(x).*cos(y/sqrt(2)))+1; % Griewank

xl = [-5; -5];
xu = [5; 5];

G = 150;
N = 50;
D = 2;

w_max = 0.8;
w_min = 0.1;
c1 = 2;
c2 = 2;

x = zeros(D,N);
v = zeros(D,N);
xb = zeros(D,N);
fitness = zeros(1,N);

f_plot = zeros(1,G);

for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    v(:,i) = 0.3*randn(D,1);

    xb(:,i) = x(:,i);

    fitness(i) = f(x(1,i), x(2,i));
end

for g=1:G
%     Plot_Contour(f,x,xl,xu); % Grafica

    for i=1:N
        fx = f(x(1,i), x(2,i));

        if fx<fitness(i)
            xb(:,i) = x(:,i);
            fitness(i) = fx;
        end
    end

    [fx_best,I_best] = min(fitness);

    w = w_max - (g/G)*(w_max-w_min);

    for i=1:N
        v(:,i) = w*v(:,i) + rand()*c1*(xb(:,i)-x(:,i)) + rand()*c2*(xb(:,I_best)-x(:,i));
        x(:,i) = x(:,i) + v(:,i);
    end

    f_plot(g) =fx_best;
end

figure 
Plot_Contour(f,x,xl,xu); % Grafica

figure
Plot_Surf(f,x(:,I_best),xl,xu) % Gráfica
disp([' mínimo global en: x=' num2str(xb(1,I_best)) ', y= ' num2str(xb(2,I_best)) ', f(x,y)=' num2str(f(xb(1,I_best),xb(2,I_best)))])

figure
plot(f_plot,'b-','LineWidth',2);
title('Gráficas de convergencia')