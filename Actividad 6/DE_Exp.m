clear all
close all
clc

f = @(x,y) (x-2).^2 + (y-2).^2; % Sphere
% f = @(x,y) -20*exp(-0.2*sqrt(0.5*(x.^2 + y.^2))) - exp(0.5*(cos(2*pi*x)+cos(2*pi*y)))+20+exp(1); % Ackley
% f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2)); % Drop-Wave Function
%f = @(x,y) 10*2 + x.^2 + y.^2 -10*cos(2*pi*x) - 10*cos(2*pi*y); %Rastrigin
% f = @(x,y) ((x.^2/4000)+(y.^2/4000))-(cos(x).*cos(y/sqrt(2)))+1; % Griewank

xl = [-5; -5];
xu = [5; 5];

G = 150;
N = 50;
D = 2;

F = 0.6;
CR = 0.9;

x = zeros(D,N);
fitness = zeros(1,N);

for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    fitness(i) = f(x(1,i), x(2,i));
end

fx_plot = zeros(1,G);

for n=1:G
     %Plot_Contour(f,x,xl,xu); % Grafica

    for i=1:N
        %Mutacion
        I = randperm(N);
        I(I == i) = [];

        r1 = I(1);
        r2 = I(2);
        r3 = I(3);
        r4 = I(4);
        r5 = I(5);

         %DE/rand/1/exp
         v = x(:,r1) + F*(x(:,r2)-x(:,r3));

%         %DE/rand/2/exp
%         v = x(:,r1) + F*(x(:,r2)-x(:,r3)) + F*(x(:,r4)-x(:,r5));

%         %DE/best/1/exp
%         [~,best] = min(fitness);
%         v = x(:,best) + F*(x(:,r1)-x(:,r2));

        %Recombinacion
        u = x(:,i); %vector de prueba
        j = randi([1,D]); %con randi entre {1,D}
        L = 1;

        while rand()<=CR && L<=D
            u(j) = v(j);
            j = 1 + mod(j,D);
            L = L + 1;
        end

        %seleccion
        fu = f(u(1),u(2));

        if fu<fitness(i)
            x(:,i) = u;
            fitness(i) = fu;
        end
    end
 
    [fx_plot(n),~] = min(fitness);
end

[~,igb] = min(fitness);

figure 
Plot_Surf(f,x,xl,xu) % Gráfica
disp([' mínimo global en: x=' num2str(x(1,igb)) ', y= ' num2str(x(2,igb)) ', f(x,y)=' num2str(fitness(igb))])

figure
grid on
hold on
plot(fx_plot,'LineWidth',2);
title('Gráficas de convergencia')
xlabel('iteraciones')
ylabel('f(x)')


figure 
Plot_Contour(f,x,xl,xu); % Grafica

