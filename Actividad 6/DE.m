clear all
close all
clc

f = @(x,y) sin(x+y)+(x-y).^2-1.5*x+2.5*y+1;

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
        r1 = i;
        while r1==i
            r1 = randi([1,N]);
        end

        r2 = r1;
        while r2==r1 || r2==i
            r2 = randi([1,N]);
        end

        r3 = r2;
        while r3==r2 || r3==r1 || r3==i
            r3 = randi([1,N]);
        end

        v = x(:,r1) + F*(x(:,r2)-x(:,r3));

        %Recombinacion
        u = zeros(D,1);

        for j=1:D
            if rand()<=CR
                u(j) = v(j);
            else
                u(j) = x(j,i);
            end
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

