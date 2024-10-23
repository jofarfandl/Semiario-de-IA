clear all
close all
clc

A = [4 2; 3 1];%valores izquierda
b = [8; 2];%valores resultados

m = numel(b);

f = @(x) (1/(2*m))*sum((b-A*x).^2);

xl = [-5 -5]';%dimensiones
xu = [5 5]';%dimensiones

G = 150;
N = 50;
D = 2;%dimensiones

F = 0.6;
CR = 0.9;

x = zeros(D,N);
fitness = zeros(1,N);

for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    fitness(i) = f(x(:,i));
end

for n=1:G

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
        fu = f(u);

        if fu<fitness(i)
            x(:,i) = u;
            fitness(i) = fu;
        end
    end
 
    [fx_plot(n),~] = min(fitness);
end

[~,igb] = min(fitness);

x(:,igb)

f(x(:,igb))

A*x(:,igb)

figure
grid on
hold on
plot(fx_plot,'LineWidth',2);
title('Gráficas de convergencia')
xlabel('iteraciones')
ylabel('f(x)')

