clear all
close all
clc


f = @(x,y) sin(x+y)+(x-y).^2-1.5*x+2.5*y+1; %McCormick Function

fp = @(x,xl,xu) f(x(1),x(2)) + 1000 * Penalty(x,xl,xu);

xl = [-1.5 -3]';
xu = [4 4]';

G = 150;
N = 50;
D = 2;

F = 0.6;
CR = 0.9;

x = zeros(D,N);
fitness = zeros(1,N);

for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    fitness(i) = fp(x(:,i),xl,xu);
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
            r = rand;

            if r<=CR
                u(j) = v(j);
            else
                u(j) = x(j,i);
            end
        end

        %seleccion
        fitness_u = fp(u,xl,xu);
        
        if fitness_u < fitness(i)
            x(:,i) = u;
            fitness(i) = fitness_u;
        end
    end
 
    fx_plot(n) = min(fitness);
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


%SEGUNDO METODO
function z = Penalty(x,xl,xu)
    z = 0;
    D = numel(xl);

    for j=1:D
        if xl(j)<x(j) 
            z = z + 0;
        else
            z = z + (xl(j)-x(j))^2;
        end

        if x(j)<xu(j) 
            z = z + 0;
        else
            z = z + (xu(j)-x(j))^2;
        end
    end
end

%PRIMER METODO
% function z = Penalty(x,xl,xu)
%     z = 0;
%     D = numel(xl);
% 
%     for j=1:D
%         if xl(j)<x(j) && x(j)<xu(j)
%             z = z + 0;
%         else
%             z = z + 1;
%         end
%     end
% end