clear all
close all
clc

f = @(x,y) (x-2).^2 + (y-2).^2; % Sphere
% f = @(x,y) -20*exp(-0.2*sqrt(0.5*(x.^2 + y.^2))) - exp(0.5*(cos(2*pi*x)+cos(2*pi*y)))+20+exp(1); % Ackley
% f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2)); % Drop-Wave Function
% f = @(x,y) 10*2 + x.^2 + y.^2 -10*cos(2*pi*x) - 10*cos(2*pi*y); %Rastrigin
% f = @(x,y) ((x.^2/4000)+(y.^2/4000))-(cos(x).*cos(y/sqrt(2)))+1; % Griewank

xl = [-5; -5]; %limite inferior 
xu = [5; 5]; %limite superior

G = 150; %iteraciones
N = 50; %poblacion
D = 2;  %dimension

L = 35; %numero de intentos maximos
Pf = 30; %abejas empleadas
Po = N-Pf; %abejas observadoras
l = zeros(1,Pf);%limite de intentos por pf/c

x = zeros(D,Pf); %soluciones
fitness = zeros(1,Pf); %variable auxiliar para la evaluacion directa de la funcion objetivo
aptitud = zeros(1,Pf); %calidad de cada fuente de alimento

f_plot = zeros(1,G);

for i=1:Pf
    x(:,i) = xl+(xu-xl).*rand(D,1);
    fitness(i) = f(x(1,i), x(2,i));
    
    fx = fitness(i);
    if fx >= 0
        aptitud(i) = 1/(1+fx);
    else
        aptitud(i) = 1+abs(fx);
    end

end

for g=1:G
    %Plot_Contour(f,x,xl,xu); % Grafica

    %abejas empleadas
    for i=1:Pf
        k=i;

        while k == i
            k = randi([1 Pf]);
        end

        j = randi([1 D]);
        phi = 2*rand()-1;

        v = x(:,i); %calculo de solucion candidata
        v(j) = x(j,i) + phi*(x(j,i)-x(j,k));

        fv = f(v(1),v(2));%evaluar la nueva solucion

        if fv < fitness(i)
            x(:,i) = v;
            fitness(i) = fv;
            l(i) = 0;
        else
            l(i) = l(i) + 1;
        end

        fx = fitness(i);
        if fx >= 0
            aptitud(i) = 1/(1+fx);
        else
            aptitud(i) = 1+abs(fx);
        end

    end

    %Abejas observadoras
    for i=1:Po
        m = Seleccion(aptitud);
        k=m;

        while k == m
            k = randi([1 Pf]);
        end

        j = randi([1 D]);
        phi = 2*rand()-1;

        v = x(:,m); %calculo de solucion candidata
        v(j) = x(j,m) + phi*(x(j,m)-x(j,k));

        fv = f(v(1),v(2));%evaluar la nueva solucion

        if fv < fitness(m)
            x(:,m) = v;
            fitness(m) = fv;
            l(m) = 0;
        else
            l(m) = l(m) + 1;
        end

        fx = fitness(m);
        if fx >= 0
            aptitud(m) = 1/(1+fx);
        else
            aptitud(m) = 1+abs(fx);
        end

    end

    %Abejas exploradoras
    for i=1:Pf
        if l(i)>L
            x(:,i) = xl+(xu-xl).*rand(D,1);
            fitness(i) = f(x(1,i), x(2,i));
            
            fx = fitness(i);
            if fx >= 0
                aptitud(i) = 1/(1+fx);
            else
                aptitud(i) = 1+abs(fx);
            end

            l(i) = 0;
            
        end
    end
    %Validar la mejor solucion por generacion...
    [fx_best,I_best] = min(fitness);
    f_plot(g) = fx_best;
end

[~,igb] = min(fitness);

figure
Plot_Surf(f,x,xl,xu) % Gráfica
disp([' mínimo global en: x=' num2str(x(1,igb)) ', y= ' num2str(x(2,igb)) ', f(x,y)=' num2str(fitness(igb))])

figure
plot(f_plot);

figure
Plot_Contour(f,x,xl,xu); % Grafica

%% Funciones
function [n] = Seleccion (aptitud)
    aptitud_total = sum(aptitud);
    N = numel(aptitud);
    
    r = rand;
    P_sum = 0;

    for i=1:N
        P_sum = P_sum + aptitud(i)/aptitud_total;

        if P_sum >= r
            n = i;
        return
        end
    end
end
