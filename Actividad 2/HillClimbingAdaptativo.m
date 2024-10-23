clear all
close all
clc

%f = @(x,y) x.*exp(-x.^2-y.^2);  
f = @(x,y) (x-2).^2 + (y-2).^2;

xl = [-5 -5]';
xu = [5 5]';

pm = 0.5;

D = 2;
G=100;%Ahora representara generaciones


x = xl + (xu-xl).*rand(D,1);

fx_plot = zeros(1,G);
for g=1:G
    
    y = x;

    for j=1:D
        r = rand();
        if r<pm
             y(j) = xl(j) + (xu(j)-xl(j)).*r; 
        end
    end

    fy = f(y(1),y(2));
    fx = f(x(1),x(2));
    if fy<fx
        x = y;
    end
    fb = f(x(1),x(2));
    fx_plot(g) = fb;

end

xb = x;
fb = f(x(1),x(2));

disp(['minimo global en: x=' num2str(xb(1))]);

x_lim = linspace(-5,5,50); % límites para eje x, -5 es inferior, 5 es superior, con 50 puntos
y_lim = linspace(-5,5,50); % límites para eje y, -5 es inferior, 5 es superior, con 50 puntos
    
[X,Y] = meshgrid(x_lim,y_lim); % creamos una rejilla de puntos (x,y) para crear el plot
Z = f(X,Y); % evaluación de cada elemento en la rejilla para crear su valor en el eje z

figure
hold on
grid on

contour(X,Y,Z,20) 
plot(xb(1),xb(2),'r*','LineWidth',2,'MarkerSize',10) %mostramos la linea con el mejor acercamiento
legend({'función','óptimo'},'FontSize',15)

figure
hold on
grid on

surf(X,Y,Z) 
plot3(xb(1),xb(2),f(xb(1),xb(2)),'r*','LineWidth',2,'MarkerSize',10) % plot de un punto cualqueira en 3D
legend({'función','óptimo'},'FontSize',15)

view([-20,60]) 

figure
grid on
hold on
plot(fx_plot,'b-','LineWidth',2)
title('Grafica de convergencia')
xlabel('iteracion')
ylabel('f(x)')