
clear all
close all
clc
%f = @(x,y) x.*exp(-x.^2-y.^2); % Primera Funcion
f = @(x,y) (x-2).^2 + (y-2).^2; % Segunda Funcion

x_lim = linspace(-10,10,50); 
y_lim = linspace(-10,10,50); 

[X,Y] = meshgrid(x_lim,y_lim); 
Z = f(X,Y); 

xi = [-1 0]';

h = 0.1;

for i=1:100
%     Gx = 1.*exp(-xi(1).^2) -2.*exp(-xi(1).^2)* xi(1).^2;% Primera funcion
%     Gy = -2.*exp(-xi(2).^2)*xi(2); % Primera funcion
    Gx = 2*xi(1) - 4; % Segunda funcion
    Gy = 2*xi(2) - 4; % Segunda funcion
    G = [Gx Gy]';

    xi = xi - h*G;
end


figure
hold on
grid on

surf(X,Y,Z) 
plot3(xi(1),xi(2),f(xi(1),xi(2)),'r*','LineWidth',2,'MarkerSize',10) % plot de un punto cualqueira en 3D
legend({'función','óptimo'},'FontSize',15)

title('Gráfica en 3D','FontSize',15)
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
zlabel('f(x,y)','FontSize',15)
view([-20,60]) 

figure
hold on
grid on

contour(X,Y,Z,20) 
plot(xi(1),xi(2),'r*','LineWidth',2,'MarkerSize',10) 
legend({'función','óptimo'},'FontSize',15)

title('Gráfica en 2D','FontSize',15)
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)