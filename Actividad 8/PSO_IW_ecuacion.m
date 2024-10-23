clear all
close all
clc

A = [4 2; 3 1];
b = [8; 2];

m = numel(b);
f = @(x) (1/(2*m))*sum((b-A*x).^2);

G = 300;
N = 20;
D = size(A,2);
  
xl = [-5 -5]';
xu = [5 5]';

w_max = 0.8;
w_min = 0.1;
c1 = 2;
c2 = 2;

x = zeros(D,N);
v = zeros(D,N);
xb = zeros(D,N);

fitness = zeros(1,N);

f_plot = zeros(1,G);

for  i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    v(:,i) = 0.5*randn(D,1);

    xb(:,i) = x(:,i);

    fitness(i) = f(x(:,i));
end    

for g=1:G
    
    for i=1:N
        fx = f(x(:,i));
        if fx < fitness(i)
            xb(:,i) = x(:,i);
            fitness(i) = fx;
        end
    end

    [fx_best,I_best] = min(fitness);

    f_plot(g) = fx_best;

    w = w_max - (g/G)*(w_max-w_min);

    for i=1:N
        v(:,i) = w*v(:,i) + rand()*c1*(xb(:,i)-x(:,i)) + rand()*c2*(xb(:,I_best)-x(:,i));
        x(:,i) = x(:,i) + v(:,i);
    end
end

disp("Valor de las incognitas: ")
disp(x(:,I_best))

disp("Error resultante: ")
disp(f(x(:,I_best)))

disp("Resultado obtenido: ")
disp(A*x(:,I_best))

figure
plot(f_plot);