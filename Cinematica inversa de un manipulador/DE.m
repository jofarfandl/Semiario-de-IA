clear all
close all
clc

a_1 = 0.35;
a_2 = 0.35;
a_3 = 0.25;

pd = [0.5; 0.5];

beta = 1000;
xl = [-160 -130 -100]'*(pi/180);
xu = [160 130 100]'*(pi/180);

p = @(theta_1,theta_2,theta_3) [a_2*cos(theta_1 + theta_2) + a_1*cos(theta_1) + a_3*cos(theta_1+theta_2+theta_3); ...
                                a_2*sin(theta_1 + theta_2) + a_1*sin(theta_1) + a_3*sin(theta_1+theta_2+theta_3)];

f = @(q) (1/4)*sum((pd-p(q(1),q(2),q(3))).^2);

fp = @(x,xl,xu) f(x) + beta*Penalty(x,xl,xu);

G = 1000;
N = 50;
D = 3;

F = 0.6;
CR = 0.9;

x = zeros(D,N);
fitness = zeros(1,N);

for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);

    q = x(:,i);

    fitness(i) = f(q);
end    

for n=1:G
    for i=1:N
        %Mutacion
        I = randperm(N,4);
        I(I == i) = [];
        r1 = I(1);
        r2 = I(2);
        r3 = I(3);

        v = x(:,r1) + F*(x(:,r2)-x(:,r3));

        %Recombinacion
        u = zeros(D,1);
        k = randi([1,D]);

        for j=1:D     
            if rand()<=CR || j == k
                u(j) = v(j);
            else
                u(j) = x(j,i);
            end
        end

        %Seleccion
        q = u;
        fu = fp(u,xl,xu);

        if fu < fitness(i)
            x(:,i) = u;
            fitness(i) = fu;
        end
    end
    [fx_best,I_best] = min(fitness);
end    

cla
hold on
grid on
Dibujar_Manipulador(q)
plot(pd(1),pd(2),'go','LineWidth',2,'MarkerSize',10)
pd

%Fuciones
function z = Penalty (x, xl, xu)
    z = 0;
    m = numel(xl);

    for i=1:m
        if xl(i) < x(i)
            z = z + 0;
        else
            z = z + (x(i)-xl(i))^2;
        end

        if x(i) < xu(i)
            z = z + 0;
        else
            z = z + (x(i)-xu(i))^2;
        end
    end
end