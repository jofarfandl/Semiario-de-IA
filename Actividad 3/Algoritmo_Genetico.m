clear all
close all
clc

% f = @(x,y) (x-2).^2+(y-2).^2;
% xl = [-5 -5]';
% xu = [5 5]';

f = @(x,y) x.*exp(-x.^2-y.^2);
xl = [-2 -2]';
xu = [2 2]';

D = 2;
G = 100;
N = 100;
pm = 0.1;

x = zeros(D, N);
fit = zeros(1, N);
aptitud = zeros(1, N);

for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
end    

f_plot = zeros(1,G);

for g=1:G
    Plot_Contour(f, x, xl, xu);

    for i=1:N
        fit(i) = f(x(1,i),x(2,i));
        if fit(i) >= 0
            aptitud(i) = 1/(1+fit(i));
        else
            aptitud(i) = 1 + abs(fit(i));
        end
    end

    f_plot(g) = min(fit);

    y = zeros(D,N);

    for i=1:2:N
        % Seleccion
        r1 = Torneo(aptitud);
        r2 = r1;

        while r1==r2
            r2 = Torneo(aptitud);
        end

        x1 = x(:,r1);
        x2 = x(:,r2);

        % Cruce
%         pc = randi([1 D]);
%         y1 = x1;
%         y2 = x2;
% 
%         y1(pc:D) = x2(pc:D);
%         y2(pc:D) = x1(pc:D);
        
        r = rand();

        y1 = r*x1 + (1-r)*x2;
        y2 = (1-r)*x1 + r*x2;

        y(:,i) = y1;
        y(:,i+1) = y2;
    end
    % Mutacion
    for i=1:N
        for j=1:D
            if rand() <pm
%                 y(j,i) = xl(j) + (xu(j)-xl(j))*rand();
                y(j,i) = y(j,i) + normrnd(0,1);
            end
        end
    end

    x = y;
end

for i=1:N
    fit(i) = f(x(1,i),x(2,i));
    if fit(i) >= 0
        aptitud(i) = 1/(1+fit(i));
    else
        aptitud(i) = 1 + abs(fit(i));
    end
end

[~, i_mejor] = max(aptitud);

figure
Plot_Surf(f,x,xl,xu);

disp(['Minimo global en: x=' num2str(x(1,i_mejor)) ', y=' num2str(x(2,i_mejor)) ', f(x, y)=' num2str(f(x(1,i_mejor), x(2,i_mejor)))]);

figure
plot(f_plot);