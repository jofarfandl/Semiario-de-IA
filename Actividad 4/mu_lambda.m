clear all
close all
clc

f = @(x,y) (x-2).^2+(y-2).^2;
xl = [-5 -5]';
xu = [5 5]';

% f = @(x,y) x.*exp(-x.^2-y.^2);
% xl = [-2 -2]';
% xu = [2 2]';

G = 500; %iteraciones
D = 2; %dimension del problema
mu = 30;%padres
lambda = 20;

x = zeros(D,mu+lambda);%inicializa valores de las soluiones
sigma = zeros(D,mu+lambda);%valores con la misma dimension de los padres
fitnes = zeros(1,mu+lambda);

for i=1:mu
    x(:,i) = xl+(xu-xl).*rand(D,1);
    sigma(:,i) =0.1*rand(D,1);
end

f_plot = zeros(1,G);

for t=1:G
    for k=1:lambda
        %seleccion
        r1 = randi([1 mu]);
        r2 = r1;
    
        while r1==r2
            r2 = randi([1 mu]);
        end
    
        %recombinacion    
        for j=1:D
            if randi([0 1])
                x(j,mu+k) = x(j,r1);
                sigma(j,mu+k) = sigma(j,r1);
            else
                x(j,mu+k) = x(j,r2);
                sigma(j,mu+k) = sigma(j,r2);
            end
        end
    
        %mutacion
        r = normrnd(0,sigma(:,mu+k));
        x(:,mu+k) = x(:,mu+k) + r;

    end

    %el peor se elimina
    for i=1:mu+lambda
        fitnes(i) = f(x(1,i),x(2,i));
    end

    [~,I] = sort(fitnes);

    x = x(:,I);
    sigma = sigma(:,I);
    fitnes = fitnes(I);
    
    f_plot(t) = f(x(1,1),x(2,1));
    %fx = f(x(1,1),x(2,1));
    %

end

xb = x(:,1);

figure 
Plot_Contour(f,x(:,1:mu), xl, xu);

figure
Plot_Surf(f, xb, xl, xu);

disp(['Minimo global en: x=' num2str(xb(1)) ', y=' num2str(xb(2)) ', f(x, y)=' num2str(f(xb(1), xb(2)))]);

figure
plot(f_plot);