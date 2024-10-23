clear all
close all
clc





f = @(x,y) (x-2).^2 + (y-2).^2;
xl = [-5 -5]';
xu = [5 5]';


D = 2;
G = 100;%Ahora representara generaciones
N = 10;

x = zeros(D,N);
fitness = zeros(1,N);
aptitud = zeros(1,N);

for i=1:N
    x(:,i) = xl + (xu-xl).*rand(D,1);%a todos los elementos de la primer columna asignarles el valor

    fitness(i) = f(x(1,i),x(2,i));

    if fitness(i) >=0
        aptitud(i) = 1/(1+fitness(i));
    else
        aptitud(i) = 1 + abs(fitness(i));
    end
end

Plot_Contour (f,x,xl,xu);



aptitud = [10 100 5 10 6 50];

p1 = Ruleta(aptitud);
p2 = p1;

while p1 == p2
    p2 = Ruleta(aptitud);
end

x1 = x(:,p1);
x2 = x(:,p2);

%cruce

y1 = x1;
y2 = x2;

pc = randi([1 D]);

y1(pc:D) = x1(pc:D);
y2(pc:D) = x2(pc:D);

return








for g=1:G
    Plot_Contour (f,x,xl,xu);

end