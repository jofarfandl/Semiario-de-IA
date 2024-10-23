clear all
close all
clc

%% 
img_ref = imread("ref_3.png");
[~,~,P] = readBarcode(img_ref,"QR-CODE");
[N,M,~] = size(img_ref);

img_des = imread('des.png');
[n,m,~] = size(img_des);

%%
X1 = P(1,:)';
X2 = P(2,:)';
X3 = P(3,:)';

x1 = [1 n]';
x2 = [1 1]';
x3 = [m 1]';

xl = [1; 1; -pi; 0];
xu = [N; M; pi; 10];

G= 150;
N = 50;
D = 4;

F = 0.6;
CR = 0.9;

x = zeros(D,N);
fitness = zeros(1,N);

for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    
    q = x(:,i);

    xp1 = Transformacion_Similitud(q,x1);
    xp2 = Transformacion_Similitud(q,x2);
    xp3 = Transformacion_Similitud(q,x3);
    
    e1 = Distancia_Euclidiana(X1,xp1);
    e2 = Distancia_Euclidiana(X2,xp2);
    e3 = Distancia_Euclidiana(X3,xp3);
    
    f = (1/6)*(e1^2+e2^2+e3^2);

    fitness(i) = f;
end

fx_plot = zeros(1,G);

for t=1:G  
    for i=1:N
        % Mutacion
        I = randperm(N);
        I(I == i) = [];

        r1 = I(1);
        r2 = I(2);
        r3 = I(3);
        r4 = I(4);
        r5 = I(5);

        % DE/rand/2/bin
        v = x(:,r1) + F*(x(:,r2)-x(:,r3)) + F*(x(:,r4)-x(:,r5));
        
        % Recombinacion
        u = x(:,i); % vector de prueba
        j = randi([1,D]); % con randi entre {1,D}
        L = 1;
        
        while rand()<=CR && L<=D
            u(j) = v(j);
            j = 1 + mod(j,D);
            L = L + 1;
        end
        
        q = u;

        xp1 = Transformacion_Similitud(q,x1);
        xp2 = Transformacion_Similitud(q,x2);
        xp3 = Transformacion_Similitud(q,x3);
        
        e1 = Distancia_Euclidiana(X1,xp1);
        e2 = Distancia_Euclidiana(X2,xp2);
        e3 = Distancia_Euclidiana(X3,xp3);
        
        f = (1/6)*(e1^2+e2^2+e3^2);
    
        fu = f + 10000*Penalty(u,xl,xu);

        if fu<fitness(i)
            x(:,i) = u;
            fitness(i) = fu;
        end
    end

    [fx_plot(t),~] = min(fitness);
end

[~,igb] = min(fitness);

fitness(igb)
q= x(:,igb)

Imprimir_Imagenes(q,img_des,img_ref)

%% Funciones
function xp = Transformacion_Similitud (qi,xi)
    dx = qi(1);
    dy = qi(2);
    theta = qi(3);
    s = qi(4);
    
    xp = [s*cos(theta) -s*sin(theta); s*sin(theta) s*cos(theta)]*xi + [dx dy]';
end

function e = Distancia_Euclidiana (X,x)
    e = sqrt((X(1)-x(1))^2+(X(2)-x(2))^2);
end

function Imprimir_Imagenes (q,img_des,img_ref)
    dx = q(1);
    dy = q(2);
    theta = q(3);
    s = q(4);
    
    T = [s*cos(theta) -s*sin(theta) dx; s*sin(theta) s*cos(theta) dy; 0 0 1];
    Tp = projective2d(T');
    
    [N,M,~] = size(img_ref);
    [n,m,~] = size(img_des);

    panoramaView = imref2d([N M]);
    Iwarp = imwarp(img_des,Tp,'OutputView',panoramaView);
    Imask = imwarp(true(n,m),Tp,'OutputView',panoramaView);
    
    blender = vision.AlphaBlender('Operation','Binary mask','MaskSource','Input port');
    panorama = step(blender,img_ref,Iwarp,Imask);
    
    imshow(panorama)
end

function z = Penalty (x,xl,xu)
    z = 0;
    m = numel(xl);

    for i=1:m
        if xl(i)<x(i)
           z = z + 0;
        else
            z = z + (x(i)-xl(i))^2;
        end

        if x(i)<xu(i)
            z = z + 0;
        else
            z = z + (x(i)-xu(i))^2;
        end
    end
end