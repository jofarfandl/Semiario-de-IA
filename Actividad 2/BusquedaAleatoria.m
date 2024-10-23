clear all
close all
clc

img = imread('Image_2.bmp');
temp = imread('Template.bmp');

img_g = rgb2gray(img);
temp_g = rgb2gray(temp);

[img_H,img_W] = size(img_g);
[temp_H,temp_W] = size(temp_g);

f = @(x,y) (x-2).^2 + (y-2).^2;

xl = [1 1];
xu = [img_W-temp_W img_H-temp_H];

G=500;

fb = 9999999;
xb =[0; 0];

fx_plot = zeros(1,G);

max = -1;
xp = 0;
yp = 0;

for g=1:G
    x = xl + (xu-xl).*rand(2,1);

    fx = f(x(1),x(2));

    if fx<fb
        fb = fx;
        xb = x;
    end

    fx_plot(g) = fb;

end

x_lim = linspace(-5,5,50); % límites para eje x, -5 es inferior, 5 es superior, con 50 puntos
y_lim = linspace(-5,5,50); % límites para eje y, -5 es inferior, 5 es superior, con 50 puntos
    
[X,Y] = meshgrid(x_lim,y_lim); % creamos una rejilla de puntos (x,y) para crear el plot
Z = f(X,Y); % evaluación de cada elemento en la rejilla para crear su valor en el eje z

