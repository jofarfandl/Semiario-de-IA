clear all
% close all
clc

f = @(x) (20-2*x).*(20-2*x).*x;
fp = @(x) 12*x.^2-160*x+400;
fpp = @(x) 24*x-160;

x = 0:0.1:10;%%%<------


xr = 0.5;
N = 10;

for i=1:N
    xr = xr - fp(xr)/fpp(xr);
end


if fpp(xr)>=0
    disp(['minimo en x = ' num2str(xr)])
else
    disp(['maximo en x = ' num2str(xr)])
end

% figure
cla
grid on
hold on
plot(x,f(x),'b-','LineWidth',2)
plot(x,fp(x),'g--','LineWidth',2)
plot(x,fpp(x),'r-.','LineWidth',2)

plot(xr,fp(xr),'ro','LineWidth',2)
