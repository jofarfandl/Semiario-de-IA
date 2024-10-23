function Plot_Contour (f,x,xl,xu)
    cla
    hold on
    grid on

    x_lim = linspace(xl(1),xu(1),50);
    y_lim = linspace(xl(2),xu(2),50);
    [X,Y] = meshgrid(x_lim,y_lim);
    Z = f(X,Y);

    contour(X,Y,Z,20);
    plot(x(1,:),x(2,:),'*r','LineWidth',2,'MarkerSize',10)

    xlabel('x','FontSize',15)
    ylabel('y','FontSize',15)

    legend({'función','óptimo'},'FontSize',15)
    title('Gráfica en 2D','FontSize',15)
    xlabel('eje x','FontSize',15) % etiqueta en el eje x
    ylabel('eje y','FontSize',15) % etiqueta en el eje y

    axis([xl(1) xu(1) xl(2) xu(2)])
    pause(0.1)