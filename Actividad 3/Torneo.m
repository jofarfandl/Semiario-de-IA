function k = Torneo (aptitud)
    N = numel(aptitud);
    tao = round(N*0.3);

    I = randi(N, [1 tao]);
    [~,i] = max(aptitud(I));

    k = I(i);