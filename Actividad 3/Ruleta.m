function k = Ruleta (aptitud)
    p = aptitud/sum(aptitud);
    r = rand();

    p_sum = 0;
    N = numel(aptitud);
    k = N;

    for i=1:N
        p_sum = p_sum + p(i);

        if p_sum >= r
            k = i;
        return
        end
    end