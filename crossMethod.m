function [O1, O2] = crossMethod (P1, P2)

    N = size(P1, 2);

    O1 = zeros(1, N);
    O2 = zeros(1, N);

    low = 0;
    high = 0;
    while ( low == high )
        crosspoints = sort(randi([2, N-1], 1, 2));
        low = crosspoints(1);
        high = crosspoints(2);
    end

    
    for i = low : high
        O1(i) = P2(i);
    end

    x = [1: low-1, high+1: N];

    for x = x
        G = find(O1 == P1(x));
        if isempty(G)
            O1(x) = P1(x);
        else
            y = P1(G);
            G2 = find( O1 == y);
            while ~ isempty(G2)
                y = P1(G2);
                G2 = find (O1 == y);
            end
            O1(x) = y;
        end
    end


    for i = low : high
        O2(i) = P1(i);
    end

    x = [1: low-1, high+1: N];

    for x = x
        G = find(O2 == P2(x));
        if isempty(G)
            O2(x) = P2(x);
        else
            y = P2(G);
            G2 = find( O2 == y);
            while ~ isempty(G2)
                y = P2(G2);
                G2 = find (O2 == y);
            end
            O2(x) = y;
        end
    end

end
