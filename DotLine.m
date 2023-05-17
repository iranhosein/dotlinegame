clc, clear, close all;

%% load data
load('Sample1.mat');

%% find color for data
redPoints = find((data(:,3) == 0));
bluePoints = find(data(:,3) == 1);

%% plot data
% plot(data(redPoints,1), data(redPoints,2), 'r.', 'MarkerSize', 15);
% hold on;
% plot(data(bluePoints,1), data(bluePoints,2), 'b.', 'MarkerSize', 15);


%% Problem definition
[N, ~] = size(redPoints);
[M, ~] = size(bluePoints);

redPointsCity = data(redPoints, 1:2);
bluePointsCity = data(bluePoints, 1:2);

fullNull = zeros(N+1, M+1);
for i = 1 : N
    fullNull(redPointsCity(i, 1)+1, redPointsCity(i, 2)+1) = 1;
end
for i = 1 : M
    fullNull(bluePointsCity(i, 1)+1, bluePointsCity(i, 2)+1) = 1;
end


%% Initial Population
popSize = 256;
Population = zeros(popSize, M);
for i = 1 : popSize
    A2 = 1:M;
    sizeA2 = M;
    A3 = zeros(1, M);
        
    for j = 1 : M
        T1 = find(bluePointsCity(:, 1) == redPointsCity(j, 1));
        T2 = find(bluePointsCity(:, 2) == redPointsCity(j, 2));
        T = [T1; T2];
        
        [sizeT, ~] = size(T);
                
        TT = zeros(1, sizeT);
        k = 0;
        for jj = 1 : sizeT
            ll = find ( A2 == T(jj));
            if ~isempty(ll)
                TT(k+1) = T(jj);
                k = k + 1;
            end
        end
        TT = TT(1:k);
        
        F = 0;
        if k == 0
            F = randi(sizeA2);
            A3(j) = A2(F);
        else
            for t = 1: k
                F = t;
                numeric = A2(F);
                if bluePointsCity(numeric, 1) == redPointsCity(j, 1)
                    numx = redPointsCity(j, 1);
                    numMin = min(bluePointsCity(numeric, 2), redPointsCity(j, 2));
                    numMax = max(bluePointsCity(numeric, 2), redPointsCity(j, 2));
                    numMin = numMin + 1;
                    while numMin ~= numMax
                        if fullNull(numx+1, numMin+1) == 1
                            break;
                        end
                        numMin = numMin + 1;
                    end
                    if numMin == numMax
                        break;
                    end
                else
                    numy = redPointsCity(j, 2);
                    numMin = min(bluePointsCity(numeric, 1), redPointsCity(j, 1));
                    numMax = max(bluePointsCity(numeric, 1), redPointsCity(j, 1));
                    numMin = numMin + 1;
                    while numMin ~= numMax
                        if fullNull(numMin+1, numy+1) == 1
                            break;
                        end
                        numMin = numMin + 1;
                    end
                    if numMin == numMax
                        break;
                    end
                end
            end
            A3(j) = A2(F);
        end
        A2([F, sizeA2]) = A2([sizeA2, F]);
        sizeA2 = sizeA2 - 1;
        A2 = A2(1: sizeA2);
    end
    
    Population(i, :) = A3;
end




NoIteration = 1000;
popCost = zeros(1, popSize);
history = zeros(1, NoIteration);
globalMin = Inf;
popMutation = zeros(popSize/2, M);
popCrossover = zeros(popSize/2, M);

% Main Loop
for i = 1:NoIteration
    i
    for j= 1: popSize
        Fn = fullNull;
        cost = 0;
        for jj = 1: M
            A = redPointsCity(jj, 1:2);
            B = bluePointsCity(Population(j, jj), 1:2);
            if A(1) ~= B(1) && A(2) ~= B(2)
                cost = cost + 1;
            elseif A(1) == B(1)
                Nx = A(1);
                Nmin = min(A(2), B(2));
                Nmax = max(A(2), B(2));
                Nmin = Nmin + 1;
                
                Fn1 = Fn;
                while Nmin ~= Nmax
                    if Fn(Nx+1, Nmin+1) == 1
                        cost = cost + 1;
                        break;
                    else
                        Fn1(Nx+1, Nmin+1) = 1;
                    end
                    Nmin = Nmin + 1;
                end
                if Nmin == Nmax
                    Fn = Fn1;
                end
            elseif A(2) == B(2)
                Ny = A(2);
                Nmin = min(A(1), B(1));
                Nmax = max(A(1), B(1));
                Nmin = Nmin + 1;
                
                Fn2 = Fn;
                while Nmin ~= Nmax
                    if Fn(Nmin+1, Ny+1) == 1
                        cost = cost + 1;
                        break;
                    else
                        Fn2(Nmin+1, Ny+1) = 1;
                    end
                    Nmin = Nmin + 1;
                end
                if Nmin == Nmax
                    Fn = Fn2;
                end
            end
        end
        popCost (j) = cost;
    end
    [minCost, minIndex] = min(popCost);
    history(i) = minCost;
    
    if minCost < globalMin
        globalMin = minCost;
        bestRoute = Population (minIndex, :);
       
%         for iii = 1 : M
%             A = zeros(1, 2);
%             B = zeros(1, 2);
%             A(1) = redPointsCity(iii, 1);
%             A(2) = bluePointsCity(Population(minIndex, iii), 1);
%             B(1) = redPointsCity(iii, 2);
%             B(2) = bluePointsCity(Population(minIndex, iii), 2);
%             
%             if A(1) == A(2) || B(1) == B(2)
%                 plot(A, B, 'k-', 'MarkerSize', 10);
%             end
% 
%         end         
%         drawnow;
        
    end
    
    
    %% GA operators
    %CrossOver
    
    for ii = 1:2:popSize/2
        r1 = randi(popSize);
        r2 = randi(popSize);
        
        [popCrossover(ii, :), popCrossover(ii+1, :)] = ...
                    crossMethod(Population(r1, :), Population(r2, :));
    end
    
    
    
    %% Mutation
    [~, index] = sort(popCost);
    sortedPopulation = Population(index, :);
   
    best = rouletteWheel(sortedPopulation, popSize, M);
    
    k = 1;
    
    for ii = 1:size(best, 1)
        crosspoints = sort(randi(M, 1, 2));
        low = crosspoints(1);
        high = crosspoints(2);
        
        popMutation(k, :) = best(ii, :);
        
        popMutation(k+1, :) = best(ii, :);
        popMutation(k+1, [low high]) = popMutation(k+1, [high low]);
        
        popMutation(k+2, :) = best(ii, :);  
        popMutation(k+2, [low:high]) = popMutation(k+2, [high:-1:low]);
        
        popMutation(k+3, :) = best(ii, :);
        popMutation(k+3, [low:high]) = popMutation(k+3, [low+1 : high low]);
        
        k = k + 4;
        
    end
    
    Population = [popMutation; popCrossover;];
    
end



%% Plotting
figure
subplot(2, 2, 1)
    plot(data(redPoints,1), data(redPoints,2), 'r.', 'MarkerSize', 15);
    hold on;
    plot(data(bluePoints,1), data(bluePoints,2), 'b.', 'MarkerSize', 15);

subplot(2, 2, 2)
    plot(data(redPoints,1), data(redPoints,2), 'r.', 'MarkerSize', 15);
    hold on;
    plot(data(bluePoints,1), data(bluePoints,2), 'b.', 'MarkerSize', 15);
    for iii = 1 : M
            A = zeros(1, 2);
            B = zeros(1, 2);
            A(1) = redPointsCity(iii, 1);
            A(2) = bluePointsCity(Population(minIndex, iii), 1);
            B(1) = redPointsCity(iii, 2);
            B(2) = bluePointsCity(Population(minIndex, iii), 2);
            
            if A(1) == A(2) || B(1) == B(2)
                plot(A, B, 'k-', 'MarkerSize', 10);
            end

     end         

subplot(2, 2, 3:4)
    plot(1:NoIteration, history);





