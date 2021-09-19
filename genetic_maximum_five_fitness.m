function cost = genetic_maximum_five_fitness(in)
    
    % Load data
    x = load('DataNew.mat');
    x = x.DataNew;
    
    % Number of nodes in each cluster
    N = zeros(1, 5);
    % Centers
    C = zeros(2, 5);
    
    % Find centers
    c = 0;
    for i= 1:1000 
        c = round(in(i));
        C(:,c) = C(:,c) + x(:,i);
        N(1, c) = N(1, c) + 1;
    end
    
    C = C ./ N;
    
    % Calculate cost function
    cost = 0;
    for i= 1:1000 
        c = round(in(i));
        cost = cost + norm(C(:,c)-x(:,i));
    end
    
end

