function cost = genetic_exactly_five_fitness(in)
    
    BIG_F = 10^15;
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
    
    available = zeros(1,5);
    % Calculate cost function
    cost = 0;
    for i= 1:1000 
        c = round(in(i));
        available(c) = 1;
        cost = cost + norm(C(:,c)-x(:,i));
    end
    
    for i = 1:5
        if(available(i) == 0)
            cost = BIG_F;
        end
    end
end
