function cost = genetic_exactly_five_fitness_1(in)
   
    % Load data
    x = load('DataNew.mat');
    x = x.DataNew;
    
    % Check if 5 clusters available
    available = ones(1,5);
    
    % Calculate cost function
    cost = 0;
    minn = 0;
    minnj = 1;
    for i = 1:1000 
        minn = norm([in(1), in(2)] - [x(1,i), x(2,i)]);
        for j = 2:5
            if(norm([in(2*j-1), in(2*j)] - [x(1,i), x(2,i)]) < minn)
                minn = min(minn, norm([in(2*j-1), in(2*j)] - [x(1,i), x(2,i)]));
                minnj = j;
            end
        end
        cost = cost + minn;
        available(minnj) = 0;
    end
    
    for i = 1:5
        if(available(i))
            cost = 200000;
        end
    end
    
end