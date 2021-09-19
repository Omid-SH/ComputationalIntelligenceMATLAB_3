function cost = genetic_maximum_five_fitness_1(in)
    
    % Load data
    x = load('DataNew.mat');
    x = x.DataNew;
    
    % Calculate cost function
    cost = 0;
    minn = 0;
    for i = 1:1000 
        minn = norm([in(1), in(2)] - [x(1,i), x(2,i)]);
        for j = 2:5
            minn = min(minn, norm([in(2*j-1), in(2*j)] - [x(1,i), x(2,i)]));
        end
        cost = cost + minn;
    end
    
end

