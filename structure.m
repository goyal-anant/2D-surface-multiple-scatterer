%Function to define the scattering object
%it takes the discretization level 'N', lambda, radius of the disc,
%plot_flag as input arguments and returns the coordinates of the disc(or
%any structure).

function [X, Y] = structure(N,lambda,radius1,radius2,plot_flag)
    %generating the square space
    x = linspace(-2*lambda,2*lambda,N);
    y = x;
    
    [X1, Y1] = meshgrid(x,y);
    X1       = reshape(X1,[N^2,1]);
    Y1       = reshape(Y1,[N^2,1]);
    
    %for selecting out disk in the space
    t = 0;
    for i = 1:length(X1)
        if ((X1(i) + lambda)^2 + Y1(i)^2 <= radius1^2) || ((X1(i) - lambda)^2 + Y1(i)^2 <= radius2^2)
            t    = t + 1;
            X(t) = X1(i);
            Y(t) = Y1(i);
        end
    end
    
    if plot_flag == 1
        scatter(X1,Y1,'black','filled');
        hold on; grid on; axis('equal','tight');
        scatter(X,Y,'red','filled');
    end
end
