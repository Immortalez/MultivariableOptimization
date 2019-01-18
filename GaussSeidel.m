clear
clc
format long

% ############## Providing necessary information ##############
% Go to the very bottom of the script.


% ## Prints additional info with each iteration (1 ON, 0 OFF)
verbose = 1;

% ## Set the initial point, bounds for x and y, directions & precision
x0 = [-1, -1];
x_range = [-1, 1];
y_range = [-1, 1];
e1 = [1 0]; % Moving along X axis -- in odd iterations
e2 = [0 1]; % Moving along Y axis -- in even iterations
prec = 0.01;

% ##################### Informational box #####################
fprintf('\n/////////////////////////////////////////////////////////////////////');
fprintf('\n//                                                                 //');
fprintf('\n//    OaDM Laboratory                                17.01.2019    //');
fprintf('\n//                                                                 //');
fprintf('\n//                   Multivariable Optimization                    //');
fprintf('\n//                      Gauss-Seidel method                        //');
fprintf('\n//                                                                 //');
fprintf('\n/////////////////////////////////////////////////////////////////////\n\n');
fprintf('  Initial point x0: (%d, %d)\n\n', x0(1), x0(2));

% #################### Gauss Seidel Method ####################
x = x0; % Stores current point coordinates (x, y)
e = []; % Stores current direction of moving
ax = -1; % Stores calculated alpha for X (displacement along X axis)
ay = -1; % Stores calculated alpha for Y (displacement along Y axis)
a_range = []; % Stores calculated alpha interval
xn = []; % Stores the calculations for next x
syms alpha; % Symbolic representation of alpha
syms minEq; % Symbolic representation of obtained equation to be minimized
iterationNumber = 1;
amountOfIterations = iterationNumber - 1;

while(ax ~= 0 || ay ~= 0)
    goingX = 1; % Going along X axis - 1; Along Y axis - 0
    % ## Choosing the proper direction
    if(mod(iterationNumber, 2) == 1)
        % Odd iteration -- moving along X axis
        e = e1;
        goingX = 1;
    else
        % Even iteration -- moving along Y axis
        e = e2;
        goingX = 0;
    end
    
    % ## Calculating next x
    xn = x + alpha*e; % Works OK
    minEq = expand(F(xn)); % Works OK --> function to be minimized in interval
    
    % ## Calculating the acceptable interval for alpha
    if(goingX == 1)
        a_range(1) = x_range(1) - x(1);
        a_range(2) = x_range(2) - x(1);
    else
        a_range(1) = x_range(1) - x(2);
        a_range(2) = x_range(2) - x(2);
    end
    
    % ## Calculating the alpha and its value (minimizing minEq)
    a_domain = a_range(1):prec:a_range(2);
    a_values = double(subs(minEq, a_domain));
    minA_val = min(a_values);
    minA = a_domain(find(minA_val == a_values));
    
    if(numel(minA) > 1)
        % If it happens that there are two equal minA_val, take the first
        minA = minA(1);
    end
    
    % ## Assigning the displacement and calculating new x
    if(goingX == 1)
        ax = minA;
        x = double(subs(xn, ax));
    else
        ay = minA;
        x = double(subs(xn, ay));
    end
    
    % ## Printing additional information if 'verbose' enabled
    if(verbose == 1)
        fprintf('  Iteration [%d]\n', iterationNumber)
        if(goingX == 1)
            fprintf('      a%d = %g\n', iterationNumber, ax);
            fprintf('      Moving along X axis by [%g].\n', ax);
        else
            fprintf('      a%d = %g\n', iterationNumber, ay);
            fprintf('      Moving along Y axis by [%g].\n', ay);
        end
        fprintf('      x%d = (%g, %g)\n', iterationNumber, x(1), x(2));
        fprintf('      F(x%d) = %g\n\n', iterationNumber, F(x));
    end
    % ## Incrementing the total iterations amount and number
    iterationNumber = iterationNumber + 1;
    amountOfIterations = amountOfIterations + 1;
end

% ## Summary of what we obtained as result
fprintf('\n##############################################################\n\n');
fprintf('    Displacement in both directions is equal to 0. Algorithm stops here.\n');
fprintf('      Minimum found at (%g, %g) with value: %g\n\n\n', x(1), x(2), F(x));


% ############## Providing necessary information ##############
% Provide the formula of the function
function value = F(args)
x = args(1);
y = args(2);
% Functions from the lab materials (for all x,y E [-1, 1]):
value = (x-0.5)^2 + (y-0.5)^2;
%value = (sin(x-0.5))^2 + (atan(y-0.5))^2;
%value = (sin(x-0.5))^2 + (sin(y-0.5))^2;
%value = acot(x-0.5)^2 + (atan(y-0.5))^2;
%value = exp(cos(x-0.5)) + exp(cos(y-0.5));
%value = (y*(1-x^2))/(1.2 - sin(2*y));
%value = (2-3*x^4)*tan(sin(2*y+1));
%value = (2-y^2)*exp(-sin(2*x-1));
end
