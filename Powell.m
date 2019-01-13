clear
clc
format long

% OaDM Laboratory
% MULTIVARIABLE OPTIMIZATION
% 17.01.2019

% Section:
%   Tomasz GOROL
%   Magdalena LOREK
%   Sebastian CHLODEK

% ############## Providing necessary information ##############
% Go to the very bottom of the script.


% ## Prints additional info with each iteration (1 ON, 0 OFF)
verbose = 1;

% ## Set the initial point, bounds for x and y, directions & precision
x0 = [-1, -1];
x_range = [-1, 1];
y_range = [-1, 1];
d1 = [1 0]; % Moving along X axis -- in odd iterations
d2 = [0 1]; % Moving along Y axis -- in even iterations
prec = 0.01;

% ##################### Informational box #####################
fprintf('\n/////////////////////////////////////////////////////////////////////');
fprintf('\n//                                                                 //');
fprintf('\n//    OaDM Laboratory                                17.01.2019    //');
fprintf('\n//                                                                 //');
fprintf('\n//                   Multivariable Optimization                    //');
fprintf('\n//                          Powell method                          //');
fprintf('\n//                                                                 //');
fprintf('\n/////////////////////////////////////////////////////////////////////\n\n');
fprintf('  Initial point x0: (%d, %d)\n\n', x0(1), x0(2));

% ####################### Powell Method #######################
% ###### The Gauss-Seidel Part ######
% Full code comments can be found at GaussSeidel.m

x = x0; % Stores current point coordinates (x, y)
d = []; % Stores current direction of moving
ax = -1; % Stores calculated alpha for X (displacement along X axis)
ay = -1; % Stores calculated alpha for Y (displacement along Y axis)
a_range = []; % Stores calculated alpha interval
syms alpha; % Symbolic representation of alpha
syms minEq; % Symbolic representation of obtained equation to be minimized
iterationNumber = 1; % Number of iteration, counting as Gauss-Seidel
amountOfIterations = iterationNumber - 1;
x_list = zeros(15, 2); % Stores all chosen X's in format [x, y; ...] 

if(verbose == 1)
    fprintf('  Iteration [%d]\n', iterationNumber)
    fprintf('    Gauss-Seidel part:\n');
end
while((ax ~= 0 || ay ~= 0) && amountOfIterations < 2)
    goingX = 1; % Going along X axis - 1; Along Y axis - 0
    if(mod(iterationNumber, 2) == 1)
        d = d1;
        goingX = 1;
    else
        d = d2;
        goingX = 0;
    end
    xn = x + alpha*d;
    minEq = expand(F(xn));
    if(goingX == 1)
        a_range(1) = x_range(1) - x(1);
        a_range(2) = x_range(2) - x(1);
    else
        a_range(1) = x_range(1) - x(2);
        a_range(2) = x_range(2) - x(2);
    end
    a_domain = a_range(1):prec:a_range(2);
    a_values = double(subs(minEq, a_domain));
    minA_val = min(a_values);
    minA = a_domain(find(minA_val == a_values));
    if(numel(minA) > 1)
        minA = minA(1);
    end
    if(goingX == 1)
        ax = minA;
        x = double(subs(xn, ax));
    else
        ay = minA;
        x = double(subs(xn, ay));
    end
    x_list(iterationNumber, 1) = x(1); % Saving x
    x_list(iterationNumber, 2) = x(2); % Saving y
    if(verbose == 1)
        if(goingX == 1)
            fprintf('        Moving along X axis by [%g].\n', ax);
        else
            fprintf('        Moving along Y axis by [%g].\n', ay);
        end
        fprintf('        x%d = (%g, %g)\n', iterationNumber, x(1), x(2));
    end
    iterationNumber = iterationNumber + 1;
    amountOfIterations = amountOfIterations + 1;
end

if(amountOfIterations < 2)
    % Minimum already found
    fprintf('\n##############################################################\n\n');
    fprintf('    The minimum was found during Gauss-Seidel part. Algorithm stops here.\n');
    fprintf('      Minimum found at (%g, %g) with value: %g\n\n\n', x(1), x(2), F(x));
    return
end

% ###### Powell Part ######
ax_range = []; % Stores calculated alpha interval concerning X
ay_range = []; % Stores calculated alpha interval concerning Y
a_range = []; % The intersection of intervals ax_range and ay_range
a = -1; % Calculated alpha
powellIterNum = 1; % Number of iteration, counting as Powell
powellIterAmount = powellIterNum - 1; % Amount of iterations, counting as Powell

while(a ~= 0)
    % ## Calculating the direction of movement
    if(powellIterNum == 1)
        d(1) = x_list(iterationNumber-1, 1) - x0(1);
        d(2) = x_list(iterationNumber-1, 2) - x0(2);
    else
        d(1) = x_list(iterationNumber-1, 1) - x_list(iterationNumber-3, 1);
        d(2) = x_list(iterationNumber-1, 2) - x_list(iterationNumber-3, 2);
    end
    
    % ## Calculating acceptable alpha interval
    ax_range(1) = x_range(1) - x(1);
    ax_range(2) = x_range(2) - x(1);
    ay_range(1) = y_range(1) - x(2);
    ay_range(2) = y_range(2) - x(2);
    a_range(1) = max(ax_range(1), ay_range(1));
    a_range(2) = min(ax_range(2), ay_range(2));
    
    % ## Calculating the displacement coefficient 'a'
    xn = x + alpha*d;
    minEq = expand(F(xn));
    a_domain = a_range(1):prec:a_range(2);
    a_values = double(subs(minEq, a_domain));
    minA_val = min(a_values);
    minA = a_domain(find(minA_val == a_values));
    if(numel(minA) > 1)
        % If there's more than one minimum, take the first from the left
        minA = minA(1);
    end
    a = minA;
    x = double(subs(xn, a));
    
    % ## Displaying additional information if enabled
    if(verbose == 1)
        if(powellIterNum ~= 1)
            fprintf('  Iteration [%d]\n', powellIterNum);
        else
            fprintf('\n    Powell part:\n');
        end
        fprintf('        Moving in direction [%g, %g] by [%g].\n', d(1), d(2), a);
        fprintf('        x%d = (%g, %g)\n', iterationNumber, x(1), x(2));
    end
    
    % ## Incerementation of counters
    iterationNumber = iterationNumber + 1;
    powellIterNum = powellIterNum + 1;
    powellIterAmount = powellIterAmount + 1;
end

% ## Summary of what we obtained as result
fprintf('\n##############################################################\n\n');
fprintf('    Displacement in chosen direction is equal to 0. Algorithm stops here.\n');
fprintf('      Minimum found at (%g, %g) with value: %g\n\n\n', x(1), x(2), F(x));

% ############## Providing necessary information ##############
% Provide the formula of the function
function value = F(args)
x = args(1);
y = args(2);
% Functions from the lab materials (for all x,y E [-1, 1]):
%value = (x-0.5)^2 + (y-0.5)^2;
%value = (sin(x-0.5))^2 + (atan(y-0.5))^2;
value = (sin(x-0.5))^2 + (sin(y-0.5))^2;
%value = acot(x-0.5)^2 + (atan(y-0.5))^2;
%value = exp(cos(x-0.5)) + exp(cos(y-0.5));
%value = (y*(1-x^2))/(1.2 - sin(2*y));
%value = (2-3*x^4)*tan(sin(2*y+1));
%value = (2-y^2)*exp(-sin(2*x-1));
%value = (x-0.5)^2 + (y-1)^2;
end




