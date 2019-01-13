clear
clc
syms f(x, y);
% OaDM Laboratory
% MULTIVARIABLE OPTIMIZATION
% 17.01.2019

% Section:
%   Tomasz GOROL
%   Magdalena LOREK
%   Sebastian CHLODEK

% ############## Providing necessary information ##############
%f(x, y) = (x-0.5)^2 + (y-0.5)^2;
%f(x, y) = (sin(x-0.5))^2 + (atan(y-0.5))^2;
%f(x, y) = (sin(x-0.5))^2 + (sin(y-0.5))^2;
f(x, y) = acot(x-0.5)^2 + (atan(y-0.5))^2;
%f(x, y) = exp(cos(x-0.5)) + exp(cos(y-0.5));
%f(x, y) = (y*(1-x^2))/(1.2 - sin(2*y));
%f(x, y) = (2-3*x^4)*tan(sin(2*y+1));
%f(x, y) = (2-y^2)*exp(-sin(2*x-1));

% ## Set the initial point, variables bounds and precision
X0 = [-1, -1];
x_range = [-1, 1];
y_range = [-1, 1];
prec = 0.01;

% ## Prints additional info with each iteration (1 ON, 0 OFF)
verbose = 1;

% ##################### Informational box #####################
fprintf('\n/////////////////////////////////////////////////////////////////////');
fprintf('\n//                                                                 //');
fprintf('\n//    OaDM Laboratory                                17.01.2019    //');
fprintf('\n//                                                                 //');
fprintf('\n//                   Multivariable Optimization                    //');
fprintf('\n//                    Steepest Descent Method                      //');
fprintf('\n//                                                                 //');
fprintf('\n/////////////////////////////////////////////////////////////////////\n\n');
fprintf('  Initial point x0: (%d, %d)\n\n', X0(1), X0(2));

% ################## Steepest Descent Method ##################
gradFormula = gradient(f)'; % Gradient of function f, symbolics
gradF = double(subs(gradFormula, {x, y}, {X0(1), X0(2)})); % Stores the gradient at x
X = X0; % Stores current point coordinates (x, y)
Xn = []; % Stores formula for the next x
syms lambda; % Symbolic representation of lambda
syms minEq; % Symbolic representation of obtained equation to be minimized
L = -1; % Lambda - displacement in gradient (gradF) direction
iterNum = 1; % Iteration number
iterAmount = iterNum - 1; % Amount of iterations
stopCalculs = 0;
prevFx = []; % Stores previous function value at X

while(~(L == 0 || (gradF(1) == 0 && gradF(2) == 0) || stopCalculs == 1))
    % ## Calculating the gradient vector
    gradF = double(subs(gradFormula, {x, y}, {X(1), X(2)}));
    
    % ## Calculating acceptable alpha interval
    ax_range(1) = x_range(1) - X(1);
    ax_range(2) = x_range(2) - X(1);
    ay_range(1) = y_range(1) - X(2);
    ay_range(2) = y_range(2) - X(2);
    a_range(1) = max(ax_range(1), ay_range(1));
    a_range(2) = min(ax_range(2), ay_range(2));
    
    % ## Calculating the displacement coefficient 'a'
    Xn = X - lambda*gradF;
    minEq = expand(f(Xn(1), Xn(2)));
    a_domain = a_range(1):prec:a_range(2);
    a_values = double(subs(minEq, a_domain));
    minA_val = min(a_values);
    minA = a_domain(find(minA_val == a_values));
    if(numel(minA) > 1)
        % If there's more than one minimum, take the first from the left
        minA = minA(1);
    end
    L = minA;
    prevFx = double(f(X(1), X(2)));
    X = double(subs(Xn, L));
    
    % ## Some problem with memory?
    if(numel(X) == 0)
        % Probably problem with asymptotes?
        error('Encountered a problem with calculations on enormously huge numbers.');
    end
    
    % ## Avoiding too much precision
    currFx = double(f(X(1), X(2)));
    if(prevFx * prec * 0.999999999 <= currFx && iterAmount >= 10)
        stopCalculs = 1;
    end
    
    % ## Displaying additional information if enabled
    if(verbose == 1)
        fprintf('  Iteration [%d]\n', iterNum);
        fprintf('        Moving in direction [%g, %g] by [%g].\n', gradF(1), gradF(2), L);
        fprintf('        x%d = (%g, %g)\n', iterNum, X(1), X(2));
        fprintf('        F(x) = %g\n\n', f(X(1), X(2)));
    end
    
    % ## Incrementing the total iterations amount and number
    iterNum = iterNum + 1;
    iterAmount = iterAmount +1;
end

% ## Summary of what we obtained as result
fprintf('\n##############################################################\n\n');
if(stopCalculs == 1)
    fprintf('    Stopping the calculations because of too little differences in F(x).\n');
else
    fprintf('    Displacement in chosen direction is equal to 0. Algorithm stops here.\n');
end
fprintf('      Minimum found at (%g, %g) with value: %g\n\n\n', X(1), X(2), f(X(1), X(2)));





