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
%f(x, y) = x - y + 2*x^2 + 2*x*y + y^2; % OK
f(x, y) = (x-0.5)^2 + (y-0.5)^2; % OK
%f(x, y) = (sin(x-0.5))^2 + (atan(y-0.5))^2; % OK
%f(x, y) = (sin(x-0.5))^2 + (sin(y-0.5))^2; % OK
%f(x, y) = acot(x-0.5)^2 + (atan(y-0.5))^2; % NOT YET
%f(x, y) = exp(cos(x-0.5)) + exp(cos(y-0.5)); % OK
%f(x, y) = (y*(1-x^2))/(1.2 - sin(2*y)); % NOT YET
%f(x, y) = (2-3*x^4)*tan(sin(2*y+1)); % NOT YET
%f(x, y) = (2-y^2)*exp(-sin(2*x-1)); % OK

% ## Set the initial point, variables bounds and precision
X0 = [0, 0];
x_range = [-2, 2];
y_range = [-2, 2];
prec = 0.01;
sameValueThreshold = 0.99999; % If F(x_i-1) * {this var} <= F(X_i), STOP

% ## Prints additional info with each iteration (1 ON, 0 OFF)
verbose = 1;

% ##################### Informational box #####################
fprintf('\n/////////////////////////////////////////////////////////////////////');
fprintf('\n//                                                                 //');
fprintf('\n//    OaDM Laboratory                                17.01.2019    //');
fprintf('\n//                                                                 //');
fprintf('\n//                   Multivariable Optimization                    //');
fprintf('\n//                   Conjugate Gradient Method                     //');
fprintf('\n//                                                                 //');
fprintf('\n/////////////////////////////////////////////////////////////////////\n\n');
fprintf('  Function: ');
disp(formula(f));
fprintf('  Initial point x0: (%d, %d)\n\n\n', X0(1), X0(2));

% ################# Conjugate Gradient Method #################
gradFormula = gradient(f)'; % Gradient of function f, symbolics
gradF = double(subs(gradFormula, {x, y}, {X0(1), X0(2)})); % Stores the gradient at x
prevGradF = []; % Stores previous gradient at x
S = -gradF; % Stores current direction of displacement
prevS = []; % Stores previous direction of displacement
B = -1; % Stores the beta coefficient
X = X0; % Stores current point coordinates (x, y)
prevX = []; % Stores previous point
Xn = []; % Stores formula for the next x
syms lambda; % Symbolic representation of lambda
syms minEq; % Symbolic representation of obtained equation to be minimized
L = -1; % Lambda - displacement in gradient (gradF) direction
iterNum = 1; % Iteration number
iterAmount = iterNum - 1; % Amount of iterations
stopCalculs = 0; % Breaks the loop if 1
prevFx = []; % Stores previous function value at X

while(~(L == 0 || (S(1) == 0 && S(2) == 0) || stopCalculs == 1))
    % ## Calculating the S vector
    if(iterNum == 1)
        gradF = double(subs(gradFormula, {x, y}, {X(1), X(2)}));
        S = -gradF;
    else
        prevS = S;
        prevGradF = gradF;
        gradF = double(subs(gradFormula, {x, y}, {X(1), X(2)}));
        %B = norm(gradF)^2 / norm(prevGradF)^2
        subGradK = subs(gradFormula, [x, y], [prevX(1), prevX(2)])';
        subGradKp1 = subs(gradFormula, [x, y], [X(1), X(2)])';
        B = double((subGradKp1' * (subGradKp1-subGradK))/(norm(subGradK)^2));
        S = -gradF + B * prevS;
    end
    
    
    
    % ## Calculating acceptable lambda interval
    ax_range(1) = x_range(1) - X(1);
    ax_range(2) = x_range(2) - X(1);
    ay_range(1) = y_range(1) - X(2);
    ay_range(2) = y_range(2) - X(2);
    a_range(1) = max(ax_range(1), ay_range(1));
    a_range(2) = min(ax_range(2), ay_range(2));
    
    % ## Calculating the displacement coefficient 'a'
    Xn = X + lambda*S;
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
    prevX = X;
    X = double(subs(Xn, L));
    
    % ## Some problem with memory?
    if(numel(X) == 0)
        % Probably problem with asymptotes?
        error('Encountered a problem with calculations on enormously huge numbers.');
    end
    
     % ## Avoiding too much precision
    currFx = double(f(X(1), X(2)));
    if(prevFx * sameValueThreshold <= currFx && iterAmount >= 10)
        stopCalculs = 1;
    end
    
    % ## Displaying additional information if enabled
    if(verbose == 1)
        fprintf('  Iteration [%d]\n', iterNum);
        fprintf('        Moving in direction [%g, %g] by [%g].\n', S(1), S(2), L);
        fprintf('        x%d = (%g, %g)\n', iterNum, X(1), X(2));
        fprintf('        F(x) = %g\n\n', f(X(1), X(2)));
    end
    
    % ## Incrementing the total iterations amount and number
    iterNum = iterNum + 1;
    iterAmount = iterAmount + 1;
end

% ## Summary of what we obtained as result
fprintf('\n##############################################################\n\n');
if(stopCalculs == 1)
    fprintf('    Stopping the calculations because of too little differences in F(x).\n');
else
    fprintf('    Displacement in chosen direction is equal to 0. Algorithm stops here.\n');
end
fprintf('      Minimum found at (%g, %g) with value: %g\n\n\n', X(1), X(2), f(X(1), X(2)));





