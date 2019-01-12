clear all
clc
pkg load symbolic % Octave GNU way of loading packages -- syms, subs
setenv PYTHON C:\ProgramData\Anaconda3\python.exe % No comment...
sympref reset % No comment...
warning('off')
% OaDM Laboratory
% MULTIVARIABLE OPTIMIZATION
% 17.01.2019

% Section:
%   Tomasz GOROL
%   Magdalena LOREK
%   Sebastian CHLODEK

############## Providing necessary information ##############
% Provide the formula of the function
function value = F(args)
  x = args(1);
  y = args(2);
  %value = (x-0.5)^2 + (y-0.5)^2; % Works OK
  value = (sin(x-0.5))^2 + (atan(y-0.5))^2;
end

## Prints additional info with each iteration (1 ON, 0 OFF)
verbose = 1;

## Set the initial point, bounds for x and y & directions
x0 = [-1, -1];
x_range = [-1, 1];
y_range = [-1, 1];
e1 = [1 0]; % Moving along X axis -- in odd iterations
e2 = [0 1]; % Moving along Y axis -- in even iterations

##################### Informational box #####################
printf('\n/////////////////////////////////////////////////////////////////////');
printf('\n//                                                                 //');
printf('\n//                                                   10.01.2019    //');
printf('\n//    Tomasz Gorol                                                 //');
printf('\n//    Magdalena Lorek                                              //');
printf('\n//    Sebastian Chlodek                                            //');
printf('\n//                                                                 //');
printf('\n//                   Multivariable Optimization                    //');
printf('\n//                      Gauss-Seidel method                        //');
printf('\n//                                                                 //');
printf('\n/////////////////////////////////////////////////////////////////////\n\n');
printf('  Initial point x0: (%d, %d)\n\n', x0(1), x0(2));

#################### Gauss Seidel Method ####################
x = x0; % Stores current point coordinates (x, y)
e = []; % Stores current direction of moving
a = -1; % Stores calculated alpha (displacement coefficient in 'e' direction]
xn = []; % Stores the calculations for next x
syms alpha; % Symbolic representation of alpha
syms minEq; % Symbolic representation of obtained equation to be minimized
coeffs = []; % Stores the coefficients of the symbolic quadratic equation
iterationNumber = 1;

while(a != 0)
  ## Choosing the proper direction
  if(mod(iterationNumber, 2) == 1)
    % Odd iteration -- moving along X axis
    e = e1;
  else
    % Even iteration -- moving along Y axis
    e = e2;
  end
  
  
  ## Calculating next x
  xn = x + alpha*e % Works OK
  minEq = expand(F(xn)) % Works OK
  
  #######----> FROM HERE
  ####### Program works only if obtained equation 'minEq' is quadraic
  coeffs = fliplr(sym2poly(minEq)) %x^k is at k+1; [numel(coeffs)-1] is the equation degree
  if(coeffs(3) > 0)
    % The parabola is facing up, minimum at the peak
    a = -coeffs(2)/(2*coeffs(3)); % Works OK
  else
    % The parabola is facing down, minimum at one of the interval bounds
    if(e == e1)
      % Taking into account X bounds
      val(1) = double(subs(minEq), x_range(1));
      val(2) = double(subs(minEq), x_range(2));
      a = min(val);
    else
      % Taking into account Y bounds
      val(1) = double(subs(minEq), y_range(1));
      val(2) = double(subs(minEq), y_range(2));
      a = min(val);
    end
  end
  x = double(subs(xn, a));
  #######----> UPTO HERE
  
  ## Printing additional information if 'verbose' enabled
  if(verbose == 1)
    printf('  Iteration [%d]\n', iterationNumber);
    printf('    a = %d\n', a);
    if(e == e1)
      printf('    Moving along X axis by [%d].\n', a);
    else
      printf('    Moving along Y axis by [%d].\n', a);
    end
    printf('    x%d = (%d, %d)\n\n', iterationNumber, x(1), x(2));
  end
  
  ## Incrementing the total iterations amount  
  iterationNumber = iterationNumber + 1;
end

## Summary of what we obtained as result
printf('\n##############################################################\n\n');
printf('    Displacement is equal to 0. Algorithm stops here.\n');
printf('      Minimum found at (%d, %d) with value: %d\n\n\n', x(1), x(2), F(x));
