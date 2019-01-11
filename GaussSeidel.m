clear all
clc
pkg load symbolic % Octave GNU way of loading packages -- syms, subs
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
  value = (x-0.5)^2 + (y-0.5)^2;
end

## Prints additional info with each iteration (1 ON, 0 OFF)
verbose = 1;

## Set the initial point, bounds for x and y & directions
x0 = [-1, -1];
x_range = [-1, 1];
y_range = [-1, 1];
e1 = [1 0]; % Moving along X axis -- in odd iterations
e2 = [0 1]; % Moving along Y axis -- in even iterations


#################### Gauss Seidel Method ####################
x = x0; % Stores current point coordinates (x, y)
e = []; % Stores current direction of moving
a = []; % Stores calculated alpha (displacement coefficient in 'e' direction]
syms alpha; % Symbolic representation of alpha
iterationNumber = 1;

while(iterationNumber < 30)
  ## Choosing the proper direction
  if(mod(iterationNumber, 2) == 1)
    % Odd iteration -- moving along X axis
    e = e1;
  else
    % Even iteration -- moving along Y axis
    e = e2;
  end
  
  
  ##
  
  
  ## Printing additional information if 'verbose' enabled
  if(verbose == 1)
    printf('  Iteration [%d]\n', iterationNumber);
    printf('    x%d = (%d, %d)\n', iterationNumber, x(1), x(2));
    printf('    a = %d\n', a);
    if(e == e1)
      printf('    Moving along X axis by [%d].\n\n', a);
    else
      printf('    Moving along Y axis by [%d].\n\n', a);
    end
  end
  
  ## Incrementing the total iterations amount  
  iterationNumber = iterationNumber + 1;
end

s = F(x);