% Contour

syms f(x, y);
f(x, y) = (sin(x-0.9))^2 + (atan(y-0.2))^2 + 0.1 * x^2;
fcontour(f, [-1, 2])
