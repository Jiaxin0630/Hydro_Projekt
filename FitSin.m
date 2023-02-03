function [x_dach,y_dach,Ax] = FitSin(t,y,f)
t = seconds(t - t(1));
% asin(wt)+bcos(wt)
% y0 = a0*sin(2*pi*f*t + phi0);
A = [];
for i = 1:length(f)
    A = [A sin(2*pi*f(i)*t) cos(2*pi*f(i)*t)];    
end
x_dach = inv(A'*A)*A'*y;
y_dach = y - A*x_dach;
Ax = A*x_dach;
end