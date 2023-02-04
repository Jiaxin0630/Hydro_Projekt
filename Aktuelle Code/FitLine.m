function par = FitLine(x,y,P)
A = x;
par = inv(A'*P*A)*A'*P*y;
end