function [ x ] = fixediter( phi,x0,tol,maxit )
x1 = phi(x0);
it = 1;
while norm(x1 - x0) > tol && it < maxit
    x0 = x1;
    x1 = phi(x0);
    it = it + 1;
end
x = x0;
end