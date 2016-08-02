function [ x,y ] = ITR( f,x0,xN,y0,N,tol,maxit )
    h = ( xN - x0 )/N;
    x = x0:h:xN;
    y = zeros(length(y0), N + 1);
    y(:,1) = y0;
    for k = 1:N
        phi = @(y1) y(:,k) + h/2 * (f(x(k),y(:,k)) + f(x(k + 1), y1));
        y(:,k + 1) = fixediter( phi,y(:,k),tol );
    end
end