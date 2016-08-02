function [x, y] = DOPRI5( f,x0,xN,y0,h0 )
abstol = repmat(1e-6, length(y0), 1);
reltol = repmat(1e-6, length(y0), 1);
A = [
		[0,0,0,0,0,0,0],
		[1./5,0,0,0,0,0,0],
		[3./40,9./40,0,0,0,0,0],
		[44./45,-56./15,32./9,0,0,0,0],
		[19372./6561,-25360./2187,64448./6561,-212./729,0,0,0],
		[9017./3168,-355./33,46732./5247,49./176,-5103./18656,0,0],
		[35./384,0,500./1113,125./192,-2187./6784,11./84,0]
];
b = [
        35./384, 
        0., 
        500./1113, 
        125./192, 
        -2187./6784, 
        11./84,
        0.
];
B = [
        5179./57600, 
        0., 
        7571./16695, 
        393./640, 
        -92097./339200, 
        187./2100,
        1./40
];
facmax = 2;
power = 5;
fac = (.25)^(1./power);
x = x0;
y = y0;
n = length(y0);
while x(end) < xN
    if x(end) + h0 >= xN
        h0 = xN - x(end);
    end    
    c = sum(A, 2);
    k = zeros(length(y(:,end)), length(A));
    for i = 1:length(A)
        k(:,i) = f(x(end) + h0 * c(i), y(:,end) + h0 * k(:,1:i) * A(i,1:i)');
    end
    y1 = y(:,end) + h0 * k * b;
    Y1 = y(:,end) + h0 * k * B;                       
    sc = abstol;                     
    for i = 1:n
        sc(i) = sc(i) + (max(abs(y(i,end)) , abs(y1(i))) * reltol(i));
    end   
    err = sqrt( 1/n * sum(((y1 - Y1)./sc).^2) );     
    if err >= realmin
        r = min(facmax, max(0.1, fac * (1/err)^(1/power)) );
    else
        r = facmax;
    end                                         
    if err <= 1                   
        x(end + 1) = x(end) + h0;
        y(:,end + 1) = Y1;
        h0 = h0 * r;
        facmax = 5;
    else
        h0 = h0 * r;
        facmax = 1;
    end  
end
end