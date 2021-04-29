%Calculate 1D Gaussian quadrature on interval [-1,1]

function [points,weights] = gauss_1d(n)

eps = 3.0E-15;

m = floor((n+1)/2);
z1 = 3;

for i=1:m
    z = cos(pi*(i-0.25)/(n+0.5));
    while(abs(z-z1)>eps)
        p1 = 1;
        p2 = 0;
        for j=1:n
            p3=p2;
            p2=p1;
            p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
        end
        pp = n*(z*p1-p2)/(z*z-1.0);
        z1 = z;
        z = z1-p1/pp;
    end
    
    points(i) = -z;
    points(n+1-i) = z;
    weights(i) = 2.0/((1.0-z*z)*pp*pp);
    weights(n+1-i) = weights(i);
end

end
        
        