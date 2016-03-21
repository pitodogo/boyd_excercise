AXITER = 200;
alpha = 0.01;
beta = 0.5;
tolerance = 1e-8;
mu = 20;
duality_gap_tolerance = 1e-4;
t = 1;
m = 30;
n = 100;
l = rand(n,1);
x = l + rand(n,1);
u = x + rand(n,1);

A = rand(m,n);
s = 1.1*maximum([maximum(A*x),maximum(1./(A*x))]);




        

function Value(t,s,A,x)
   return   t*s - sum(log(s - A*x)) - sum(log(u-x)) - sum(log(x-l)) - sum(log(s*(A*x) - 1));
end


function Grad(t,s,A,x)
    
return [t - sum(1./(s - A*x)) - sum( (A*x).* (1./((s*(A*x)) - 1)) ); 
        A'*(1./(s - A*x)) +  1./(u-x) - 1./(x-l) + s*A'*(1./(s*(A*x) -1))]
end


function Hess(t,s,A,x)
    y = A*x


    z =    [  sum( (s - y).^-2 +  y./(s*y - 1).^(2) )   (-(s - y).^(-2) + (s*y - 1).^(-2) )' * A;
              A'*( -(s - y).^(-2)  - (s*y - 1).^(-2) )  A'*( diagm(vec( (s-y).^(-2) )) + diagm(vec((s./(s*y-1)).^(-2))) )*A + diagm(vec( -(u-x).^(-2) + (x-l).^(-2) ))]
    println(size(z))


    
end


println(size(Hess(t,s,A,x)))


