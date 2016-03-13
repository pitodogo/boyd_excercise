n = 100
p = 30
A = 0
alpha = 0.01
beta = 0.5
tolerance = 10e-8

while true
    A = rand(p,n);
    if( rank(A) == p)
        break
    end
        
end

    x = rand(n,1);
    b = A*x

    

    
    function Value(x)
        return sum(x.*log(x))
    end

    function Grad(x)
        return log(x) + 1
    end

    function Hess(x)


        return diagm(vec(1./x))
    end
    






