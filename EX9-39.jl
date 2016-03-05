function GradientDescentMethod(A,x,alpha,beta,tolerance)
    MAXITER = 1000
    val = Value(A,x)
    for i = 1:MAXITER
    ### DIRECTION
    delta_x = -Gradient(A,x)

 
 
    ### LINE SEARCH        
    ### PRACITICALLY, A FEASIBLE t must be bound first before line search
        t = 1
        while((maximum(A*(x+t * delta_x)) >= 1 || maximum(abs(x + t*delta_x)) >= 1))
            t = beta*t
        end


        while( Value(A,x + t * delta_x) > val - alpha * t * (delta_x'*delta_x)[1])
            t = beta * t
        end

    ## UPDATE
    x = x + t * delta_x

    ### CHECK EXIT CRITERIA
        if(norm(delta_x) < tolerance)
            return x,val
        end                
        println(x)
    end

    return x,val
end


   function Value(A,x)    
    return -1*( sum(log(1-A*x)) + sum(log(1-map(x->x^2,x))))

end

function Gradient(A,x)
   return A'*(1./(1-A*x)) - 1./(1+x) + 1./(1-x)
end

function main()
    n = 3
    m = 10
    x = zeros(n,1)
    A = rand(m,n)
    alpha = 0.01
    beta = 0.5
    tolerance = 0.001
    GradientDescentMethod(A,x,alpha,beta,tolerance)
end

main()
