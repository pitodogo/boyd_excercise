using PyPlot
function GradientDescentMethod(A,x,alpha,beta,tolerance)

    MAXITER = 1000

    history_val = []
    
    for i = 1:MAXITER
    val = Value(A,x)
    ### DIRECTION
    delta_x = -Gradient(A,x)

 
    ### LINE SEARCH        
    ### PRACITICALLY, A FEASIBLE t must be bound first before line search
        t = 1
        while((maximum(A*(x+t * delta_x)) >= 1 || maximum(abs(x + t*delta_x)) >= 1))
            t = beta * t
        end

        while( Value(A,x + t * delta_x) > val - alpha * t * ((delta_x'*delta_x)[1]))
            t = beta * t
        end

    ## UPDATE
            x = x + t * delta_x

    push!(history_val,Value(A,x))

    ### CHECK EXIT CRITERIA
        if(norm(-delta_x) < tolerance)
            return x,val,history_val
        end                

    end

    return x,val,history_val
end


   function Value(A,x)    
    return -1*( sum(log(1-A*x)) + sum(log(1+x)) + sum(log(1-x)))
   end

function Gradient(A,x)
   return A'*(1./(1-A*x)) - 1./(1+x) + 1./(1-x)
end

           function main()

    n = 100
    m = 200
    x = zeros(n,1)
    A = rand(m,n)
    alpha = 0.01
    beta = 0.5
    tolerance = 0.001
    answer,value,history = GradientDescentMethod(A,x,alpha,beta,tolerance)

    x_axis = 1:length(history)
               plot(x_axis,history-value)
               show()
end

main()
