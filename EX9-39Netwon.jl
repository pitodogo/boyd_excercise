using PyPlot
function NewtonMethod(A,x,alpha,beta,tolerance)

    MAXITER = 1000

    history_val = []
    
    for i = 1:MAXITER
    val = Value(A,x)
    ### DIRECTION
    grad = Gradient(A,x)
    hess = Hessian(A,x)
        
    x_nt = -inv(hess)*grad
    lambda2 = grad'*inv(hess)*grad

 
    ### LINE SEARCH        
    ### PRACITICALLY, A FEASIBLE t must be bound first before line search
        t = 1
        while((maximum(A*(x+t * x_nt)) >= 1 || maximum(abs(x + t * x_nt)) >= 1))
            t = beta * t
        end

        while( Value(A,x + t * x_nt) > val - alpha * t * ((x_nt'*x_nt)[1]))
            t = beta * t
        end

    ## UPDATE
            x = x + t * x_nt

    push!(history_val,Value(A,x))

    ### CHECK EXIT CRITERIA
        if(abs(lambda2[1]/2)  <= tolerance)
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

function Hessian(A,x)
    d = 1./(1-A*x)
    m = length(x)
    n = length(d.^2)

   
    return A'*Diagonal(reshape(d.^2,n,))*A + Diagonal(reshape(1./(1+x).^2 + 1./(1-x).^2,m,))
end
            
function main()

    n = 100
    m = 200
    x = zeros(n,1)
    A = rand(m,n)
    alpha = 0.01
    beta = 0.5
    tolerance = 0.001
    answer,value,history = NewtonMethod(A,x,alpha,beta,tolerance)

    x_axis = 1:length(history)
               plot(x_axis,history-value)
               show()
end

main()
