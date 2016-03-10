using PyPlot
function NewtonMethod(A,x,alpha,beta,tolerance)

    MAXITER = 1000

    history_val = []
    val = 0
    for i = 1:MAXITER
    val = Value(A,x)
    ### DIRECTION
        grad = Gradient(A,x)
        
    hess = Diag_Hessian(A,x)
        println(i)
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


function Diag_Hessian(A,x)
    ppp = Hessian(A,x)
    dp = diag(ppp)
    return Diagonal(dp)
end
#function Diag_Hessian(A,x)
#    d = 1./(1-A*x)
#    p = A'*d
#    d = diag(p*p')
#    m = size(A)[1]
#    n = size(A)[2]
#    return -1*Diagonal(d) + Diagonal(reshape(1./(1+x).^2 - 1./(1-x).^2,n,))
#end


function Hessian(A,x)

        d = 1./(1-A*x)
        m = length(x)
        n = length(d.^2)
        return A'*Diagonal(reshape(d.^2,n,))*A + Diagonal(reshape(1./(1+x).^2 + 1./(1-x).^2,m,))

end
            
            
function main()

    n = 3
    m = 3
    x = zeros(n,1)
    A = rand(m,n)
    alpha = 0.01
    beta = 0.5
    tolerance = 10e-8
#println( Hessian(A,x))
#println( Diag_Hessian(A,x))
    answer,value,history = NewtonMethod(A,x,alpha,beta,tolerance)

    x_axis = 1:length(history)
               plot(x_axis,history-value)
    show()

end

main()
