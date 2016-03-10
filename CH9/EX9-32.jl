using PyPlot

function generate_psd_matrix(n)
    A = rand(n,n)
    A = A + A' #Gaurantee Symmetric
    A = A + n*eye(n)
    return A
end


function main()
    As = []
    bs = []
    n = 50
    m = 100
    x = zeros(n,1)
    for i = 1:m
    A = generate_psd_matrix(n)
    b = rand(n,1)
        push!(As,A)
        push!(bs,b)
    end
    GNMethod(As,bs,m,x)
end

function Second(As,bs,x,m)
    ans = 0
    for i = 1:m
        ans = ans + Value(As[i],bs[i],x)[1]*Grad(As[i],bs[i],x)
    end
    return ans
end

function First_Inverse(As,bs,x,m)
    ans = 0
    for i = 1:m
        grad = Grad(As[i],bs[i],x)

        ans = ans + grad*grad'
    end

    return inv(ans)
end

function Sum_Grad(As,bs,x,m)
    ans = 0
    for i = 1:m
        grad = Grad(As[i],bs[i],x)
        ans = ans + grad
    end
    return ans
end


function Sum_Value(As,bs,x,m)
    ans = 0
    for i = 1:m
        val = Value(As[i],bs[i],x)
        ans = ans + val
    end
    return ans
end
function GNMethod(As,bs,m,x)
    alpha = 0.01
    beta = 0.5
    tolerance = 10e-8
    MAXITER = 1000
    history_val = []
    hess = 0
    for i = 1:MAXITER
        println(i)
        val = Sum_Value(As,bs,x,m) 
        x_nt = First_Inverse(As,bs,x,m)*Second(As,bs,x,m)
        lambda2 = Sum_Grad(As,bs,x,m)'*x_nt        
        ### LINE SEARCH
        ### PRACITICALLY, A FEASIBLE t must be bound first before line search
        t = 1

            while( Sum_Value(As,bs,x + t * x_nt,m)[1] > val[1] - alpha * t * ((x_nt'*x_nt)[1]))
                t = beta * t
            end
        ## UPDATE
           x = x + t * x_nt

                ### CHECK EXIT CRITERIA
                if(abs(lambda2[1]/2)  <= tolerance)
                    return x,val,history_val
                end
                    

        
    end
end

function Hess(A,b,x)
    return A
end

function Grad(A,b,x)

return A*x + b
end

function Value(A,b,x)
    return 0.5*x'*A*x+b'*x + 1
end





main()
