scenario = [1/3,1/6,1/3,1/6]
P = [1 1.3 1;2 0.5 1;0.5 1.3 1;0.5 0.5 1]
x = [0 0 1]'
A = [1 1 1]
b = 1
alpha = 0.01
beta = 0.5
tolerance = 10e-8


function Compute_Xnt()
    hess = Hess()
    AA = [hess A';A 0]
    bb = [-1*Grad();0]
    return (inv(AA)*bb)[1:3]
end
function Grad()
    return P'*(scenario.*(1./(P*x)))
end

function Hess()
    d = 1./(P*x)
    d = d.^2
    d = d.*scenario
    return P'* Diagonal(vec(d))*P  
end
function Value(x)
    ans = 0
    for i = 1:4
        ans = ans + scenario[i]*P[i,:]*x
    end
    return -1*ans[1]
end



for i = 1:200
    println(i)
    val = Value(x)
    x_nt = Compute_Xnt()
    println(val)
    hess = Hess()
    lambda = (x_nt'*hess*x_nt)[1]/2

    if(lambda < tolerance)
        break
    end
        ## LINE SEARCH
        t = 1
        while( Value(x + t*x_nt) > val + alpha*t*(Grad()'*x_nt)[1] )
            t = beta*t
        end
    x = x + t*x_nt

end






