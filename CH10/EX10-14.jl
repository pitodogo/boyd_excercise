scenario = [1/3,1/6,1/3,1/6]
P = [2 1.3 1;2 0.5 1;0.5 1.3 1;0.5 0.5 1]
x = [0,0,1]
A = [1 1 1]
b = 1
alpha = 0.01
beta = 0.5
tolerance = 10e-8


function Compute_Xnt()
    hess = Hess()
    AA = [hess A';A 0]
    bb = [-1*Grad();0]
    return -1*(AA\bb)[1:3]
end

function Grad()
    return -1*P'*(scenario.*(1./(P*x)))
end

function Hess()
    d = 1./(P*x)
    d = d.^2
    d = d.*scenario
    return P'* diagm(d)*P  
end
function Value(x)


    ans = sum(scenario.*P*x)

    return ans*-1
end


function all_nonzero(x)
    for i = 1:length(x)
        if(x[i] < 0)
            return false
        end
        end
        return true
end


for i = 1:200
    println(i)
    val = Value(x)
    println(x)
    x_nt = Compute_Xnt()
    println(val)
    hess = Hess()
    lambda = (x_nt'*hess*x_nt)[1]/2

    if(lambda < tolerance)
        break
    end
        ## LINE SEARCH
        t = 1
        while( Value(x + t*x_nt) > val + alpha*t*(Grad()'*x_nt)[1] && all_nonzero(x+t*x_nt) )
            t = beta*t
        end
    x = x + t*x_nt

end






