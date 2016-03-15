include("parameter_generator.jl")
iteration = 0;

for iteration = 1:1000

    val = (b'* v)[1] +sum(exp(-A'*v - 1))
    grad = b - A*exp(-A'v - 1)

    hess = A * Diagonal(vec(exp(-A'*v - 1))) * A'
    v_nt = -hess\grad
    fprime = grad'*v_nt

    if(abs(fprime[1]) <= tolerance)
        break
    end
        
        t = 1




        while( (b'*(v + t * v_nt))[1] + sum( exp(-A'*(v + t * v_nt) - 1 )) > val + t * alpha * (fprime[1]))
            t = beta * t
        end
    v = v + t * v_nt
        end
        println(iteration)
       
