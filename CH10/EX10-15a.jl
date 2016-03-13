include("parameter_generator.jl")
using PyPlot
### Calculate X_nt
val_hist = []
for iteration = 1:1000
    val = Value(x)
    push!(val_hist,val)
    hess = Hess(x)
    grad = Grad(x)

    AA = [hess A';A 0*eye(p)]
    bb = [-grad;zeros(p,1)]
    s = AA\bb
    x_nt = s[1:n]
    w = s[n+1:end]

    lambda = (x_nt'*hess*x_nt)[1]

    if(lambda / 2 < tolerance)
        break
    end

        t = 1

        ######### CAUTION #########
        while(minimum(x + t*x_nt) <=0)
            t = beta*t
        end
            

        while(Value(x + t * x_nt) > Value(x) + t * alpha * (grad'*x_nt)[1] )
            t = beta*t
        end
    x = x + x_nt
end


x_axis = 1:length(val_hist)
plot(x_axis,val_hist-val_hist[end])
show()
        

