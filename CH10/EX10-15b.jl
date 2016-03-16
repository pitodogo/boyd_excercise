include("parameter_generator.jl")
x = rand(n,1);
v = zeros(p,1)
b = A*x



for iteration = 1:1000
println(iteration)
grad = Grad(x)
hess = Hess(x)

AA = [hess A';A 0*eye(p)]



bb = -1*[Grad(x);A*x-b]
s = AA\bb
x_nt = s[1:n]
v_nt = s[n+1:end] - v
r = [Grad(x) + A'*v;A*x-b]    


t = 1
while(minimum(x+t*x_nt) <= 0)
     t = beta*t;
end

    #### LINE SEARCH
    
    while(norm([Grad(x+t*x_nt) + A'*(v+t*v_nt);A*(x + t*x_nt)-b]) > (1-alpha*t)*norm(r))
        t = beta*t
    end

        x = x + t * x_nt
        v = v + t * v_nt
        if (norm(A*x - b) < tolerance && norm(r) <= tolerance)
            break
        end
        
end

println(x)    
