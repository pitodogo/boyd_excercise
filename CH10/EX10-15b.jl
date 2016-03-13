include("EX10-15a.jl")


for iteration = 1:1000
grad = Grad(x)
hess = Hess(x)

    AA = [hess A';A 0*eye(p)]
    bb = [-grad();Ax-b]
    s = AA\bb
    x_nt = s[1:n]
    v_nt = s[n+1:end]
    t = 1
    while(minimum(x+t*x_nt) <= 0)
        t = beta*t;
    end

        while(
end
