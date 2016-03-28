include("parameter_generator_11-21.jl")


for i = 1:MAXITER

### Compute optimal z t*f0 + phi



    #### NEWTON METHOD ####

        val = Value(t,s,A,x);
        grad = Grad(t,s,A,x);
        hess = Hess(t,s,A,x);
        step = -1*(hess\grad);
        lambda = grad'*step;
    println(s)
        if(abs(lambda[1]) < tolerance) ### SATISFY NEWTON METHOD
            ###  !!!!!!!!!!!!!!!! DESERVE NOTICING !!!!!!!!!!!!!!!! ###
            ### m in the textbook refers to the number of inequality constraints which is not m in this case but 3m + 2n
            if( (3m + 2n)/t < duality_gap_tolerance )
                break;
            end
            t = mu*t;
        else

            t_ls = 1;
            ds = step[1];
            dx = step[2:end];
            dy = A*dx;

            new_s = s + t_ls * ds;
            new_x = x + t_ls * dx;
            new_y = A * new_x;
            ### CHECK DOMAIN ###
                while( minimum([new_s - new_y;u - new_x;new_x - l;new_s - 1./new_y;new_y]) <= 0)
                    t_ls = beta * t_ls;
                    new_s = s + t_ls * ds;
                    new_x = x + t_ls * dx;
                    new_y = A * new_x
                end
                println(t_ls)
                while(Value(t,new_s,A,new_x) > val + alpha*t_ls*lambda[1])
                    t_ls = beta * t_ls;
                    new_s = s + t_ls * ds;
                    new_x = x + t_ls * dx;
                end

                s = s + t_ls * ds;
                x = x + t_ls * dx;                        
        end
    end
    
    





println(A*x)
