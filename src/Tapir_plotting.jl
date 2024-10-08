using Plots

function plot_x(model::Tapir_model)
    p=Plots.plot();
    X=model.solution.X;
    T=model.solution.T;
    for i=1:model.mesh.N
        x=reshape(X[model.mesh.st_len_x*(i-1)+1:model.mesh.st_len_x*(i-1)+model.nx*(model.mesh.Px+1)],(model.nx,model.mesh.Px+1))

        if model.settings.flex_mesh
            t1=T[i];
            t2=T[i+1];
        else
            t1=model.mesh.outer_mesh[i]*(T[2]-T[1])/2+(T[2]+T[1])/2;
            t2=model.mesh.outer_mesh[i+1]*(T[2]-T[1])/2+(T[2]+T[1])/2;
        end
        tp=model.mesh.inner_x*0.5*(t2-t1).+0.5*(t2+t1);

        for d=1:model.nx
            scatter!(p,tp,x[d,:],color=:red)
            xf(t) = (x[d,:]'*L(2*(t)/(t2-t1)-(t2+t1)/(t2-t1), model.mesh.inner_x))[1]
            plot!(p,xf,t1,t2,legend=false,color=:red);
        end
    end
    return p
end

function plot_u(model::Tapir_model)
    p=Plots.plot();
    U=model.solution.U;
    T=model.solution.T;
    for i=1:model.mesh.N
        u=reshape(U[model.mesh.st_len_u*(i-1)+1:model.mesh.st_len_u*(i-1)+model.nu*(model.mesh.Pu+1)],(model.nu,model.mesh.Pu+1))

        if model.settings.flex_mesh
            t1=T[i];
            t2=T[i+1];
        else
            t1=model.mesh.outer_mesh[i]*(T[2]-T[1])/2+(T[2]+T[1])/2;
            t2=model.mesh.outer_mesh[i+1]*(T[2]-T[1])/2+(T[2]+T[1])/2;
        end
        tp=model.mesh.inner_u*0.5*(t2-t1).+0.5*(t2+t1);
        
        for d=1:model.nu
            scatter!(p,tp,u[d,:],color=:black)
            uf(t) = (u[d,:]'*L(2*(t)/(t2-t1)-(t2+t1)/(t2-t1), model.mesh.inner_u))[1]
            plot!(p,uf,t1,t2,legend=false,color=:black);
        end
    end
    return p
end

function plot(model::Tapir_model)
    p=Plots.plot(plot_x(model),plot_u(model), layout = (2, 1));
    return p
end

function eval_X(model::Tapir_model,t)
    X=model.solution.X;
    T=model.solution.T;
    if t<model.solution.T[1] #reorder this
        return X[1:model.nx]
    elseif t>model.solution.T[end]
        return X[end-model.nx+1:end]
    else
        for i=1:model.mesh.N
            if model.settings.flex_mesh
                t1=T[i];
                t2=T[i+1];
            else
                t1=model.mesh.outer_mesh[i]*(T[2]-T[1])/2+(T[2]+T[1])/2;
                t2=model.mesh.outer_mesh[i+1]*(T[2]-T[1])/2+(T[2]+T[1])/2;
            end
            if t<=t2 || i==model.mesh.N
                x=reshape(X[model.mesh.st_len_x*(i-1)+1:model.mesh.st_len_x*(i-1)+model.nx*(model.mesh.Px+1)],(model.nx,model.mesh.Px+1));
                x_t=x*L(2*(t)/(t2-t1)-(t2+t1)/(t2-t1), model.mesh.inner_x);
                return x_t;
            end
        end   
    end     
end

function eval_U(model::Tapir_model,t)
    U=model.solution.U;
    T=model.solution.T;
    if t<model.solution.T[1] #reorder this
        return U[1:model.nu]
    elseif t>model.solution.T[end]
        return U[end-model.nu+1:end]
    else
        for i=1:model.mesh.N
            if model.settings.flex_mesh
                t1=T[i];
                t2=T[i+1];
            else
                t1=model.mesh.outer_mesh[i]*(T[2]-T[1])/2+(T[2]+T[1])/2;
                t2=model.mesh.outer_mesh[i+1]*(T[2]-T[1])/2+(T[2]+T[1])/2;
            end
            if t<=t2 || i==model.mesh.N
                u=reshape(U[model.mesh.st_len_u*(i-1)+1:model.mesh.st_len_u*(i-1)+model.nu*(model.mesh.Pu+1)],(model.nu,model.mesh.Pu+1))
                u_t=u*L(2*(t)/(t2-t1)-(t2+t1)/(t2-t1), model.mesh.inner_u);
                return u_t;
            end
        end 
    end
end



