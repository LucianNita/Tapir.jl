function feasibility_optimize!(model::Tapir_model)

    model_JuMP=Model(Ipopt.Optimizer);
    @variable(model_JuMP, X[i=1:length(model.Xguess)], start=model.Xguess[i]);
    @variable(model_JuMP, U[i=1:length(model.Uguess)], start=model.Uguess[i]);
    @variable(model_JuMP, T[i=1:length(model.Tguess)], start=model.Tguess[i]);

    @constraint(model_JuMP,model.boundary_eq_const(X[1:model.nx],X[end-model.nx+1:end],T[1],T[end]).=0) #compare to vectorwise
    @constraint(model_JuMP,model.boundary_ineq_const(X[1:model.nx],X[end-model.nx+1:end],T[1],T[end]).<0)

    for i=1:model.mesh.N
        u=reshape(U[model.mesh.st_len_u*(i-1)+1:model.mesh.st_len_u*(i-1)+model.nu*(model.mesh.Pu+1)],(model.nu,model.mesh.Pu+1))
        u=u*model.mesh.evalXU;
        x=reshape(X[model.mesh.st_len_x*(i-1)+1:model.mesh.st_len_x*(i-1)+model.nx*(model.mesh.Px+1)],(model.nx,model.mesh.Px+1))
        
        if model.settings.flex_mesh
            t1=T[i];
            t2=T[i+1];
        else
            t1=model.mesh.outer_mesh[i]*(T[2]-T[1])/2+(T[2]+T[1])/2;
            t2=model.mesh.outer_mesh[i+1]*(T[2]-T[1])/2+(T[2]+T[1])/2;
        end
        
        dx=x*model.mesh.Dx*2.0./(t2-t1);
        tp=model.mesh.inner_x*0.5*(t2-t1).+0.5*(t2+t1);

        for j=1:model.mesh.Px+1
            if i!=1 && j==1 && model.settings.cont_type_x=="samevar"
                continue
            end
            @constraint(model_JuMP,model.path_ineq_const(dx[:,j],x[:,j],u[:,j],tp[j]).<=0);
        end
    end

    
    if model.settings.cont_type_x=="constr"
        #ADD constraints
    end

    if model.settings.cont_type_u=="constr"
        #ADD constraints
    end
    

    @objective(model_JuMP,Min,sum(sum(residual(X,U,T,model))))

    optimize!(model_JuMP)

    solution_update!(model, value.(X), value.(U), value.(T), status)

    return value.(X), value.(U), value.(T), objective_value(model_JuMP)
end

function residual(X,U,T,model::Tapir_model)
    ir = Matrix{Float64}(undef, model.ne, model.mesh.N);

    for i=1:model.mesh.N
        x=reshape(X[model.mesh.st_len_x*(i-1)+1:model.mesh.st_len_x*(i-1)+model.nx*(model.mesh.Px+1)],(model.nx,model.mesh.Px+1))
        u=reshape(U[model.mesh.st_len_u*(i-1)+1:model.mesh.st_len_u*(i-1)+model.nu*(model.mesh.Pu+1)],(model.nu,model.mesh.Pu+1))

        if model.settings.flex_mesh
            t1=T[i];
            t2=T[i+1];
        else
            t1=model.mesh.outer_mesh[i]*(T[2]-T[1])/2+(T[2]+T[1])/2;
            t2=model.mesh.outer_mesh[i+1]*(T[2]-T[1])/2+(T[2]+T[1])/2;
        end

        dx=x*model.mesh.Dx*model.mesh.evalQX*2.0./(t2-t1);
        x=x*model.mesh.evalQX;
        u=u*model.mesh.evalQU;

        qp=model.mesh.quad_pts*0.5*(t2-t1).+0.5*(t2+t1);
        qw=model.mesh.quad_weights*0.5*(t2-t1);

        ir[:,i]=sum(qw[k]*model.path_eq_const(dx[:,k],x[:,k],u[:,k],qp[k]).^2 for k in eachindex(qp))
    end
    return ir;
end






