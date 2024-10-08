function accuracy_check(model::Tapir_model, ref_model::Tapir_model; dt=0.001)
    accuracy=0.0;
    
    t0=min(model.solution.T[1],ref_model.solution.T[1]);
    tf=max(model.solution.T[end],ref_model.solution.T[end]);

    for t=t0:dt:tf
        u=eval_U(model,t);
        u_ref=eval_U(ref_model,t);
        if length(u)!=length(u_ref)
            error("Control inputs for model and reference should be the same size")
        end

        for i in eachindex(u)
            accuracy+=(u[i]-u_ref[i])^2;
        end
    end
    return accuracy*dt
end
