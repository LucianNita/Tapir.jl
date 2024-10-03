function accuracy_check(model::Tapir_model, ref::Function, t0,tf)
    accuracy=0.0;
    dt=0.001
    #make sure t0 and tf match 

    for t=t0:dt:tf
        accuracy+=sum((eval_U(model,t)-ref).^2);
    end
    return accuracy*dt
end
