mutable struct Tapir_solution
    X::Vector{Float64}
    U::Vector{Float64}
    T::Vector{Float64}

    residual::Matrix{Float64}
    objective::Float64

    status::MOI.TerminationStatusCode
    iterations::Int
    time::Float64

    Tapir_solution(X,U,T,resid,obj,stat,iter,time) = new(X,U,T,resid,obj,stat,iter,time);
end
    
    
function Tapir_solution!(sol,X,U,T,resid,obj,stat,iter,time) 
    sol.X=X;
    sol.U=U;
    sol.T=T;
    sol.residual=resid;
    sol.objective=obj;
    sol.status=stat;
    sol.iterations=iter;
    sol.time=time;
    return sol
end

function solution_update!(model,X,U,T,resid,obj,stat,iter,time)
    if isdefined(model,:solution)
        model.solution=Tapir_solution!(model.solution,X,U,T,resid,obj,stat,iter,time)
        push!(model.solution_hist,model.solution)
    else
        model.solution=Tapir_solution(X,U,T,resid,obj,stat,iter,time);
        model.solution_hist=[model.solution]
    end
end


