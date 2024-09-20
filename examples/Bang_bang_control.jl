using Tapir 

function F(x,dx,u,t)

    eq1=dx[1]-x[2];
    eq2=dx[2]-u[1];

    return [eq1,eq2]
end

function Lgr_cost(x,dx,u,t)
    return 0;
end

function Mayer_cost(x0,xf,t0,tf)
    return tf;
end

function boundary_eq(x0,xf,t0,tf)
    c1=x0-[0.0,0.0];
    c2=xf-[300.0,0.0];
    c3=t0;
    c4=tf-30;

    return [c1; c2; c3; c4]
end

function boundary_ineq(x0,xf,t0,tf)
    return [];
end

function path_ineq(x,dx,u,t)
    c1=u-1;
    c2=-2-u;
    return [c1; c2]
end

model=Tapir_model(2,1);

add_pathcost!(model, Lgr_cost)
add_boundarycost!(model, Mayer_cost)

add_path_eq_constr!(model, F)
add_path_ineq_constr!(model, path_ineq)

add_boundary_eq_constr!(model, boundary_eq)
add_boundary_ineq_constr!(model, boundary_ineq)

#update_settings!()
feasibility_optimize!(model)

plot_solution(model);