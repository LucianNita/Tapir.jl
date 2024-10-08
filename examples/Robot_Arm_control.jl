using Tapir 

function F(dx,x,u,t)
    eq1= dx[1] - ((sin(x[3])*(9/4*cos(x[3])*x[1]^2)+2*x[2]^2 + 4/3*(u[1]-u[2]) - 3/2*cos(x[3])*u[2] )/ (31/36 + 9/4*sin(x[3])^2));
    eq2= dx[2] + ((sin(x[3])*(9/4*cos(x[3])*x[2]^2)+7/2*x[1]^2 - 7/3*u[2] + 3/2*cos(x[3])*(u[1]-u[2]) )/ (31/36 + 9/4*sin(x[3])^2));
    eq3= dx[3] - (x[2]-x[1]);
    eq4= dx[4] - x[1];
    eq5= dx[5] - 0.01*(u[1]*u[1]+u[2]*u[2])
    
    return [eq1,eq2,eq3,eq4,eq5]
end

function Lgr_cost(x,dx,u,t)
    return 0.01*(u[1]*u[1]+u[2]*u[2]);
end

function Mayer_cost(x0,xf,t0,tf)
    return tf+xf[5];
end

function boundary_eq(x0,xf,t0,tf)
    c1=x0-[0.0,0.0,0.0,0.0,0.0];
    c2=xf[1:4]-[0.0,0.0,0.5,0.522];
    c3=t0;
    #c4=tf-30.;

    return [c1; c2; c3]
end

function boundary_ineq(x0,xf,t0,tf)
    return [];
end

function path_ineq(x,dx,u,t)
    c1=u[1]-1;
    c2=-1-u[1];
    c3=u[2]-1;
    c4=-1-u[2];
    return [c1; c2; c3; c4]
end

#=
model=Tapir_model(5,2);
model.mesh=Tapir_mesh(200,3,2,7)
model.settings=Tapir_settings("feasibility","Ipopt",discretization="Inte_Res",resid_tol=5.0*10^(-8));

add_pathcost!(model, Lgr_cost);
add_boundarycost!(model, Mayer_cost);

add_path_eq_constr!(model, F);
add_path_ineq_constr!(model, path_ineq);

add_boundary_eq_constr!(model, boundary_eq);
add_boundary_ineq_constr!(model, boundary_ineq);

update_settings!(model);
#feasibility_optimize!(model);
control_optimize!(model);

plot(model);
=#

#Plots.scatter()
for N=5:10
model2=Tapir_model(5,2);
model2.mesh=Tapir_mesh(N,3,2,7)
model2.settings=Tapir_settings("feasibility","Ipopt",discretization="Inte_Res",flex_mesh=false, resid_tol=10^(-4)/N);

add_pathcost!(model2, Lgr_cost);
add_boundarycost!(model2, Mayer_cost);

add_path_eq_constr!(model2, F);
add_path_ineq_constr!(model2, path_ineq);

add_boundary_eq_constr!(model2, boundary_eq);
add_boundary_ineq_constr!(model2, boundary_ineq);

update_settings!(model2);
#feasibility_optimize!(model);
control_optimize!(model2);

if N==5
    lbl="Integrated residuals & Fixed mesh"#Integrated residuals & Flexible mesh"#Integrated residuals
else
    lbl=""
end
Plots.scatter!([model2.solution.time],[accuracy_check(model2,model)],color=:blue,marker=:rect,label=lbl);#marker=, legend=false);,marker=:rect
end
#Tapir.plot(model2);,,1/N,model2.solution.objective
#xlabel!("1/N");
#ylabel!("Optimal cost")