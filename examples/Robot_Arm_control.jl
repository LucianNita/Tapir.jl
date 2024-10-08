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
function run()
    #Plots.scatter()
    Nruns=1;
    Nstart=8;
    Nend=10;
    mk=:dtriangle;
    clr=:orange;
    lbl=""#Integrated residual flexible mesh"
    disc="Collocation"
    fm=true;
    wst=false;
    xtot=zeros(Nend-Nstart+1);
    ytot=zeros(Nend-Nstart+1);

    for run=1:Nruns
        xplot=[];
        yplot=[];

        N=Nstart;
        model2=Tapir_model(5,2);
        model2.mesh=Tapir_mesh(N,3,2,7)
        model2.settings=Tapir_settings("feasibility","Ipopt",discretization=disc,flex_mesh=fm, resid_tol=10^(-4)/N);

        add_pathcost!(model2, Lgr_cost);
        add_boundarycost!(model2, Mayer_cost);

        add_path_eq_constr!(model2, F);
        add_path_ineq_constr!(model2, path_ineq);

        add_boundary_eq_constr!(model2, boundary_eq);
        add_boundary_ineq_constr!(model2, boundary_ineq);

        update_settings!(model2);

        #feasibility_optimize!(model);
        control_optimize!(model2);

        push!(xplot,1/N);
        push!(yplot,model2.solution.objective);
        #Plots.scatter!([model2.solution.time],[accuracy_check(model2,model)],color=clr,marker=mk,label=lbl);#marker=, legend=false);,marker=:rect

        for N=Nstart+1:Nend
            model2.mesh=Tapir_mesh(N,3,2,7)
            model2.settings=Tapir_settings("feasibility","Ipopt",discretization=disc,flex_mesh=fm, resid_tol=10^(-4)/N);

            update_settings!(model2,ws=wst);

            #feasibility_optimize!(model);
            control_optimize!(model2);

            push!(xplot,1/N);
            push!(yplot,model2.solution.objective);
            #Plots.scatter!([model2.solution.time],[accuracy_check(model2,model)],color=clr,marker=mk,label="");#marker=, legend=false);,marker=:rect
        end
        xtot+=xplot;
        ytot+=yplot;
    end

    scatter!(xtot/Nruns,ytot/Nruns,color=clr,marker=mk,label=lbl);

    #Tapir.plot(model2);,,1/N,model2.solution.objective,accuracy_check(model2,model)
    xlabel!("1/N");
    ylabel!("Optimal cost");#Total absolute error ϵₜ");
end