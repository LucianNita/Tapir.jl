mutable struct Tapir_model
    nx::Int
    nu::Int

    path_cost::Function
    terminal_cost::Function

    path_eq_const::Function
    path_ineq_const::Function

    boundary_eq_const::Function
    boundary_ineq_const::Function

    settings::Tapir_settings;
    mesh::Tapir_mesh;
    mesh_hist::Vector{Tapir_mesh}
    solution::Tapir_solution;
    solution_hist::Vector{Tapir_solution}

    Xguess::Vector{Float64}
    Uguess::Vector{Float64}
    Tguess::Vector{Float64}
    
    ne::Int
    Tapir_model(nx,nu) = new(nx,nu)
end

function add_pathcost!(model, Lgr)
    model.path_cost=Lgr;
end

function add_boundarycost!(model, Mayer)
    model.terminal_cost=Mayer;
end

function add_path_eq_constr!(model, eq_con)
    model.path_eq_const=eq_con;
    model.ne=length(eq_con(zeros(model.nx),zeros(model.nx),zeros(model.nu),0.0));
end

function add_path_ineq_constr!(model, ineq_con)
    model.path_ineq_const=ineq_con;
end

function add_boundary_eq_constr!(model, b_eq)
    model.boundary_eq_const=b_eq;
end

function add_boundary_ineq_constr!(model, b_ineq)
    model.boundary_ineq_const=b_ineq;
end

function incremented_T_guess(old_mesh,new_N)
    T_new=zeros(new_N+1);
    old_N=length(old_mesh)-1;
    if old_N+1!=new_N
        warning("step larger than 1, redirecting to the main function")
        return get_T_guess(old_mesh,new_N)
    end
    interval_lengths=diff(old_mesh);
    ~, index = findmax(interval_lengths);
    T_new[1:index]=copy(old_mesh[1:index]); #deep?
    T_new[index+1]=0.5*(old_mesh[index]+old_mesh[index+1]);
    T_new[index+2:end]=copy(old_mesh[index+1:end]);
    return T_new;
end

function update_settings!(model::Tapir_model;ws=false)
    if !isdefined(model,:settings)
        model.settings=Tapir_settings("feasibility","Ipopt");
    end
    if !isdefined(model,:mesh)
        model.mesh=Tapir_mesh(10,3,2,7);
    end
        set_length!(model)
        generate_quad!(model.mesh, model.settings.quad_type)
        generate_outer!(model.mesh, model.settings.outer_mesh_start)
        generate_interpolation_mesh!(model.mesh,model.settings.inner_mesh_type)
        generate_interpolation!(model.mesh,model.settings.interpolation_type)
        model.mesh_hist=[model.mesh];

        if model.settings.cont_type_x=="same_var"
            model.Xguess=zeros(model.mesh.st_len_x*model.mesh.N+model.nx);
        else
            model.Xguess=zeros(model.mesh.st_len_x*model.mesh.N);
        end
        if model.settings.cont_type_u=="same_var"
            model.Uguess=zeros(model.mesh.st_len_u*model.mesh.N+model.nu);
        else
            model.Uguess=zeros(model.mesh.st_len_u*model.mesh.N);
        end

        if ws && isdefined(model,:solution)
            if model.settings.flex_mesh
                model.Tguess=incremented_T_guess(model.solution.T,model.mesh.N);#[3.1*i/model.mesh.N for i=0:model.mesh.N];
            else
                model.Tguess=model.solution.T;
            end
            
            crtx=0;
            crtu=0;
            for i=1:model.mesh.N
                t1=model.Tguess[i];
                t2=model.Tguess[i+1];
                ipx=model.mesh.inner_x*0.5*(t2-t1).+0.5*(t2+t1);
                ipu=model.mesh.inner_u*0.5*(t2-t1).+0.5*(t2+t1);
                for j in eachindex(ipx)
                    if model.settings.cont_type_x=="same_var" && i!=1 
                        continue
                    end
                    model.Xguess[crtx+1:crtx+model.nx]=eval_X(model,ipx[j])
                    crtx+=model.nx;
                end
                for j in eachindex(ipu)
                    if model.settings.cont_type_u=="same_var" && i!=1 
                        continue
                    end
                    model.Uguess[crtu+1:crtu+model.nu]=eval_U(model,ipu[j])
                    crtu+=model.nu;
                end
            end
        else
            if model.settings.flex_mesh
                model.Tguess=[3.1*i/model.mesh.N for i=0:model.mesh.N];
            else
                model.Tguess=[0.0,3.1];
            end
        end

    return model;
end

