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
        if ws
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
            if model.settings.flex_mesh==false
                model.Tguess=[0.0,3.1];
            else
                model.Tguess=[3.1*i/model.mesh.N for i=0:model.mesh.N];
            end
        else
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
            if model.settings.flex_mesh==false
                model.Tguess=[0.0,3.1];
            else
                model.Tguess=[3.1*i/model.mesh.N for i=0:model.mesh.N];
            end
        end
    return model;
end

