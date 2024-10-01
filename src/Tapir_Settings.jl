mutable struct Tapir_settings
    discretization_method::String #IR/colloc/
    flex_mesh::Bool          #Fixed/flex
    cont_type_x::String       #constr/samevar/discont
    cont_type_u::String       #constr/samevar/discont
    cost_type::String         #path_cost/aug_state
    
    outer_mesh_start::String  #Equi
    quad_type::String   #lgr
    inner_mesh_type::String   #Cheb2
    interpolation_type::String#Lagr_Bary
    
    mode::String #feas/OC/feas+OC
    solver::String #Ipopt/other

    resid_tol::Float64  #dyn constraints enforcement
    solver_tol::Float64 #solver accuracy
    quad_tol::Float64   #quadrature tolerance 
    flex_tol::Float64   #flexible mesh tolerance

    Tapir_settings(mode, solver; discretization="Collocation", flex_mesh=true, cont_type_x="same_var", cont_type_u="discontinuous", cost_type="path_cost", outer_mesh_start="equidistant", quad_type="gausslegendre", inner_mesh_type="cheb2", interpol="lagrange_bary", resid_tol=10^-6, solver_tol=10^-8, quad_tol=10^-10, flex_tol=0.01) = new(discretization, flex_mesh, cont_type_x, cont_type_u, cost_type, outer_mesh_start, quad_type, inner_mesh_type, interpol, mode, solver, resid_tol, solver_tol, quad_tol, flex_tol) #Basic constructor function with defaults
end

