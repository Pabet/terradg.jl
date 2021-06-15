struct Advection <: Equation end
@declare_dofs Advection [:ρ1, :ρ2, :ρ3]

struct PlanarWaves <: Scenario
end

function is_periodic_boundary(equation::Advection, scenario::PlanarWaves)
    true
end

function get_initial_values(eq::Advection, scenario::PlanarWaves, global_position; t=0.0)
    x, y = global_position
    [sinpi(2(x+y-2t)), sinpi(2(y-t)), 1.0]
end

function is_analytical_solution(equation::Advection, scenario::PlanarWaves)
    true
end


function evaluate_flux(eq::Advection, scenario::Scenario, celldofs, cellflux)
    cellflux .= [celldofs; celldofs]
end

function max_eigenval(eq::Advection, celldata, normalidx)
    # Is actually correct!
    1.0
end