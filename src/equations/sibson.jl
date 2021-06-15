struct Sibson <: Equation end
@declare_dofs Sibson [:ρ1]

struct Worksheet2 <: Scenario 
end

function is_periodic_boundary(eq::Sibson, scenario::Worksheet2)
    true
end

function get_initial_values(eq::Sibson, scenario::Worksheet2, global_position; t=0.0)
    x, y = global_position
    [cospi(4*sqrt((x - .25)^2 + (y - .25)^2))]
end

function is_analytical_solution(eq::Sibson, scenario::Worksheet2)
    true
end

function evaluate_flux(eq::Sibson, scenariop::Scenario, celldofs, cellflux)
    cellflux .= [celldofs; celldofs]
end

function max_eigenval(eq::Sibson, celldata, normalidx)
    1.0
end