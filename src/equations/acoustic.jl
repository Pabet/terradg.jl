struct Acoustic <: Equation end
@declare_dofs Acoustic [:u, :v, :pressure, :rho, :lambda] 2

struct GaussianWave <: Scenario end

function is_periodic_boundary(eq::Acoustic, scenario::GaussianWave)
    false
end

function evaluate_boundary(eq::Acoustic, scenario::GaussianWave, face, normalidx, dofsface, dofsfaceneigh)
    s = AcousticShortcuts()
    dofsfaceneigh .= dofsface
    if normalidx == 1
        dofsfaceneigh[:,s.u] .= -dofsface[:,s.u]
    else
        dofsfaceneigh[:,s.v] .= -dofsface[:,s.v]
    end
end

function get_initial_values(eq::Acoustic, scenario::GaussianWave, global_pos; t=0.0)
    x, y = global_pos
    [0, 0, exp(-100*(x-0.5)^2 - 100*(y-0.5)^2), 1, x <= 0.5 ? 0.2 : 1]
end

function is_analytical_solution(eq::Acoustic, scenario::GaussianWave)
    false
end

function evaluate_flux(eq::Acoustic, scenario::Scenario, celldofs, cellflux)
    s = AcousticShortcuts()
    x_range = 1:size(celldofs, 1)
    y_range = x_range .+ size(celldofs, 1)
    v_flux = celldofs[:, s.pressure] ./ celldofs[:, s.rho]

    cellflux[x_range, s.u] .= v_flux
    cellflux[x_range, s.v] .= 0
    cellflux[x_range, s.pressure] .= celldofs[:, s.lambda] .* celldofs[:, s.u]
    
    cellflux[y_range, s.u] .= 0
    cellflux[y_range, s.v] .= v_flux
    cellflux[y_range, s.pressure] .= celldofs[:, s.lambda] .* celldofs[:, s.v]

    cellflux[:, s.rho:s.lambda] .= 0
end

function max_eigenval(eq::Acoustic, celldata, normalidx)
    c = celldata[:, 5] ./ celldata[:,4]
    c = sqrt(maximum(c))
end