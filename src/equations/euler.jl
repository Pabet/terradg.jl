struct Euler <: Equation end
@declare_dofs Euler [:rhou, :rhov, :rho, :rhoE] 

struct ShockTube <: Scenario end

function is_periodic_boundary(eq::Euler, scenario::GaussianWave)
    true
end

function is_periodic_boundary(eq::Euler, scenario::ShockTube)
    false
end

function evaluate_boundary(eq::Euler, scenario::ShockTube, face, normalidx, dofsface, dofsfaceneigh)
    s = EulerShortcuts()
    dofsfaceneigh .= dofsface
    if normalidx == 1
        dofsfaceneigh[:,s.rhou] .= -dofsface[:,s.rhou]
    else
        dofsfaceneigh[:,s.rhov] .= -dofsface[:,s.rhov]
    end
end

function get_initial_values(eq::Euler, scenario::GaussianWave, global_pos; t=0.0)
    pxg, pyg = global_pos
    p = exp(-100 * (pxg - 0.5)^2 - 100 *(pyg - 0.5)^2) + 1
    rho = 1.0
    rhoE = evaluate_energy(eq, 0,0,rho,p)
    [0.0, 0.0, rho, rhoE]
end

function get_initial_values(eq::Euler, scenario::ShockTube, global_pos; t=0.0)
    pxg, pyg = global_pos
    if pxg < 0.5
        rho = 0.125
        p = 0.1
    else
        rho = 1.0
        p = 1.0
    end
    rhoE = evaluate_energy(eq, 0,0,rho,p)
    [0.0, 0.0, rho, rhoE]
end

function is_analytical_solution(eq::Euler, scenario::Scenario)
    false
end

function evaluate_flux(eq::Euler, scenario::Scenario, celldofs, cellflux)
    s = EulerShortcuts()
    x_range = 1:size(celldofs, 1)
    y_range = x_range .+ size(celldofs, 1)
    eval_pressure = (rho, rhoE, rhou, rhov) -> evaluate_pressure(eq, rho, rhoE, rhou, rhov)

    # TODO compare performance between calculating this each time it is required vs allocating here
    p = eval_pressure.(celldofs[:,s.rho], celldofs[:,s.rhoE], celldofs[:,s.rhou], celldofs[:,s.rhov])

    cellflux[x_range, s.rhou] .= celldofs[:,s.rhou] .* celldofs[:,s.rhou] ./ celldofs[:,s.rho] .+ p
    cellflux[x_range, s.rhov] .= celldofs[:,s.rhou] .* celldofs[:,s.rhov] ./ celldofs[:,s.rho]
    cellflux[x_range, s.rho] .= celldofs[:,s.rhou]
    cellflux[x_range, s.rhoE] .= celldofs[:,s.rhou] .* (celldofs[:,s.rhoE] .+ p) ./ celldofs[:,s.rho]
    
    cellflux[y_range, s.rhou] .= celldofs[:,s.rhou] .* celldofs[:,s.rhov] ./ celldofs[:,s.rho]
    cellflux[y_range, s.rhov] .= celldofs[:,s.rhov] .* celldofs[:,s.rhov] ./ celldofs[:,s.rho] .+ p
    cellflux[y_range, s.rho] .= celldofs[:,s.rhov]
    cellflux[y_range, s.rhoE] .= celldofs[:,s.rhov] .* (celldofs[:, s.rhoE] .+ p) ./ celldofs[:,s.rho]
end

function max_eigenval(eq::Euler, celldata, normalidx)
    s = EulerShortcuts()

    v_idx = normalidx == 1 ? s.rhou : s.rhov
    max_eigenval = 0
    for i = 1:(size(celldata, 1))
        inv_rho = inv(celldata[i, s.rho])
        cs = sqrt(1.4 * evaluate_pressure(eq, celldata[i, s.rho], celldata[i, s.rhoE], celldata[i, s.rhou], celldata[i, s.rhov]) * inv_rho)
        vₙ = celldata[i, v_idx] * inv_rho
        max_eigenval = max(max_eigenval, vₙ + cs)
    end

    max_eigenval
end

function evaluate_energy(eq::Euler, imx, imy, rho, p)
    #0.4 for γ=1.4
    p / 0.4 + 0.5 / rho * (imx^2 + imy^2)
end

function evaluate_pressure(eq::Euler, rho, E, imx, imy)
    #0.4 for γ=1.4
    0.4 * (E - 0.5 / rho * (imx^2 + imy^2))
end