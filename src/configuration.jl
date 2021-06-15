"""
    Configuration(configfile::String)

Parses configuration file `configfile`.
"""
struct Configuration
    equation_name::String
    scenario_name::String

    end_time::Float64
    order::Int64
    timeintegrator_name::String
    courant::Float64

    filter_name::String
    filter_order::Int64
    plot_start::Float64
    plot_step::Float64
    grid_elements::Int64
    physicalsize::Array{Float64,1}
    cellsize::Array{Float64,1}

    function Configuration(configfile::String)
        config = YAML.load(open(configfile))
        equation_name = config["equation"]["equation"]
        scenario_name = config["equation"]["scenario"]

        end_time = config["simulation"]["end_time"]
        order = config["solver"]["order"]
        timeintegrator_name = config["solver"]["timeintegrator"]
        courant = config["solver"]["courant"]

        filter_name = get(config["solver"], "filter", "identity")
        filter_order = get(config["solver"], "filter_order", order)
        plot_start = config["output"]["start"]
        plot_step = config["output"]["step"]
        grid_elements = config["simulation"]["grid_elements"]
        physicalsize_1d = config["simulation"]["grid_size"]
        physicalsize = [physicalsize_1d, physicalsize_1d]
        cellsize = physicalsize ./ grid_elements

        new(
            equation_name,
            scenario_name,
            end_time,
            order,
            timeintegrator_name,
            courant,
            filter_name,
            filter_order,
            plot_start,
            plot_step,
            grid_elements,
            physicalsize,
            cellsize
        )
    end
end