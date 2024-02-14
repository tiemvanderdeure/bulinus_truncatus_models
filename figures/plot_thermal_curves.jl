###########################
### Plotting trait fits ###
###########################

# Import packages
using CairoMakie, StatsBase

# Import trait data
include("../growth_rate_model/traits.jl")

# Include the model for the final pop_growth_rate_normalized
include("../growth_rate_model/model.jl")

##########################
## First calculate life history traits at each temperature
##########################

# Which temperatures
temperatures = collect(5:0.01:40)

# Basic functions
quadratic(t, T0, Tm, q) = q * (t - T0) * (Tm - t)
gdd_function(t, T0, GDD) = ifelse(t > T0, (t-T0)/GDD, Float32(0))

# Egg-laying rate
ELR = reduce(hcat,
    [quadratic.(temperatures, draw["egglaying_T0"], draw["egglaying_Tm"], draw["egglaying_q"]) for draw in trait_draws]
)
ELR = max.(ELR, Ref(0))

# Maturation time
mat_time = reduce(hcat,
 [1 ./ gdd_function.(temperatures, draw["maturation_T0"], draw["maturation_GDD"]) for draw in trait_draws]
)

# Lifespan
lifespan = reduce(hcat,
    [exp.(quadratic.(temperatures, draw["lifespan_T0"], draw["lifespan_Tm"], draw["lifespan_q"])) for draw in trait_draws]
)

# Hatching success
hatching_success = [draw["hatchingsuccess_q"] for t in temperatures, draw in trait_draws]

# Hatching time
hatching_time = reduce(hcat,
    [Float32(1) ./ gdd_function.(temperatures, draw["hatchingtime_T0"], draw["hatchingtime_GDD"]) for draw in trait_draws]
)

# Population growth rate
pop_growth_rate = mapreduce(traits -> get_pop_growth_rate(temperatures, traits), hcat, trait_draws)
pop_growth_rate_normalized = max.(pop_growth_rate ./ maximum(pop_growth_rate, dims = 1), Ref(0)) # Normalize to 0-1

# Get 2.5%, 50%, and 97.5% quantile for plots
posteriors_to_plot = [hatching_success, hatching_time, mat_time, lifespan, ELR, pop_growth_rate_normalized]
quantiles_to_plot = [hcat(quantile.(eachrow(posterior), [[0.025, 0.5, 0.975]])...) for posterior in posteriors_to_plot]

########################
####### Plotting #######
########################

# parameters for each plot
plot_titles = ["Hatching success", "Hatching time", "Maturation time", "Lifespan", "Egg-laying rate", "Population growth"]
ylabels = ["Rate", "Time (days)", "Time (weeks)", "Time (weeks)", "Eggs/week", "Relative rate"]
ylimits = [1, 40, 25, nothing, nothing, 1] # upper limits for y axis

# Import data points for each plot
data_points = CSV.read("growth_rate_model/trait_fits/data/data_points_plot.csv", DataFrame) # these are generated in R
data_to_plot = [filter(x -> x.parameter == par, data_points) for par in ["hatching_success", "hatching_time", "maturation_time", "lifespan", "eggs"]]

# Make the plot
f = Figure(size = (1000, 600), fontsize = 16, grid = false, backgroundcolor = :transparent) 

axes = Axis[]
indices = [(row,col) for col in 1:3, row in 1:2]

# Loop through plots
for (i, data) in (enumerate(quantiles_to_plot))
    ax = Axis(f[indices[i]...],
        title = plot_titles[i],
        xlabel = "Temperature (°C)",
        ylabel = ylabels[i],
        limits = (5, 40, 0, ylimits[i]),
        ygridvisible = false, xgridvisible = false, 
        backgroundcolor = :transparent
    )

    lines!(ax, temperatures, data[2,:], color = :black)

    for j in [1,3] lines!(ax, temperatures, data[j,:], color = :black, linestyle=:dash) end

    if i < 6
        scatter!(ax, 
            data_to_plot[i].temperature, data_to_plot[i].trait, 
            color = :black,
            marker = :circle,
        )
    end
    append!(axes, [ax])
end

# Lifespan on pseudolog scale
axes[4].yscale = Makie.pseudolog10
axes[4].yticks = [0, 1, 3, 10, 30, 100, 300]

## Save the output
save("figures/figure_03.pdf", f)
save("figures/figure_03.png", f, px_per_unit = 10)

##########################
## Find values mentioned in text
##########################

# Optimum temperature
T_opt_pgr = map(x -> temperatures[findmax(x)[2]], eachcol(pop_growth_rate))
Statistics.median(T_opt_pgr)
quantile(T_opt_pgr, [0.025, 0.975])

# Minimum suitable temperature
T_min_pgr = map(x -> temperatures[findfirst(x .> 0)], eachcol(pop_growth_rate))
Statistics.median(T_min_pgr)
quantile(T_min_pgr, [0.025, 0.975])

# Maximum suitable temperature
T_max_pgr = map(x -> temperatures[findlast(x .> 0)], eachcol(pop_growth_rate))
Statistics.median(T_max_pgr)
quantile(T_max_pgr, [0.025, 0.975])

# Quantiles for egg laying max range
quantile(map(ls -> temperatures[findlast(ls .> 0)], eachcol(ELR)), quantiles)

# Quantiles for peak lifespan
quantile(map(x -> temperatures[findmax(x)[2]], eachcol(lifespan)), quantiles)

# Quantiles for hatching success
quantile(map(x -> temperatures[findmax(x)[2]], eachcol(lifespan)), quantiles)