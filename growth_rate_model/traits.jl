using CSV, DataFrames, Distributions

# Get parameter estimates
traits = [
    "egglaying", 
    "maturation", 
    "hatchingtime", 
    "hatchingsuccess", 
    "lifespan"
]
trait_values = Float64[]
trait_names = String[]

for trait in traits
    fit = CSV.read("$(@__DIR__)/trait_fits/fits/$trait.csv", DataFrame)
    fit.name .= trait .* "_" .* replace.(fit.Column1, "_mu" => "")

    append!(trait_names, fit[!, "name"])
    append!(trait_values, fit[!, "mean"])
end

traits_values = Dict{String, Float64}(zip(trait_names, trait_values))

traits_values["burrow_temp"] = 28.0
traits_values["burrow_surv_rate"] = 0.05
traits_values["juvenile_death_factor"] = 1
traits_values["max_eggs_to_cc_ratio"] = 3.0
traits_values["max_DDR"] = 1.0
traits_values["egg_death_rate"] = 0.5

#another version as Float32
traits_values_f32 = Dict{String, Float32}(zip(keys(traits_values), Float32.(values(traits_values))))

######### Get trait values
trait_posterior_draws = CSV.read("$(@__DIR__)/trait_fits/fits/draws_from_posterior.csv", DataFrame)

ndraws = nrow(trait_posterior_draws)

trait_posterior_draws.egg_death_rate = rand(Uniform(0.1, 0.3), ndraws)

trait_names = DataFrames.names(trait_posterior_draws)

trait_draws = [Dict{String, Float32}(zip(trait_names, Vector(trait_posterior_draws[i,:]))) for i in 1:ndraws]
