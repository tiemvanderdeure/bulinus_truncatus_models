#####################################
# Performance for mechanistic model #
#####################################

# This code runs 5000 iterations of the mechanistic model for different water temperatures 
# at presence and background points and calculates performance statistics
# Code by: Tiem van der Deure, tvd@sund.ku.dk

# Load packages
using Rasters, Statistics, NCDatasets, ArchGDAL

# Load global variables
include("../globals.jl")
# Load posterior trait data
include("traits.jl")
# Load the model itself
include("model.jl")
# Number of replicates
nruns = 5000

####################
#### Load data #####
####################

# Presences and background points
presences = CSV.read("occurrence_data/thinned_presences.csv", DataFrame)
background10k = CSV.read("occurrence_data/background_points.csv", DataFrame)
npresences = size(presences, 1); nbackground = size(background10k, 1)

presence_points = (GeoInterface.Point(p.x, p.y) for p in eachrow(presences))
background_points = (GeoInterface.Point(p.x, p.y) for p in eachrow(background10k))
all_points = Iterators.flatten((eachrow(presences), eachrow(background10k)))

# Load temperature data
T_air = Raster("$datadir/climate_data/worldclim/air_temp.tif")
T_air_smoothed = Raster("$datadir/climate_data/worldclim/smoothed_air_temp.tif")
T_water_Paaijmans = Raster("$datadir/climate_data/worldclim/T_water_body_model.tif")

# Extract temperature data at presence points
T_air_points = reduce(vcat, [collect(T_air[X = Near(p.x), Y = Near(p.y)]) for p in all_points]')
T_air_smoothed_points = reduce(vcat, [collect(T_air_smoothed[X = Near(p.x), Y = Near(p.y)]) for p in all_points]')
T_water_Paaijmans_points = reduce(vcat, [collect(T_water_Paaijmans[X = Near(p.x), Y = Near(p.y)]) for p in all_points]')

####################
#### Run model #####
####################

# Initialise lists to save outcomes
AUC_list = [zeros(nruns) for i in 1:9]
sensitivity_list = [zeros(nruns) for i in 1:9]
specificity_list = [zeros(nruns) for i in 1:9]

# Initialise matrices to save output
water_temperature = copy(T_air_points)
output_data = copy(T_air_points)
mean_output_data = copy(output_data[:,1])
boolean_matrix_AUC = Matrix{Bool}(undef, npresences, nbackground)

n_exp = 1 # keep track of which experiment 

function run_and_get_metrics(
    i, 
    u_output_data, u_water_temperature, 
    output_data, water_temperature, mean_output_data, boolean_matrix_AUC;
    npresences = npresences, nbackground = nbackground
)
    get_pop_growth_rate!(u_output_data, u_water_temperature, trait_draws[i])
    d = Dict(w => o for (w, o) in zip(u_water_temperature, u_output_data))

    output_data .= getindex.(Ref(d), water_temperature)
    mean!(mean_output_data, output_data)

    presences_output = view(mean_output_data, 1:npresences)
    background_output = view(mean_output_data, npresences+1:npresences+nbackground)

    boolean_matrix_AUC .= [p > b for p in presences_output, b in background_output]
    AUC = Statistics.mean(boolean_matrix_AUC)
    sens = Statistics.mean(presences_output .>= 0)
    sel = Statistics.mean(background_output .< 0)

    return (AUC, sens, sel)
end

# calculate AUC for air temperature and smoothed air temperature with an offset
for air_temp in [T_air_points, T_air_smoothed_points]
    u_air_temp = unique(air_temp)
    u_output_data = similar(u_air_temp)

    for T_offset in 0.0:1:3.0
        water_temperature .= air_temp .- T_offset
        u_water_temperature = u_air_temp .- T_offset

        for i in 1:nruns
            # Generate and save AUC, sensitivity, 
            AUC_list[n_exp][i], sensitivity_list[n_exp][i], specificity_list[n_exp][i] = run_and_get_metrics(
                i,
                u_output_data, u_water_temperature, 
                output_data, water_temperature, mean_output_data, boolean_matrix_AUC
            )
        end
        n_exp += 1
    end
end

## Do the same for shallow water body output
u_water = unique(T_water_Paaijmans_points)
u_output_data = similar(u_water)

for i in 1:nruns
    AUC_list[n_exp][i], sensitivity_list[n_exp][i], specificity_list[n_exp][i] = run_and_get_metrics(
        i,
        u_output_data, u_water, 
        output_data, T_water_Paaijmans_points, mean_output_data, boolean_matrix_AUC
    )
end

### summarize, format and write outcomes
TSS_list = [sens .+ spec .- 1 for (sens, spec) in zip(sensitivity_list, specificity_list)]

AUC_quantiles = [quantile(AUC, quantiles) for AUC in AUC_list]
TSS_quantiles = [quantile(TSS, quantiles) for TSS in TSS_list]
sensitivity_quantiles = [quantile(sens, quantiles) for sens in sensitivity_list]
specificity_quantiles = [quantile(spec, quantiles) for spec in specificity_list]

AUC_df = DataFrame(reduce(vcat, AUC_quantiles'), string.(quantiles))
TSS_df = DataFrame(reduce(vcat, TSS_quantiles'), string.(quantiles))
sensitivity_df = DataFrame(reduce(vcat, sensitivity_quantiles'), string.(quantiles))
specificity_df = DataFrame(reduce(vcat, specificity_quantiles'), string.(quantiles))

performance_df = DataFrame(
    reduce(hcat, [["$(r[2]) (95% CI $(r[1])-$(r[3]))" for r in eachrow(round.(df; sigdigits = 2))] for df in [AUC_df, TSS_df, sensitivity_df, specificity_df]]),
    ["AUC", "TSS", "sensitivity", "specificity"]
)

CSV.write("growth_rate_model/outputs/AUC.csv", AUC_df)
CSV.write("growth_rate_model/outputs/AUC.csv", TSS_df)
CSV.write("growth_rate_model/outputs/sensitivity.csv", sensitivity_df)
CSV.write("growth_rate_model/outputs/specificity.csv", specificity_df)
        
CSV.write("growth_rate_model/outputs/performance.csv", performance_df) ## Table in supplementaries