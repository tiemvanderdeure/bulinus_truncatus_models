#####################################
####### Run mechanistic model #######
#####################################

## This code runs 5000 iterations of the mechanistic model for the entire raster
## Running this can take hours to days, depending the number of threads and processor speed
## Code by: Tiem van der Deure, tvd@sund.ku.dk


# Load packages
using Rasters, Statistics, NCDatasets, ArchGDAL, Optim

# Load global variables
include("../globals.jl")
# Load posterior trait data
include("traits.jl")
# Load the model itself
include("model.jl")
# Number of replicates
nruns = 5000

### Get peak growth rate for each draw for normalizations
peak_growth_rates = Vector{Float32}(undef, nruns)

for i in 1:nruns
    function growth_rate_Optim(x)
        return -get_pop_growth_rate(x, trait_draws[i])[1]
    end
    peak_growth_rates[i] = -optimize(growth_rate_Optim, [25.2], g_tol = 1e-12).minimum
end

### Read in data
watertemp = Float32.(Raster("$datadir/climate_data/worldclim/air_temp.tif", missingval = NaN32) .- 3.0)

s = size(watertemp)

# Preallocate
med_raster = copy(watertemp[Band = 1])
med_raster_gcm = fill(NaN32, (dims(watertemp[Band = 1])..., Rasters.Band(eachindex(models))), missingval = NaN32)
std_raster = copy(med_raster)
std_raster_gcm = copy(med_raster_gcm)
mean_raster = copy(med_raster)
mean_raster_gcm = copy(med_raster_gcm)
presence_raster = copy(med_raster)
presence_raster_gcm = copy(med_raster_gcm)
nas = isnan.(watertemp[Band = 1])

overwrite = true
nmodels = length(models)

chunksize = 100 # higher is faster but uses more RAM
output_data = zeros(Float32, s[1], chunksize, 12, Threads.nthreads())
future_data = zeros(Float32, s[1], chunksize, nmodels, nruns)
current_data = zeros(Float32, s[1], chunksize, nruns)

if ~isdir(mech_dir) mkdir(mech_dir) end

function growth_rate_and_mean!(output_data, future_data, future_climate_data, m, j)
    monthdata = view(output_data,:,:,:,Threads.threadid())
    get_pop_growth_rate!(monthdata, future_climate_data, trait_draws[j])
    mean!(view(future_data, :,:,m, j), monthdata)
    view(future_data, :,:, m, j) ./= peak_growth_rates[j]
end


for (ssp, time_period) in Iterators.product(SSPs, time_periods)
    println("starting on ssp $ssp, period $time_period")

    outdir = "$mech_dir" * join([ssp, time_period], '_')
    if ~isdir(outdir)
        mkdir(outdir)
        newdir = true
    else
        newdir = false
    end

    if overwrite | newdir
        # find files and read in
        climate_rasters = [Raster("$datadir/climate_data/worldclim/future_temp/$(join([ssp, time_period, model], '_')).tif") .- 3.0 for model in models]

        for i in 1:chunksize:s[2]
            e = i+chunksize-1

            for m in eachindex(models)
                future_climate_data = climate_rasters[m][Y(i:e)]
                @inbounds Threads.@threads for j in 1:nruns
                    growth_rate_and_mean!(output_data, future_data, future_climate_data, m, j)
                end
            end    
            
            @sync begin 
                Threads.@spawn med_raster[Y(i:e)] .= median(future_data; dims = (3, 4))
                Threads.@spawn med_raster_gcm[Y(i:e)] .= median(future_data; dims = (4))

                Threads.@spawn mean_raster[Y(i:e)] .= mean(future_data, dims = (3, 4)) # mean value
                Threads.@spawn mean_raster_gcm[Y(i:e)] .= mean(future_data, dims = (4)) # mean value

                Threads.@spawn presence_raster[Y(i:e)] .= sum(future_data .> 0, dims = (3, 4)) # how often bigger than 0?
                Threads.@spawn presence_raster_gcm[Y(i:e)] .= sum(future_data .> 0, dims = (4)) # how often bigger than 0?

                Threads.@spawn std_raster[Y(i:e)] .= std(future_data; dims = (3, 4))
                Threads.@spawn std_raster_gcm[Y(i:e)] .= std(future_data; dims = (4))       
            end
             
        end

        presence_raster[nas] .= NaN
        presence_raster ./= (nruns*nmodels)

        write("$outdir/median.tif", med_raster, force = true)
        write("$outdir/std.tif", std_raster, force = true)
        write("$outdir/mean.tif", mean_raster, force = true)
        write("$outdir/prob_present.tif", presence_raster, force = true)

        write("$outdir/median_gcm.nc", med_raster_gcm, force = true)
        write("$outdir/std_gcm.nc", std_raster_gcm, force = true)
        write("$outdir/mean_gcm.nc", mean_raster_gcm, force = true)
        write("$outdir/prob_present_gcm.nc", presence_raster_gcm, force = true)
    end
end

# quantiles and probability of presence for current
for i in 1:chunksize:s[2]
    @show i
    e = i+chunksize-1

    current_climate_data = watertemp[Y(i:e)]

    Threads.@threads for j in 1:nruns
        monthdata = view(output_data,:,:,:,Threads.threadid())
        get_pop_growth_rate!(monthdata, current_climate_data, trait_draws[j])
        mean!(view(current_data, :,:,j), monthdata)
        view(current_data, :,:,j) ./= peak_growth_rates[j]
    end            

            
    med_raster[Y(i:e)] .= median(current_data; dims = (3))
    mean_raster[Y(i:e)] .= mean(current_data, dims = (3)) # mean value
    presence_raster[Y(i:e)] .= sum(current_data .> 0, dims = (3)) # how often bigger than 0?
    std_raster[Y(i:e)] .= std(current_data; dims = (3))
end

presence_raster[nas] .= NaN
presence_raster ./= nruns

outdir = "$mech_dir/current"
if ~isdir(outdir) mkdir(outdir) end

write("$outdir/mean.tif", mean_raster, force = true)
write("$outdir/median.tif", med_raster, force = true)
write("$outdir/prob_present.tif", presence_raster, force = true)
write("$outdir/std_present.tif", std_raster, force = true)