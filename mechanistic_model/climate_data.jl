###########################
## Download climate data ##
###########################

## This script downloads the climate data used in mechanistic modelling
## The first part generates all temperature data for which performance was calculated, including for the Paaijmans (2008) model
## The second part downloads temperature data used for future prediction, which is only air temperature
## All data is downloaded from the internet and saved in datadir (defined in globals.jl)

# Load globals
include("../globals.jl")

##################
## Current data ##
##################

# Load packages
using Rasters, ArchGDAL, Statistics, NLsolve, ZipFile, HTTP

# Variables
vars = [:prec, :srad, :vapr, :wind, :tavg]
datarasters = Raster[]

# Load data. First download and write all the files, if they are not already downloaded
for var in vars
    dir = "$datadir/climate_data/worldclim/current/wc2.1_2.5m_$var"
    if ~isdir(dir) 
        zipped_bytes = HTTP.get("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_$var.zip")
        zip = ZipFile.Reader(IOBuffer(zipped_bytes.body))

        mkdir(dir)

        for i in eachindex(zip.files)
            filename_i = string(i+100)[2:3] # add 0 where appropriate

            write("$dir/wc2.1_2.5m_$(var)_$filename_i.tif", read(zip.files[i]))
        end
    end

    # Read in rasters
    rasters = [read(view(Raster(file, lazy = true), X(-20 .. 60), Y(-40 .. 60))) 
        for file in readdir(dir, join = true)]
    rast = cat(rasters..., dims = Rasters.Band(1:12))
    rast = replace_missing(rast, NaN)
    @eval $var = $rast # Assign raster to variable name
end

##########################################
# Shallow pool model from Paaijmans 2008 #
##########################################

# Convert to the right units
vapr .*= 10 # kpa to hpa
srad .*= (1000/(3600*24)) # kJ/day to Watt
tavg_K = tavg .+ 273.15 # Temperature in K

# Physical constants
const σ = 5.67e-8 # Boltzmanns constant
const gas_constant = 287.052874 # specific gas constant for air
const c_p = 1005.0 # heat capacity of air at constant pressure in J kg-1 K-1
const C_h = 3e-3 # Heat transfer coefficient
const ϵ_w = 0.98 # emissivity of water
const P = 101325 # atmospheric pressure, assuming standard atmospheric pressure
const m = 4.81e-26 # molecular mass of air

# function for the change in water temperature, as described in Paaijmans (2008)
function delta_Tw(Tw, Ta, K_in, vapr, wind; α = 0.15, d = 30.0)
    Tw_C = Tw-273.15

    L_in = (0.56 + 0.067*sqrt(vapr)) * σ * Ta^4 + d
    L_out = ϵ_w * σ * Tw^4

    air_density = P / (gas_constant * Ta)  # Density of Air
    H = air_density * c_p * wind * C_h * (Tw - Ta)

    sat_vapor_pressure_surface = 611.2*exp(17.67*Tw_C/(Tw_C+243.5))
    mass_mixing_ratio_surface = 0.622*sat_vapor_pressure_surface/P
    specific_humidity_surface = mass_mixing_ratio_surface/(mass_mixing_ratio_surface+1.0)

    mass_mixing_ratio_air = 0.622*vapr*100/P
    specific_humidity_air = mass_mixing_ratio_air/(mass_mixing_ratio_air+1.0)

    LE = air_density * 2.260e6 * wind * C_h * (specific_humidity_surface - specific_humidity_air)

    R_net = K_in*(1-α) + L_in - L_out # net radiation
    Gs = 0.15*R_net # soil heat flux

    delta = R_net - H - LE - Gs

    return delta
end

# set up solver to solve for delta_Tw = 0
function find_equilibrium_Tw(Ta, K_in, vapr, wind)
    if isnan(Ta) return NaN
    else
        solver(F, x) = F[1] = delta_Tw(x[1], Ta, K_in, vapr, wind)
        solve = nlsolve(solver, [Ta])
        return solve.zero[1]
    end
end

# get water temperature assuming equilibrium
Tw = broadcast(tavg_K, srad, vapr, wind) do T, s, v, w
    find_equilibrium_Tw(T, s, v, w)
end

Tw_C = max.(Tw .- 273.15, Ref(0)) # convert to C and cap to 0
write("$datadir/climate_data/worldclim/T_water_body_model.tif", Tw_C) # write file

##########################
# Smooth air temperature #
##########################

tavg_smoothed = copy(tavg)
for i in 1:12
    tavg_smoothed[Band = i] .= mean(tavg[Band = mod1.([i,i-1, i-2], 12)],dims = 3)[:,:,1]
end

write("$datadir/climate_data/worldclim/smoothed_air_temp.tif", tavg_smoothed) # write file

########################
# Just air temperature #
########################

write("$datadir/climate_data/worldclim/air_temp.tif", tavg) # write file



###############################
##### Future climate data #####
###############################

## Mean is average of min and max. Downloads data directly from worldclim
for (model, ssp, time_period) in Iterators.product(models, SSPs, time_periods) 
    println("reading in for model $model, ssp $ssp, period $time_period")
    filename = "$datadir/climate_data/worldclim/future_temp/$(join([ssp, time_period, model], '_')).tif"

    if ~isfile(filename)
        rasters = [
            read(
                view(
                    Raster("https://geodata.ucdavis.edu/cmip6/2.5m/$model/$ssp/wc2.1_2.5m_$(join([var, model, ssp, time_period], '_')).tif", lazy = true),
                    X(-20 .. 60), Y(-40 .. 60)
                    )
                )
            for var in ["tmin", "tmax"]]

        tmean = broadcast(+, rasters...) ./ length(rasters)

        write(filename, tmean, force = true)
    else
        println("file exists, moving to next file")
    end
end