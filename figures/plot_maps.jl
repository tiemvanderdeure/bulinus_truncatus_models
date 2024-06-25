###########################
###### Plotting maps ######
###########################

# The code in this scripts generates figures 3, 4, 5 in the paper.

using Statistics, Rasters, CairoMakie, Colors, ArchGDAL, NCDatasets, Dates, GeoInterface, CSV, DataFrames

include("../globals.jl")

futures = Iterators.product(time_periods, SSPs)
indices = CartesianIndices(collect(futures))

limits = (-20, 60, -36, 60) # coordinates to plot

################
# Read in data #
################ 

############ mechanistic output
datasets = [mech_dir * join([ssp, time_period], '_') for (time_period, ssp) in futures]

mean_r0 = [Raster(dataset * "/mean.tif") for dataset in datasets]
median_r0 = [Raster(dataset * "/median.tif") for dataset in datasets]
prob_suitable = [Raster(dataset * "/prob_present.tif") for dataset in datasets]
std_r0 = [Raster(dataset * "/std.tif") for dataset in datasets]

median_r0_gcm = [Raster(dataset * "/median_gcm.nc") for dataset in datasets]
r0_std_gcm = [std(max.(m, Ref(0)); dims = 3)[:, :, 1] for m in median_r0_gcm]

median_r0_current = Raster(mech_dir * "current/median.tif")

mech_cutoff = 0. # cutoff value for presence/absence maps

########## Correlative output
cor_mean_set1 = Raster("$cor_dir/current_set1.tif")[Band(1)] ./ 1000
cor_mean_set2 = Raster("$cor_dir/current_set2.tif") ./ 1000

cor_futures = [Raster("$cor_dir/average_$(yr)$(ssp).tif") ./ 1000 for yr in [50, 90], ssp in [126, 370]]

cor_gcm_files = [readdir("$cor_dir/by_gcm/$(yr)$(ssp)"; join = true) for yr in [50, 90], ssp in [126, 370]]
cor_by_gcm = broadcast(cor_gcm_files) do files
    rs = RasterStack(files)[Rasters.Band(1)] # read in all files
    r = Raster(rs) # convert to a single raster, take mean (not CA)
    replace_missing(r, NaN) ./ 1000 # convert to 0-1
end

# Calculate std between GCMs
cor_std_gcm = [std(c, dims = 3)[:, :, 1] for c in cor_by_gcm]

# calculate cut-off value
presences = CSV.read("occurrence_data/thinned_presences.csv", DataFrame)
suitability_presences = reduce(vcat, [collect(cor_mean_set2[X = Near(p.x), Y = Near(p.y)]) for p in eachrow(presences)])
cor_cutoff = quantile(suitability_presences, 0.1) # cutoff value for presence/absence maps is 10th percentile

##### haematobium data
haem = Raster("$datadir/haem_2000_2010.tif")
haem_med = median(haem; dims = 3)[:, :, 1]

###################
# Main text plots #
###################
cmap = cgrad(:Spectral, rev = true)

##### Figure 3
begin
    fsize = 100
    fsize_labels = fsize * 0.8

    f = Figure(size = (6000, 2200), backgroundcolor = :transparent, fontsize = fsize)

    axes = [Axis(f[1, i], limits = limits, backgroundcolor = :transparent) for i in 1:3]
    for axis in axes hidedecorations!(axis); hidespines!(axis) end

    plot!(axes[1], median_r0_current, colormap = cmap, colorrange = (0,1))
    plot!(axes[2], cor_mean_set2, colormap = cmap, colorrange = (0,1))
    plot!(axes[3], median_r0_current, colormap = cgrad([:lightgrey, :lightgrey])) # grey background where there is no data
    plot!(axes[3], haem_med, colormap = cmap, nan_color = :transparent, colorrange = (0,0.5), highclip = :darkred)

    for (i, text) in enumerate(["Mechanistic model", "Correlative model", "S. haematobium prevalence"])
        Label(f[1,i, Top()], text =text)
    end
    for (i, lab) in enumerate(["A", "B", "C"])
        Label(f[1,i], text =lab, font = :bold, tellheight = false, tellwidth = false, halign = :left, valign = :top)
    end

    for i in 1:2
        Colorbar(f[1,i], vertical = true,
            colorrange = (0,1), colormap = cmap,
            ticksvisible = false, ticklabelsize = fsize_labels, labelsize = fsize_labels,
            tellheight = false, tellwidth = false, halign = 0.15, valign = 0.1, 
            height = 520, size = 64, ticklabelpad = 16,
            label = "Suitability")
    end

    Colorbar(f[1,3], vertical = true,
        colorrange = (0,0.5), colormap = cmap, highclip = :darkred,
        ticksvisible = false, ticklabelsize = fsize_labels, labelsize = fsize_labels,
        tellheight = false, tellwidth = false, halign = 0.15, valign = 0.14, 
        height = 460, size = 64, ticklabelpad = 16,
        label = "Prevalence")

    nodata = Legend(f[1,3],
        [PolyElement(color = :lightgrey)],
        ["No data"],
        framevisible = false,
        tellheight = false, tellwidth = false,
        labelsize = fsize_labels, halign = 0.18, valign = 0.03, patchsize = (64,64), patchlabelgap = 16,
    )

    colgap!(f.layout, 0)

    # Add a title
    Label(f[0,1:3], text = "B. truncatus suitability and S. haematobium prevalence", tellwidth = false, font = :bold)
    save("figures/figure_03.png", f)
end

####### Figure 4 - Future suitability for mechanistic and correlative model in columns
begin
    f = Figure(
        colormap = cmap, 
        fontsize = 25,
        size = (800, 1200),
        backgroundcolor = :transparent)
    maps = f[1,1] = GridLayout()
    axes = [Axis(maps[r,c], limits = limits, backgroundcolor = :transparent) for r in 1:4, c in 1:2]
    for ax in axes hidedecorations!(ax); hidespines!(ax) end

    for t in 1:2, s in 1:2
        row = t+(s-1)*2
        plot!(axes[row, 1], cor_futures[t, s])
        plot!(axes[row, 2], median_r0[t, s], colorrange = (0,1))

        Box(maps[row,0], color = (:lightgrey, 0.5), strokevisible = false)

        Label(
            maps[row,0], time_periods[t]; 
            rotation = pi/2, font = :bold, padding = (0,0,0,0), tellheight = false)
    end

    Label(maps[1,1,Top()], "Correlative"; font = :bold, padding = (0,0,5,0))
    Label(maps[1,2,Top()], "Mechanistic"; font = :bold, padding = (0,0,5,0))

    for i in 1:2
        r = i*2-1:i*2
        Box(maps[r,-1], color = (:lightgrey, 0.5), strokevisible = false)
        Label(maps[r,-1], SSPs[i]; rotation = pi/2, font = :bold)
    end

    colgap!(maps, 10); rowgap!(maps, 10)

    cb = Colorbar(f[1,2], colorrange = (0,1), label = "Suitability")

    save("figures/figure_04.png", f)
end

##### Figure 5 - correlative and mechanistic over each other - current and future
begin
    # Current
    mech_binary = median_r0_current .> mech_cutoff
    cor_binary = cor_mean_set2 .> cor_cutoff
    haem_binary = haem_med .> 0.1
    nas = isnan.(median_r0_current)

    # Future
    cor_future_binary = [c .> (cor_cutoff) for c in cor_futures]
    mech_future_binary = [c .> 0 for c in median_r0]

    #### Plot!!
    ## Set up colours
    neither_col = RGBA(0.75, 0.75, 0.75, 1.0)
    mech_col = RGBA(0.02, 0.58, 0.93, 1.0)
    cor_col = RGBA(0.99, 0.91, 0.06, 1.0)
    both_col = RGBA(0.5, 0.01, 0.3, 1.0)
    nan_col = RGBA(1,1,1,0)
    colors = [neither_col, mech_col, cor_col, both_col]

    haem_col = RGBA(0,0,0,0.5)# Makie.LinePattern(direction= [1, 1]; width = 1, tilesize=(5,5), linecolor = :black)

    function pick_col(cor, mech, nas, colors = colors, nan_col = nan_col) # function to go from binary rasters to RGB raster
        map = mech .+ cor .* 2 .+ 1 # assign 1 if both are 0, 2 if only mech, 3 if only cor, 4 if both
        colmap = getindex.([colors], map)
        colmap[nas] .= nan_col
        return colmap
    end

    # Maps from thresholded data
    current_map = pick_col(cor_binary, mech_binary, nas)
    future_maps = [pick_col(cor, mech, nas) for (mech, cor) in zip(mech_future_binary, cor_future_binary)]

    #### Suitability area calculations
    suitable_binary = mech_binary .&& cor_binary
    future_suitable = [mech .&& cor for (mech, cor) in zip(mech_future_binary, cor_future_binary)]
    #haem_binary = replace_missing(rasterize(haem_10; to = mech_binary, fill = 1), 0)
    cs = cellsize(suitable_binary)

    haem_rs = resample(haem_med, to = suitable_binary)
    haem_rs_binary = haem_rs .> 0.1
    suitable_size = sum(suitable_binary .* cs)
    round(suitable_size / 1e6; digits = 1)

    f_suitable_size = [sum(f_s .* cs) for f_s in future_suitable]
    round(f_suitable_size[2,2] / 1e6; digits = 1)

    pct_change_in_suitable_area = round.(Int64, (f_suitable_size ./ suitable_size .- 1) .* 100)


    sum((future_suitable[2,2] .&& suitable_binary) .* cs)
    sum((future_suitable[2,2] .&& .~suitable_binary) .* cs)
    sum((.~future_suitable[2,2] .&& suitable_binary) .* cs)

    sum((suitable_binary .&& haem_rs_binary) .* cs) # currently suitable with high haem
    # share of currently suitable with high haem that won't be suitable anymore under SSP370
    sum((.~future_suitable[2,2] .&& suitable_binary .&& haem_rs_binary) .* cs) / sum((suitable_binary .&& haem_rs_binary) .* cs)

    #### Make figure
    f = Figure(size =  (4800, 3200), fontsize = 88, backgroundcolor = :transparent)
    f_current = f[1,1]
    f_future = f[1,2] = GridLayout()

    ax_current = Axis(f_current[1,1], limits = limits, backgroundcolor = :transparent)
    axes_future = [Axis(f_future[i, j], limits = limits, backgroundcolor = :transparent) for i in 1:2, j in 1:2]

    heatmap!(ax_current, lookup(dims(current_map, (X, Y)))..., current_map)
    heatmap!(ax_current, haem_binary, colormap = cgrad([:transparent, haem_col]), colorrange = (0,1))

    for i in 1:2, j in 1:2
        heatmap!(axes_future[j,i], lookup(dims(future_maps[i,j], (X, Y)))..., future_maps[i,j])
        heatmap!(axes_future[j,i], haem_binary, colormap = cgrad([:transparent, haem_col]), colorrange = (0,1))
    #  Makie.poly!(axes_future[j,i], borders_shape, color = :transparent, strokecolor = :black, strokewidth = 0.1, overdraw = true)
        text!(axes_future[j,i], -10, -10; text = "+$(pct_change_in_suitable_area[i,j])%", color = both_col, font = :bold)
    end

    # Legend
    l = Legend(f_current[1,1],
        [[PolyElement(color = col) for col in colors], [PolyElement(color = haem_col)]],
        [["Unsuitable", "Mechanistic only", "Correlative only", "Suitable"], ["High prevalence"]],
        ["B. truncatus", "S. haematobium"],
        framevisible = false, patchsize = (80,80), patchlabelgap = 20,
        tellheight = false, tellwidth = false, halign = 0.05, valign = 0.1,
        gridshalign = :left, titlehalign = :left
    )  

    ## Labels
    Label(f_current[1, 1, Top()], text = "Current", font = :bold)

    for (i, tp) in enumerate(time_periods)
        Label(f_future[1, i, Top()], text = tp, font = :bold)
    end

    for (i, ssp) in enumerate(SSPs)
        Label(f_future[i, 2, Right()], text = ssp, font = :bold, rotation = -pi/2)
    end

    for ax in [ax_current; axes_future] hidedecorations!(ax); hidespines!(ax) end

    colgap!(f.layout, 0), rowgap!(f.layout, 10)
    colgap!(f_future, 0), rowgap!(f_future, 10)

    save("figures/figure_05.png", f)
end


#########################
# Supplementary figures #
#########################

####### Figure S8 - EMca scores correlative
begin
    EMcas = [Raster("$cor_dir/EMca_set$i.tif") for i in 1:2]

    f = Figure(size = (1000, 600), fontsize = 20, backgroundcolor = :transparent)

    label = ["Set 1", "Set 2"]
    ab = ["A", "B"]

    for i in 1:2
        ax = Axis(f[1,i], limits = limits, backgroundcolor = :transparent)
        hidedecorations!(ax); hidespines!(ax)
        plot!(ax, EMcas[i], colormap = cmap, colorrange = (0,1000))

        Colorbar(f[1,i], vertical = true,
            colorrange = (0,1), colormap = cmap,
            ticksvisible = false,
            tellheight = false, tellwidth = false, halign = 0.15, valign = 0.1, height = 130,
            label = "EMca score")

        Label(f[1, i, Top()], text = label[i], padding = (0,0,5,0))
        Label(f[1, i, TopLeft()], text = ab[i], font = :bold)
    end

    rowgap!(f.layout, 0)

    save("figures/figure_S08.png", f)
end


#### Figure S10 - std between GCMs 
begin
    m_range = max(
        Statistics.maximum(Statistics.maximum.(skipmissing.(cor_std_gcm))),
        Statistics.maximum(Statistics.maximum.(skipmissing.(r0_std_gcm))),
    )

    f = Figure(
        colormap = cmap, 
        fontsize = 25,
        size = (800, 1200), 
        backgroundcolor = :transparent)
    maps = f[1,1] = GridLayout()
    axes = [Axis(maps[r,c], limits = limits, backgroundcolor = :transparent) for r in 1:4, c in 1:2]
    for ax in axes hidedecorations!(ax); hidespines!(ax) end

    for t in 1:2, s in 1:2
        row = t+(s-1)*2
        plot!(axes[row, 1], cor_std_gcm[t, s], colorrange = (0,m_range))
        plot!(axes[row, 2], r0_std_gcm[t, s], colorrange = (0,m_range))

        Label(
            maps[row,1,Left()], time_periods[t]; 
            rotation = pi/2, font = :bold, padding = (0,-40,0,0))
    end

    Label(maps[1,1,Top()], "Correlative"; font = :bold, padding = (0,0,5,0))
    Label(maps[1,2,Top()], "Mechanistic"; font = :bold, padding = (0,0,5,0))

    [Label(maps[(i*2-1):(i*2),1,Left()], SSPs[i]; rotation = pi/2, font = :bold, padding = (0,30,0,0)) for i in 1:2]

    colgap!(maps, 10); rowgap!(maps, 10)

    cb = Colorbar(f[2,1], colorrange = (0,0.32), vertical = false, label = "Standard Deviation")

    save("figures/figure_S10.png", f)
end



#### Figure S11 - Output per GCM mechanistic
begin
    f = Figure(colormap = cmap, fontsize = 20, size = (1200, 800), backgroundcolor = :transparent)
    maps = f[1,1] = GridLayout()
    axes = [Axis(maps[r, c], backgroundcolor = :transparent) for r in 1:4, c in 1:7]
    hidedecorations!.(axes); hidespines!.(axes)

    for c in 1:7, r in 1:4
        plot!(axes[r,c], median_r0_gcm[r][Band = c], colorrange = (0,1))
    end

    colgap!(maps, 0); rowgap!(maps, 0)

    for c in 1:7 Label(maps[1,c,Top()], models[c]; font = :bold, padding = (0,0,5,0)) end

    [Label(maps[(i*2-1):(i*2),1,Left()], SSPs[i]; rotation = pi/2, font = :bold, padding = (0,30,0,0)) for i in 1:2]
    [Label(maps[i,1,Left()], time_periods[mod1(i, 2)]; rotation = pi/2, font = :bold, padding = (0,-30,0,0)) for i in 1:4]

    cb = Colorbar(f[1,2], colorrange = (0,1), label = "Suitability", labelpadding = 2.0)

    save("figures/figure_S11.png", f)
end


#### Figure S12 - Output per GCM Correlative
begin
    f = Figure(colormap = cmap, fontsize = 20, size = (1200, 800), backgroundcolor = :transparent)
    maps = f[1,1] = GridLayout()
    axes = [Axis(maps[r, c], backgroundcolor = :transparent) for r in 1:4, c in 1:7]
    hidedecorations!.(axes); hidespines!.(axes)

    for c in 1:7, r in 1:4
        plot!(axes[r,c], cor_by_gcm[r][Band = c], colorrange = (0,1))
    end

    colgap!(maps, 0); rowgap!(maps, 0)

    for c in 1:7 Label(maps[1,c,Top()], models[c]; font = :bold, padding = (0,0,5,0)) end

    [Label(maps[(i*2-1):(i*2),1,Left()], SSPs[i]; rotation = pi/2, font = :bold, padding = (0,30,0,0)) for i in 1:2]
    [Label(maps[i,1,Left()], time_periods[mod1(i, 2)]; rotation = pi/2, font = :bold, padding = (0,-30,0,0)) for i in 1:4]

    cb = Colorbar(f[1,2], colorrange = (0,1), label = "Suitability", labelpadding = 2.0)
    
    save("figures/figure_S12.png", f)
end
