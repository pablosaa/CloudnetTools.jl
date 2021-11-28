# Main CloudNet.jl module.
"""
A set of tools to process and analyze data outputs from CloudNet classification algorithm.

(c) 2020, Pablo Saavedra Garfias
University of Leipzig
Faculty of Physics and Geosciences
LIM

See LICENSE
"""


module CloudnetTools

using NCDatasets, DataStructures
using Dates
using Interpolations
using Printf
using Statistics
using ImageFiltering
using UUIDs
using ARMtools
using Plots

# Including transformation files:
include("arm_mwr2nc.jl")
include("arm_lidar2nc.jl")
include("arm_kazr2nc.jl")
include("arm_hsrl2nc.jl")


"""
Read time from netCDF file as fraction of the day and
convert it to Julia DateTime format.
USAGE:
> nctime = convert_time2DateTime(year, month, day, hour_of_day)
WHERE:
* year, month and day -> Float32
* hour_of_day -> Vector{Float32} with fraction of day [0 to 24]

or
> nctime = convert_time2DateTime(nc::NCDataset)
WHERE:
* nc -> NCDataset identifier for an opened netCDF file,
* nctime -> Vector{DataTime}

"""
function convert_time2DateTime(yy::Int64, mm::Int64, dd::Int64, hr_time::AbstractVector)::Vector{DateTime}
    hh = floor.(Int64, hr_time)
    mi = @. mod(hr_time*60, 60)
    ss = @. mod(mi*60, 60)
    ms = @. mod(ss, 1)*1000
        
    mi = floor.(Int64, mi)
    ss = floor.(Int64, ss)
    ms = floor.(Int64, ms)
        
    return @. DateTime(yy, mm, dd, hh, mi, ss, ms)
end
function convert_time2DateTime(nc; time_var="time")::Vector{DateTime}

    eltype(nc[time_var]) <: DateTime && (return nc[time_var])
    yy = Int64(nc.attrib["year"])
    mm = Int64(nc.attrib["month"])
    dd = Int64(nc.attrib["day"])
    hr_time = float.(nc[time_var])

    return convert_time2DateTime(yy, mm, dd, hr_time)
end
# ----/

# ***************************************************************
# Read Cloudnet Ice water content IWC file
"""
# Read Cloudnet Ice water content IWC file

> iwc = readIWCFile(nfile::String; modelreso=false)

"""
function readIWCFile(nfile::String; modelreso=false)

    @assert isfile(nfile) error("$nfile cannot be found!")

    vars_categorize = Dict(
        :height => "height",
        :iwc => "iwc_inc_rain",
        :flag => "iwc_retrieval_status",
    )

    # Defining Output varaible:
    var_output = Dict{Symbol, Any}()

    NCDataset(nfile; format=:netcdf4_classic) do nc
        var_output[:time] = convert_time2DateTime(nc)

        for inkey ∈ keys(vars_categorize)
            x = vars_categorize[inkey]
            println(x)
            tmp = nc[x][:,:]
            if haskey(nc[x].attrib, "missing_value")
                miss_val = nc[x].attrib["missing_value"]
            elseif haskey(nc[x].attrib, "_FillValue")
                miss_val = nc[x].attrib["_FillValue"]
            else
                miss_val = 9.96921f36
            end

            # Cleaning missing values from variables :
            eltype(tmp) <: AbstractFloat && (tmp[tmp .≈ miss_val] .= NaN)
            
            var_output[inkey] = tmp
        end

    end

    # Adding computed variables:
    # integration of iwc only for pixels with flag=1 or 2:
    # Cloudnetpy units iwc [kg m⁻³] and height [m]
    var_output[:IWP] = let tmp = 1f3var_output[:iwc]

        @. tmp[!(0 < var_output[:flag] < 3)] = NaN
        # IWP [g m⁻²]
        ∫fdh(tmp, var_output[:height])
    end
    
    return var_output
end
# ----/

# ***************************************************************
# Read Liquid water content LWC file
"""
# Read Liquid water content LWC file

> lwc = readLWCFile(nfile::String; modelreso=false)

"""
function readLWCFile(nfile::String; modelreso=false)

    @assert isfile(nfile) error("$nfile cannot be found!")

    vars_categorize = Dict(
        :height => "height",
        :lwc => "lwc",
        :flag => "lwc_retrieval_status",
        :LWP => "lwp",
    )

    # Defining Output varaible:
    var_output = Dict{Symbol, Any}()

    NCDataset(nfile; format=:netcdf4_classic) do nc
        var_output[:time] = convert_time2DateTime(nc)

        for inkey ∈ keys(vars_categorize)
            x = vars_categorize[inkey]
            println(x)
            tmp = nc[x][:,:]
            if haskey(nc[x].attrib, "missing_value")
                miss_val = nc[x].attrib["missing_value"]
            elseif haskey(nc[x].attrib, "_FillValue")
                miss_val = nc[x].attrib["_FillValue"]
            else
                miss_val = 9.96921f36
            end

            # Cleaning missing values from variables :
            eltype(tmp) <: AbstractFloat && (tmp[tmp .≈ miss_val] .= NaN)
            
            var_output[inkey] = tmp
        end

    end

    # Additional computation:
    # integration of iwc only for pixels with flag=1 or 2:
    @. var_output[:lwc][!(0 < var_output[:flag] < 3)] = NaN
    
    return var_output
end
# ----/

# ******************************************************************
# Reading Classification & Categorize files both at once:
"""
# Fuction to read Cloudnet classification or category files:

> classi = readCLNFile(nfile::String; modelreso=false)

if optional parameter modelreso=true, then the model variables are
interpolated from hourly resolution to Cloudnet resolution.
"""
function readCLNFile(nfile::String; modelreso=false)
    @assert isfile(nfile) error("$nfile cannot be found!")
    if contains(nfile, "categorize")
        println("reading Categorize file")
    elseif contains(nfile, "classific")
        println("reading Classification file")
    else
        error("$nfile does not apear to be a Cloudnet file!")
    end

    # Categorize variables to read:
    vars_categorize = Dict(
        # RADAR
        :Z => "Z",
        :V => "v",
        :σV => "v_sigma",
        :ωV => "width",
        # LIDAR
        :β => "beta",
        :δ => "depol",
        #MWR
        :LWP => "lwp",
        #CLOUDNETpy
        :P_INSECT => "insect_prob",
        :QUALBITS => "quality_bits",
        :CATEBITS => "category_bits",
        #MODEL
        :model_height => "model_height",
        :T => "temperature",
        :Tw => "Tw",
        :Pa => "pressure",
        :QV => "q",
        :UWIND => "uwind",
        :VWIND => "vwind",
    )


    # Defining Output varaible:
    var_output = Dict{Symbol, Any}()

    # Starting reading CloudNet files:
    NCDataset(nfile; format=:netcdf4_classic) do nc

        var_output[:time] = convert_time2DateTime(nc)

        var_output[:height] = nc["height"][:]

      
        for (inkey, x) ∈ vars_categorize
            #x = vars_categorize[inkey]
            !haskey(nc, x) && continue

            println(x)
            tmp = nc[x][:,:]
            if haskey(nc[x].attrib, "missing_value")
                miss_val = nc[x].attrib["missing_value"]
            elseif haskey(nc[x].attrib, "_FillValue")
                miss_val = nc[x].attrib["_FillValue"]
            else
                miss_val = 9.96921f36
            end

            # Cleaning missing values from variables :
            eltype(tmp) <: AbstractFloat && (tmp[tmp .≈ miss_val] .= NaN)

            var_output[inkey] = tmp
        end

        # If modelreso = true, interpolate model data to cloudnet resolution:

        haskey(nc, "model_time") && (var_output[:model_time] = convert_time2DateTime(nc, time_var="model_time"))

        #(model_time = float.(nc["model_time"]))
        #model_time = float.(nc["model_time"])
        #model_height = nc["model_height"][:]
        if !modelreso
            
            #var_output[:model_height] = model_height
        else
            # To convert model_time to cloudnet time, model_time needs to be
            # a Vector with elements containing the fraction of day:
            model_time = let x = haskey(var_output, :model_time)
                x ? var_output[:model_time] : var_output[:time]
            end
            
            ConvertModelResolution(var_output,
                                   model_time,
                                   var_output[:model_height])
            #    cln_time=tmp_time)
        end

    end

        
    # For Classification dataset:
    # Classification variables to read
    vars_classific = Dict(
        :CLASSIFY => "target_classification",
        :DETECTST =>"detection_status",
    )

    classfile = replace(nfile, "categorize" => "classific")

    if isfile(classfile)
    
        NCDataset(classfile; format=:netcdf4_classic) do nc
            [var_output[k] = nc[v][:,:] for (k,v) ∈ vars_classific if haskey(nc, v)]
            #[var_output[x] = nc[vars_classific[x]][:,:] for x ∈ keys(vars_classific)];
        
        end
    else
        @warn("$(classfile) cannot be found and not loaded!")
    end
    
    return var_output;

end  # end of function
# ----/

# *********************************************************
# Interpolate Meteo data from Model to CloudNet resolution
"""
# Interpolate Meteo data from Model to CloudNet resolution.

> var_out = ConvertModelResolution(cln_in::Dict{Symbol, Any},
                                   model_time::Vector{Float32},
                                   model_height::Vector{Float32};
                                   cln_time=nothing,
                                   cln_height=nothing)

the output containg the model variables :T, :Pa, :UWIND, :VWIND, :QV
but at the same resolution as cloudnet.

"""
function ConvertModelResolution(cln_in::Dict{Symbol, Any},
                                model_time::Vector{Any},
                                model_height::Vector{Real};
                                cln_time=nothing,
                                cln_height=nothing)


    # creating modes for interpolation depending on typeof model_time:
    if eltype(model_time) <: DateTime
        model_ts = let ts = model_time
            @. hours(ts) + minute(ts)/60 + seconds(ts)/3600
        end
        nodel = (model_height, model_ts)
    else
        nodes = (model_height, model_time)
    end
        
    if isnothing(cln_time)
        let ts = cln_in[:time]
            cln_time = @. hour(ts) + minute(ts)/60 + seconds(ts)/3600
        end
    end
    
    if isnothing(cln_height)
        cln_height = cln_in[:height]
    end

    MODELVAR = (:T, :Pa, :UWIND, :VWIND, :QV)

    map(MODELVAR) do var
        itp = interpolate(nodes, cln_in[var], Gridded(Linear()))
        outvar = [itp(i,j) for i ∈ cln_height, j ∈ cln_time]
        cln_in[var] = outvar
    end

    return cln_in
end
# ----/

# ------------------------------------------------------
# AUX FUNCTIONS
# ------------------------------------------------------

# *******************************************************
# Integration of variable (height, time) over altitudes (dims=1)
"""
Integration of matrix x over dims=1 skipping NaN, optionally
the integration is perform with respect of vector dh.
USAGE:

Xt = ∫fdh(xi::Matrix)
Xt = ∫fdh(xi::Matrix; dh::Vector)

"""
function ∫fdh(x::AbstractArray, dx::AbstractVector)
    δx = Vector{eltype(x)}(undef, length(dx)) .= 0
    δx[2:end] = dx |> diff
    δx[1] =  δx[2:end] |> minimum

    x[isnan.(x)] .= 0.0
    
    return sum(x.*δx, dims=1) |> vec
end
#----/

# ***************************************************
# Function to normalize matrix to [0,1]
# in case the input variable has log units (e.g. dBz),
# then the optional parameter could be isdB=true to treat
# input array as a linear variable, the output is based on
# the non-logarithm variable.
#
function TransformZeroOne(X::T; isdB=false) where T<:AbstractArray
    
    nonans = .!isnan.(X)
    if isdB
        X = @. 10^(X/10)
    end
    x₀, x₁ = extrema(X[nonans])
    Xout = (X .- x₀)./(x₁ - x₀);

    return Xout
end  # end of function
# ----/

# ***************************************************
# Function to mimic meshgrid MATLAB
function mymeshgrid(x::Vector, y::Vector)
    nx = length(x)
    ny = length(y)

    return ones(ny).*x', y.*ones(nx)'
end
# ----/

end  # end of Module CloudnetTools

# ******************************************************************
# +++++++++++++++ MODULE CloudnetViz  ++++++++++++++++++++++++++++++
# ******************************************************************

module CloudnetViz

using Plots
using Dates
using Printf
using LaTeXStrings

function show_LWP_IWP(iwc::Dict, liq::Dict)
    p1=plot(iwc[:time], iwc[:height], log10.(1f3iwc[:iwc]), st=:heatmap, color=:lapaz);
    p2=plot(liq[:time], liq[:height], log10.(liq[:lwc]), st=:heatmap, color=:tokyo);
    plot(p1,p2, layout=(2,1));
end

# *********************************************************************
# Function to plot the Classification data
"""
# Function to plot the Classification data

> show_classific(cnt::Dict; 

WHERE:
* cnt::Dict dictionary ouput from read_CNTfile(cloudnet_file),
* SITENANE::String (optional) string with name of site,
* maxhgt::Number (optional) indicating the maximum height in km, default=8,
* showlegend::Bool (optional) show Cloudnet legend colors, default=true
* addons::Bool to add meteo data to the plot, "wind" & "temp", default=true
Output:
* plt::Plot output plot object.
"""
function show_classific(cnt::Dict; SITENAME="", maxhgt=8, showlegend=true)

    # defining time axis ticks:
    tm_tick = cnt[:time][1]:Minute(90):cnt[:time][end];

    Xstrname = "TIME UTC from "*Dates.format(cnt[:time][2], "dd.u.yyyy")

    strtitle = (tm_tick[1], 7.5, text("CloudNet Target Classification "*SITENAME, 11, halign=:left))
    strtitle = "CloudNet Target Classification "*SITENAME
    
    cldnet = CloudNetPalette("classific")
    l = grid(2,1, heights=(0.2,0.8)) #[a{0.25h}; b];
    Y_LIM = (0, maxhgt)
    
    classplt = plot(cnt[:time], 1f-3cnt[:height], cnt[:CLASSIFY],
                    st=:heatmap, colorbar = false, framestyle = :box,
                    color=palette(cldnet, 11), ylim=Y_LIM, clim=(0,10),
                    ylabel="Height A.G.L. [km]", ytickfontsize=11, minorticks=true,
                    xlabel= Xstrname, tick_direction=:out, 
                    xticks=(tm_tick, Dates.format.(tm_tick, "H:MM")), xrot=45,
                    xtickfontsize=11, xguidefont=font(12), title=strtitle);

    # Preparing to add graphs for isotherms and wind vectors:
    ihmax = findlast(1f-3cnt[:model_height] .≤ maxhgt)

    
    Xin = haskey(cnt, :model_time) ? Vector(1:length(cnt[:model_time])) : round.(Int64, range(1, stop=length(cnt[:time]), length=25))
    
    Yin = 1f-3cnt[:model_height]
    #TLEV = [5, 0, -5, -10, -15, -20, -25, -30]
    BB = bbox(0,0,1,1)
    
    Attach_Isotherms(classplt, Xin, Yin[1:ihmax], cnt[:T][1:ihmax, Xin] .- 273.15,
                     (1, BB), 2, TLEV = [5, 0, -5, -10, -15, -20, -25, -30])
    

    Attach_Windvector(classplt, Xin[1:2:end], Yin[2:4:ihmax],
                      cnt[:UWIND][2:4:ihmax, Xin[1:2:end]],
                      cnt[:VWIND][2:4:ihmax, Xin[1:2:end]],
                      (1, BB), 3)

    #
    #θv = @. cnt[:T][1:ihmax,2:2:end]*(1024f0/cnt[:Pa][1:ihmax,2:2:end])^0.286;
    
    #Attach_Profile_Cascate(classplt, Xin[2:2:end], Yin[1:ihmax], θv, (1,BB), 4)
    
    if showlegend
        classleg = ShowLegendCloudNetClassification("classific")

        pltout = plot(classleg, classplt, layout=l, size=(800,700))
    else
        pltout = plot(classplt, size=(800,600))
    end

    return pltout
end

function CloudNetPalette(ColorType::String)
    if ColorType == "classific"
        tmp = [
            (1.00,  1.0 ,  1.00);
            (0.44,  1.0 ,  0.92);
            (0.17,  0.62,  0.95);
            (0.75,  0.6 ,  1.0);
            (0.9 ,  0.9 ,  0.92);
            (0.28,  0.27,  0.72);
            (0.99,  0.65,  0.0);
            (0.78,  0.98,  0.19);
            (0.8 ,  0.73,  0.53);
            (0.89,  0.29,  0.13);
            (0.7 ,  0.21,  0.34)
        ];
    elseif ColorType == "detection"
        tmp = [
            (0.247,  0.704,  0.43);
            (0.44 ,  0.926,  0.34); 
            (0.996,  0.91 ,  0.24);
            (0.8  ,  0.96 ,  0.96); 
            (0.455,  0.51 ,  0.41); 
            (0.89 ,  0.292,  0.13);
        ];
    else
        error("$ColorType not supported!")
    end
    cldnet = map(x->RGB(x...), tmp);
    return cldnet
end

function ShowLegendCloudNetClassification(LegendType::String)
    if LegendType == "classific"
        txt_labels = [
            #X, Y, "short_name", "long_name"
            (40, 1,"Clear sky", "Clear sky"),
            (0, 1,"Liquid drops", "Cloud liquid droplets only"), 
            (20, 1,"Drizzle | Rain", "Drizzle or rain"), 
            (20, 2,"Drizzle & Cloud", "Drizzle or rain coexisting with cloud liquid droplets"),
            (0, 3,"Ice", "Ice particles"), 
            (0, 2, "SLC & Ice", "Ice coexisting with supercooled liquid droplets"),
            (20, 3,"Melting", "Melting ice particle"),
            (20, 4,"Melting & cloud", "Melting ice particles coexisting with cloud liquid droplets"), 
            (40, 2,"Aerosol", "Aerosol particles/no cloud or precipitation"), 
            (40, 3,"Insects", "Insects/no cloud or precipitation"), 
            (40, 4,"Aerosol & Insects", "Aerosol coexisting with insects/no cloud or precipitation")
        ];
    elseif LegendType == "detection"
        txt_labels = [
            (1, 1, "Ra + Li", "Good Radar & Lidar");
            (1, 2, "Ra", "Good Radar Echo"); 
            (4, 1, "Li", "Lidar Echo Only");
            (4, 2, "Ra + Att", "Radar corrected \nfor liquid attenuation");
            (7, 1, "Ra - Att", "Radar uncorrected \nfor liquid attenuation");
            (7, 2, "Ra clutter", "Radar ground clutter")
        ];
    else
        error("$LegendType not supported!")
    end 
    cldnet = CloudNetPalette(LegendType)

    scatter(map(x->(x[1], 1.5x[2]), txt_labels),
            marker=:square, grid=:false, markerstroke=:none,
            framestyle=:none, widen=true,
            markersize=8, markercolor=cldnet, legend=false, axis=false,
            clipping=:false, xlim=(-0.1, 60), ylim=(0, 8)); #1+maximum(txt_labels)[2]));
    
    tmp = map(x->(x[1]+1.5, 1.5x[2], text(x[3], 11, :left)), txt_labels) |> annotate!;

    return tmp
end
# ----/

function Attach_Isotherms(pltin, Xin, Yin, Zvar, BB, sp; maxhgt=(0,8), TLEV=10)
    #TLEV = [5, 0, -5, -10, -15, -20, -25, -30]
    
    Plots.contour!(pltin, Xin, Yin, Zvar,
                   levels = TLEV, ylim=maxhgt, xlim = extrema(Xin), ticks=:none,
                   linecolor=:black, alpha=:.5, lw=1, contour_labels = true,
                   axis=false, colorbar=:none, subplot=sp, inset=BB, 
                   background_color_subplot=:transparent);

    return pltin
end
# ----/

# **********************************************************
# Subroutine to attach the wind vector to plot
function Attach_Windvector(pltin, Xin, Yin, Uin, Vin, BB, sp; maxhgt=(0,8))

    TT, HH, UU, VV = prepare_quiver(Xin, Yin, Uin, Vin)
    
    Plots.quiver!(pltin, TT[:], HH[:], quiver=(UU[:],VV[:]),
                  ylim=maxhgt, axis=false, color=:gray, alpha=0.3,
                  ticks=:none, shape=:none, inset=BB, subplot=sp,
                  background_color_subplot=:transparent);

    return pltin
end
# ----/                                
                                
# ***********************************************************
# Subroutine to create cascate plot of profile variables
function Attach_Profile_Cascate(pltin, ts, hkm, Z, BB, sp; maxhgt=(0, 8), δts=1.5)

    Znn = normalize_2Dvar(Z, δx=δts);
    Plots.plot!(pltin, Znn .+ ts', hkm, xlim=extrema(ts), ylim=maxhgt,
                color=:lightblue, legend=false, axis=:none, ticks=:none,
                inset=BB, background_color_subplot=:transparent, subplot=sp)

    return pltin
end
# ----/

# ************************************************************
# ++++++++++ AUXILIARY FUNCTIONS FOR VISUALIZATION +++++++++++
# ************************************************************

# ************************************************************
# Prepare vectors to be used with quiver
using ..CloudnetTools

function prepare_quiver(x::Vector, y::Vector, u::Matrix, v::Matrix; factor=2)

    # 2D meshing of coordinates:
    X, Y = CloudnetTools.mymeshgrid(x, y)

    # normalizing the wind variables
    UU, VV = let WS = @. factor*sqrt(u^2 + v^2)
    
        u./WS, v./WS
    end

    return X, Y, UU, VV
    
end
# ----/

# ****************
function normalize_2Dvar(X::Matrix; δx=1.5)

    δx < 1 && error("Optional input δx needs to be ≥ 1")
    T01 = [];
    for T ∈ eachcol(X)
        let lims = extrema(T)
            (2(T .- lims[1])./(lims[2] - lims[1]) .- 1)δx |> x->push!(T01, x)
            #push!(T01, (T .- lims[1])./(lims[2] - lims[1]) |> v-> (2δx.*v) .- δx)
        end
    end
    return hcat(T01...)
end
# ----/


end  # end of Module CloudnetViz

# --end of script
