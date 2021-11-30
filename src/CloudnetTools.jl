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
using Plots


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

# ****************************************************
# Including further files:
# * For plotting CloudNet products:
Base.include(CloudnetTools, "CloudnetVisualization.jl")
Base.include(CloudnetTools, "Cloudnet_for_ARM.jl")

# --end of script
