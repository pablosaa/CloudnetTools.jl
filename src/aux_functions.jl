#=
********************************************************************
File contains auxiliary functions for the CloudnetTools.jl package
********************************************************************
=#

# ****
"""
Function to return a String containing the cloudnet site information
based on information from the data file or extra given to the function:
USAGE:
> get_SITE(data, extra_params)

returns: "mosaic-arm"
if the data or extra_params contain the keys: [:site] = "arm", [:campaign]="mosaic"

"""
function get_SITE(data::Dict, extra_params::Dict; inkeys=(:site, :campaing))
    # checking if inkeys exist in extra_params or data:
    tmpkeys = filter(k->haskey(extra_params, k), inkeys)
    
    buffer = isempty(tmpkeys) ? data : extra_params

    isempty(tmpkeys) && (tmpkeys = filter(k->haskey(buffer, k), inkeys))
        
    return getindex.((buffer,), tmpkeys) |> x->join( filter(!isempty, x), "-")
end
# /----

# ****
"""
Returm given parameter from either data input Dict or extra_params input Dict,
if not present in any of them, then return a default
USAGE:
> get_parameter(lidar, :λ_nm, input_params, default=510, lims=(450, 1005))
> 910
if lidar[:λ_nm] has the value of 910 and that value belongs to lim. In case lidar[:λ_nm] does not exist, then the fiven default value of 510 is returned.

The function priority return values are from: input_params (if exist), data Dict (if exist), or default (if given), otherwise returns nothing.
"""
function get_parameter(data::Dict, param::Symbol, extra_params::Dict{Symbol, Any};
                       default=nothing, lims=())

    out_value = if haskey(extra_params, param)
        extra_params[param]
    elseif haskey(data, param)
        data[param]
    else
        default
    end

    if !isempty(lims) && !isnothing(out_value)
        return check_limits(out_value, lims) ? out_value : default
    else
        return out_value
    end
    
end
# ----/

# ******************************************
"""
Function to check whether the data is between the given limits:
USAGE:
> check_limits(mwr[:lat], (-90, 90))
> true
If latitude data in mwr is correct, otherwise returns false.

"""
function check_limits(data, limits::Tuple{Any, Any})
    return all(limits[1] .≤ extrema(data) .≤ limits[2])
end
# ----/

# ****************************************
"""
Function to convert a vector of DateTime to the fraction of day as Real numbers:
USAGE:
> hourofday = datatime24hours(nc[:time])
WHERE:
* nc[:time] is a Vector{Datetime}
* hourofday is the fraciton of the 24hr day, e.g. 12:30:00 is 12.5

"""
function datetime24hours(time_in::Vector{DateTime})
    return hour.(time_in) + minute.(time_in)/60. + second.(time_in)/3600.0;
end
# ----/

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

    hr_time = nc[time_var]
    typeof(hr_time[1]) <: DateTime && (return hr_time[:])
    yy = Int64(nc.attrib["year"])
    mm = Int64(nc.attrib["month"])
    dd = Int64(nc.attrib["day"])
    hr_time = float.(nc[time_var])

    return convert_time2DateTime(yy, mm, dd, hr_time)
end
# ----/

#=
**** MATHEMATICAL FUNCTIONS AND HELPERS *******
=#

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

# *********************************************************
# Interpolate Meteo data from Model to CloudNet resolution
"""
# Interpolate Meteorological data from Model to CloudNet resolution.

> var_out = ConvertModelResolution(cln_in::Dict{Symbol, Any},
                                   model_time::Vector{Float32},
                                   model_height::Vector{Float32};
                                   cln_time=nothing,
                                   cln_height=nothing)

the output containg the model variables :T, :Pa, :UWIND, :VWIND, :QV
but at the same resolution as cloudnet.

"""
function ConvertModelResolution(cln_in::Dict{Symbol, Any},
                                model_time::Vector{<:Any},
                                model_height::Vector{<:Real};
                                cln_time=nothing,
                                cln_height=nothing)


    # creating modes for interpolation depending on typeof model_time:
    function day_fraction(time_in::Vector{DateTime})
        time_out = let ts = time_in
            tmp = datetime24hours(ts)
            tmp ./= 24
            tmp .+= day.(ts)
        end
        return time_out
    end

    if eltype(model_time) <: DateTime
        model_ts = day_fraction(model_time)
        
        nodes = (model_height, model_ts)
    else
        nodes = (model_height, model_time)
    end
        
    if isnothing(cln_time)
        cln_time = day_fraction(cln_in[:time])
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


# end of file
