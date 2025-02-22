#=
********************************************************************
File contains auxiliary functions for the CloudnetTools.jl package

(c) 2022 P. Saavedra Garfias
University of Leipzig

********************************************************************
=#
using ImageFiltering
using Statistics

# ****
"""
Function to return a String containing the cloudnet site information
based on information from the data file or extra given to the function:
USAGE:
```julia-repl
julia> extra_params = Dict(:site => "arm", :campaign => "mosaic", :location=>"arctic");
julia> get_SITE(data, extra_params)
julia> "mosaic-arm"
julia> get_SITE(data, extra_params; inkeys=(:location, :site))
julia> "arctic-arm"
```
WHERE:
* ```data::Dict```
* ```extra_params::Dict``` contain the keys: ```:site = "arm", [:campaign]="mosaic"```
* ```inkeys::Tuple(Symbol,)``` Optional keys to use, default ```(:site, :campaign)```

the function has priority on ```extra_params```, if empty then try with ```data```.

Part of ```CloudnetTools.jl```. See LICENSE.TXT
"""
function get_SITE(data::Dict, extra_params::Dict; inkeys=(:site, :campaign))
    # checking if inkeys exist in extra_params or data:

    thesite = ""
    for indat ∈ (extra_params, data)
        tmpkeys = filter(k->haskey(indat, k), inkeys)
        isempty(tmpkeys) && continue
        theinputs = getindex.((indat,), tmpkeys)
        thesite = join( filter(!isempty, theinputs), "-")
        thesite = replace(thesite, ","=>"_", " "=>"")
        break
    end
            
    return thesite
end
# /----

# ****
"""
Returm given parameter from either data input Dict or extra_params input Dict,
if not present in any of them, then return a default
USAGE:
```julia-repl
julia> get_parameter(lidar, :λ_nm, input_params, default=510, lims=(450, 1005))
julia> 910
```
if lidar[:λ_nm] has the value of 910 and that value belongs to lim. In case lidar[:λ_nm] does not exist, then the fiven default value of 510 is returned.

The function priority return values are from: input_params (if exist), data Dict (if exist), or default (if given), otherwise returns nothing.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function get_parameter(data::Dict, param::Symbol, extra_params::Dict;
                       default=nothing, lims=(), silent=true)

    out_value = if haskey(extra_params, param)
        extra_params[param]
    elseif haskey(data, param)
        data[param]
    else
        !silent && @warn("$(param) not found in data/extra_params, returning default=$(default)")
        default
    end

    if !isempty(lims) && !isnothing(out_value)
        
        if check_limits(out_value, lims)
            return out_value
        elseif isnan(out_value)
            !silent && @warn("get_parameter found $(out_value) to be $(param)")
            return out_value
        else
            !silent && @warn("get_parameter found $(param) out of lims=$(lims)")
            return out_value
        end
    else
        return out_value
    end
    
end
# ----/

# ******************************************
"""
Function to check whether the data is between the given limits:
USAGE:
```julia-repl
julia> check_limits(mwr[:lat], (-90, 90))
julia> true
```
If latitude data in mwr is correct, otherwise returns false.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function check_limits(data, limits::Tuple{Any, Any})
    return all(limits[1] .≤ extrema(data) .≤ limits[2])
end
# ----/

# ****************************************
"""
Function to convert a vector of DateTime to the fraction of day as Real numbers:

USAGE:
```julia-repl
julia> hourofday = datatime24hours(nc[:time])
```
WHERE:
* nc[:time]::Vector{DateTime} or
* nc[:time]::Datetime

OUTPUT:
* hourofday is the fraciton of the 24hr day, e.g. 12:30:00 is 12.5

NOTE: If the input DateTime vector contains multiple days, 24 is added to the
subsequent days, e.g. 01:30:00 of next days will output 25.5

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function datetime24hours(time_in::DateTime)
    return hour.(time_in) + minute.(time_in)/60. + (second.(time_in) .+ 1f-3millisecond.(time_in))/3600.0;
end
function datetime24hours(time_in::Vector{DateTime})
    days_in_date = day.(time_in)
    unique(days_in_date) |> x-> length(x)>1 && @warn "Multiple days in DateTime vector to convert as a fraction of 24 hours!"
    return 24(day.(time_in) .- day.(time_in[1])) .+ datetime24hours.(time_in)
end
# ----/

"""
Convert given year, month, day, hour_day to Julia DateTime format.

When a netCDF NCDataset is given, the variable time is converted it to Julia DateTime format.
This is useful when the variable time in the netCDF is given as a fraction of day.
USAGE:
```julia-repl
julia> nctime = convert_time2DateTime(year, month, day, hour_of_day)
julia> nctime = convert_time2DateTime(nc; time_var="other_time_variable")
```
WHERE:
* year, month and day -> ```Float32```
* hour_of_day -> ```Vector{Float32}``` with fraction of day [0 to 24]
* ```nc::NCDataset``` identifier for the netCDF file to read.
* ```time_var::String``` (Optional) given the netCDF variable to read instead.
RETURN:
* ```nctime::DateTime``` for ```::Vector{DateTime}```

EXAMPLE:
```julia-repl
julia> nc = NCDataset("/data/mwr.nc", "r");
julia> nctime = convert_time2DateTime(nc; time_var="time_hr");
```
Part of ```CloudnetTools.jl```, see LICENSE.TXT
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

    hr_time = nc[time_var][:]
    typeof(hr_time) <: Vector{DateTime} && (return hr_time)

    yy = Int64(nc.attrib["year"])
    mm = Int64(nc.attrib["month"])
    dd = Int64(nc.attrib["day"])
    hr_time = float.(hr_time)

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
```julia-repl
julia> Xt = ∫fdh(Xi::Matrix, H::Vector)
julia> Xt = ∫fdh(Xi::Matrix, H::Vector; h₀=CBH, hₜ=CBT)
```
WHERE:
* ```Xi::Matrix``` - to integrate over 1st dimension,
* ```H::Vector```  - with the integrating variable, same length as Xi 1st dimension,
* ```h₀::Vector``` - (Optional) with low limit height to integrate from,
* ```hₜ::Vector``` - (Opitonal) with top limit height to integrate to.

If neither h₀ nor hₜ are provided, then the integral is performed over whole Xi 1st
dimension. Note that h₀ and hₜ, if provided, need to have same length as Xi 2nd dimension.

OUTPUT:
* ```Xt::Vector``` - with the integrated value, same length as Xi 2nd dimension.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function ∫fdh(x::AbstractArray, H::AbstractVector; h₀=Real[], hₜ=Real[])
    # getting dimensions:
    nheight, ntime = size(x)

    # checking dimensions:
    nheight != length(H) && @error "Heights muss have same length as matrix 1st dim."
    !isempty(h₀) && ntime != length(h₀) && @error "Optional variable h₀ must have length same as matrix 2nd dim."
    !isempty(hₜ) && ntime != length(hₜ) && @error "Optional variable hₜ must have length same as matrix 2nd dim."
	
    # defining the derivate of H
    δh = Vector{eltype(x)}(undef, length(H)) .= 0
    δh[2:end] = diff(H)
    δh[1] =  δh[2:end] |> minimum

    # integrating variable within limits:
    𝐼₀ₜ = fill(NaN32, ntime)
    for (i, X) ∈ enumerate(eachcol(x))
        # finding the limits for integration:
        i0 = isempty(h₀) ? 1 : !isnan(h₀[i]) ? argmin(abs.(H .- h₀[i])) : nothing
        it = isempty(hₜ) ? nheight : !isnan(hₜ[i]) ? argmin(abs.(H .- hₜ[i])) : nothing

        isnothing(i0) && continue
        isnothing(it) && continue

	lims = i0:it

        Xh = let Xh=X[lims] 
            tmp = ismissing.(Xh) .|| isnan.(Xh)
            all(tmp) && continue
            Xh[tmp] .= 0
            Xh
        end

        #X = (collect∘skipmissing)(X)
        
	𝐼₀ₜ[i] = Xh'*vec(δh[lims]) 
        
    end
    return 𝐼₀ₜ
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
```julia-repl
julia> var_out = ConvertModelResolution(cln_in::Dict{Symbol, Any},
                                   model_time::Vector{Float32},
                                   model_height::Vector{Float32};
                                   cln_time=nothing,
                                   cln_height=nothing)
```
the output containg the model variables :T, :Pa, :UWIND, :VWIND, :QV
but at the same resolution as cloudnet.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
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
        if haskey(cln_in, var)
            itp = interpolate(nodes, cln_in[var], Gridded(Linear()))
            outvar = [itp(i,j) for i ∈ cln_height, j ∈ cln_time]
            cln_in[var] = outvar
        end
    end

    return cln_in
end
# ----/

#=
**** ATMOSPHERIC RELATED FUNCTIONS  *******
=#
"""
Estimate base and top height of up to 3 cloud layers base on Cloudnet classification.

USAGE:
```julia-repl
julia> CBH, CTH = estimate_cloud_layers(clnet)
julia> CBH, CTH = estimate_cloud_layers(clnet; lidar=ceilometer, nlayers=2)
julia> CBH, CTH, CLB = estimate_cloud_layers(clnet; nlayers=2, liquid_base=true)
```
WHERE:
* ```clnet::Dict``` is the Cloudnet data from categorize and classification output files,
* ```lidar::Dict``` (Optional) the lidar data with keys :time and :CBH [m]
* ```nlayers::Int``` (Optional) number of cloud layer to try to detect, default 3
* ```alttime::Vector{DateTime}``` (Optional) to which the ouput will be fit, defualt none
* ```smooth_classify::Bool``` (Optional) to smooth the 2D CLASSIFICATION to avoid spikes.
* ```liquid_base::Bool``` (Optional) if true then add extra output with liquid cloud base. 

OUTPUT:
* ```CBH::Matrix(ntime, nlayers)``` with cloud base height in m
* ```CTH::Matrix(ntime, nlayers)``` with cloud top height in m
* ```CLB::Matrix(ntime, nlayers)``` Optional if input flag 'liquid_base' is true.

NOTE: be sure clnet[:height] and lidar[:CBH] have the same units, e.g. m

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function estimate_cloud_layers(clnet::Dict; lidar=nothing, nlayers=3, alttime=nothing, smooth_classify=false, liquid_base=false)

    # Defining constants:
    ntime = length(clnet[:time])
    cloud_flags = (1,3,4,5,7)
    hydro_flags = (1:7)
    liquid_flags = (1,3,5,7)

    # Smoothing Cloudnet classification array to minimize noise:
    fmin(x) = filter(>(0), x) |> z-> isempty(z) ? 0 : maximum(z)
    
    CLASSIFY = let tmp=clnet[:CLASSIFY]
        smooth_classify ? round.(mapwindow(minimum, clnet[:CLASSIFY], (1,5))) : tmp
    end
    #CATEBITS = round.(imfilter(clnet[:CATEBITS], ker2d));
    
    # creating output arrays:
    CBH = fill(NaN32, ntime, nlayers)
    CTH = fill(NaN32, ntime, nlayers)
    CLB = fill(NaN32, ntime, nlayers)
    
    # cheking if optional lidar data is provided:
    if !isnothing(lidar) && isa(lidar, Dict)
        CBH_lidar = Interpolate2Cloudnet(clnet, lidar[:time], lidar[:CBH])
    end
    
    # starting iteration over time dimension:
    foreach(1:ntime) do k
        CT = 1
	
        # assigning true/false pixels corresponding to cloud_flags:
        tmp = map(j->any(j .∈ cloud_flags), CLASSIFY[:, k])
        liq = map(j->any(j .∈ liquid_flags), CLASSIFY[:, k])
        
        for ih ∈ (1:nlayers)
            CB = if isnothing(lidar)
                findfirst(tmp[CT:end])
            elseif isnan(CBH_lidar[k])
                nothing
            else
                abs.(clnet[:height] .- CBH_lidar[k]) |> argmin
            end
            isnothing(CB) && continue
            CB += CT-1
            
            # cloud liquid base:
            CL = liq[CB:end] |> findfirst

            # cloud top:
            CT = findfirst(j->all(j .∉ hydro_flags), CLASSIFY[CB:end, k])
            isnothing(CT) && (@warn "Cloud base found without cloud top at layer $(ih)! $k $(CB)"; break)
            
            CT += CB - 1
            
            CBH[k, ih] = clnet[:height][CB]
            CTH[k, ih] = clnet[:height][CT]
            CLB[k, ih] = !isnothing(CL) ? clnet[:height][CL+CB-1] : NaN32
        end
        
    end

    if !isnothing(alttime) && typeof(alttime)<:Vector{DateTime}
	ntime = length(alttime)
	CBH_out = fill(NaN32, ntime, nlayers)
	CTH_out = fill(NaN32, ntime, nlayers)
        CLB_out = fill(NaN32, ntime, nlayers)
        
	for (idxT, T) ∈ enumerate(alttime)
	    δt, idxin = findmin(abs.(T .- clnet[:time]))
	    δt > Minute(1) && continue
	    CBH_out[idxT,:] = CBH[idxin, :]
	    CTH_out[idxT,:] = CTH[idxin, :]
            CLB_out[idxT,:] = CLB[idxin, :]
	end
	#idxin = [findmin(abs.(T .- clnet[:time])) |> x->x[] for T ∈ alttime]
    else
	CBH_out = CBH
	CTH_out = CTH
        CLB_out = CLB
    end

    if liquid_base
        return CBH_out, CTH_out, CLB_out
    end
    
    return CBH_out, CTH_out
end
# ----/

# ****************************************************
"""
Function to interpolate variables to the cloudnet time (and height)

USAGE:
```julia-repl
julia> var_out = Interpolate2Cloudnet(clnet, time_in, var_in)
```
WHERE:
* ```clnet::Dict``` with cloudnet output,
* ```time_in::Vector{DateTime}``` with the time of variable to interpolate,
* ```var_in::Vector{Any}``` with the variable to interpolate.

OUTPUT:
* ```var_out::Vector{Any}``` with the interpolated variable at cloudnet time.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function Interpolate2Cloudnet(clnet::Dict, time_in::Vector, var_in::Vector)
    cnt_hr = datetime24hours(clnet[:time])
    var_hr = datetime24hours(time_in)
    itp = LinearInterpolation(var_hr, var_in, extrapolation_bc=Linear())

    return itp(cnt_hr)
end
function Interpolate2Cloudnet(clnet::Dict, data_in::Dict, var::Symbol)
    cnt_hr = datetime24hours(clnet[:time])
    var_hr = datetime24hours(data_in[:time])
    
    nodes = (data_in[:height], var_hr)

    itp = interpolate(nodes, data_in[var], Gridded(Linear()))
    extp = extrapolate(itp, Flat())

    outvar = [extp(i,j) for i ∈ clnet[:height], j ∈ cnt_hr]    

    return outvar
end
# ----/

# **************************************************
"""
Function to estimate the temperature of given altitudes,
For instance the temperature at cloud base & top heights.

```julia-repl
julia> T_h = cloud_temperature(data, H)
julia> T_h = cloud_temperature(data, H, clnet=clnet, var=:Temp)
```
WHERE:
* ```data::Dict()``` data with variable keys :time and :T for temperature, (can be clnet),
* ```H::Array{T}(time, layers)``` with the altitudes at which the temperature is estimated,
* ```clnet::Dict``` (Optional) when data is not same resolution, interpolate to clnet[:time],
* ```var::Symbol``` (Optional) variable from data to be used, default :T

Input Dict data must have variables with key :T for temperature and :height.
The output will have the same units as ```data[:T]``` and same dimentions as H.

NOTE: be sure H and data[:height] have the same units, e.g. m

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function cloud_temperature(data::Dict, H_in::Array; clnet=nothing, var=:T)
    
    !haskey(data, :height) && @error "data has not key :height"
    !haskey(data, var) && @error "data has not key $(var)"
    
    nt = size(H_in,1)
    nlyr = size(H_in, 2)
        
    T_h = fill(NaN32, size(H_in)...)

    T_var = if !isnothing(clnet) && isa(clnet, Dict)
        Interpolate2Cloudnet(clnet, data, var)
    elseif length(data[:time]) == nt
        data[var]
    else
        @error "data :time length not compatible with length of input H!"    
    end

    for h ∈ 1:nlyr
	for i ∈ 1:nt
	    isnan(H_in[i, h]) && continue
	    idx_h = abs.(data[:height] .- H_in[i, h]) |> argmin
	    T_h[i, h] = T_var[idx_h, i]
	end
    end

    return T_h
end
# ----/

# end of file
