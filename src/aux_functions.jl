# *************************************************

# ****
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
# /----

# ****
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
# /----

# ****************************************
"""
"""
function datetime24hours(time_in::Vector{DateTime})
    return hour.(time_in) + minute.(time_in)/60. + second.(time_in)/3600.0;
end

# end of file
