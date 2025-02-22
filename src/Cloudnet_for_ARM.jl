# Part of CloudnetTools.jl
# Module including tools and functions to work with ARM data
# ******************************************************************
# +++++++++++++++ MODULE ARM for working with arm data files +++++++
# ******************************************************************

module ARM

using ARMtools
using NCDatasets, DataStructures
using Dates
using Interpolations
using Printf
using Statistics
using ImageFiltering
using UUIDs

# Including transformation functions: from ARM files to Cloudnet input files.
Base.include(ARM, "aux_functions.jl")
Base.include(ARM, "arm_mwr2nc.jl")
Base.include(ARM, "arm_lidar2nc.jl")
Base.include(ARM, "arm_kazr2nc.jl")
Base.include(ARM, "arm_hsrl2nc.jl")
Base.include(ARM, "arm_rsonde2nc.jl")


# *****************************************************************
"""
Function to convert model input data.

USAGE:
```julia-repl
julia> result = ARM.model2nc("/data/ECMWF/20180121_arm-nsa_ecmwf.nc", "/data/output")
```
NOTE:
This function is still dummy. When ECMWF input file is provides, it only copies
that file to the ouput_path where other Cloudnet input files might be located.
In the future, other alternative MWP models will be processed and converted
to Cloudnet input.
    
Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function model2nc(infile::String, output_path::String)

    !contains(infile, "ecmwf") && (return missing)
    CMD = `cp -up $infile $output_path`
    run(CMD)

    return true
end
# ----/

# *****************************************************************
"""
Function to obtain the instrument type from netcdf attributes:

USAGE:
```julia-repl
julia> get_instrument_type("/data/KAZR/ARSCL/nsaarsclkazr1kolliasC1.c0.20171109.000000.nc")
julia> "arsclkazr1kollias"
```
Returning the value of parameters "platform_id" or "datastream".
In case the NetCDF file does not include any of those paramters, `nothing` is
returned. An alternative output can be used by given an optional paramter
`default_type="kazr"`.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function get_instrument_type(ns::NCDataset; default_type=nothing)
    list_attrib = ("platform_id", "datastream")

    tmp = [ns.attrib[var] for var ∈ list_attrib if haskey(ns.attrib, var)]
        
    type = if !isempty(tmp)
        tmp[1]
    else
        @warn "Instrument type not found! return default $(default_type)"
        default_type
    end
    return type
end
function get_instrument_type(fn::String; default_type=nothing)
    
    type = NCDataset(fn) do ns
        get_instrument_type(ns; default_type=default_type)
    end
    return type
end


# *****************************************************************
"""
Function to invoque the respective ARM data converter given a ARM NetCDF:

USAGE:
```julia-repl
julia> ARM.converter(keys, input_file, output_dir)
julia> ARM.converter(list_of_files, output_dir)
julia> ARM.converter(yyyy, mm, dd, keys, input_params)
```
WHERE:
* ```keys::Symbol``` - Indicator of instrument type, it can be one of (:radar, :lidar, :mwr, :ceilometer, :radiosonde),
* ```input_file::String``` - full path to ARM input file to convert,
* ```output_dir::String``` - path to the ouptut directory to place the output files,
* ```list_of_files::Dict{Symbol, String}``` - with keys to input files e.g. ```:radar=>"/data/KAZR/armradar.20180121.nc"```,
* ```yyyy::Int``` - year
* ```mm::Int``` - month
* ```dd::Int``` - day
* ```extra_params::Dict{Symbol, Any}``` - Optional Dictionary with extra arguments, default empty.
* ```input_params::Dict{Symbol, Any}``` - must include keys :data_path, :output_path, :products
* * ```:data_path => ::String``` - indicating the path for the ARM input data,
* * ```:output_path => ::String``` - indicating the path where to put the converted data,
* * ```:products => ::Dict{Symbol, Any}``` - with ```keys``` indicating the ARM product to use.


EXAMPLE:
To convert the file "/data/KAZR/armradar.20180121.nc" into CloudNet input file, with extra parameters given
by the dictionary `extras`:
```julia-repl
julia> extras = Dict(:snr_filter=>nothing, :site=>"arctic", :campaign=>"mosaic", :lat=>82, :lon=>105);
julia> result = ARM.converter(:radar, "/data/KAZR/armradar.20180121.nc", "/output"; extra_params=extras)
julia> "20180121_arctic-mosaic_radar.nc"
```
1. Alternative:
```julia-repl
julia> list_of_files=Dict(:radar=>"/data/KAZR/armkazr.20180121.nc",
                          :lidar=>"/data/CEIL10m/armceil10m.20180121.nc",
                          :mwr=>"/data/MWR/armret1.20180121.nc",
                          :radiosonde=>"/data/INTERPOLATESONDE/arminterpolrs.20180121.nc");
julia> ARM.converter(list_of_files, "/output"; extra_params=extras)
```
2. Alternative:
```julia-repl
julia> NSAproduct = Dict(
                      :radar => ("KAZR/ARSCL", "KAZR/CORGE"),
                      :lidar => "CEIL10m",
                      :mwr => ("MWR/RET", "MWR/LOS"),
                 );

julia> input_params = Dict(
    :site => "utqiagvik", 
    :campaign => "nsa",
    :products => NSAproduct,
    :data_path => "/data/",
    :output_path => "/output")

julia> ARM.converter(2019,01,19, :radar, input_params; owndir=false)
julia> "20190119_nsa-utqiagvik_radar.nc"
```
The function looks into a data path as ```/data/KAZR/ARSCL/2019/``` to find a data file for ```20190119```,
in case the product ```KAZR/ARSCL``` does not exist, alternatively look at ```KAZR/CORGE``` product.
If the optional argument ```owndir=false``` , then all converted data files
will be placed at ```/output```. Conversely ```owndir=true``` (default) then the converted files are
placed at ```/output/KAZR/ARSCL/2019``` or ```/output/KAZR/CORGE/2019```, etc.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
(c) Pablo Saavedra Garfias
"""
function converter(the_key::Symbol, arm_filenc::String, out_path::String; extra_params=Dict{Symbol, Any}())

    func_for_lidar = let tmp=get_instrument_type(arm_filenc)
        if tmp=="hsrl"
            ARM.hsrl2nc
        elseif tmp=="ceil10m"
            ARM.lidar2nc
        elseif the_key ∈ (:lidar, :ceilometer)
            @info "Type of $(tmp) is neither CEIL10m nor HSRL! using default as CEIL10m. Output file can be compromised!"
            ARM.lidar2nc
        end
    end
    
    # Here the values of dictionary keys needs to be replaces by ARM data converter function:
    list_of_func = Dict(
        :radar => ARM.kazr2nc,
        :lidar => func_for_lidar,
        :mwr => ARM.mwr2nc,
        :model => ARM.model2nc,
        :ceilometer => ARM.lidar2nc,
        :radiosonde => ARM.rsonde2nc
    );

    !isdir(out_path) && mkpath(out_path);

    ex = :($(list_of_func[the_key])($arm_filenc, $out_path, extra_params=$extra_params))
    try
        return eval(ex)
        
    catch e
        @warn "$(the_key) cannot be converted! $(e)"
        return nothing
    end
end

# -- OR --:
function converter(list_of_data::Dict, out_path::String; extra_params=Dict{Symbol, Any}())

    list_of_func = Dict(
        :radar => ARM.kazr2nc,
        :lidar => ARM.hsrl2nc,
        :radiometer => ARM.mwr2nc,
        :model => ARM.model2nc,
        :ceilometer => ARM.lidar2nc
    );

    output_files = String[]
    foreach(list_of_data) do (the_key, arm_data)
        output_ex = converter(the_key, arm_data, out_path, extra_params=extra_params)
        #ex = :($(list_of_func[the_key])($arm_data, $out_path, extra_params=$extra_params))
        isnothing(output_ex) && @warn "$(the_key) cannot be converted!. Return Nothing instead."
        !isnothing(output_ex) && push!(output_files, output_ex )
    end
    
    return output_files
end

# -- OR --:
function converter(yy, mm, dd, thekey::Symbol, input_params::Dict; owndir=true)

    product = input_params[:products][thekey]

    DATA_PATH = input_params[:data_path]

    !isdir(DATA_PATH) && (@warn("bad $DATA_PATH"); return missing)

    # creating the file name for every listed product
    tmp = ARMtools.getFilePattern(DATA_PATH, product, yy, mm, dd)

    isnothing(tmp) && (@warn("NO $product found!"); return missing)

    arm_filenc = typeof(tmp) <: Vector ? tmp[1] : tmp

    # create output path according to optional input "owndir"
    OUT_PATH = let tmp_path=input_params[:output_path]
        owndir ? joinpath(tmp_path, product, "$(yy)") : tmp_path
    end

    # if output path does not exist, then create it
    !isdir(OUT_PATH) && mkpath(OUT_PATH)
    
    return converter(thekey, arm_filenc, OUT_PATH; extra_params=input_params)
end
# ----/

end  # end of module ARM
# ----/

# end of file
