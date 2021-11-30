# Part of CloudnetTools.jl
# Module including tools and functions to work with ARM data

module Cloudnet4ARM

using ARMtools
using NCDatasets, DataStructures
using Dates
using Interpolations
using Printf
using Statistics
using ImageFiltering
using UUIDs

# Including transformation functions: from ARM files to Cloudnet input files.
Base.include(Cloudnet4ARM, "arm_mwr2nc.jl")
Base.include(Cloudnet4ARM, "arm_lidar2nc.jl")
Base.include(Cloudnet4ARM, "arm_kazr2nc.jl")
Base.include(Cloudnet4ARM, "arm_hsrl2nc.jl")

function model2nc(infile::String, output_path::String)

    !contains(infile, "ecmwf") && (return missing)
    CMD = `cp -up $infile $output_path`
    run(CMD)

    return true
end
# ----/


# *****************************************************************
"""
Function to invoque the respective ARM data converter function:
To convert the file _armradar.20180121.nc_ into CloudNet input file:

USAGE:

>result = ARMconverter(:radar, "/data/KAZR/armradar.20180121.nc", "/data/output")

    
    (c) Pablo Saavedra G.
"""
function ARMconverter(the_key::Symbol, arm_filenc::String, out_path::String)
    # Here the values of dictionary keys needs to be replaces by ARM data converter function:
    list_of_func = Dict(
        :radar => Cloudnet4ARM.kazr2nc,
        :lidar => Cloudnet4ARM.hsrl2nc,
        :radiometer => Cloudnet4ARM.mwr2nc,
        :model => Cloudnet4ARM.model2nc,
        :ceilometer => Cloudnet4ARM.lidar2nc
    );

    !isdir(out_path) && mkpath(out_path);
 
    ex = :($(list_of_func[the_key])($arm_filenc, $out_path))

    return eval(ex)
end
# -- OR --:
function ARMconverter(yy, mm, dd, thekey::Symbol, input_params::Dict; owndir=true)

    product = input_params[:products][thekey]

    DATA_PATH = input_params[:data_path]

    !isdir(DATA_PATH) && (print("bad $DATA_PATH"); return missing)

    # creating the file name for every listed product
    tmp = ARMtools.getFilePattern(DATA_PATH, product, yy, mm, dd)

    isnothing(tmp) && (print("NO $product found!"); return missing)

    arm_filenc = typeof(tmp) <: Vector ? tmp[1] : tmp

    # create output path according to optional input "owndir"
    OUT_PATH = let tmp_path=input_params[:output_path]
        owndir ? joinpath(tmp_path, product, "$(yy)") : tmp_path
    end

    # if output path does not exist, then create it
    !isdir(OUT_PATH) && mkpath(OUT_PATH)
    
    return ARMconverter(thekey, arm_filenc, OUT_PATH)
end
# ----/

end  # end of module
# ----/

# end of file
