#!/opt/julia-1.6.0/bin/julia

# Main script to run the conversion of ARM data files to netCDF input CloudNetpy files.

#using NCDatasets, DataStructures
#using Dates
using Printf
#using UUIDs
#using Statistics
#using ImageFiltering

using CloudnetTools

# Declaration of parameter to be used to run
# Possible input parameters are:
# * lidar => "CEIL10m" or "HSRL"
# * radar => "KAZR/ARSCL" or "KAZR/xxx"
# * radiometer => "MWR/RET" or "MWR/LOS"
# * model => "ECMWF" or "ICON"
input_params = Dict(
    "site" => "utqiagvik", 
    "campaign" => "nsa",
    "products" => Dict("radar" => "KAZR/ARSCL",
                "ceilometer" => "CEIL10m", #
                "lidar" => "HSRL", #"", #
                "model" => "",
                "radiometer" => "MWR/RET"), #),
    "years" => (2018),
    "months" => (1),
    "days" => (31),
    "data_path" => "/home/psgarfias/LIM/data/",
    "output_path" => "/home/psgarfias/LIM/data/CloudNet/input");


# *****************************************************************
#include("./arm_mwr2nc.jl")
#include("./arm_lidar2nc.jl")
#include("./arm_hsrl2nc.jl")
#include("./arm_kazr2nc.jl")
# *****************************************************************
"""
Function to invoque the respective ARM data converter function:
    result = ARMconverter("radar", "/data/KAZR/armradar.20180121.nc")
    will convert the file _armradar.20180121.nc_ into CloudNet input file.
    
    (c) Pablo Saavedra G.
"""
function ARMconverter(the_key, arm_filenc, out_path)
    # Here the values of dictionary keys needs to be replaces by ARM data converter function:
    list_of_func = Dict(
        "radar" => CloudnetTools.kazr2nc,
        "lidar" => CloudnetTools.hsrl2nc,
        "radiometer" => CloudnetTools.mwr2nc,
        "model" => isempty,
        "ceilometer" => CloudnetTools.lidar2nc
    );

    !isdir(out_path) && mkpath(out_path);
 
    ex = :($list_of_func[$the_key]($arm_filenc, $out_path))

    return eval(ex)
end

# Running over all specified products:
for year ∈ input_params["years"]
    for month ∈ input_params["months"]
        for day ∈ input_params["days"]
            date_pattern = @sprintf("%04d%02d%02d", year, month, day);
            for (thekey, product) ∈ filter(x->!isempty(x[2]), input_params["products"] |> collect)
            
                DATA_PATH = @sprintf("%s-%s/%s/%04d", input_params["site"], input_params["campaign"], product, year);
                DATA_PATH = joinpath(input_params["data_path"], DATA_PATH)
                
                @assert isdir(DATA_PATH) "$DATA_PATH does not exist!"
                # creating the file names for every product:
                # RADAR:
                # LIDAR:
                # MWR:
                
                file_pattern = lowercase(input_params["campaign"]);
                files_in_dir = readdir(DATA_PATH, join=true)
                arm_filenc = filter(x->all(occursin.([file_pattern, date_pattern], x)), files_in_dir);
                
                isempty(arm_filenc) && continue
                #@assert !isempty(arm_filenc) "Data not found for $product at $DATA_PATH !!!"
                
                arm_filenc = arm_filenc[1]
                out_path = joinpath(input_params["output_path"], product, "$(year)")
                #isfile(arm_filenc) ? println("Workin on $arm_filenc") : error("MWR file does not exist!")
                a = ARMconverter(thekey, arm_filenc, out_path)
                println(thekey," ",a)
                # MODEL:

            end
        end
    end
end


# --- end of script.

