# Main CloudNet.jl module.
"""
A set of tools to process and analyze data outputs from CloudNet classification algorithm.

(c) 2020, Pablo Saavedra Garfias
University of Leipzig
Faculty of Physics and Geosciences
LIM

See LICENSE
"""


module CloudNet

using NCDatasets
using Dates
using Interpolations

function readCLNFile(nfile::String; modelreso=false)
    @assert isfile(nfile) error("$nfile cannot be found!")
    if contains(nfile, "categorize")
        println("reading Categorize file")
    elseif contains(nfile, "classific")
        println("readin Classification file")
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
        #MWR
        :LWP => "lwp",
        #CLOUDNETpy
        :P_INSECT => "insect_prob",
        :QUALBITS => "quality_bits",
        :CATEBITS => "category_bits",
        #MODEL
        :T => "temperature",
        :Tw => "Tw",
        :Pa => "pressure",
        :QV => "q",
        :UWIND => "uwind",
        :VWIND => "vwind",
    )

    # Classification variables to read:
    vars_classific = Dict(
        :CLASSIFY => "target_classification",
        :DETECTST =>"detection_status",
    )

    # Defining Output varaible:
    var_output = Dict{Symbol, Any}()

    # Starting reading CloudNet files:
    NCDataset(nfile; format=:netcdf4_classic) do nc
        yy = Int64(nc.attrib["year"])
        mm = Int64(nc.attrib["month"])
        dd = Int64(nc.attrib["day"])
        tmp_time = float.(nc["time"])
        hh = floor.(Int64, tmp_time)
        mi = @. mod(tmp_time*60, 60)
        ss = @. mod(mi*60, 60)
        ms = @. mod(ss, 1)*1000
        
        mi = floor.(Int64, mi)
        ss = floor.(Int64, ss)
        ms = floor.(Int64, ms)
        
        tit = @. DateTime(yy, mm, dd, hh, mi, ss, ms)

        var_output[:time] = tit

        var_output[:height] = nc["height"][:]

      
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

        # If modelreso = true, interpolate model data to cloudnet resolution:

        model_time = float.(nc["model_time"])
        model_height = nc["model_height"][:]
        if !modelreso
            var_output[:model_time] = model_time
            var_output[:model_height] = model_height
        else
            var_output = ConvertModelResolution(var_output,
                                                model_time,
                                                model_height;
                                                cln_time=tmp_time)
        end

    end

    
    
    # For Classification dataset:
    classfile = replace(nfile, "categorize" => "classific")
    @assert isfile(classfile) errro("$classfile cannot be found!")
    NCDataset(classfile; format=:netcdf4_classic) do nc
        
        [var_output[x] = nc[vars_classific[x]][:,:] for x ∈ keys(vars_classific)];
        
    end

    
    return var_output;

end  # end of function
# ----/

# ***************************************************
# Function to normalize matrix to [0,1]
function TransformZeroOne(X::T) where T<:AbstractArray
    nonans = .!isnan.(X)
    x₀, x₁ = extrema(X[nonans])
    return (X .- x₀)./(x₁ - x₀);
end  # end of function
# ----/

# *********************************************************
# Interpolate Meteo data from Model to CloudNet resolution
function ConvertModelResolution(cln_in::Dict{Symbol, Any}, model_time, model_height;
                                cln_time=nothing, cln_height=nothing)

    nodes = (model_height, model_time)
    if isnothing(cln_time)
        cln_time = @. hour(cln_in[:time]) + minute(cln_in[:time])/60 + seconds(cln_in[:time])/3600
    end
    
    if isnothing(cln_height)
        cln_height = cln_in[:height]
    end

    MODELVAR = (:T, :Pa, :UWIND, :VWIND, :QV)

    for var ∈ MODELVAR
        itp = interpolate(nodes, cln_in[var], Gridded(Linear()))
        outvar = [itp(i,j) for i ∈ cln_height, j ∈ cln_time]
        cln_in[var] = outvar
    end

    return cln_in
end
# ----/


end  # end of Module
# --end of script
