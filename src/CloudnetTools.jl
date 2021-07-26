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

function readCLNFile(nfile::String)
    @assert isfile(nfile) error("$nfile cannot be found!")
    if contains(nfile, "categorize")
        println("reading Categorize file")
    elseif contains(nfile, "classific")
        println("readin Classification file")
    else
        error("$nfile does not apear to be a Cloudnet file!")
    end
    NCDataset(nfile; format=:netcdf4_classic) do nc
        yy = Int64(nc.attrib["year"])
        mm = Int64(nc.attrib["month"])
        dd = Int64(nc.attrib["day"])
        tmp = nc["time"]
        hh = floor.(Int64, tmp)
        mi = @. mod(tmp*60, 60)
        ss = @. mod(mi*60, 60)
        ms = @. mod(ss, 1)*1000
        
        mi = floor.(Int64, mi)
        ss = floor.(Int64, ss)
        ms = floor.(Int64, ms)
        
        global tit = @. DateTime(yy, mm, dd, hh, mi, ss, ms)

        global height = nc["height"][:]
        # RADAR
        global Ze = nc["Z"][:,:]
        global V = nc["v"][:,:]
        global σV = nc["v_sigma"][:,:]
        global ωV = nc["width"][:,:]
        # LIDAR
        global β = nc["beta"][:,:]
        # MWR
        global LWP = nc["lwp"][:,:]
        # CLOUDNETpy
        global P_insect = nc["insect_prob"][:,:]
        global Q_bits = nc["quality_bits"][:,:]
        global CAT_bits = nc["category_bits"][:,:]
        # MODEL
        global T = nc["temperature"][:,:]
        global Tw= nc["Tw"][:,:]
        global Pa = nc["pressure"][:,:]
        global QV = nc["q"][:,:]
        global UWIND = nc["uwind"][:,:]
        global VWIND = nc["vwind"][:,:]
        
    end

    # For Classification dataset:
    classfile = replace(nfile, "categorize" => "classific")
    @assert isfile(classfile) errro("$classfile cannot be found!")
    NCDataset(classfile; format=:netcdf4_classic) do nc
        global CLASS = nc["target_classification"][:,:]
        global detection = nc["detection_status"][:,:]
    end

    # Cleaning variables:
    Ze[Ze .> 100] .= NaN;
    
    return Dict(:time=>tit,
                :height=>height,
                :Z => Ze,
                :V => V,
                :σV => σV,
                :ωV => ωV,
                :β => β,
                :LWP => LWP,
                :T => T,
                :Tw => Tw,
                :Pa => Pa,
                :Qv => QV,
                :UWIND => UWIND,
                :VWIND => VWIND,
                :QUBITS => Q_bits,
                :P_INSECT => P_insect,
                :CATEBITS => CAT_bits,
                :CLASSIFY => CLASS,
                :DETECTST => detection,
                )
end

end
# --end of script
