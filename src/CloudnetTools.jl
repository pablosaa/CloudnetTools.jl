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
        ss = @. floor(mod(mi*60, 60))
        mi = floor.(mi)
        global tit = @. DateTime(yy, mm, dd, hh, mi, ss)

        global Ze = nc["Z"][:,:]
    end
    Ze[Ze .> 100] .= NaN;
    return Dict(:time=>tit, :Z=>Ze)
end

end
# --end of script
