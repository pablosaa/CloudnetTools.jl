# Part of CloudnetTools.jl
# Module including to check Quality Control rounitnes of Cloudnetpy
# ******************************************************************
# ++++ MODULE QC for checking generated arm input data files +++++++
# ******************************************************************

module QC

using DataFrames
using PyCall

cloudnet_qc = pyimport("cloudnetpy_qc")

function get_cloudnet_qc(data_file::Dict)

    meta_df = DataFrame()
    data_df = DataFrame()
    foreach(data_file) do (k, fn)
        println("***** ", k, ":")
        tmp = get_cloudnet_qc(fn)
        #meta_df = vcat(meta_df, hcat(DataFrame(:type=>k), tmp[:metadata_qc]))
        #data_df = vcat(data_df, hcat(DataFrame(:type=>k), tmp[:data_qc]))
    end
    return #meta_df, data_df
end
function get_cloudnet_qc(fn::String)
    
    quality = cloudnet_qc.Quality(fn)
    metadata_results = quality.check_metadata()
    data_results = quality.check_data()
    println("METADATA_QC:")
    println(Dict2DF(metadata_results))
    println("DATA_QC:")
    println(Dict2DF(data_results))
    return Dict(:metadata_qc => Dict2DF(metadata_results),
                :data_qc => Dict2DF(data_results))
end

function Dict2DF(dict_in::Dict)
    return DataFrame(; [Symbol(k)=>v for (k,v) ∈ dict_in if !isempty(v)]...)
end
              
end
# /----
# end of module
