# Part of CloudnetTools.jl

#=
 Module including function to generate categorizing, classifing
 and the all available products delivered by cloudnetpy
 ****************************************************************
 +++ MODULE FMI for interfacing cloudnetpy to ARM data files ++++
 ****************************************************************
=#
module ACTRIS

using PyCall
using DataFrames

const catcnet = PyNULL()
const procnet = PyNULL()
const cloudnetpy_qc = PyNULL()

function __init__()
    copy!(catcnet, pyimport("cloudnetpy.categorize"))
    copy!(procnet, pyimport("cloudnetpy.products"))
    copy!(cloudnetpy_qc, pyimport("cloudnetpy_qc"))
end

# ***************************************************************
"""
Function to generate the categorization file from cloudnetpy
> using CloudnetTools.ACTRIS
> uuid = ACTRIS.categorize_it(input_files, "/tmp/categorize.nc")

INPUT:
* input_files::Dict{
    :radar => "arm_radar.nc",
    :lidar => "arm_lidar.nc",
    :mwr   => "arm_mwr.nc",
    :model => "ecmwf_model.nc"}
* output_file::String with the file name of the output file.

RETURN:
If successful returns the UUID of the categorization file,
if failed returns a warning message and nothing

"""
function categorize_it(data_files::Dict, categorize_file::String)

    try
	uuid = catcnet.generate_categorize(data_files, categorize_file)
    catch e
	@warn "No categorization possible: $(categorize_file)\n $e"
	uuid = nothing
    end
    
    return uuid
end
# ----/


# *****************************************************************
"""
Function to generate cloudnetpy products i.e. iwc, lwc, der, ier, etc.
USAGE:
> uuid = ACTRIS.generate_products(PROD, categorize_file, output_file)

INPUTs:
* PROD::Symbol with the key for the product to produce, e.g. :lwc,
* categorize_file::String full path to the categorize file created by cloudnetpy,
* output_file:Stirng full path and file name for the generated product

OR alternatively the function can be used as:

> uuid = ACTRIS.generate_products(products_files, cloudnet_products)

INPUTS: dictionaries containing informations as follow,
* product_files::Dict(
    :categorize => "/data/arm/input/arm_categorize.nc",
    :iwc => "/data/arm/cloud_product/arm_iwc.nc",
    :lwc => "/data/arm/cloud_product/arm_lwc.nc",
    :der => "/data/arm/effradius/arm_der.nc")
* cloudnet_products::Dict(
    :classification => false,
    :iwc => false,
    :lwc => true,
    :der => true,
    :drizzle => false)

RETURNS:
If successful the UUID of the lates product generated, otherwise type nothing

"""
function generate_products(K::Symbol, fn_categ::String, output_file::String)

    # checking if categorize file exist:
    if !isfile(fn_categ)
        @warn "$(fn_categ) cannot be found!"
        return nothing
    end

    # checking if Symbol K correspont to valid cloudnetpy product:
    if !any(==(K), (:lwc, :iwc, :drizzle, :der, :ier, :classification))
        @warn "$K is not a valid/supported cloudnetpy product"
        return nothing
    end
    
    # creating expression for the corresponding cloudnetpy product:
    ex = Meta.parse("uuid=procnet.generate_$K(\"$(fn_categ)\", \"$(output_file)\")")
    try
        eval(ex)
    catch e
        println("\e[1m\e[38;2;230;30;30;249m","* Error trying to run $K from $(fn_categ)");
	println(e)
	uuid = nothing
    end
    
    return uuid
end
function generate_products(clnet_file::Dict, CLNTprod::Dict)

    if !haskey(clnet_file, :categorize)
	K = findfirst(==(1), CLNTprod)

	clnet_file[:categorize] = try
	    replace(clnet_file[K], String(K)=>"categorize")
	catch e
	    @warn "cannot generage categorize file name: $e"
	    "nothing"
	end
    end

    if !isfile(clnet_file[:categorize])
	@warn "categorize file not found: $(clnet_file[:categorize]) "
    end
		    
    uuids = Any[]
    for (K, V) in CLNTprod
	K == :categorize && continue
	!V && continue

        uuid = generate_products(K, clnet_file[:categorize], clnet_file[K])
        push!(uuids, uuid)
    end

    return uuids
end
# ----/


# ================================================================

#=
 Module including to check Quality Control rounitnes of Cloudnetpy
 ******************************************************************
 ++++ MODULE QC for checking generated arm input data files +++++++
 ******************************************************************
=#

function cloudnet_qc(data_file::Dict)

    meta_df = DataFrame()
    data_df = DataFrame()
    foreach(data_file) do (k, fn)
        println("***** ", k, ":")
        tmp = cloudnet_qc(fn)
        #meta_df = vcat(meta_df, hcat(DataFrame(:type=>k), tmp[:metadata_qc]))
        #data_df = vcat(data_df, hcat(DataFrame(:type=>k), tmp[:data_qc]))
    end
    return #meta_df, data_df
end
function cloudnet_qc(fn::String)

    function Dict2DF(dict_in::Dict)
        return DataFrame(; [Symbol(k)=>v for (k,v) ∈ dict_in if !isempty(v)]...)
    end
    
    quality = cloudnetpy_qc.Quality(fn)
    
    metadata_results = quality.check_metadata()
    data_results = quality.check_data()
    println("METADATA_QC:")
    println(Dict2DF(metadata_results))
    println("DATA_QC:")
    println(Dict2DF(data_results))
    return Dict(:metadata_qc => Dict2DF(metadata_results),
                :data_qc => Dict2DF(data_results))
end
# ----/

end # --- end of module

# end of file
