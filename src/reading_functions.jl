#=
********************************************************************
File contains reading functions for the CloudnetTools.jl package
********************************************************************
=#

# ***************************************************************
# Read Cloudnet Ice water content IWC file
"""
# Read Cloudnet Ice water content IWC file

> iwc = readIWCFile(nfile::String; modelreso=false)

"""
function readIWCFile(nfile::String; modelreso=false)

    @assert isfile(nfile) error("$nfile cannot be found!")

    vars_categorize = Dict(
        :height => "height",
        :iwc => "iwc_inc_rain",
        :flag => "iwc_retrieval_status",
    )

    # Defining Output varaible:
    var_output = Dict{Symbol, Any}()

    NCDataset(nfile; format=:netcdf4_classic) do nc
        var_output[:time] = convert_time2DateTime(nc)

        for inkey ∈ keys(vars_categorize)
            x = vars_categorize[inkey]

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

    end

    # Adding computed variables:
    # integration of iwc only for pixels with flag=1 or 2:
    # Cloudnetpy units iwc [kg m⁻³] and height [m]
    ##var_output[:IWP] = let tmp = 1f3var_output[:iwc]
    ##
    ##    @. tmp[!(0 < var_output[:flag] < 3)] = NaN
    ##    # IWP [g m⁻²]
    ##    ∫fdh(tmp, var_output[:height])
    ##end
    
    return var_output
end
# ----/

# ***************************************************************
# Read Liquid water content LWC file
"""
# Read Liquid water content LWC file

> lwc = readLWCFile(nfile::String; modelreso=false)

"""
function readLWCFile(nfile::String; modelreso=false)

    @assert isfile(nfile) error("$nfile cannot be found!")

    vars_categorize = Dict(
        :height => "height",
        :lwc => "lwc",
        :flag => "lwc_retrieval_status",
        :LWP => "lwp",
    )

    # Defining Output varaible:
    var_output = Dict{Symbol, Any}()

    NCDataset(nfile; format=:netcdf4_classic) do nc
        var_output[:time] = convert_time2DateTime(nc)

        for inkey ∈ keys(vars_categorize)
            x = vars_categorize[inkey]

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

    end

    # Additional computation:
    # integration of iwc only for pixels with flag=1 or 2:
    @. var_output[:lwc][!(0 < var_output[:flag] < 3)] = NaN
    
    return var_output
end
# ----/

# ******************************************************************
# Reading Classification & Categorize files both at once:
"""
# Fuction to read Cloudnet classification or category files:

> classi = readCLNFile(nfile::String; modelreso=false)

if optional parameter modelreso=true, then the model variables are
interpolated from hourly resolution to Cloudnet resolution.
"""
function readCLNFile(nfile::Vector{String}; modelreso=false, altfile=nothing)
    cln_out = Dict{Symbol, Any}()
    catvar = Dict{Symbol, Union{Nothing, Int}}()
    ntime = -1

    getdim(x,n) = findall(==(n), size(x))

    foreach(nfile) do fn
        # reading single file
        cln = readCLNFile(fn, modelreso=modelreso, altfile=altfile)

        if isempty(cln_out)
            cln_out = cln
            catvar = let ntime=length(cln[:time])
                tmp = [k=>getdim(v, ntime) for (k,v) ∈ cln if typeof(v)<:Array]
                filter(p->!isnothing(p.second), tmp) |> Dict
            end
        else
            # merging every variable following time dimension v
            [cln_out[k] = cat(cln_out[k], cln[k]; dims=v) for (k,v) ∈ catvar if !isempty(v)]
        end
    end
    return cln_out
end
function readCLNFile(nfile::String; modelreso=false, altfile=nothing)
    @assert isfile(nfile) error("$nfile cannot be found!")
    if contains(nfile, "categorize")
        println("reading Categorize file")
    elseif contains(nfile, "classific")
        println("reading Classification file")
    else
        error("$nfile does not apear to be a Cloudnet file!")
    end

    # Categorize variables to read:
    vars_categorize = Dict(
        # General data
        :alt => "altitude",
        :lat => "latitude",
        :lon => "longitude",
        # RADAR
        :Ze => "Z",
        :V => "v",
        :σV => "v_sigma",
        :ωV => "width",
        # LIDAR
        :β => "beta",
        :δ => "depol",
        #MWR
        :LWP => "lwp",
        #CLOUDNETpy
        :P_INSECT => "insect_prob",
        :QUALBITS => "quality_bits",
        :CATEBITS => "category_bits",
        #MODEL
        :model_height => "model_height",
        :T => "temperature",
        :Tw => "Tw",
        :Pa => "pressure",
        :QV => "q",
        :UWIND => "uwind",
        :VWIND => "vwind",
        :gas_atten => "radar_gas_atten",
        :liq_atten => "radar_liquid_atten",
    )


    # Defining Output varaible:
    var_output = Dict{Symbol, Any}()

    # Starting reading CloudNet files:
    NCDataset(nfile; format=:netcdf4_classic) do nc

        var_output[:time] = convert_time2DateTime(nc)

        var_output[:height] = nc["height"][:]

        # Reading some global attributes:
        if haskey(nc.attrib, "software_version")
            var_output[:version] = nc.attrib["software_version"]
            var_output[:algorithm] = "tropos"
            vars_categorize[:QV] = "specific_humidity"
        end
    
        if haskey(nc.attrib, "cloudnetpy_version")
            var_output[:version] = nc.attrib["cloudnetpy_version"]
            var_output[:algorithm] = "cloudnetpy"
        end
      
        for (inkey, x) ∈ vars_categorize
            #x = vars_categorize[inkey]
            !haskey(nc, x) && continue

            tmp = nc[x][:,:]
            
            if haskey(nc[x].attrib, "missing_value")
                miss_val = nc[x].attrib["missing_value"]
            elseif haskey(nc[x].attrib, "_FillValue")
                miss_val = nc[x].attrib["_FillValue"]
            else
                miss_val = 9.96921f36
            end

            # Cleaning missing values from variables :
            
            varout = let dd = size(tmp)
                isempty(dd) ? NaN : fill(NaN, dd)
            end
            if typeof(tmp) <: Number
                varout = tmp
                
            elseif eltype(tmp) <: Union{Missing, AbstractFloat}
                idxnan = .!ismissing.(tmp)
                varout[idxnan] .= tmp[idxnan]
                varout[varout .≈ miss_val] .= NaN
                
            elseif eltype(tmp) <: AbstractFloat
                idxnan = tmp .≈ miss_val
                varout[.!idxnan] .= tmp[.!idxnan]
            else
                # in case of variables with strings
                varout = tmp
            end
            
            var_output[inkey] = varout
        end

        # If modelreso = true, interpolate model data to cloudnet resolution:

        haskey(nc, "model_time") && (var_output[:model_time] = convert_time2DateTime(nc, time_var="model_time"))

        #(model_time = float.(nc["model_time"]))
        #model_time = float.(nc["model_time"])
        #model_height = nc["model_height"][:]
        if !modelreso
            
            #var_output[:model_height] = model_height
        else
            # To convert model_time to cloudnet time, model_time needs to be
            # a Vector with elements containing the fraction of day:
            model_time = let x = haskey(var_output, :model_time)
                x ? var_output[:model_time] : var_output[:time]
            end
            
            ConvertModelResolution(var_output,
                                   model_time,
                                   var_output[:model_height])
            #    cln_time=tmp_time)
        end

    end

        
    # For Classification dataset:
    # Classification variables to read
    vars_classific = Dict(
        :CLASSIFY => "target_classification",
        :DETECTST => "detection_status",
    )

    classfile = ifelse(isnothing(altfile),
                       replace(nfile, "categorize" => "classification"),
                       altfile)

    if isfile(classfile)
    
        NCDataset(classfile; format=:netcdf4_classic) do nc
            [var_output[k] = nc[v][:,:] for (k,v) ∈ vars_classific if haskey(nc, v)]
            #[var_output[x] = nc[vars_classific[x]][:,:] for x ∈ keys(vars_classific)];
        
        end
    else
        @warn("Classification $(classfile) cannot be found and not loaded!")
    end
    
    return var_output;

end  # end of function
# ----/

# ***************************************************************
# Read Cloudnet Ice and Liquid Effective radius files
"""
 Read Cloudnet Effective Radius for Ice and Water files
The function reads 'infile' assuming it is a Der file and it assumes
that the Ier file is in the same directory. Otherwise Ier file needs to
be explicitly given as optional argument 'altfile=""' as follow:

> reff = readReffFile(infile)
> reff = readReffFile(infile; altfile="ier_file.nc")

WHERE:
infile::String Input file for Der or Ier cloudnet file
altfile::String (Optional) if String given should be a Der file of Ier file

For legacy Cloudnet files, the reff_Fischer and reff_ice files needs to be
explicitly assigned to 'infile' and 'altfile', respectively.

OUTPUT:
reff::Dict Dictionary contaning the data collected from both files Der and Ier.
"""
function readReffFile(nfile::String; altfile::String="")

    der_file, ier_file = if contains(nfile, "der") && isempty(altfile)
        nfile, replace(nfile, "der"=>"ier")
    elseif contains(nfile, "ier") && isempty(altfile)
        replace(nfile, "ier"=>"der"), nfile
    elseif isfile(nfile) && isfile(altfile) 
        nfile, altfile
    elseif isfile(nfile) && isempty(altfile)
        nfile, nfile
    else
        @warn "Are you sure $nfile is a Der or Ier cloudnet file?"
        
    end
    
    @assert isfile(der_file) error("$(der_file) cannot be found!")


    # Defining Output varaible:
    var_output = Dict{Symbol, Any}()

    # Reading first D_er data files:
    NCDataset(der_file; format=:netcdf4_classic) do nc
        var_output[:time] = convert_time2DateTime(nc)

        isfmi = haskey(nc.attrib, "cloudnetpy_version")
        !isfmi && @warn("It seems not a Der Cloudnetpy file! ...assuming legacy Cloudnet!")
        
        vars_cloudnet = Dict(
            :height => "height",
            :Dₑᵣ => isfmi ? "der_scaled" : "r_eff_Frisch_scaled",
            :δDₑᵣ => isfmi ? "der_scaled_error" : "r_eff_Frisch_scaled_error",
            :dₑᵣ => isfmi ? "der" : "r_eff_Frisch",
            :δdₑᵣ => isfmi ? "der_error" : "r_eff_Frisch_error",
            :Nₗ => "N_scaled",
            :Dflag => isfmi ? "der_retrieval_status" : "retrieval_status"
        )

        for (inkey, x) ∈ vars_cloudnet
            
            !haskey(nc, x) && continue
            
            tmp = if eltype(nc[x]) <: Integer
                nc[x][:,:]
            else
                nomissing(nc[x][:,:], NaN32)
            end
            
            if haskey(nc[x].attrib, "missing_value")
                miss_val = nc[x].attrib["missing_value"]
            elseif haskey(nc[x].attrib, "_FillValue")
                miss_val = nc[x].attrib["_FillValue"]
            else
                miss_val = 9.96921f36
            end
            
            # Cleaning missing values from variables :
            #eltype(tmp) <: Union{Missing, AbstractFloat} && (tmp[tmp .≈ miss_val] .= NaN32)
            #eltype(tmp) <: AbstractFloat && (tmp[tmp .≈ miss_val] .= NaN32)
            
            var_output[inkey] = tmp
        end

    end

    # For I_er data files:

    @assert isfile(ier_file) @warn("$(ier_file) cannot be found!")
    
    # Reading for i_er data files:
    NCDataset(ier_file; format=:netcdf4_classic) do nc

        isfmi = haskey(nc.attrib, "cloudnetpy_version")

        !isfmi && @warn("It seems not a Ier Cloudnetpy file! ...assuming legacy Cloudnet!")
        
        #var_output[:time] = convert_time2DateTime(nc)
        vars_cloudnet = Dict(
            :iₑᵣ  => isfmi ? "ier_inc_rain" : "reff_ice_inc_rain",
            :Iₑᵣ  => isfmi ? "ier" : "reff_ice",
            :δIₑᵣ => isfmi ? "ier_error" : "reff_ice_error",
            :Iflag=> isfmi ? "ier_retrieval_status" : "reff_ice_retrieval_status",
        )

        for (inkey, x) ∈ vars_cloudnet

            !haskey(nc, x) && continue
            
            tmp = if eltype(nc[x]) <: Integer
                nc[x][:,:]
            else
                nomissing(nc[x][:,:], NaN32)
            end
            
            var_output[inkey] = tmp
        end

    end

    
    return var_output
end
# ----/


# end of file
