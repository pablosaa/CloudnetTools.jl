#=
********************************************************************
File contains reading functions for the CloudnetTools.jl package
********************************************************************
=#

# **************************************************************
# Auxiliary function to retrieve NCDataset variables based on dimension.
"""
Function to retrieve netCDF variable considering variable dimensions.
The function will return then either Matrix, Vector, or scalar.

```julia-repl
julia> var = getNCvariable(nc::Dataset, varname::String)
```
WHERE:
* ```nc::Dataset``` the netCDF to read pointer,
* ```varname::String``` the variable name to be read.
OUTPUT:
* ```var::Union{Matrix, Vector, Number}``` the retrieved Matrix, Vector, or Scalar from netCDF.

"""
function getNCvariable(nc::NCDataset, var::String)::Union{Matrix, Vector, Number}
    ndx = ndims(nc[var])
    data = if ndx==2
        nc[var][:,:]
    elseif ndx==1
        nc[var][:]
    else
        nc[var][]
    end

    # Cleaning the data from missing_values, _FillValues, NaN:
    miss_val = nothing
    if eltype(data) <: Union{Missing, AbstractFloat} && isa(data, AbstractArray)
        if haskey(nc[var].attrib, "missing_value") || haskey(nc[var].attrib, "_FillValue")
            data = nomissing(data, NaN32) #miss_val = nc[var].attrib["missing_value"]

            #elseif haskey(nc[var].attrib, "_FillValue")
            #    miss_val = nc[var].attrib["_FillValue"]

        else
            miss_val = ifelse(eltype(data) <: AbstractFloat, 9.96921e36, nothing)
        end
    end

    # Cleaning missing values from variables :
    if !isnothing(miss_val) && isa(data, AbstractArray)
        data[data .≈ miss_val] .= NaN32
        !contains(var, "height") && @warn("No 'missing_value' nor '_FillValue' attribute found in variable '$(var)'. Default used: $(miss_val)")                
    end
    
    #if eltype(data) <: Union{Missing, AbstractFloat} && isa(data, AbstractArray)
    #    data = nomissing(data, NaN32)

    #elseif !isnothing(miss_val) && isa(data, AbstractArray)
        
    #end

    return data
end
# ----/

# **************************************************************
"""
Function the retrieve the units of netCDF variable:
```julia-repl
julia> units = getNCunits(nc::Dataset, varname::String)
```
WHERE:
* ```nc::Dataset``` the netCDF to read pointer,
* ```varname::String``` the variable name to be read.
OUTPUT:
* ```units::String``` the netCDF variable's unit, if not existing returns "-".

"""
function getNCunits(nc::NCDataset, var::String)::String

    # Reading the Units of variables:
    einheit = if haskey(nc[var].attrib, "units")
        nc[var].attrib["units"] |> x-> x=="1" ? "-" : x
    else
        "-"
    end

    einheit = replace(einheit, "-1"=>"⁻¹", "-2"=>"⁻²", "-3"=>"⁻³")
    
    return einheit
end
# ----/

# ***************************************************************
# Read Cloudnet Ice water content IWC file
"""
# Read Cloudnet Ice water content IWC file

```julia-repl
julia> iwc = readIWCFile(nfile::String; inc_rain=false, apply_flag=[1:3...])
```
WHERE:
* ```nfile::String``` Full paht to Cloudnet IWC file,
* ```iwc::Dict``` Data containing :time, :height, :iwc, :flag

OPTIONAL:
* ```inc_rain::Bool``` Whether or not use the IWC Cloudnet variable ```iwc_inc_rain``` (Default true),
* ```apply_flag::Vector{Int}``` This vector indicates flagged values to be included (Default include all)

To have the definition of the flag values, see Metadata in IWC Cloudnet file.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function readIWCFile(nfile::String; inc_rain=true, apply_flag::Vector{Int}=Int[])

    @assert isfile(nfile) error("$nfile cannot be found!")

    vars_categorize = Dict(
        :height => "height",
        :iwc => inc_rain ? "iwc_inc_rain" : "iwc",
        :flag => "iwc_retrieval_status",
    )

    # Defining Output varaible:
    var_output = Dict{Symbol, Any}()

    NCDataset(nfile; format=:netcdf4_classic) do nc
        var_output[:time] = convert_time2DateTime(nc)
        var_output[:units] = (time="Time UTC [hour]",)
        
        for (inkey, x) ∈ vars_categorize

            var_output[inkey] = getNCvariable(nc, x)
            # Adding NamedTuple :units e.g. (time="UTC", (:lwc=>"kg m-3",)...) --> (time="UTC", lwc="kg m-3")
            var_output[:units] = (; var_output[:units]..., (inkey=>getNCunits(nc, x),)... )
        end

    end

    # Checking variable flag:
    # Only pass values that are indicated by apply_flag, otherwise values are converted to NaN32.
    if !isempty(apply_flag)
        idxnan = findall(x->x ∉ apply_flag, var_output[:flag])
        var_output[:iwc][idxnan] .= NaN32
    end
    
    return var_output
end
# ----/

# ***************************************************************
# Read Liquid water content LWC file
"""
# Read Liquid water content LWC file
```julia-repl
julia> lwc = readLWCFile(nfile::String; apply_flag=[1,3,4])
```
WHERE:
* nfile::String Full paht to Cloudnet LWC file,
* lwc::Dict Data containing :time, :height, :lwc, :LWP, :flag

OPTIONAL:
* ```apply_flag::Vector{Int}``` This vector indicates flagged values to be included (Default include all)

To have the definition of the flag values, see Metadata in LWC Cloudnet file.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function readLWCFile(nfile::String; apply_flag::Vector{Int}=Int[])

    @assert isfile(nfile) error("$nfile cannot be found!")

    vars_categorize = Dict(
        :height => "height",
        :lwc => "lwc",
        :flag => "lwc_retrieval_status",
        :LWP => "lwp",
    )

    # Defining Output varaible:
    var_output = Dict{Symbol, Any}()

    NCDataset(nfile) do nc  # ; format=:netcdf4_classic
        # reading first the dimensions:
        var_output[:time] = convert_time2DateTime(nc)
        var_output[:units] = (time="Time UTC [hour]",)
        
        for (inkey, x) ∈ vars_categorize

            var_output[inkey] = getNCvariable(nc, x)
            
            var_output[:units] = (; var_output[:units]..., (inkey=>getNCunits(nc, x),)... )
        end

    end

    # Additional computation:
    # Only pass values that are indicated by apply_flag, otherwise values are converted to NaN32.
    if !isempty(apply_flag)
        idxnan = findall(x->x ∉ apply_flag, var_output[:flag])
        var_output[:lwc][idxnan] .= NaN32
    end

    return var_output
end
# ----/

# ******************************************************************
# Reading Classification & Categorize files both at once:
"""
# Fuction to read Cloudnet classification or category files:

USAGE:
```julia-repl
julia> classi = readCLNFile(catfile; modelreso=false)
julia> classi = readCLNFile(catfile; altfile="cloudnet_classify.nc")
```
WHERE:
* ```catfile::String``` full path to cloudnet categorization file,
* ```modelreso::Bool``` (Optional) if true, the model variables are
interpolated from hourly resolution to Cloudnet resolution.
* ```altfile::String``` (Optional) alternative classification file.

NOTE: ```catfile``` needs to have the "categorize" in string, then that is
automatically converted to classification file by replacing "categorize"->"classification".
If Classification file has a different name, then use ```altfile``` string.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
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
        :height => "height",
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
        :IWV => "iwv",
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
        var_output[:units] = (time="Time UTC [hour]",)
        
        #var_output[:height] = getNCvariable(nc, "height") #nc["height"][:]

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
            
            !haskey(nc, x) && continue

            var_output[inkey] = getNCvariable(nc, x)
            var_output[:units] = (; var_output[:units]..., (inkey=>getNCunits(nc, x),)... )            
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
            [var_output[k] = getNCvariable(nc, v) for (k,v) ∈ vars_classific if haskey(nc, v)]
        
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
```julia-repl
julia> reff = readReffFile(infile)
julia> reff = readReffFile(infile; altfile="ier_file.nc")
```
WHERE:
* ```infile::String``` Input file for Der or Ier cloudnet file
* ```altfile::String``` (Optional) if String given should be a Der file of Ier file

For legacy Cloudnet files, the reff_Fischer and reff_ice files needs to be
explicitly assigned to 'infile' and 'altfile', respectively.

OUTPUT:
* ```reff::Dict``` Dictionary contaning the data collected from both files Der and Ier.

Tha variables loaded from the Cloudnet NetCDF files by default are:
```julia-repl
:height => "height",
:Dₑᵣ => "der_scaled" or "r_eff_Frisch_scaled"  #(legacy)
:δDₑᵣ => "der_scaled_error" or "r_eff_Frisch_scaled_error"  #(legacy)
:dₑᵣ => "der" or "r_eff_Frisch"  #(legacy)
:δdₑᵣ => "der_error" or "r_eff_Frisch_error",  #(legacy)
:Nₗ => "N_scaled",
:Dflag => "der_retrieval_status" or "retrieval_status" #(legacy)
```
The variable names marked as (legacy) are refered for the old MATLAB Cloudnet version.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
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
        var_output[:units] = (time="Time UTC [hour]",)
        
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
            
            var_output[inkey] = getNCvariable(nc, x)
            var_output[:units] = (; var_output[:units]..., (inkey=>getNCunits(nc, x),)... )
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
            
            var_output[inkey] = getNCvariable(nc, x)
            var_output[:units] = (; var_output[:units]..., (inkey=>getNCunits(nc, x),)... )
        end

    end

    
    return var_output
end
# ----/


# end of file
