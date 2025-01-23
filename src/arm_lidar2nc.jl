"""
Function to convert ARM lidar/ceilometer file to Cloudnet lidar input file.
USAGE:
```julia-repl
julia> lidar2nc(lidar_file, output_path)
julia> lidar2nc(lidar_file, output_path; extra_params=extras)
julia> lidar2nc(list_files, output_path)
julia> lidar2nc(data, output_path)
```
WHERE:
* ```lidar_file::String``` full path to ARM lidar netCDF file,
* ```output_path::String``` path to put the converted file,
* ```list_files::Vector{String}``` several ARM files to be concatenated,
* ```data::Dict``` dataset readed by ARMtools.getLidarData(armfile) or getCeil10mData()
* ```extra_params::Dict``` (Optional) dictionary with alternative parameter to pass.

The ```extra_params::Dict``` could be for example values not present in the netCDF file or values present but need to be adjusted like Latitude and Longitude for measurements in a Research Vessel, e.g. ```Dict(:SITE=>"Polarstern", :altitude_m=>10f0, :tilt_angle=>5f0, :λ_nm=>910)```, etc.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""

function lidar2nc(lidar_file::String, output_path::String; extra_params=Dict{Symbol,Any}())
    if !isfile(lidar_file)
        error("$lidar_file does not exist!")
    end

    # Reading input ARM netCDF file:
    lidar = ARMtools.getLidarData(lidar_file)

    return lidar2nc(lidar, output_path, extra_params=extra_params)
end
function lidar2nc(lidar_file::Vector{String}, output_path::String; extra_params=Dict{Symbol,Any}())

    # Reading input ARM netCDF file:
    lidar = ARMtools.getLidarData(lidar_file)

    return lidar2nc(lidar, output_path, extra_params=extra_params)
end
function lidar2nc(lidar::Dict, output_path::String; extra_params=Dict{Symbol,Any}())
                  

    
    time = lidar[:time];
    file_time = datetime24hours(time)
    #        hour.(time) + minute.(time)/60 + second.(time)/3600;

    arm_year = year(time[1])
    arm_month = month(time[1])
    arm_day = day(time[1])

    range = lidar[:height];
    beta = lidar[:β_raw]; # converted from 1/(sr km 10000) to 1/(sr m)

    # Default wavelength=910nm obtained from the ARM ceilometer handbook [nm]
    λ_nm = get_parameter(lidar, :λ_nm, extra_params, default=910, lims=(300, 2000))
    #haskey(lidar, :λ_nm) ? lidar[:λ_nm] : λ_nm
    
    tilt_angle = get_parameter(lidar, :TILT, extra_params, default=0)

    altitude_m = get_parameter(lidar, :alt, extra_params, default=0)

    SITE = get_SITE(lidar, extra_params, inkeys=(:site, :campaign, :location))

    # finding DOI from data file or as extra parameter?
    doi = get_parameter(lidar, :doi, extra_params, default="none") ;

    # title from the input data file?
    arm_title = get_parameter(lidar, :title, extra_params, default="")

    # generating UUID for the file:
    file_uuid = string(uuid1());

    # generating file history:
    file_history = get_parameter(lidar, :history, extra_params,
                                 default="Created by Julia Lang (CloudnetTools.jl)"*string(", ", today() ));


    lidar_source = join(lidar[:instrumentmodel]);


    # Adaptation for CloudNetpy 
    noise_params = (100, 1e-12, 2e-7, (1.1e-8, 2.9e-8)) #(100, 1e-12, 3e-6, (1.1e-8, 2.9e-8))
    range_square = (range*1e-3).^2;  # [km²]
    beta_new = Array(1e-7*beta./range_square);  # converting to [sr⁻¹ m⁻¹] 
    signal_var = var(beta_new[end-noise_params[1]:end,:], dims=1);
    is_saturation = findall(signal_var[:] .< noise_params[2]);
    noise = std(beta_new[end-noise_params[1]:end,:], dims=1);
    noise_min = noise_params[3];
    noise[noise .< noise_min] .= noise_min;
    
    function reset_low_values_above_saturation(beta_in, is_saturation, saturation_noise)
        beta_new = Array(beta_in);
        for sat_prof in is_saturation
        
            profile = beta_in[:, sat_prof]
            peak_ind = argmax(profile)
            alt_ind = findall(profile[peak_ind:end] .< saturation_noise) .+ peak_ind .- 1
            #beta_new = Array(beta_in);
            beta_new[alt_ind, sat_prof] .= NaN
            #return beta_new
        end 
        return beta_new
    end
    beta_new = reset_low_values_above_saturation(beta_new, is_saturation, noise_params[4][1]);
    snr = beta_new./noise;
    snr_limit = 5;
    β = beta_new; #1e-7*beta;
    ind_snr_limit = findall(snr .< snr_limit);
    β[ind_snr_limit] .= NaN;
    β .*= range_square;
    β[isnan.(β)] .= NCDatasets.fillvalue(eltype(β))

    # Creating cloudnetpy input netCDF file:
    # Creating output file for CloudNetpy
    ARM_OUTFILE = @sprintf("%04d%02d%02d_%s_lidar.nc", arm_year, arm_month, arm_day,SITE);
    outputfile = joinpath(output_path, ARM_OUTFILE);
    
    ds = NCDataset(outputfile, "c", attrib = OrderedDict(
        "Conventions"               => "CF-1.7",
        "file_uuid"                 => file_uuid,
        "cloudnet_file_type"        => "lidar",
        "title"                     => arm_title,
        "year"                      => Int16(arm_year),
        "month"                     => Int16(arm_month),
        "day"                       => Int16(arm_day),
        "location"                  => SITE,
        "history"                   => file_history,
        "source"                    => lidar_source,
        "references"                => doi,
    ))

    # Dimensions

    ds.dim["time"] = length(file_time); # unlimited dimension 2804
    ds.dim["range"] = length(range); # 1024

    # Declare variables

    ncbeta_raw = defVar(ds,"beta_raw", Float32, ("range", "time"), attrib = OrderedDict(
        "units"                     => "sr-1 m-1",
        "long_name"                 => "Attenuated backscatter coefficient",
        "comment"                   => "Range corrected, attenuated backscatter.",
        "missing_value"             => NCDatasets.fillvalue(eltype(β)),
    ))

    ncbeta = defVar(ds,"beta", Float32, ("range", "time"), attrib = OrderedDict(
        "units"                     => "sr-1 m-1",
        "long_name"                 => "Attenuated backscatter coefficient",
        "comment"                   => "Range corrected, SNR screened, attenuated backscatter.",
        "missing_value"             => NCDatasets.fillvalue(eltype(β)),
    ))

    ncbeta_smooth = defVar(ds,"beta_smooth", Float32, ("range", "time"), attrib = OrderedDict(
        "units"                     => "sr-1 m-1",
        "long_name"                 => "Attenuated backscatter coefficient",
        "comment"                   => "Range corrected, SNR screened backscatter coefficient.\n Weak background is smoothed using Gaussian 2D-kernel.",
        "missing_value"             => NCDatasets.fillvalue(eltype(β)),
    ))

    ncrange = defVar(ds,"range", Float32, ("range",), attrib = OrderedDict(
        "units"                     => "m",
        "long_name"                 => "Range from instrument",
        "comment"                   => "Height grid from the instrument towards the line of sight.",
    ))

    nctime = defVar(ds,"time", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "decimal hours since midnight",
        "long_name"                 => "Time UTC",
        "standard_name"             => "time",
    ))

    nctilt_angle = defVar(ds,"zenith_angle", Float32, (), attrib = OrderedDict(
        "units"                     => "degree",
        "long_name"                 => "Zenith angle",
        "standard_name"             => "zenith_angle",
    ))

    ncheight = defVar(ds,"height", Float32, ("range",), attrib = OrderedDict(
        "units"                     => "m",
        "long_name"                 => "Height above mean sea level",
        "standard_name"             => "height_above_mean_sea_level",
        "comment"                   => "Height grid from the mean sea level towards zenith.",
    ))

    ncwavelength = defVar(ds,"wavelength", Float32, (), attrib = OrderedDict(
        "units"                     => "nm",
        "long_name"                 => "Laser wavelength",
    ))

    ncaltitude = defVar(ds,"altitude", Float32, (), attrib = OrderedDict(
        "units"                     => "m",
        "long_name"                 => "Altitude of site",
    ))


    # Define variables

    ncbeta_raw[:] = 1e-7*beta; #...
    ncbeta[:] = β; #...
    ncbeta_smooth[:] = β; # β_smooth isn't used by cloudnetpy!
    ncrange[:] = range; #...
    nctime[:] = file_time; #...
    nctilt_angle[:] = tilt_angle[1]; #...
    ncheight[:] = altitude_m .+ range; #...
    ncwavelength[:] = λ_nm; #...
    ncaltitude[:] = altitude_m;
    close(ds)

    return ARM_OUTFILE
end
# ----/
# end of function

