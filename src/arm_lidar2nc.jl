function lidar2nc(lidar_file::String, output_path::String; SITE="not-defined", altitude_m=0f0, tilt_angle=0f0, λ_nm=910)

    if !isfile(lidar_file)
        error("$lidar_file does not exist!")
    end

    # Reading input ARM netCDF file:
    lidar = ARMtools.getLidarData(lidar_file)
    
    time = lidar[:time];
    arm_year = year(time[1])
    arm_month = month(time[1])
    arm_day = day(time[1])
    range = lidar[:height];
    beta = lidar[:β_raw]; # converted from 1/(sr km 10000) to 1/(sr m)
    # Default wavelength=910nm obtained from the ARM ceilometer handbook [nm]
    λ_nm = haskey(lidar, :λ_nm) ? lidar[:λ_nm] : λ_nm
    tilt_angle = max(tilt_angle, maximum(lidar[:TILT]));
    altitude_m = max(lidar[:ALT], altitude_m);
    lidar_location = join(lidar[:location]);
    lidar_source = join(lidar[:instrumentmodel]);
    lidar_doi = join(lidar[:doi]);
    file_uuid = string(uuid1());
    file_history = string(now(), " - ceilometer file created");
    file_time = hour.(time) + minute.(time)/60 + second.(time)/3600;
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
        "title"                     => lidar_location,
        "year"                      => Int16(arm_year),
        "month"                     => Int16(arm_month),
        "day"                       => Int16(arm_day),
        "location"                  => lidar_location,
        "history"                   => file_history,
        "source"                    => lidar_source,
        "references"                => lidar_doi,
    ))

    # Dimensions

    ds.dim["time"] = length(file_time); # unlimited dimension 2804
    ds.dim["range"] = length(range); # 1024

    # Declare variables

    ncbeta_raw = defVar(ds,"beta_raw", Float32, ("range", "time"), attrib = OrderedDict(
        "units"                     => "sr-1 m-1",
        "long_name"                 => "Raw attenuated backscatter coefficient",
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
        "long_name"                 => "Smoothed attenuated backscatter coefficient",
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
    ))

    nctilt_angle = defVar(ds,"tilt_angle", Float32, (), attrib = OrderedDict(
        "units"                     => "degrees",
        "long_name"                 => "Tilt angle from vertical",
    ))

    ncheight = defVar(ds,"height", Float32, ("range",), attrib = OrderedDict(
        "units"                     => "m",
        "long_name"                 => "Height above mean sea level",
        "comment"                   => "Height grid from the mean sea level towards zenith.",
    ))

    ncwavelength = defVar(ds,"wavelength", Float32, (), attrib = OrderedDict(
        "units"                     => "nm",
        "long_name"                 => "laser wavelength",
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
    nctilt_angle[:] = tilt_angle; #...
    ncheight[:] = altitude_m .+ range; #...
    ncwavelength[:] = λ_nm; #...
    ncaltitude[:] = altitude_m;
    close(ds)

    return nothing
end
# ----/
# end of function

