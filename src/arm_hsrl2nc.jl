function hsrl2nc(lidar_file::String, output_path::String; extra_params=Dict())
                 #SITE="not-defined", altitude_m=0f0, tilt_angle=0f0, λ_nm=510)

    # script to read the ARM NSA HSRL data to convert to CloudNet lidar input
    # Part of (AC)3 B07 project.
    
    # Reading input file:
    if !isfile(lidar_file)
        error("$lidar_file does not exist!")
    end
    
    # reading data:
    lidar = ARMtools.getLidarData(lidar_file)

    arm_year = year(lidar[:time][1])
    arm_month = month(lidar[:time][1])
    arm_day = day(lidar[:time][1])
    
    time = lidar[:time];
    file_time = datetime24hours(time)
    
    range = lidar[:height]
    beta = lidar[:β_raw]  # 1/(m sr) plot scale logarithmic

    std_beta = lidar[:SNR]

    depol_c = lidar[:δ];
    
    
    # Default wavelength=510 [nm] this is from documentation, ncfile doesn't have it
    λ_nm = get_parameter(lidar, :λ_nm, extra_params, default=510, lims=(300, 2000))

    #λ_nm = haskey(lidar, :λ_nm) ? lidar[:λ_nm] : λ_nm ;
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
    ## end of reading parameters


    # DEFINING FUNCTIONS
    function interpolate_15s(tt, β_in)
        nh, nt = size(β_in)
    
        tt_out = zeros(2nt)
        tt_out[1:2:end-1] = tt
        tt_out[2:2:end-1] = @. tt[1:end-1] + 0.5*(tt[2:end] - tt[1:end-1])
        tt_out[end] = tt[end] + 0.5*(tt[end] - tt[end-1])

        β_out = zeros(nh, 2*nt)
        β_out[:, 1:2:end-1] = β_in
        β_out[:, 2:2:end-1] = 0.5*(β_in[:, 1:end-1] .+ β_in[:, 2:end])
        β_out[:, end] = β_in[:, end]

        return tt_out, β_out
    end

    function estimate_clouds_from_beta(beta)
        cloud_limit = 1e-6
        cloud_idx = findall(beta .> cloud_limit)
        return cloud_idx, beta[cloud_idx], cloud_limit
    end
    
    # Funciton definitions:
    # Calculates Gaussian peak std parameters.
    # The amount of smoothing is hard coded. This function calculates
    # how many steps in time and height corresponds to this smoothing.

    # Args:
    #    time (ndarray): 1D vector (fraction hour).
    #    range_instru (ndarray): 1D vector (m).

    # Returns:range
    #    tuple: Two element tuple containing number of steps in time and height
    #        to achieve wanted smoothing.
    function calc_sigma_units(time, range_instru)
        mdiff(x) = mean(diff(x))
        minutes_in_hour = 60
        sigma_minutes = 2
        sigma_metres = 5
        time_step = mdiff(time) * minutes_in_hour
        alt_step = mdiff(range_instru)
        y_std = sigma_minutes / time_step
        x_std = sigma_metres / alt_step
        return x_std, y_std
    end


    # signal noisefloor = 2* mean(β(:,end-50:end)) + 5* std(β(:,end-50:end))
    # below this noisefloor is β(β < noisefloor) = nan;

    function reset_low_values_above_saturation(β_in, is_saturation, saturation_noise)
        (nh, ) = size(β_in);
        Δr = ceil(Int, nh/200);
        β_out = copy(β_in);
        
        for idx in is_saturation
            profile = β_out[:, idx];
        
            for i in (nh:-1:2*Δr)
                i_0 = max(i-Δr, 1)
                i_n = min(i+Δr, nh)
                all(isnan.(profile[i_0:i_n])) ? continue : 1
            
                idx_low = findall(profile[i_0:i_n] .≤ 1e-4);
                isempty(idx_low) ? continue : 1
                non_nan = findall(.!isnan.( profile[i_0:i_n] ) );
                # percentage of non NaN low values:
                n_per = length(idx_low)/length(non_nan)
            
                weight = 1 .- abs.((i .- (i_0:i_n))./Δr)
                if n_per > 0.95
                    garbage = mean(weight .* profile[i_0:i_n])
                    β_out[i, idx] = garbage
                else
                    #β_out[i, idx] = NaN
                end
                
            end     
        end 
        return β_out
    end

    if @isdefined std_beta
        β_att, δ = ARMtools.calc_lidar_β(beta, range, std_beta = std_beta, depol_c=depol_c)
    else
        # the following parameters are empirically for HSRL
        lidar_noise_params = (100, 1e-13, 1e-9, (1.1e-9, 2.9e-8));

        β_att = ARMtools.calc_lidar_β(beta, range, noise_params = noise_params);
    end

    ###return file_time, range, β, δ
    
    # Interpolate the time dimension (also change length of file_time):
    _, β_raw = interpolate_15s(file_time, beta)
    file_time, β = interpolate_15s(file_time, β_att)

    # Calculate the smoothed backscattering:
    β_smooth = copy(β);
    cloud_idx, cloud_values, cloud_limit = estimate_clouds_from_beta(β);
    β_smooth[cloud_idx] .= cloud_limit;
    σ_d = calc_sigma_units(file_time, range)
    β_smooth = imfilter(β_smooth, Kernel.gaussian(σ_d));
    β_smooth[cloud_idx] = cloud_values;

    #β_smooth = calc_β(β_smooth, range, noise_params);

    # For storing as netCDF:
    β[isnan.(β)] .= NCDatasets.fillvalue(eltype(β));
    β_smooth[isnan.(β_smooth)] .= NCDatasets.fillvalue(eltype(β));

    # Creating cloudnetpy input netCDF file:
    
    ARM_OUTFILE = @sprintf("%04d%02d%02d_%s_lidar.nc", arm_year, arm_month, arm_day, SITE);
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

    ds.dim["time"] = length(file_time); # unlimited dimension
    ds.dim["range"] = length(range);

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

    ncdepol = defVar(ds,"depol", Float32, ("range", "time"), attrib = OrderedDict(
        "units"                     => "-",
        "long_name"                 => "lidar linear depolarization ratio",
        "comment"                   => "Only present for HSRL. No smoothing.",
        "missing_value"             => NCDatasets.fillvalue(eltype(δ)),
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
    ncbeta_raw[:] = β_raw;
    ncbeta[:] = β;
    ncbeta_smooth[:] = β_smooth;
    ncdepol[:] = δ;
    ncrange[:] = range;
    nctime[:] = file_time;
    nctilt_angle[:] = first(tilt_angle);
    ncheight[:] = altitude_m .+ range; #... PSG to be fixed
    ncwavelength[:] = λ_nm;
    ncaltitude[:] = altitude_m
    close(ds)

    return nothing
end
#----/
# end of Function


