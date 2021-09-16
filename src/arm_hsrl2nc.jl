function hsrl2nc(lidar_file::String, output_path::String)

    # script to read the ARM NSA HSRL data to convert to CloudNet lidar input
    # Part of (AC)3 B07 project.

    
    # Reading input file:
    if !isfile(lidar_file)
        error("$lidar_file does not exist!")
    end

    lidar = ARMtools.getLidarData(lidar_file)

    
    time = lidar[:time];
    #arm_year = year(time[1])
    #arm_month = month(time[1])
    #arm_day = day(time[1])
    range = lidar[:height]
    beta = lidar[:β_raw]  # 1/(m sr) plot scale logarithmic

    std_beta = lidar[:SNR]

    depol_c = lidar[:δ];
    
    file_time = hour.(time) + minute.(time)/60 + second.(time)/3600;
    file_uuid = string(uuid1());
    file_history = string(now(), " - hsrl file created");
    lidar_source = join(lidar[:location])
    lidar_location = join(lidar[:instrumentmodel])
    lidar_doi = haskey(lidar, :doi) ? join(lidar[:doi]) : "none" ;
    λ_nm = 510;  # [nm]  this is from documentation, ncfile doesn't have it
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

    function calc_β(beta, range; std_beta=nothing, noise_params=nothing, depol_c=nothing)
    
        range_square = (range*1e-3).^2;  # [km²]
        range_square[1] = range[1] ≈ 0 ? range_square[2] : range_square[1]
        beta_new = Array(beta./range_square);
        δ = copy(beta_new)
        if !isnothing(std_beta)
            SNR = beta./std_beta
            beta_new[SNR .< 5] .= NaN
            # convert circular- to linear-depolarization
            if !isnothing(depol_c)
                δ = depol_c./(2.0 .+ depol_c)
                δ[SNR .< 5] .= NaN
            end
        else
            mystd(x) = std(x[.!isnan.(x)]);
            myvar(x) = var(x[.!isnan.(x)]);
            mymean(x) = mean(x[.!isnan.(x)]);

            signal_var = mapslices(myvar, beta[end-noise_params[1]:end, :], dims=1);
            is_saturation = findall(signal_var[:] .> noise_params[2]);

    
            # calculating noise_floor for backscattering signal:
            noise_std = mapslices(mystd, beta_new[end-noise_params[1]:end, :], dims=1);
            noise_ave = mapslices(mymean, beta_new[end-noise_params[1]:end, :], dims=1);

            #noisefloor = 2*noise_ave + 9*noise_std;
            noisefloor = noise_ave + noise_std;
            beta_new[beta_new .< noisefloor] .= NaN;
            beta_new = reset_low_values_above_saturation(beta_new, is_saturation, noise_params[4][1]);
        end
        β = copy(beta_new);
        β .*= range_square;


        return isnothing(depol_c) ? β : β, δ
    end

    if @isdefined std_beta
        β_att, δ = calc_β(beta, range, std_beta = std_beta, depol_c=depol_c)
    else
        # the following parameters are empirically for HSRL
        lidar_noise_params = (100, 1e-13, 1e-9, (1.1e-9, 2.9e-8));

        β_att = calc_β(beta, range, noise_params = noise_params);
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
    ARM_OUTPATH = output_path;
    arm_year = year(lidar[:time][1])
    arm_month = month(lidar[:time][1])
    arm_day = day(lidar[:time][1])
    
    ARM_OUTFILE = @sprintf("arm_lidar_%04d%02d%02d.nc", arm_year, arm_month, arm_day);
    outputfile = joinpath(ARM_OUTPATH, ARM_OUTFILE);
    ds = NCDataset(outputfile, "c", attrib = OrderedDict(
        "Conventions"               => "CF-1.7",
        "cloudnetpy_version"        => "1.3.2",
        "file_uuid"                 => file_uuid,
        "cloudnet_file_type"        => "lidar",
        "title"                     => "HSRL made by Julia",
        "year"                      => Int16(arm_year),
        "month"                     => Int16(arm_month),
        "day"                       => Int16(arm_day),
        "location"                  => lidar_location,
        "history"                   => file_history,
        "source"                    => lidar_source,
        "references"                => lidar_doi, 
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

    # Define variables
    ncbeta_raw[:] = β_raw;
    ncbeta[:] = β;
    ncbeta_smooth[:] = β_smooth;
    ncrange[:] = range;
    nctime[:] = file_time;
    nctilt_angle[:] = 0;
    ncheight[:] = 15 .+ range; #... PSG to be fixed
    ncwavelength[:] = λ_nm;

    close(ds)

    return nothing
end
# end of Function
###


