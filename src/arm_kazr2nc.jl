function kazr2nc(radar_file::String, output_path::String; extra_params=Dict{Symbol, Any}())
    # Reading input ARM netCDF file:
    
    if !isfile(radar_file)
        error("$radar_file does not exist!")
    end 

    radar = ARMtools.getKAZRData(radar_file)

    return kazr2nc(radar, output_path, extra_params=extra_params)
end
function kazr2nc(radar_file::Vector{String}, output_path::String; extra_params=Dict{Symbol, Any}())
    # Reading input ARM netCDF file:
    
    radar = ARMtools.getKAZRData(radar_file)

    return kazr2nc(radar, output_path, extra_params=extra_params)
end

function kazr2nc(radar::Dict, output_path::String; extra_params=Dict{Symbol, Any}())

    # Script to adapt ARM KAZR radar to netCDF input for cloudnetpy
    # 

    time = radar[:time];
    file_time = datetime24hours(time)
        #hour.(time) + minute.(time)/60 + second.(time)/3600;
    arm_year  = year(time[1])
    arm_month = month(time[1])
    arm_day   = day(time[1])

    # if LDR is not included in data file, usa proxy diff Zxpol-Zcpol:
    !haskey(radar, :LDR) && (radar[:LDR] = radar[:Zxpol] .- radar[:Ze])

    # considering the global attribute radar_frequency has the form: 35.6 GHz
    # the value will be retrieved as a Vector with two elemens (35.6, "GHz")
    radar_frequency = let x=radar[:radar_frequency][1]
        typeof(x) <: Number ? x : NaN32
    end

    nzrg = length(radar[:height]);  # Number of range gates

    # the following four variables are not present in the ARSCL data files:
    # nyquist velocity [m/s] 
    nyquist_velocity = get_parameter(radar, :nyquist_velocity, extra_params, default=-999)
    #haskey(radar, :nyquist_velocity) ? radar[:nyquist_velocity][1] : -999;
    # [counts] number of fast fourier transform
    num_nfft = haskey(radar, :fft_len) ? radar[:fft_len][1] : -999;
    # [Hz] pulse repetition frequency
    num_prf  = haskey(radar, :prf) ? radar[:prf][1] : -999;
    # [counts] number of spectral average ((10-0)/0.001
    num_nave = haskey(radar, :number_spectral_ave) ? radar[:number_spectral_ave][1] : 9999.;

    # Looking for information about instrument site or location/campaign
    SITE = get_SITE(radar, extra_params, inkeys=(:site, :campaign, :location))

    # finding DOI from data file or as extra parameter?
    doi = get_parameter(radar, :doi, extra_params, default="none")

    # title from the input data file?
    arm_title = get_parameter(radar, :title, extra_params, default="")

    # generating UUID for the file:
    file_uuid = string(uuid1());

    # generating file history:
    file_history = get_parameter(radar, :history, extra_params,
                                 default="Created by Julia Lang (CloudnetTools.jl)"*string(", ", today() ));
    

    # Creating output file for CloudNetpy
    ARM_OUTFILE = @sprintf("%04d%02d%02d_%s_radar.nc", arm_year, arm_month, arm_day,SITE);
    outputfile = joinpath(output_path, ARM_OUTFILE);

    ds = NCDataset(outputfile, "c", attrib = OrderedDict(
        "Conventions"               => "CF-1.7",
        "file_uuid"                 => file_uuid,
        "cloudnet_file_type"        => "radar",
        "title"                     => arm_title,
        "year"                      => Int16(arm_year),
        "month"                     => Int16(arm_month),
        "day"                       => Int16(arm_day),
        "location"                  => SITE,
        "history"                   => file_history,
        "source"                    => haskey(radar, :instrumentmodel) ? join(radar[:instrumentmodel]) : "arm-radar",
        "references"                => doi,
    ))

    # Dimensions

    ds.dim["time"] = Inf; # unlimitted dimension 
    ds.dim["range"] = nzrg;

    # Declare variables

    ncZe = defVar(ds,"Zh", Float32, ("range", "time"), attrib = OrderedDict(
        "units"                     => "dBZ",
        "long_name"                 => "Radar reflectivity factor",
        "missing_value"             => NCDatasets.fillvalue(eltype(radar[:Ze])),
        "attenuation_corr_gas"      => 1, #true,
        "attenuation_corr_liq"      => 0, #false,
    ))
    
    ncv = defVar(ds,"v", Float32, ("range", "time"), attrib = OrderedDict(
        "units"                     => "m s-1",
        "long_name"                 => "Doppler velocity",
        "comment"                   => "This parameter is the radial component of the velocity, with positive
    velocities are away from the radar.",
        "ancillary_variables"       => "v_sigma",
        "positive"                  => "up",
        "missing_value"             => NCDatasets.fillvalue(eltype(radar[:MDV])),
    ))

    ncwidth = defVar(ds,"width", Float32, ("range", "time"), attrib = OrderedDict(
        "units"                     => "m s-1",
        "long_name"                 => "Spectral width",
        "comment"                   => "This parameter is the standard deviation of the reflectivity-weighted
    velocities in the radar pulse volume.",
        "missing_value"             => NCDatasets.fillvalue(eltype(radar[:SPW])),
))

    ncldr = defVar(ds,"ldr", Float32, ("range", "time"), attrib = OrderedDict(
        "units"                     => "dB",
        "long_name"                 => "Linear depolarisation ratio",
        "comment"                   => "This parameter is the ratio of cross-polar to co-polar reflectivity.",
        "missing_value"             => NCDatasets.fillvalue(eltype(radar[:LDR])),
))

    ncSNR = defVar(ds,"SNR", Float32, ("range", "time"), attrib = OrderedDict(
        "units"                     => "dB",
        "long_name"                 => "Signal-to-noise ratio",
        "missing_value"             => NCDatasets.fillvalue(eltype(radar[:SNR])),
    ))

    nclatitude = defVar(ds,"latitude", Float32, (), attrib = OrderedDict(
        "units"                     => "degree_north",
        "long_name"                 => "Latitude of site",
        "standard_name"             => "latitude",
    ))

    nclongitude = defVar(ds,"longitude", Float32, (), attrib = OrderedDict(
        "units"                     => "degree_east",
        "long_name"                 => "Longitude of site",
        "standard_name"             => "longitude",
    ))

    ncaltitude = defVar(ds,"altitude", Float32, (), attrib = OrderedDict(
        "units"                     => "m",
        "long_name"                 => "Altitude of site",
    ))

    nctime = defVar(ds,"time", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "decimal hours since midnight",
        "long_name"                 => "Time UTC",
        "standard_name"             => "time",
    ))

    ncrange = defVar(ds,"range", Float32, ("range",), attrib = OrderedDict(
        "units"                     => "m",
        "long_name"                 => "Range from instrument",
        "comment"                   => "Height grid from the instrument towards the line of sight.",
    ))

    ncradar_frequency = defVar(ds,"radar_frequency", Float32, (), attrib = OrderedDict(
        "units"                     => "GHz",
        "long_name"                 => "Radar transmit frequency",
    ))

    ncnyquist_velocity = defVar(ds,"nyquist_velocity", Float32, (), attrib = OrderedDict(
        "units"                     => "m s-1",
        "long_name"                 => "Nyquist velocity",
    ))

    ncnfft = defVar(ds,"nfft", Int32, (), attrib = OrderedDict(
        "long_name"                 => "Number of FFT Points",
        "units"                     => "1", #"count",
    ))

    ncprf = defVar(ds,"prf", Int32, (), attrib = OrderedDict(
        "long_name"                 => "Pulse Repetition Frequency",
        "units"                     => "Hz",
    ))

    ncnave = defVar(ds,"nave", Int32, (), attrib = OrderedDict(
        "long_name"                 => "Number of Spectral Avreages",
        "units"                     => "1", #"count",
    ))

    nczrg = defVar(ds,"zrg", Int32, (), attrib = OrderedDict(
        "long_name"                 => "Number of Range Gates",
        "units"                     => "count",
    ))

    ncrg0 = defVar(ds,"rg0", Int32, (), attrib = OrderedDict(
        "long_name"                 => "Number of Lowest Range Gates",
        "units"                     => "1", #"count",
    ))

    ncdrg = defVar(ds,"drg", Float32, (), attrib = OrderedDict(
        "long_name"                 => "Range Resolution",
        "unit"                      => "m",
    ))


    # Define variables

    ncZe[:] = radar[:Ze];
    ncv[:] = radar[:MDV];
    ncwidth[:] = radar[:SPW];
    ncldr[:] = radar[:LDR];
    ncSNR[:] = radar[:SNR];
    nclatitude[:] = get_parameter(radar, :lat, extra_params, lims=(-90, 90))
    nclongitude[:] = get_parameter(radar, :lon, extra_params, lims=(-180, 360))
    ncaltitude[:] = get_parameter(radar, :alt, extra_params, default=0)
    nctime[:] = file_time;
    ncrange[:] = radar[:height];
    ncradar_frequency[:] = radar_frequency;
    ncnyquist_velocity[:] = nyquist_velocity;
    ncnfft[:] = num_nfft;
    ncprf[:] = num_prf;
    ncnave[:] = num_nave;
    nczrg[:] = nzrg;
    ncrg0[:] = 5;
    ncdrg[:] = haskey(radar,:drg) ? radar[:drg][1] : 30;

    close(ds)

    return ARM_OUTFILE
end
#----/
#end of function

