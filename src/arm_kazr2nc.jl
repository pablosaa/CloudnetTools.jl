"""
Function to convert ARM radar file to Cloudnet radar input file.
USAGE:
```julia-repl
julia> kazr2nc(radar_file, output_path)
julia> kazr2nc(radar_file, output_path; extra_params=extras)
julia> kazr2nc(list_files, output_path)
julia> kazr2nc(data, output_path)
```
WHERE:
* ```radar_file::String``` full path to ARM radar netCDF file,
* ```output_path::String``` path to put the converted file,
* ```list_files::Vector{String}``` several ARM files to be concatenated,
* ```data::Dict``` dataset readed by ARMtools.getKAZRData(armfile)
* ```extra_params::Dict``` (Optional) dictionary with alternative parameter to pass.
* ```rng0::Union{Bool, Float}``` (Optional) extrapolate radar data to rng0, default false.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function kazr2nc(radar_file::String, output_path::String; extra_params=Dict())
    # Reading input ARM netCDF file:
    
    if !isfile(radar_file)
        error("$radar_file does not exist!")
    end 

    radar = ARMtools.getKAZRData(radar_file; extra_params...)

    return kazr2nc(radar, output_path, extra_params=extra_params)
end
function kazr2nc(radar_file::Vector{String}, output_path::String; extra_params=Dict())
    # Reading input ARM netCDF file:
    
    radar = ARMtools.getKAZRData(radar_file; extra_params...)

    return kazr2nc(radar, output_path, extra_params=extra_params)
end

function kazr2nc(radar::Dict, output_path::String; extra_params=Dict(), rng0=false)

    # Script to adapt ARM KAZR radar to netCDF input for cloudnetpy
    # 

    time = radar[:time];
    file_time = datetime24hours(time)
    ntime = length(file_time)
    
    arm_year  = year(time[1])
    arm_month = month(time[1])
    arm_day   = day(time[1])

    # if LDR is not included in data file, usa proxy diff Zxpol-Zcpol:
    !haskey(radar, :LDR) && (radar[:LDR] = radar[:Zxpol] .- radar[:Ze])

    # considering the global attribute radar_frequency has the form: "35.6 GHz"
    # the value will be retrieved as a Vector with two elements e.g. [35.6, "GHz"]
    # The KAZR by default is assumed the frequency 34.89 GHz
    radar_frequency = let x = get_parameter(radar, :radar_frequency, extra_params, default=[34.89, "GHz"])[1]
        typeof(x) <: Number ? x : NaN32
    end

    # Getting the range resolution
    δrg = let xrg=radar[:height][1:10] |> diff |> first
        get_parameter(radar, :drg, extra_params, default=xrg) |> first
    end

    #=
    Checking whether the optional parameter rng0 has been given:
    By default rng0=false thus no extrapolation performed,
    if rng0=true then radar range is extrapolated to δrng/2,
    otherwise if rng0 is a Float32, then is extrapolated to rng0
    =#
    if typeof(rng0)<:Union{Bool, AbstractFloat}
        let rng_first = rng0 && δrg/2
            vh = radar[:height]
            vt = file_time
            ext_h = range(radar[:height][1]-δrg, step=-δrg, stop=rng_first) |> reverse
            ext_t = file_time
            radar[:height] = vcat(ext_h, vh)
            
            for (k,v) ∈ radar
                !(typeof(v)<:Matrix) && continue

                # Interpolating variable v
                ipt = interpolate((vh,vt), v, Gridded(Linear()) )
                ept = extrapolate(ipt, Flat() )
                ext_v = [ept(h, t) for h in ext_h, t in ext_t]
                radar[k] = vcat(ext_v, v)
            end
        end
        
    end
    
    nzrg = length(radar[:height]);  # Number of range gates

    
    # the following four variables are not present in the ARSCL data files:
    # nyquist velocity [m/s] 
    nyquist_velocity = get_parameter(radar, :nyquist_velocity, extra_params, default=-999)

    # [counts] number of fast fourier transform
    num_nfft = haskey(radar, :fft_len) ? radar[:fft_len][1] : -999;
    # [Hz] pulse repetition frequency
    num_prf  = haskey(radar, :prf) ? radar[:prf][1] : -999;
    # [counts] number of spectral average ((10-0)/0.001
    num_nave = haskey(radar, :number_spectral_ave) ? radar[:number_spectral_ave][1] : 9999;

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
        "units"                     => "Time UTC",
        "long_name"                 => "decimal hours since midnight",
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

    ncprf = defVar(ds,"prf", Float32, (), attrib = OrderedDict(
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

    ncrainrate = defVar(ds,"rainfall_rate", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "m s-1",
        "long_name"                 => "Rainfall rate",
        "standard_name"             => "rainfall_rate",
        "comment"                   => "Fill values denote rain with undefined intensity.",
    ))
    
    # Define variables
    nctime[1:ntime] = file_time;
    ncZe[:] = radar[:Ze];
    ncv[:] = radar[:MDV];
    ncwidth[:] = radar[:SPW];
    ncldr[:] = radar[:LDR];
    ncSNR[:] = radar[:SNR];
    nclatitude[:] = get_parameter(radar, :lat, extra_params, lims=(-90, 90))
    nclongitude[:] = get_parameter(radar, :lon, extra_params, lims=(-180, 360))
    ncaltitude[:] = get_parameter(radar, :alt, extra_params, default=0)

    ncrange[:] = radar[:height];
    ncradar_frequency[:] = radar_frequency;
    ncnyquist_velocity[:] = nyquist_velocity[1];
    ncnfft[:] = num_nfft;
    ncprf[:] = num_prf;
    ncnave[:] = num_nave;
    nczrg[:] = nzrg;
    ncrg0[:] = 5;
    ncdrg[:] = δrg;

    # filling rainrate_fall and converting from mm h⁻¹ to m s⁻¹:
    if haskey(radar, :RR)
        ncrainrate[:] = radar[:RR]/3.6f6
    else
        ncrainrate[:] .= NaN32
    end

    
    close(ds)

    return ARM_OUTFILE
end
#----/
#end of function

