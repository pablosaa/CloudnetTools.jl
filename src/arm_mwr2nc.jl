function mwr2nc(mwr_file::String, output_path::String)

    mwr = ARMtools.getMWRData(mwr_file)

    SITE = "arm-nsa"
    # Aux variables:
    file_time = @. Second(mwr[:time] - DateTime(2001,1,1,0,0,0));

    # Creating output file for CloudNetpy
    arm_year = year(mwr[:time][1])
    arm_month = month(mwr[:time][1])
    arm_day = day(mwr[:time][1])
    ARM_OUTFILE = @sprintf("%04d%02d%02d_%s_mwr.nc", arm_year, arm_month, arm_day, SITE);
    outputfile = joinpath(output_path, ARM_OUTFILE);
    ds = NCDataset(outputfile, "c", attrib = OrderedDict(
        "netCDF_convention"         => " CF-1.0",
        "radiometer_location"       => "Utqiagvik",
        "radiometer_system"         => "Radiometrics",
        "serial_number"             => " R-xx/xxx",
        "station_altitude"          => "$(mwr[:alt])", #" 21",
        "station_longitude"         => "$(mwr[:lon])", #"9.9�' West",
        "station_latitude"          => "$(mwr[:lat])", #"53.33�' North",
        "comments"                  => "Made it by Julia",
        "radiometer_software_version" => "none",
        "host_PC_software_version"  => "none",
    ))

    # Dimensions

    ds.dim["time"] = Inf # unlimited dimension

    # Declare variables

    nctime_reference = defVar(ds,"time_reference", Int32, (), attrib = OrderedDict(
        "long_name"                 => "flag indicating the time zone reference",
        "units"                     => "unitless",
        "comment"                   => "0 = local time, 1 = UTC",
    ))

    ncnumber_integrated_samples = defVar(ds,"number_integrated_samples", Int32, ())

    ncminimum = defVar(ds,"minimum", Float32, (), attrib = OrderedDict(
        "units"                     => "g / m^2",
    ))

    ncmaximum = defVar(ds,"maximum", Float32, (), attrib = OrderedDict(
        "units"                     => "g / m^2",
    ))

    nctime = defVar(ds,"time", Int32, ("time",), attrib = OrderedDict(
        "long_name"                 => "sample time",
        "units"                     => "seconds since 1.1.2001, 00:00:00",
        "comment"                   => "reference time zone indicated in field time_reference",
    ))

    ncrain_flag = defVar(ds,"rain_flag", Int32, ("time",), attrib = OrderedDict(
        "Info"                      => "0 = No Rain, 1 = Raining",
    ))

    ncelevation_angle = defVar(ds,"elevation_angle", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "viewing_elevation angle",
        "comment"                   => "-90 is blackbody view, 0 is horizontal view (red arrow), 90 is zenith view, 180 is horizontal view (2nd quadrant)",
        "units"                     => "degrees (-90 - 180)",
    ))

    ncazimuth_angle = defVar(ds,"azimuth_angle", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "DEG (0�-360�)",
    ))

    ncretrieval = defVar(ds,"retrieval", Int32, (), attrib = OrderedDict(
        "Info"                      => "0 = Linear Regr., 1 = Quadr. Regr., 2 = Neural Network",
    ))

    ncLWP_data = defVar(ds,"LWP_data", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "g / m^2",
    ))

    ncIWV_data = defVar(ds,"IWV_data", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "kg / m^2",
    ))
    
    # Define variables

    nctime_reference[:] = 1   # dummy value
    ncnumber_integrated_samples[:] = 1   # dummy value
    ncminimum[:] = minimum(mwr[:LWP])
    ncmaximum[:] = maximum(mwr[:LWP])
    nctime[:] = map(x->x.value, file_time);
    ncrain_flag[:] = mwr[:wet];
    ncelevation_angle[:] = mwr[:elevation];
    ncazimuth_angle[:] = mwr[:azimuth];
    ncretrieval[:] = 1;   # dummy value
    ncLWP_data[:] = mwr[:LWP];
    ncIWV_data[:] = mwr[:IWV];

    close(ds)
    
end

# Starting to run over input parameters:
#N_product = input_params["products"] |> !isempty |> length;
#N_years   = input_params["year"] |> !isempty |> length;
#N_months  = input_params["month"] |> !isempty |> length;
#N_days    = input_params["day"] |> !isempty |> length;

#arm_year = input_params["year"]; #2017;
#arm_month = input_params["month"]; #11;
#arm_day = input_params["day"]; #9;

#ARM_PRODUCT = input_params["products"];  # can be "LOS", "HF"
#ARM_SITE = input_params["campaign"];  # e.g. NSA, SGP
#ARM_INSTRUMENT = input_params["instrument"];

# Reading input netCDF file:
#DATA_PATH = @sprintf("/home/psgarfias/LIM/data/utqiagvik-nsa/MWR/%s/%04d", ARM_PRODUCT, arm_year);
#DATA_NAME = @sprintf("nsamwrlosC1.b1.%04d%02d%02d.000030.cdf", arm_year, arm_month, arm_day);
#file_pattern = @sprintf("nsamwrret1liljclouC1.c1.%04d%02d%02d.", arm_year, arm_month, arm_day);
#file_pattern = string(lowercase(ARM_SITE), lowercase(ARM_INSTRUMENT), lowercase(ARM_PRODUCT));
#date_pattern = @sprintf(".%04d%02d%02d.", arm_year, arm_month, arm_day);
#file_pattern = joinpath(DATA_PATH, ARM_PRODUCT, DATA_NAME);
#mwr_filenc = filter(x->occursin(file_pattern, x) & occursin(date_pattern, x), readdir(DATA_PATH, join=true))[1];


##    ARM_PRODUCT = "RET"
### Reading input ARM netCDF file:
##
##if !isfile(mwr_file)
##    error("$mwr_file does not exist!")
##end 
##ncin = NCDataset(mwr_file, "r");
##lat = ncin["lat"]
##lon = ncin["lon"]
##alt = ncin["alt"]
###time = datetime2unix.(ncin["time"]);
##time = ncin["time"];
##arm_year = year(time[1])
##arm_month = month(time[1])
##arm_day = day(time[1])
##file_time = Second.(time .- DateTime(2001,1,1,0,0,0)); #hour.(time) + minute.(time)/60 + second.(time)/3600;
##if ARM_PRODUCT === "RET"
##    liq = ncin["phys_lwp"];   # [g/m²]
##    iwv = ncin["phys_pwv"];   # [cm]
##    elevation = Float32(90) #90*ones(Float32, length(time))
##    azimuth = Float32(0.0) #zeros(Float32, length(time))
##    wet_window = Int32(0) #zeros(Int32, length(time))
##    factor_lwp = 1;
##    factor_iwv = 997E-2;  # [cm] -> [kg m⁻²]
##elseif ARM_PRODUCT === "LOS"
##    elevation = ncin.attrib["elevation"]
##    elevation = parse(Float32, elevation[1:6]);
##    azimuth = ncin.attrib["azimuth"]
##    azimuth = parse(Float32, azimuth[1:6])
##    wet_window = ncin["wet_window"];    
##    liq = ncin["liq"];   # [cm]
##    iwv = ncin["vap"];   # [cm]
##    factor_lwp = 1E4;
##    factor_iwv = 997E-2;  # [cm] -> [kg m⁻²]
##else
##    liq = ncin["lwp"];
##    factor_lwp = 1;
##end
### Converting units for LWP. e.g. liq [cm] into LWP [g/m²] = 1E4*liq;
##lwp = factor_lwp*Array(liq);  #[g m⁻²]
##iwv = factor_iwv*Array(iwv);  #[kg m⁻²]
### Closing input netCDF file:
##close(ncin)
##
