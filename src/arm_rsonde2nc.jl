using ATMOStools
using CloudnetTools


"""
Function to convert the ARM radiosonde data to Cloudnet input file.
USAGE:
```julia-repl
julia> output_file = ARM.rsonde2nc(rs_file, output_path)
julia> output_file = ARM.rsonde2nc(rs_file, output_path; extra_params=extra_params)
julia> output_file = ARM.rsonde2nc(rs, output_path)
julia> output_file = ARM.rsonde2nc(data_path, YY, MM, DD, output_path)
```
WHERE:
```julia-repl
* rs_file::String Full path to ARM file to be converted,
* output_path::String Path to store the converted file,
* extra_params::Dict{Symbol, Any} Extra parameters passed as dictionary :key=>value pairs,
* rs::Dict{Symbol, Any} dictionary containing radiosonde data, see ARMtools.getSondeData()
* data_path::String Base path where to look for radiosonde data, by default looks at /data_path/INTERPOLATEDSONDE/YY/
* YY::Int, MM::Int, DD::Int year, month, and day to look data file,
```
OPTIONAL PARAMETERS:
Extra parameters can be passed to the function by using the extra\\_params optional variable:
```julia-repl
*  extra_params = Dict{Symbol, Any}(:site=>"mosaic", :freq=>35.5f0, :lat=>83.5, :lon=>101.0)
*  extra_params = Dict{Symbol, Any}(:product_name=>"MERGEDSONDE")
```
For the latest case, extra\\_params passes the key :product\\_name which will change the directory
where the function look for ARM radiosonde data, e.g. /data_path/MERGEDSONDE/YY

EXAMPLE:
```julia-repl
julia> ARM.rsonde2nc("/data/nsa/nsainterpolatedsondeC1.c1.20190101.000030.nc", "/cloudnet/inputs")
# To look for radiosonde data product MERGEDSONDE files for the year 2019, january 1st, then convert it:
julia> ARM.rsonde2nc("/data/nsa", 2019, 1, 1, "/cloudnet/inputs", extra_params=Dict{Symbol, Any}(:MERGEDSONDE) )
# To read first ARM radiosonde data with the ARMtools package, then convert it to Cloudnet input file:
julia> rs = ARMtools.getSondeData("/data/nsa/nsainterpolatedsondeC1.c1.20190101.000030.nc");
julia> ARM.rsonde2nc(rs, "/cloudnet/inputs", extra_params = Dict{Symbol, Any}(:site=>"mosaic") )
```

NOTE:
By default the radar attenuation coefficient is assumed Ka-band (35.5 GHz),
however this will be filled by zero since it is not yet supported.

    (c) Pablo Saavedra Garfias
"""
function rsonde2nc(rs_file::String, output_path::String; extra_params=Dict{Symbol, Any}() )
    data = ARMtools.getSondeData(rs_file, addvars=["lat", "lon", "alt"]);
    return rsonde2nc(data, output_path; extra_params=extra_params)
end
# -- OR --
function rsonde2nc(PATH_DATA::String, yy::Int, mm::Int, dd::Int, output_path::String; freq=35.5f0, extra_params=Dict{Symbol, Any}() )
    # RADIOSONDE input files:
    product_name = haskey(extra_params, :product_name) ? extra_params[:product_name] : "INTERPOLATEDSONDE"
    rs_file = ARMtools.getFilePattern(PATH_DATA, product_name, yy,mm,dd); # /2019/mosinterpolatedsondeM1.c1.20191230.000030.nc";
    isnothing(rs_file) && @error "No radiosonde file could be found under $(PATH_DATA) and $(product_name)"

    return rsonde2nc(rs_file, output_path; extra_params=extra_params)
end
# -- OR --
function rsonde2nc(data::Dict, output_path::String; extra_params=Dict{Symbol, Any}() )


    # Checking default values:
    freq = haskey(extra_params, :freq) ? extra_params[:freq] : 35.5f0
    δtime = haskey(extra_params, :δtime) ? extra_params[:δtime] : Hour(1)
    δheight = haskey(extra_params, :δheight) ? extra_params[:δheight] : 30
    
    # RADIOSONDE input files:
    #- rs_file = ARMtools.getFilePattern(PATH_DATA, "INTERPOLATEDSONDE", yy,mm,dd); # /2019/mosinterpolatedsondeM1.c1.20191230.000030.nc";
    #- data = ARMtools.getSondeData(rs_file, addvars=["lat", "lon", "alt"]);
    ##
    ### getting latitude, longitude and altitude:
    LAT, LON, ALT = let tmp=findfirst(>(-100), data[:LAT])
        if !isnothing(tmp) 
            data[:LAT][tmp], data[:LON][tmp], data[:ALT][tmp]
        else
            NaN32, NaN32, NaN32
        end
    end
    ##
    # TROPOS cloudnet categorization:
    #rs_file = ARMtools.getFilePattern("/home/psgarfias/LIM/data/B07/arctic-mosaic/CloudNet/output", "TROPOS/processed/categorize", yy, mm, dd, fileext="categorize.nc");

    #data = CloudnetTools.readCLNFile(rs_file);
    
    idx_rstime = let tmp=CloudnetTools.datetime24hours(data[:time])
        T_RS_launch = haskey(extra_params, :launch_time) ? extra_params[:launch_time] : (0:1:24)
        out = [findmin(abs.(tmp.-T)) |> V-> V[1]<0.5 && V[2] for T ∈ T_RS_launch]
        filter(!=(false), out)
    end
        #floor.(Int32, range(1, stop=length(data[:time]), length=24))
    idx_rslevel = let imax = findlast(<(15), data[:height])
        floor.(Int32, range(1, stop=imax, length=137)) |> unique
    end;

    # dimentions:
    Ntime = length(idx_rstime)
    Nlevel = length(idx_rslevel)
    Nfreq = length(freq)

    # Meteorological variables:
    # LAT, LON, ALT = data[:lat], data[:lon], data[:alt]
    HEIGHT = repeat(1f3data[:height][idx_rslevel], 1, length(idx_rstime) );

    LEVELS = (length(idx_rslevel):-1:1)

    Pa = 1f3data[:Pa][idx_rslevel, idx_rstime];  # 1f3* for radiosonde

    TK = data[:T][idx_rslevel, idx_rstime] .+ 273.15;

    U = try
        data[:UWIND][idx_rslevel, idx_rstime];  # :U for radiosonde
    catch
        data[:U][idx_rslevel, idx_rstime];
    end

    V = try
        data[:VWIND][idx_rslevel, idx_rstime];  # :V for radiosonde
    catch

        data[:V][idx_rslevel, idx_rstime];
    end

    QV = try
        data[:QV][idx_rslevel, idx_rstime];   # :qv for radiosonde
    catch
        data[:qv][idx_rslevel, idx_rstime];
    end

    RH = if haskey(data, :RH)
        data[:RH][idx_rslevel, idx_rstime]
    elseif haskey(data, :rh)
        data[:rh][idx_rslevel, idx_rstime]
    else
        @info "Key :RH not found in data. Calculating relative humidity"
        ATMOStools.qv_to_rh(QV, 1f1data[:Pa][idx_rslevel, idx_rstime], 273.15 .+ data[:T][idx_rslevel, idx_rstime]);
    end
    RH ./= 100  # to convert from % to "1" units.

    gas_atten = if haskey(data, :gas_atten)
        let tmp = data[:gas_atten][idx_rslevel, idx_rstime]
            reshape(tmp, Nlevel, Ntime, Nfreq)
        end
    else
        fill(0f0, Nlevel, Ntime, Nfreq)
    end

    liq_atten = if haskey(data, :liq_atten)
        let tmp=data[:liq_atten][idx_rslevel, idx_rstime]
            reshape(tmp, Nlevel, Ntime, Nfreq)
        end
    else
        fill(0f0, Nlevel, Ntime, Nfreq)
    end

    K2 = fill(0.91, Nlevel, Ntime, Nfreq); #cat(, fill(0.91, Nlevel, Ntime), dims=3);

    # Aux variables:
    #file_time = @. Second(mwr[:time] - DateTime(2001,1,1,0,0,0));
    file_time = CloudnetTools.datetime24hours(data[:time][idx_rstime]) .|> round
    
    # Creating output file for CloudNetpy
    arm_year = year(data[:time][1])
    arm_month = month(data[:time][1])
    arm_day = day(data[:time][1])

    SITE = CloudnetTools.get_SITE(data, extra_params, inkeys=(:site, :campaign))

    # finding DOI from data file or as extra parameter?
    doi = CloudnetTools.get_parameter(data, :doi, extra_params, default="none") ;

    # title from the input data file?
    arm_title = CloudnetTools.get_parameter(data, :title, extra_params, default="")

    # generating UUID for the file:
    file_uuid = string(uuid1());

    # generating file history:
    file_history = CloudnetTools.get_parameter(data, :history, extra_params,
                                               default="Created by Julia Lang (CloudnetTools.jl)"*string(", ", today() ));
                  

    ARM_OUTFILE = @sprintf("%04d%02d%02d_%s_ecmwf.nc", arm_year, arm_month, arm_day, SITE);

    outputfile = joinpath(output_path, ARM_OUTFILE);
    ds = NCDataset(outputfile, "c", attrib = OrderedDict(
        "Conventions"               => "CF-1.7",
        "location"                  => SITE,
        "source"                    => "Adapted ARM radiosonde data to 'Fake!' ECMWF",
        "institution"               => "DOE ARM",
        "initialization_time"       => "2019-02-05 00:00:00 +00:00",
        "file_uuid"                 => file_uuid,
        "cloudnet_file_type"        => "model",
        "year"                      => arm_year,
        "month"                     => arm_month,
        "day"                       => arm_day,
        "history"                   => file_history,
        "title"                     => "Model file from $SITE",
        #"pid"                       => "https://hdl.handle.net/21.12132/1.913566d962124797",
    ))

    # Dimensions

    ds.dim["time"] = Ntime
    ds.dim["level"] = Nlevel
    ds.dim["flux_level"] = 138
    ds.dim["frequency"] = Nfreq

    # Declare variables

    nclatitude = defVar(ds,"latitude", Float32, (), attrib = OrderedDict(
        "units"                     => "degrees_N",
        "long_name"                 => "Latitude of model gridpoint",
        "standard_name"             => "latitude",
    ))

    nclongitude = defVar(ds,"longitude", Float32, (), attrib = OrderedDict(
        "units"                     => "degrees_E",
        "long_name"                 => "Longitude of model gridpoint",
        "standard_name"             => "longitude",
    ))

    nchorizontal_resolution = defVar(ds,"horizontal_resolution", Float32, (), attrib = OrderedDict(
        "units"                     => "km",
        "long_name"                 => "Horizontal resolution of model",
    ))

    nctime = defVar(ds,"time", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "hours since 2019-02-05 00:00:00 +00:00",
        "long_name"                 => "Hours UTC",
        "standard_name"             => "time",
        "axis"                      => "T",
    ))

    ncforecast_time = defVar(ds,"forecast_time", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "hours",
        "long_name"                 => "Time since initialization of forecast",
        "comments"                  => "For each profile in the file this variable contains the time elapsed since the initialization time of the forecast from which it was taken. Note that the profiles in this file may be taken from more than one forecast.",
    ))

    nclevel = defVar(ds,"level", Int16, ("level",), attrib = OrderedDict(
        "units"                     => "1",
        "long_name"                 => "Model level",
        "standard_name"             => "model_level_number",
        "axis"                      => "Z",
        "positive"                  => "down",
    ))

    ncflux_level = defVar(ds,"flux_level", Int16, ("flux_level",), attrib = OrderedDict(
        "units"                     => "1",
        "long_name"                 => "Model flux level",
        "axis"                      => "Z",
        "positive"                  => "down",
    ))

    ncpressure = defVar(ds,"pressure", Float32, ("level", "time"), attrib = OrderedDict(
        "units"                     => "Pa",
        "long_name"                 => "Pressure",
        "standard_name"             => "air_pressure",
        "C_format"                  => "%.0f",
        "missing_value"             => Float32(-999.0),
    ))

    ncuwind = defVar(ds,"uwind", Float32, ("level", "time"), attrib = OrderedDict(
        "units"                     => "m s-1",
        "long_name"                 => "Zonal wind",
        "standard_name"             => "eastward_wind",
        "C_format"                  => "%.6f",
        "missing_value"             => Float32(-999.0),
    ))

    ncvwind = defVar(ds,"vwind", Float32, ("level", "time"), attrib = OrderedDict(
        "units"                     => "m s-1",
        "long_name"                 => "Meridional wind",
        "standard_name"             => "northward_wind",
        "C_format"                  => "%.6f",
        "missing_value"             => Float32(-999.0),
    ))

    ncomega = defVar(ds,"omega", Float32, ("level", "time"), attrib = OrderedDict(
        "units"                     => "Pa s-1",
        "long_name"                 => "Vertical wind in pressure coordinates",
        "standard_name"             => "omega",
        "C_format"                  => "%.6f",
        "missing_value"             => Float32(-999.0),
    ))

    nctemperature = defVar(ds,"temperature", Float32, ("level", "time"), attrib = OrderedDict(
        "units"                     => "K",
        "long_name"                 => "Temperature",
        "standard_name"             => "air_temperature",
        "C_format"                  => "%.2f",
        "missing_value"             => Float32(-999.0),
    ))

    ncq = defVar(ds,"q", Float32, ("level", "time"), attrib = OrderedDict(
        "units"                     => "1",
        "long_name"                 => "Specific humidity",
        "standard_name"             => "specific_humidity",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncrh = defVar(ds,"rh", Float32, ("level", "time"), attrib = OrderedDict(
        "units"                     => "1",
        "long_name"                 => "Relative humidity",
        "standard_name"             => "relative_humidity",
        "comment"                   => "With respect to liquid above 0 degrees C and with respect to ice below 0 degrees C.",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncql = defVar(ds,"ql", Float32, ("level", "time"), attrib = OrderedDict(
        "units"                     => "1",
        "long_name"                 => "Gridbox-mean liquid water mixing ratio",
        "standard_name"             => "mass_fraction_of_cloud_liquid_water_in_air",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncqi = defVar(ds,"qi", Float32, ("level", "time"), attrib = OrderedDict(
        "units"                     => "1",
        "long_name"                 => "Gridbox-mean ice water mixing ratio",
        "standard_name"             => "mass_fraction_of_cloud_ice_in_air",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    nccloud_fraction = defVar(ds,"cloud_fraction", Float32, ("level", "time"), attrib = OrderedDict(
        "units"                     => "1",
        "long_name"                 => "Cloud fraction",
        "standard_name"             => "cloud_area_fraction",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncflx_net_sw = defVar(ds,"flx_net_sw", Float32, ("flux_level", "time"), attrib = OrderedDict(
        "units"                     => "W m-2",
        "long_name"                 => "Net shortwave flux",
        "standard_name"             => "net_downward_shortwave_flux_in_air",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncflx_net_lw = defVar(ds,"flx_net_lw", Float32, ("flux_level", "time"), attrib = OrderedDict(
        "units"                     => "W m-2",
        "long_name"                 => "Net longwave flux",
        "standard_name"             => "net_downward_longwave_flux_in_air",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncflx_down_sens_heat = defVar(ds,"flx_down_sens_heat", Float32, ("flux_level", "time"), attrib = OrderedDict(
        "units"                     => "W m-2",
        "long_name"                 => "Sensible heat flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncflx_turb_moist = defVar(ds,"flx_turb_moist", Float32, ("flux_level", "time"), attrib = OrderedDict(
        "units"                     => "kg m-2 s-1",
        "long_name"                 => "Turbulent moisture flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncflx_ls_rain = defVar(ds,"flx_ls_rain", Float32, ("flux_level", "time"), attrib = OrderedDict(
        "units"                     => "kg m-2 s-1",
        "long_name"                 => "Large-scale rainfall flux",
        "standard_name"             => "large_scale_rainfall_flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncflx_ls_snow = defVar(ds,"flx_ls_snow", Float32, ("flux_level", "time"), attrib = OrderedDict(
        "units"                     => "kg m-2 s-1",
        "long_name"                 => "Large-scale snowfall flux",
        "standard_name"             => "large_scale_snowfall_flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncflx_conv_rain = defVar(ds,"flx_conv_rain", Float32, ("flux_level", "time"), attrib = OrderedDict(
        "units"                     => "kg m-2 s-1",
        "long_name"                 => "Convective rainfall flux",
        "standard_name"             => "convective_rainfall_flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncflx_conv_snow = defVar(ds,"flx_conv_snow", Float32, ("flux_level", "time"), attrib = OrderedDict(
        "units"                     => "kg m-2 s-1",
        "long_name"                 => "Convective snowfall flux",
        "standard_name"             => "convective_snowfall_flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncflx_turb_mom_u = defVar(ds,"flx_turb_mom_u", Float32, ("flux_level", "time"), attrib = OrderedDict(
        "units"                     => "kg m-1 s-2",
        "long_name"                 => "Zonal turbulent momentum flux",
        "standard_name"             => "downward_eastward_momentum_flux_in_air",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncflx_turb_mom_v = defVar(ds,"flx_turb_mom_v", Float32, ("flux_level", "time"), attrib = OrderedDict(
        "units"                     => "kg m-1 s-2",
        "long_name"                 => "Meridional turbulent momentum",
        "standard_name"             => "downward_northward_momentum_flux_in_air",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_pressure = defVar(ds,"sfc_pressure", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "Pa",
        "long_name"                 => "Surface pressure",
        "standard_name"             => "surface_pressure",
        "C_format"                  => "%.0f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_net_sw = defVar(ds,"sfc_net_sw", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "W m-2",
        "long_name"                 => "Surface net downward shortwave flux",
        "standard_name"             => "surface_net_downward_shortwave_flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_net_lw = defVar(ds,"sfc_net_lw", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "W m-2",
        "long_name"                 => "Surface net downward longwave flux",
        "standard_name"             => "surface_net_downward_longwave_flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_down_sw = defVar(ds,"sfc_down_sw", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "W m-2",
        "long_name"                 => "Surface downwelling shortwave flux",
        "standard_name"             => "surface_downwelling_shortwave_flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_down_lw = defVar(ds,"sfc_down_lw", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "W m-2",
        "long_name"                 => "Surface downwelling longwave flux",
        "standard_name"             => "surface_downwelling_longwave_flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_cs_down_sw = defVar(ds,"sfc_cs_down_sw", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "W m-2",
        "long_name"                 => "Clear sky downwelling shortwave flux",
        "standard_name"             => "surface_downwelling_shortwave_flux_in_air_assuming_clear_sky",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_cs_down_lw = defVar(ds,"sfc_cs_down_lw", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "W m-2",
        "long_name"                 => "Clear sky downwelling longwave flux",
        "standard_name"             => "surface_downwelling_longwave_flux_in_air_assuming_clear_sky",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_down_lat_heat_flx = defVar(ds,"sfc_down_lat_heat_flx", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "W m-2",
        "long_name"                 => "Latent heat flux",
        "standard_name"             => "surface_downward_latent_heat_flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_down_sens_heat_flx = defVar(ds,"sfc_down_sens_heat_flx", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "W m-2",
        "long_name"                 => "Sensible heat flux",
        "standard_name"             => "surface_downward_sensible_heat_flux",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_ls_rain = defVar(ds,"sfc_ls_rain", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "kg m-2",
        "long_name"                 => "Large-scale rainfall amount",
        "standard_name"             => "large_scale_rainfall_amount",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_conv_rain = defVar(ds,"sfc_conv_rain", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "kg m-2",
        "long_name"                 => "Convective rainfall amount",
        "standard_name"             => "convective_rainfall_amount",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_ls_snow = defVar(ds,"sfc_ls_snow", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "kg m-2",
        "long_name"                 => "Large-scale snowfall amount",
        "standard_name"             => "large_scale_snowfall_amount",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_conv_snow = defVar(ds,"sfc_conv_snow", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "kg m-2",
        "long_name"                 => "Convective snowfall amount",
        "standard_name"             => "convective_snowfall_amount",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

ncsfc_ls_precip_fraction = defVar(ds,"sfc_ls_precip_fraction", Float32, ("time",), attrib = OrderedDict(
    "units"                     => "1",
        "long_name"                 => "Large-scale precipitation fraction",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_cloud_fraction = defVar(ds,"sfc_cloud_fraction", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "1",
        "long_name"                 => "Surface total cloud fraction",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_bl_height = defVar(ds,"sfc_bl_height", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "m",
        "long_name"                 => "Boundary layer height",
        "C_format"                  => "%.3f",
        "missing_value"             => Float32(-999.0),
    ))
    
    ncsfc_albedo = defVar(ds,"sfc_albedo", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "1",
        "long_name"                 => "Surface albedo",
        "standard_name"             => "surface_albedo",
        "C_format"                  => "%.6f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_temp_2m = defVar(ds,"sfc_temp_2m", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "K",
        "long_name"                 => "Temperature at 2m",
        "C_format"                  => "%.2f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_q_2m = defVar(ds,"sfc_q_2m", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "1",
        "long_name"                 => "Specific humidity at 2m",
        "C_format"                  => "%.8f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_rough_mom = defVar(ds,"sfc_rough_mom", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "m",
        "long_name"                 => "Surface roughness for momentum",
        "C_format"                  => "%.6f",
        "missing_value"             => Float32(-999.0),
    ))

    ncsfc_rough_heat = defVar(ds,"sfc_rough_heat", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "m",
        "long_name"                 => "Surface roughness for heat",
        "C_format"                  => "%.6f",
        "missing_value"             => Float32(-999.0),
    ))

ncsfc_skin_temp = defVar(ds,"sfc_skin_temp", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "K",
        "long_name"                 => "Skin temperature",
        "C_format"                  => "%.2f",
        "missing_value"             => Float32(-999.0),
    ))

ncsfc_wind_u_10m = defVar(ds,"sfc_wind_u_10m", Float32, ("time",), attrib = OrderedDict(
    "units"                     => "m s-1",
    "long_name"                 => "Zonal wind at 10m",
    "C_format"                  => "%.6f",
    "missing_value"             => Float32(-999.0),
))

ncsfc_wind_v_10m = defVar(ds,"sfc_wind_v_10m", Float32, ("time",), attrib = OrderedDict(
    "units"                     => "m s-1",
    "long_name"                 => "Meridional wind at 10m",
    "C_format"                  => "%.6f",
    "missing_value"             => Float32(-999.0),
))

ncsfc_geopotential = defVar(ds,"sfc_geopotential", Float32, ("time",), attrib = OrderedDict(
    "units"                     => "m2 s-2",
    "long_name"                 => "Geopotential",
    "standard_name"             => "geopotential",
    "C_format"                  => "%.6f",
    "missing_value"             => Float32(-999.0),
))

ncheight = defVar(ds,"height", Float32, ("level", "time"), attrib = OrderedDict(
    "units"                     => "m",
    "long_name"                 => "Height above ground",
    "standard_name"             => "height",
    "comment"                   => "The heights have been calculated using pressure, temperature and specific humidity.",
    "C_format"                  => "%.3f",
    "missing_value"             => Float32(-999.0),
))

ncsfc_height_amsl = defVar(ds,"sfc_height_amsl", Float32, ("time",), attrib = OrderedDict(
    "units"                     => "m",
    "long_name"                 => "Surface height above mean sea level",
    "comment"                   => "The heights have been calculated from the geopotential.",
    "C_format"                  => "%.3f",
    "missing_value"             => Float32(-999.0),
))

ncflx_height = defVar(ds,"flx_height", Float32, ("flux_level", "time"), attrib = OrderedDict(
    "units"                     => "m",
    "long_name"                 => "Height above ground",
    "C_format"                  => "%.3f",
    "missing_value"             => Float32(-999.0),
))

ncwwind = defVar(ds,"wwind", Float32, ("level", "time"), attrib = OrderedDict(
    "units"                     => "m s-1",
    "long_name"                 => "Vertical wind",
    "standard_name"             => "upward_wind",
    "comment"                   => "The vertical wind has been calculated from omega (Pa s-1), height and pressure using: w=omega*dz/dp",
    "C_format"                  => "%.6f",
    "missing_value"             => Float32(-999.0),
))

ncfrequency = defVar(ds,"frequency", Float32, ("frequency",), attrib = OrderedDict(
    "units"                     => "GHz",
    "long_name"                 => "Microwave frequency",
    "C_format"                  => "%.6f",
    "missing_value"             => Float32(-999.0),
))

ncgas_atten = defVar(ds,"gas_atten", Float32, ("level", "time", "frequency"), attrib = OrderedDict(
    "units"                     => "dB",
    "long_name"                 => "Two-way attenuation from the ground due to atmospheric gases",
    "C_format"                  => "%.8f",
    "missing_value"             => Float32(-999.0),
))

ncspecific_gas_atten = defVar(ds,"specific_gas_atten", Float32, ("level", "time", "frequency"), attrib = OrderedDict(
    "units"                     => "dB km-1",
    "long_name"                 => "Specific one-way attenuation due to atmospheric gases",
    "C_format"                  => "%.8f",
    "missing_value"             => Float32(-999.0),
))

ncspecific_saturated_gas_atten = defVar(ds,"specific_saturated_gas_atten", Float32, ("level", "time", "frequency"), attrib = OrderedDict(
    "units"                     => "dB km-1",
    "long_name"                 => "Specific one-way attenuation due to atmospheric gases for saturated air (saturated with respect to ice below 0 degrees C)",
    "C_format"                  => "%.8f",
    "missing_value"             => Float32(-999.0),
))

ncspecific_dry_gas_atten = defVar(ds,"specific_dry_gas_atten", Float32, ("level", "time", "frequency"), attrib = OrderedDict(
    "units"                     => "dB km-1",
    "long_name"                 => "Specific one-way attenuation due to atmospheric gases for dry air (no water vapour)",
    "C_format"                  => "%.8f",
    "missing_value"             => Float32(-999.0),
))

ncK2 = defVar(ds,"K2", Float32, ("level", "time", "frequency"), attrib = OrderedDict(
    "units"                     => "dB km-1",
    "long_name"                 => "Dielectric parameter (|K|^2) of liquid water",
    "C_format"                  => "%.8f",
    "missing_value"             => Float32(-999.0),
))

ncspecific_liquid_atten = defVar(ds,"specific_liquid_atten", Float32, ("level", "time", "frequency"), attrib = OrderedDict(
    "units"                     => "(dB km-1)/(g m-3)",
    "long_name"                 => "Specific one-way attenuation due to liquid water, per unit liquid water content",
    "C_format"                  => "%.8f",
    "missing_value"             => Float32(-999.0),
))


# Define variables

nclatitude[:] = LAT #...**
 nclongitude[:] = LON #...**
 nchorizontal_resolution[:] = 100 #...
 nctime[:] = file_time #...**[idx_rstime]
 ncforecast_time[:] = file_time #...[idx_rstime]
 nclevel[:] = LEVELS #...**
# ncflux_level[:] = ...
 ncpressure[:] = Pa #...**
 ncuwind[:] = U #...**
 ncvwind[:] = V #...**
# ncomega[:] = ...
 nctemperature[:] = TK #...**
 ncq[:] = QV #...**
 ncrh[:] = RH #...**
# ncql[:] = ...**
# ncqi[:] = ...
# nccloud_fraction[:] = ...
# ncflx_net_sw[:] = ...
# ncflx_net_lw[:] = ...
# ncflx_down_sens_heat[:] = ...
# ncflx_turb_moist[:] = ...
# ncflx_ls_rain[:] = ...
# ncflx_ls_snow[:] = ...
# ncflx_conv_rain[:] = ...
# ncflx_conv_snow[:] = ...
# ncflx_turb_mom_u[:] = ...
# ncflx_turb_mom_v[:] = ...
# ncsfc_pressure[:] = ...
# ncsfc_net_sw[:] = ...
# ncsfc_net_lw[:] = ...
# ncsfc_down_sw[:] = ...
# ncsfc_down_lw[:] = ...
# ncsfc_cs_down_sw[:] = ...
# ncsfc_cs_down_lw[:] = ...
# ncsfc_down_lat_heat_flx[:] = ...
# ncsfc_down_sens_heat_flx[:] = ...
# ncsfc_ls_rain[:] = ...
# ncsfc_conv_rain[:] = ...
# ncsfc_ls_snow[:] = ...
# ncsfc_conv_snow[:] = ...
# ncsfc_ls_precip_fraction[:] = ...
# ncsfc_cloud_fraction[:] = ...
# ncsfc_bl_height[:] = ...
# ncsfc_albedo[:] = ...
# ncsfc_temp_2m[:] = ...
# ncsfc_q_2m[:] = ...
# ncsfc_rough_mom[:] = ...
# ncsfc_rough_heat[:] = ...
# ncsfc_skin_temp[:] = ...
# ncsfc_wind_u_10m[:] = ...
# ncsfc_wind_v_10m[:] = ...
# ncsfc_geopotential[:] = ...
 ncheight[:] = HEIGHT .- fill(ALT, 1, Ntime) #...**
 ncsfc_height_amsl[:] = ALT #...
# ncflx_height[:] = ...
# ncwwind[:] = ...
 ncfrequency[:] = 35.5f0 #...**
 ncgas_atten[:] = gas_atten #...**
 ncspecific_gas_atten[:] = 0.001f0 #...**
 ncspecific_saturated_gas_atten[:] = 0.0001f0 #...**
 ncspecific_dry_gas_atten[:] = 0.0001f0 #...**
 ncK2[:] = 0f0 #K2 #...**
 ncspecific_liquid_atten[:] = 0f0 #liq_atten #...**

close(ds)

end # end of function
# ---/
