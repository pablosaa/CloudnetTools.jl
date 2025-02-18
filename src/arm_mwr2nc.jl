"""
Function to convert ARM microwave radiometer file to Cloudnet MWR input file.
USAGE:
```julia-repl
julia> mwr2nc(mwr_file, output_path)
julia> mwr2nc(mwr_file, output_path; extra_params=extras)
julia> mwr2nc(list_files, output_path)
julia> mwr2nc(data, output_path)
```
WHERE:
* ```mwr_file::String``` full path to ARM MWR/LOS or MWR/RET netCDF file,
* ```output_path::String``` path to put the converted file,
* ```list_files::Vector{String}``` several ARM files to be concatenated,
* ```data::Dict``` dataset readed by ARMtools.getMWRData(armfile)
* ```extra_params::Dict``` (Optional) dictionary with alternative parameter to pass.

The ```extra_params::Dict``` could be for example values not present in the netCDF file or values present but need to be adjusted like Latitude and Longitude for measurements in a Research Vessel, e.g. ```Dict(:SITE=>"Polarstern", :altitude_m=>10f0, :wet=>rain_flag, :lat=>78.1, :lon=>10.5)```, etc.

Part of ```CloudnetTools.jl```, see LICENSE.TXT
"""
function mwr2nc(mwr_file::String, output_path::String; extra_params=Dict{Symbol, Any}())
    mwr = ARMtools.getMWRData(mwr_file)
    
    return mwr2nc(mwr::Dict, output_path; extra_params=Dict(:source=>mwr_file, extra_params...) )
end
function mwr2nc(mwr_file::Vector{String}, output_path::String; extra_params=Dict{Symbol, Any}())
    mwr = ARMtools.getMWRData(mwr_file)
    
    return mwr2nc(mwr::Dict, output_path; extra_params=Dict(:source=>mwr_file, extra_params...) )
end
function mwr2nc(mwr::Dict, output_path::String; extra_params=Dict{Symbol, Any}())

    file_time = datetime24hours(mwr[:time])
    ntime = length(file_time)
    
    # Creating output file for CloudNetpy
    arm_year = year(mwr[:time][1])
    arm_month = month(mwr[:time][1])
    arm_day = day(mwr[:time][1])

    SITE = get_SITE(mwr, extra_params, inkeys=(:site, :campaign))

    # Checking for LWP conversion to older cloudnetpy versions (default false):
    old_ver = get_parameter(mwr, :oldversion, extra_params, default=false)

    # Converting LWP to kg m⁻² if old_version is false:
    !old_ver && (mwr[:LWP] .*= 1f-3 )
    
    # finding DOI from data file or as extra parameter?
    doi = get_parameter(mwr, :doi, extra_params, default="none") ;

    # title from the input data file?
    arm_title = get_parameter(mwr, :title, extra_params, default="")

    # generating UUID for the file:
    file_uuid = string(uuid1());

    # generating file history:
    file_history = get_parameter(mwr, :history, extra_params,
                                 default="Created by Julia Lang (CloudnetTools.jl)"*string(", ", today() ));
                  
    # WET flag
    wet_flag = if isa(mwr[:wet], Vector)
        mwr[:wet]
    elseif haskey(extra_params, :wet)
        extra_params[:wet]
    else
        fill(NaN32, ntime)
    end
    
    ARM_OUTFILE = @sprintf("%04d%02d%02d_%s_mwr.nc", arm_year, arm_month, arm_day, SITE);

    outputfile = joinpath(output_path, ARM_OUTFILE);
    ds = NCDataset(outputfile, "c", attrib = OrderedDict(
        "Conventions"               => "CF-1.7",
        "file_uuid"                 => file_uuid,
        "cloudnet_file_type"        => "mwr",
        "title"                     => arm_title,
        "location"                  => SITE,
        "radiometer_system"         => get_parameter(mwr, :radiometer_system, extra_params, default="Radiometrics"),
        "serial_number"             => get_parameter(mwr, :serial_number, extra_params, default="not defined"),
        "year"                      => arm_year,
        "month"                     => arm_month,
        "day"                       => arm_day,
        "source"                    => get_parameter(mwr, :source, extra_params, default="arm.gov"),
        "history"                   => file_history,
        "radiometer_software_version" => get_parameter(mwr, :radiometer_software_version, extra_params, default="none"),
        "host_PC_software_version"  => get_parameter(mwr, :host_PC_software_version, extra_params, default="none"),
        "references"                => doi,
    ))

    # Dimensions

    #ds.dim["time"] = Inf # unlimited dimension
    defDim(ds, "time", ntime)
    
    # Declare variables

    nctime_reference = defVar(ds,"time_reference", Float32, (), attrib = OrderedDict(
        "long_name"                 => "flag indicating the time zone reference",
        "units"                     => "unitless",
        "comment"                   => "0 = local time, 1 = UTC",
    ))

    ncnumber_integrated_samples = defVar(ds,"number_integrated_samples", Float32, ())

    ncminimum = defVar(ds,"minimum", Float32, (), attrib = OrderedDict(
        "units"                     => "g / m^2",
    ))

    ncmaximum = defVar(ds,"maximum", Float32, (), attrib = OrderedDict(
        "units"                     => "g / m^2",
    ))

    nctime = defVar(ds,"time", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "Time UTC",
        "standard_name"             => "time",
        "units"                     => "seconds since 1.1.2001, 00:00:00",
        "comment"                   => "reference time zone indicated in field time_reference",
    ))

    ncrain_flag = defVar(ds,"rain_flag", Float32, ("time",), attrib = OrderedDict(
        "Info"                      => "0 = No Rain, 1 = Raining",
    ))

    ncelevation_angle = defVar(ds,"elevation_angle", Float32, ("time",), attrib = OrderedDict(
        "long_name"                 => "viewing_elevation angle",
        "comment"                   => "-90 is blackbody view, 0 is horizontal view (red arrow), 90 is zenith view, 180 is horizontal view (2nd quadrant)",
        "units"                     => "degrees (-90 - 180)",
    ))

    ncazimuth_angle = defVar(ds,"azimuth_angle", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "DEG (0�-360�)",
        "long_name"                 => "Azimuth angle",
        "standard_name"             => "solar_azimuth_angle",
    ))

    ncretrieval = defVar(ds,"retrieval", Float32, (), attrib = OrderedDict(
        "Info"                      => "0 = Linear Regr., 1 = Quadr. Regr., 2 = Neural Network",
    ))

    ncLWP_data = defVar(ds,"lwp", Float32, ("time",), attrib = OrderedDict(
        "units"                     => ifelse(old_ver, "g m-2", "kg m-2"),
        "long_name"                 => "Liquid water path",
        "standard_name"             => "atmosphere_cloud_liquid_water_content",
    ))

    ncIWV_data = defVar(ds,"iwv", Float32, ("time",), attrib = OrderedDict(
        "units"                     => "kg m-2",
        "long_name"                 => "Integrated water vapour",
        "standard_name"             => "atmosphere_mass_content_of_water_vapour",
    ))

    ncaltitude = defVar(ds,"altitude", Float32, (), attrib = OrderedDict(
        "units"                     => "m",
        "long_name"                 => "Altitude of site",
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
    
    
    # Assigning variable values:

    nctime[1:ntime] .= file_time
    nctime_reference[:] = get_parameter(mwr, :rerefence, extra_params, default=1f0);
    ncnumber_integrated_samples[:] = get_parameter(mwr, :integrated_samples, extra_params, default=1f0);
    ncminimum[:] = minimum(mwr[:LWP])
    ncmaximum[:] = maximum(mwr[:LWP])
    ncrain_flag[:] .= wet_flag;
    ncelevation_angle[:] .= mwr[:elevation];
    ncazimuth_angle[:] .= mwr[:azimuth];
    ncretrieval[:] = get_parameter(mwr, :retrieval, extra_params, default=1);
    ncLWP_data[:] = mwr[:LWP];
    ncIWV_data[:] = mwr[:IWV];
    ncaltitude[:] = get_parameter(mwr, :alt, extra_params, default=0, lims=(0,8000));
    nclatitude[:] = get_parameter(mwr, :lat, extra_params, lims=(-90, 90));
    nclongitude[:] = get_parameter(mwr, :lon, extra_params, lims=(-180, 360));
    
    close(ds)

    return ARM_OUTFILE
end
# ----/
# end of function
