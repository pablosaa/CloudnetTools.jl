function lidar2nc(lidar_file::String, output_path::String)

if !isfile(lidar_file)
    error("$lidar_file does not exist!")
end

# Reading input ARM netCDF file:
ncin = NCDataset(lidar_file, "r");
time = ncin["time"];
arm_year = year(time[1])
arm_month = month(time[1])
arm_day = day(time[1])
range = ncin["range"];
beta = ncin["backscatter"]; # converted from 1/(sr km 10000) to 1/(sr m)
wavelength = 910;  # Obtained from the ARM ceilometer handbook [nm]
tilt_angle = maximum(ncin["tilt_angle"]);
height = ncin["alt"];
lidar_location = ncin.attrib["location_description"];
lidar_source = ncin.attrib["ceilometer_model"];
lidar_doi = ncin.attrib["doi"];
file_uuid = string(uuid1());
file_history = string(now(), " - ceilometer file created");
file_time = hour.(time) + minute.(time)/60 + second.(time)/3600;
# Adaptation from CloudNetpy vaisala.py script:
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
β[isnan.(β)] .= 9.96921f36

# Creating cloudnetpy input netCDF file:
# Creating output file for CloudNetpy
ARM_OUTPATH = output_path;
ARM_OUTFILE = @sprintf("arm_lidar_%04d%02d%02d.nc", arm_year, arm_month, arm_day);
outputfile = joinpath(ARM_OUTPATH, ARM_OUTFILE);
ds = NCDataset(outputfile, "c", attrib = OrderedDict(
    "Conventions"               => "CF-1.7",
    "cloudnetpy_version"        => "1.3.2",
    "file_uuid"                 => file_uuid,
    "cloudnet_file_type"        => "lidar",
    "title"                     => ncin.attrib["location_description"], #"Ceilometer file from Mace-Head",
    "year"                      => Int16(arm_year),
    "month"                     => Int16(arm_month),
    "day"                       => Int16(arm_day),
    "location"                  => lidar_location, #"Mace-Head",
    "history"                   => file_history, #"2020-09-29 14:50:54 - ceilometer file created",
    "source"                    => lidar_source, #"Jenoptik CHM15k",
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
    #"missing_value"             => beta.attrib["missing_value"],
))

ncbeta = defVar(ds,"beta", Float32, ("range", "time"), attrib = OrderedDict(
    "units"                     => "sr-1 m-1",
    "long_name"                 => "Attenuated backscatter coefficient",
    "comment"                   => "Range corrected, SNR screened, attenuated backscatter.",
    #"missing_value"             => beta.attrib["missing_value"],
))

ncbeta_smooth = defVar(ds,"beta_smooth", Float32, ("range", "time"), attrib = OrderedDict(
    "units"                     => "sr-1 m-1",
    "long_name"                 => "Smoothed attenuated backscatter coefficient",
    "comment"                   => "Range corrected, SNR screened backscatter coefficient.\n Weak background is smoothed using Gaussian 2D-kernel.",
    #"missing_value"             => beta.attrib["missing_value"],
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

ncbeta_raw[:] = 1e-7*beta; #...
ncbeta[:] = β; #...
ncbeta_smooth[:] = β; #...
ncrange[:] = range; #...
nctime[:] = file_time; #...
nctilt_angle[:] = tilt_angle; #...
ncheight[:] = height[1] .+ range; #...
ncwavelength[:] = wavelength; #...

close(ds)
close(ncin)

return nothing
end

#using NCDatasets, DataStructures
#using UUIDs
#using Dates
#using Printf, Statistics
#
#
#arm_year = 2017;
#arm_month = 11;
#arm_day = 9;
#DATA_PATH = "/home/psgarfias/LIM/data/utqiagvik-nsa/CEIL10m/";
#DATA_NAME = @sprintf("nsaceil10mC1.b1.%04d%02d%02d.000008.nc", arm_year, #arm_month, arm_day);
#lidar_file = joinpath(DATA_PATH, DATA_NAME);
