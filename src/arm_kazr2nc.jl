function kazr2nc(radar_file::String, output_path::String)

# Script to adapt ARM KAZR radar to netCDF input for cloudnetpy
# 

# Reading input ARM netCDF file:

if !isfile(radar_file)
    error("$radar_file does not exist!")
end 
ncin = NCDataset(radar_file, "r");
time = ncin["time"];
arm_year = year(time[1])
arm_month = month(time[1])
arm_day = day(time[1])

range = ncin["height"];  # [m]
Ze = ncin["reflectivity"];  # [dBz] reflectivity best mode
mdv = ncin["mean_doppler_velocity"]; # [m/s] Mean Doppler velocity
spw = ncin["spectral_width"];  # [m/s] Spectral Width
LDR = ncin["linear_depolarization_ratio"]; # [dB]
SNR = ncin["signal_to_noise_ratio"];  # [dB]
latitude = ncin["lat"];
longitude = ncin["lon"];
altitude = ncin["alt"];
radar_frequency = parse(Float32, replace(ncin.attrib["radar_operating_frequency_chirp"], "GHz"=>""));  # [GHz]
nzrg = length(range);  # Number of range gates
nyquist_velocity = -999;  # [m/s] (to be confirmed) 
num_nfft = -999;  # [counts] number of fast fourier transform (to be confirmed)
num_prf  = -999;  # [Hz] pulse repetition frequency (to be confirmed)
num_nave = 9999.;  # [counts] number of spectral average ((10-0)/0.001 ? to be confirmed)
location_radar = ncin.attrib["location_description"];
radar_source = ncin.attrib["process_version"];
file_uuid = string(uuid1());
file_history = string(now(), " - radar file created");
file_time = hour.(time) + minute.(time)/60 + second.(time)/3600;

# Creating output file for CloudNetpy
ARM_OUTPATH = output_path;
ARM_OUTFILE = @sprintf("arm_radar_%04d%02d%02d.nc", arm_year, arm_month, arm_day);
outputfile = joinpath(ARM_OUTPATH, ARM_OUTFILE);
ds = NCDataset(outputfile, "c", attrib = OrderedDict(
    "Conventions"               => "CF-1.7",
    "cloudnetpy_version"        => "1.3.2",
    "file_uuid"                 => file_uuid, #"58bccd4053f74be385d527a4cb335d9b",
    "cloudnet_file_type"        => "radar",
    "title"                     => "ARSCL Radar made by Julia",
    "year"                      => Int16(arm_year),
    "month"                     => Int16(arm_month),
    "day"                       => Int16(arm_day),
    "location"                  => location_radar,
    "history"                   => file_history, #"2020-09-29 13:40:35 - radar file created",
    "source"                    => radar_source, #"METEK MIRA-35",
    "references"                => "https://doi.org/10.21105/joss.02123",
))

# Dimensions

ds.dim["time"] = Inf; # unlimitted dimension e.g. 8432
ds.dim["range"] = length(range); # 498

# Declare variables

ncZe = defVar(ds,"Ze", Float32, ("range", "time"), attrib = OrderedDict(
    "units"                     => "dBZ",
    "long_name"                 => "Radar reflectivity factor (uncorrected), vertical polarization",
    "missing_value"             => Ze.attrib["missing_value"],
    "attenuation_corr_gas"      => true,
    "attenuation_corr_liq"      => false,
))

ncv = defVar(ds,"v", Float32, ("range", "time"), attrib = OrderedDict(
    "units"                     => "m s-1",
    "long_name"                 => "Doppler velocity",
    "comment"                   => "This parameter is the radial component of the velocity, with positive
velocities are away from the radar.",
    "ancillary_variables"       => "v_sigma",
    "positive"                  => "up",
    "missing_value"             => mdv.attrib["missing_value"],
))

ncwidth = defVar(ds,"width", Float32, ("range", "time"), attrib = OrderedDict(
    "units"                     => "m s-1",
    "long_name"                 => "Spectral width",
    "comment"                   => "This parameter is the standard deviation of the reflectivity-weighted
velocities in the radar pulse volume.",
    "missing_value"             => spw.attrib["missing_value"],
))

ncldr = defVar(ds,"ldr", Float32, ("range", "time"), attrib = OrderedDict(
    "units"                     => "dB",
    "long_name"                 => "Linear depolarisation ratio",
    "comment"                   => "This parameter is the ratio of cross-polar to co-polar reflectivity.",
    "missing_value"             => LDR.attrib["missing_value"],
))

ncSNR = defVar(ds,"SNR", Float32, ("range", "time"), attrib = OrderedDict(
    "units"                     => "dB",
    "long_name"                 => "Signal-to-noise ratio",
    "missing_value"             => SNR.attrib["missing_value"],
))

nclatitude = defVar(ds,"latitude", Float32, (), attrib = OrderedDict(
    "units"                     => "degrees_north",
    "long_name"                 => "Latitude of site",
))

nclongitude = defVar(ds,"longitude", Float32, (), attrib = OrderedDict(
    "units"                     => "degrees_east",
    "long_name"                 => "Longitude of site",
))

ncaltitude = defVar(ds,"altitude", Float32, (), attrib = OrderedDict(
    "units"                     => "m",
    "long_name"                 => "Altitude of site",
))

nctime = defVar(ds,"time", Float32, ("time",), attrib = OrderedDict(
    "units"                     => "decimal hours since midnight",
    "long_name"                 => "Time UTC",
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
    "units"                     => "count",
))

ncprf = defVar(ds,"prf", Int32, (), attrib = OrderedDict(
    "long_name"                 => "Pulse Repetition Frequency",
    "units"                     => "Hz",
))

ncnave = defVar(ds,"nave", Int32, (), attrib = OrderedDict(
    "long_name"                 => "Number of Spectral Avreages",
    "units"                     => "count",
))

nczrg = defVar(ds,"zrg", Int32, (), attrib = OrderedDict(
    "long_name"                 => "Number of Range Gates",
    "units"                     => "count",
))

ncrg0 = defVar(ds,"rg0", Int32, (), attrib = OrderedDict(
    "long_name"                 => "Number of Lowest Range Gates",
    "units"                     => "count",
))

ncdrg = defVar(ds,"drg", Float32, (), attrib = OrderedDict(
    "long_name"                 => "Range Resolution",
    "unit"                      => "m",
))


# Define variables

ncZe[:] = Ze;
ncv[:] = mdv;
ncwidth[:] = spw;
ncldr[:] = LDR;
ncSNR[:] = SNR;
nclatitude[:] = latitude[1];
nclongitude[:] = longitude[1];
ncaltitude[:] = altitude[1];
nctime[:] = file_time;
ncrange[:] = range;
ncradar_frequency[:] = radar_frequency;
ncnyquist_velocity[:] = nyquist_velocity;
ncnfft[:] = num_nfft;
ncprf[:] = num_prf;
ncnave[:] = num_nave;
nczrg[:] = nzrg;
ncrg0[:] = 5;
ncdrg[:] = 30;

close(ds)
close(ncin)

end
# end of function

#arm_year = 2017;
#arm_month = 11;
#arm_day = 9;
#DATA_PATH = "/home/psgarfias/LIM/data/utqiagvik-nsa/KARZ/ARSCL/";
#DATA_NAME = @sprintf("nsaarsclkazr1kolliasC1.c0.%04d%02d%02d.000000.nc", #arm_year, arm_month, arm_day);
#radar_file = joinpath(DATA_PATH, DATA_NAME);
