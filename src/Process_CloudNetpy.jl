#!/opt/julia-1.5.0/bin/julia

# Script to process ARM data with CloudNetpy
using Printf

DATA_PATH = "/home/psgarfias/LIM/remsens/CloudNet/input/"; #"/home/psgarfias/LIM/data/CloudNet/input/";
CATE_PATH = "/home/psgarfias/LIM/data/CloudNet/output/";
PLOT_PATH = "/home/psgarfias/LIM/data/CloudNet/plots/";
years = (2019)
months = (3)
days = (23)
LIDAR_TYPE = "CEIL10m"; # or"HSRL"; # 
USE_MODEL = "ECMWF";

for yy in years
    for mm in months
        for dd in days


radar_file = joinpath(DATA_PATH, @sprintf("KAZR/ARSCL/%04d/arm_radar_%04d%02d%02d.nc", yy, yy, mm, dd));

lidar_file = joinpath(DATA_PATH, @sprintf("%s/%04d/arm_lidar_%04d%02d%02d.nc", LIDAR_TYPE, yy, yy, mm, dd));

model_file = joinpath(DATA_PATH, @sprintf("%s/%04d/%04d%02d%02d_arm-nsa_ecmwf.nc", USE_MODEL, yy, yy, mm, dd));

mwr_file = joinpath(DATA_PATH, @sprintf("MWR/RET/%04d/arm_mwr_%04d%02d%02d.nc", yy, yy, mm, dd));

# Cheking availability for all files
aval_flag = false;
aval_flag  = isfile(radar_file);
aval_flag &= isfile(lidar_file);
aval_flag &= isfile(mwr_file);

aval_flag ? println("Working on $dd.$mm.$yy") : continue

# file names for CloudNet output:
OUTPUT_PATH = joinpath(CATE_PATH, @sprintf("%s/%04d/%02d", LIDAR_TYPE, yy, mm))
isdir(OUTPUT_PATH) ? println(OUTPUT_PATH) : mkdir(OUTPUT_PATH)

categorize_file = joinpath(OUTPUT_PATH, @sprintf("nsa_%s_categorize_%04d%02d%02d.nc", USE_MODEL, yy, mm, dd))

classific_file = joinpath(OUTPUT_PATH, @sprintf("nsa_%s_classific_%04d%02d%02d.nc", USE_MODEL, yy, mm, dd))

iwc_file = joinpath(OUTPUT_PATH, @sprintf("nsa_%s_iwc_%04d%02d%02d.nc", USE_MODEL, yy, mm, dd))

lwc_file = joinpath(OUTPUT_PATH, @sprintf("nsa_%s_lwc_%04d%02d%02d.nc", USE_MODEL, yy, mm, dd))

drizzle_file = joinpath(OUTPUT_PATH, @sprintf("nsa_%s_drizzle_%04d%02d%02d.nc", USE_MODEL, yy, mm, dd))

# plots name:
OUTPUT_PATH = joinpath(PLOT_PATH, @sprintf("%s/%04d/%02d", LIDAR_TYPE, yy, mm))
isdir(OUTPUT_PATH) ? println(OUTPUT_PATH) : mkdir(OUTPUT_PATH)

zbetalwp_plotf = joinpath(OUTPUT_PATH, @sprintf("nsa_%s_ZbetaLWP_%04d%02d%02d", USE_MODEL, yy, mm, dd))

classific_plotf = joinpath(OUTPUT_PATH, @sprintf("nsa_%s_classific_%04d%02d%02d", USE_MODEL, yy, mm, dd))

iwc_plotf = joinpath(OUTPUT_PATH, @sprintf("nsa_%s_iwc_%04d%02d%02d", USE_MODEL, yy, mm, dd))

lwc_plotf = joinpath(OUTPUT_PATH, @sprintf("nsa_%s_lwc_%04d%02d%02d", USE_MODEL, yy, mm, dd))

drizzle_plotf = joinpath(OUTPUT_PATH, @sprintf("nsa_%s_drizzle_%04d%02d%02d", USE_MODEL, yy, mm, dd))

data_file = Dict(
    "radar" => radar_file,
    "lidar" => lidar_file,
    "model" => model_file,
    "mwr"   => mwr_file
);

### Starting the Py Block:
#ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python")
using PyCall
pyimport("sys").executable

# **** Template to user categorize directly without the storage of temporal nc files:
# utils = pyimport("cloudnetpy.utils")
# catcnet = pyimport("cloudnetpy.categorize")
# radar = catcnet.radar.Radar(radar_template_file)
# radar.data["Z"].data = Ze  or set!(radar.data, "Z", Ze)
# radar.data["v"].data = v, etc.
# lidar = catcnet.lidar.Lidar(lidar_template_file)
# 
# time, height = utils.time_grid(), radar.height
# wl=utils.get_wl_band(radar.radar_frequency)

catcnet = pyimport("cloudnetpy.categorize")
procnet = pyimport("cloudnetpy.products")
try
    uuid = catcnet.generate_categorize(data_file, categorize_file)

    
    uuid = procnet.generate_classification(categorize_file, classific_file)
    uuid = procnet.generate_iwc(categorize_file, iwc_file)
    uuid = procnet.generate_lwc(categorize_file, lwc_file)
    uuid = procnet.generate_drizzle(categorize_file, drizzle_file)

catch e
    println("\e[1m\e[38;2;230;30;30;249m","*** Error trying to caracterize data from $dd.$mm.$yy");
    continue
end 
    #
# For Visualization:
MaxHeight = 5;
plotcnet = pyimport("cloudnetpy.plotting")
plotcnet.generate_figure(classific_file, ["target_classification", "detection_status"], show=false, max_y=MaxHeight, image_name=classific_plotf)

plotcnet.generate_figure(lwc_file, ["lwc", "lwc_error", "lwc_retrieval_status"], show=false, max_y=MaxHeight, image_name=lwc_plotf)
plotcnet.generate_figure(iwc_file, ["iwc", "iwc_error", "iwc_retrieval_status"], show=false, max_y=MaxHeight, image_name=iwc_plotf)
plotcnet.generate_figure(drizzle_file, ["Do", "mu", "S"], max_y=MaxHeight, show=false, image_name=drizzle_plotf)

        end
    end
end

#from cloudnetpy.plotting import generate_figure

#py"""
#input_files = {
#    'radar': $radar_file,
#    'lidar': $lidar_file,
#    'model': $model_file,
#    'mwr': $mwr_file
#}
#
#from cloudnetpy.categorize import generate_categorize
#uuid = generate_categorize(input_files, $categorize_file)
##
#from cloudnetpy.products import generate_classification
#uuid = generate_classification($categorize_file, $classific_file)
##
#from cloudnetpy.products import generate_iwc
#generate_iwc($categorize_file, $iwc_file)
##
#from cloudnetpy.products import generate_lwc
#generate_lwc($categorize_file, $lwc_file)
##
#from cloudnetpy.products import generate_drizzle
#generate_drizzle($categorize_file, $drizzle_file)
#
#### Visualization:
#
#from cloudnetpy.plotting import generate_figure
#generate_figure($categorize_file, ['Z','beta','lwp'], show=False, max_y=7, image_name=$zbetalwp_plotf)
#
#generate_figure($classific_file, ['target_classification', 'detection_status'], show=False, max_y=7, image_name=$classific_plotf)
#
#generate_figure($iwc_file, ['iwc', 'iwc_error', 'iwc_retrieval_status'], show=False, max_y=5, image_name=$iwc_plotf)
#
#generate_figure($lwc_file, ['lwc', 'lwc_error', 'lwc_retrieval_status'], show=False, max_y=5, image_name=$lwc_plotf)
#
#generate_figure($drizzle_file, ['Do', 'mu', 'S'], max_y=3, show=False, image_name=$drizzle_plotf)
#
#"""
#
## end of Script
