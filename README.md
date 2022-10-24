# CloudnetTools.jl
A set of Julia routines and functions to work with Cloudnet classification algorithm and its data files.

This package is mean to interface with the Cloudnetpy version developed by [Finnish Meteorological Institute](https://fmi.fi) under the ACTRIS program. There is a minor support to work with the Cloudnet output files from the legacy code (Matlab-based).

## Installation
It is assumed that you have Cloudnetpy installed in your machine. To install Cloudnetpy follow the indications by [FMI Cloudnet installation](https://cloudnetpy.readthedocs.io/en/latest/installation.html).

### Using Cloudnet in Julia
The following packages are required to use Cloudnet from Julia:
* PyCall.jl
* ARMtools.jl

and the following 3 steps are need to install CloudnetTools.jl :

1. Install PyCall.jl from the Julia package manager and make sure that Juila's ENV variable PYTHON is indicating to the python environment where Cloudnetpy is installed and operational: 
```julia
(ac3-b07) ~/LIM/cloudnetpy$ julia
julia> ]
(@v1.8) pkg> add PyCall
julia> ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python")
julia> ENV["PYTHON"]="/usr/local/bin/python3.8"
"/usr/local/bin/python3.8"
julia> Pkg.build("PyCall")
```

2. Install ARMtools.jl which is a collection of tools to work with ARM data files. This package can be cloned to your local machine or installed from Github:
```julia
julia> ]
(@v1.8) pkg> add git@github.com:pablosaa/ARMtools.jl
julia> using ARMtools
```

3. Finally intall CloudnetTools.jl same way as with ARMtools.jl :
```julia
julia> ]
(@v1.8) pkg> add git@github.com:pablosaa/CloudnetTools.jl
julia> using CloudnetTools
``` 

# Usage
To use Cloudnet with ARM data you have to convert first the ARM data files into Cloudnet supported input files, those are data files for cloud radar, lidar, and microwave radiometer.
## Convert ARM data to Cloudnet compatible input files
### Cloud Radar:
```julia
julia> using CloudnetTools.ARM
julia> output_path = "/data/cloudnet/input"   # or whereever you want to store the files
julia> kazr_file = "/data/KAZR/ARSCL/2019/nsaarsclkazr1kolliasC1.c0.20190119.000000.nc"
julia> ARM.converter(:radar, kazr_file, output_path)
```
the last command will create a cloudnet compatible file ``` /data/cloudnet/input/2019019_nsa_radar.nc```.
Supported ARM cloud radars are KAZR GE (KAZR MD maybe, not tested), KAZR ARSCL, and MWACR.

### Lidar:
```julia
julia> ceilometer_file = "/data/CEIL10m/2019/nsaceil10mC1.b1.20190119.000008.nc"
julia> ARM.converter(:ceilometer, ceilometer_file, output_path)
```
or
```julia
julia> hsrl_file = "/data/HSRL/2019/nsahsrlC1.a1.20190119.000000.nc"
julia> ARM.converter(:lidar, hsrl_file, output_path)
```
to create the cloudnet compatible file ``` /data/cloudnet/input/2019019_nsa_lidar.nc```.
Supported ARM lidars are CEILOMETER 10m, and HSRL.

### Microwave Radiometer:
```julia
julia> mwr_file = "/data/CEIL10m/2019/nsaceil10mC1.b1.20190119.000008.nc"
julia> ARM.converter(:mwr, mwr_file, output_path)
```
to create the cloudnet compatible file ``` /data/cloudnet/input/2019019_nsa_mwr.nc```.
Supported ARM radiometer data is MWR/RET.

## Run Cloudnet
Once the ARM data for cloud radar, lidar and microwave radiometer are converted, Clounet required in addition numerical weather data for the specific site. FMI provides in its public repository those model data ready to use for cloudnet by using the [API model data](https://docs.cloudnet.fmi.fi/api/data-portal.html#get-apimodels--model). There you can get the model file like ```20190119-arm_nsa_ecmwf.nc``` which is need to work with the examples above.

The two basic output files from Cloudnet are the categorization and classification products. You need to create the first in order to obtain the second, follwing the example:
```julia
julia> using CloudnetTools.ACTRIS
julia> data_files = Dict(
        :radar  => "/data/cloudnet/input/20190119_nsa_radar.nc",
        :lidar  => "/data/cloudnet/input/20190119_nsa_lidar.nc",
        :mwr    => "/data/cloudnet/input/20190119_nsa_mwr.nc",
        :model  => /data/cloudnet/input/20190119_arm-nsa_ecmwf.nc",
        )
julia> ACTRIS.categorize_it(data_files, "/data/cloudnet/output/20190119_nsa_categorization.nc")
```
the last command will create the categorization file ```20190119_nsa_categorization.nc``` which Cloudnet needs to create the classification output file as following:
```julia
julia> categorize_files = "/data/cloudnet/output/20190119_nsa_categorization.nc"
julia> CLNTprod = Dict(
        :categorize => true,
        :classification => true,
        :lwc => true,
        :iwc => true,
        :drizzle => false,
        :der => true,
        :ier => true,
      )
julia> ACTRIS.generate_product(categorize_files, CLNTprod)
```
The last call will generate Cloudnet products and retrievals indicated as ```true``` in the dictionary of products ```CLNTprod```. The output files then are called after the product like ```20190119_nsa_classification.nc```, ```20190119_nsa_lwc.nc``` for cloud liquid water content, ```20190119_nsa_iwc.nc``` for ice water content, etc.

## Visualization
The Cloudnet output files can be visualized, for example for categorization, as follow:
```julia
julia> using CloudnetTools
julia> clnet = CloudnetTools.readCLNFile("/data/cloudnet/output/20190119_nsa_categorization.nc");
julia> CloudnetTools.Vis.show_measurements(clnet; SITENAME="NSA", maxhgt=7.5, savefig="20190119_measurements_cloudnet.png")
```
<img src="https://user-images.githubusercontent.com/24556480/197576435-cc6ee6b8-f1a2-4670-ad9a-e76e2ccf2220.png" width="700" />

the only input argument needed is the dictionary containing all Cloudnet information ```clnet```, the other arguments are optional which indicates ```SITENAME::String``` with the name of the ARM site, ```maxhgt::Real``` maximum altitude to show in plot in km, and ```savefig::String``` the name of a file for the plot to be stored.

```julia
julia> CloudnetTools.Vis.show_classific(clnet; SITENAME="NSA", maxhgt=7.5, savefig="20190119_classification_cloudnet.png")
```
<img src="https://user-images.githubusercontent.com/24556480/197576487-cd4acde4-9a6c-4785-a1c1-ce7085ce51aa.png" width="700" />

The wind vectors and iso-therm lines can be excluded from the plot if the following optional argument is pass: ```showatm=Dict(:wind=>false, :isoT=>true, :procas=>false)```, where the wind vectors are skipped and only the iso-therms will be displayes.

# Further estimations
Besides the standard Cloudnet classification and retrievals, ```CloudnetTools.jl``` enables to estimate the cloud base and top height for multiple cloud layers
(up to three layers for now) and their corresponding cloud base and top temperatures for the same layers.

# Disclamer
This package is independent to the development of ACTRIS' Cloudnetpy, and it is recomended to constantly update Cloudnetpy since this in not done by the Julia package. Moreover this package has been developed for application related to measuremnts in the Arctic, mainly for the ARM site in Utquiaǵvik and the MOSAiC expeditions. Thus the ARM data available for those sites have been tested and adapted for Cloudnet, it could happen that ARM data for other sites are different and that might produce errors, in that case we recomend to open an issue and we will correct or adapt to new data types.

### Acknowledges
Thank to FMI and ACTRIS for shareing the open source version of Cloudnet and providing the auxiliary model data.
The work in this repository has been funded by the _Deutsche Forschungsgesellschaft_ (DFG) project TR-172 [(AC)3 sub-project B07](https://ac3-tr.de/projects/cluster-b/b07).

--<br>
(c) 2021, Pablo Saavedra Garfias<br>
pablo.saavedra@uni-leipzig.de<br>
University of Leipzig<br>
Faculty of Physics and Geosciences<br>
LIM<br>

See LICENSE
