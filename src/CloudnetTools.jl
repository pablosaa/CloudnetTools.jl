# Main CloudNet.jl module.
"""
A set of tools to process and analyze data outputs from CloudNet classification algorithm.

(c) 2020, Pablo Saavedra Garfias
University of Leipzig
Faculty of Physics and Geosciences
LIM

See LICENSE
"""


module CloudnetTools

using NCDatasets, DataStructures
using Dates
using Interpolations
using Printf

include("reading_functions.jl")

# ------------------------------------------------------
# AUX FUNCTIONS
# ------------------------------------------------------
include("aux_functions.jl")

# ****************************************************
# Including further modules:
# * For plotting CloudNet products:
Base.include(CloudnetTools, "CloudnetVisualization.jl")

# * For dealing with ARM data for cloudnet:
Base.include(CloudnetTools, "Cloudnet_for_ARM.jl")

# * For interfacing with cloudnetpy from FMI's ACTRIS:
Base.include(CloudnetTools, "Cloudnet_ACTRIS.jl")

end  # end of Module CloudnetTools


# --end of script
