# Part of CloudnetTools.jl
# ******************************************************************
# +++++++++++++++ MODULE Vis  for visulaization ++++++++++++++++++++
# ******************************************************************

module Vis

using Plots
using Dates
using Printf

function show_LWP_IWP(iwc::Dict, liq::Dict)
    p1=plot(iwc[:time], iwc[:height], log10.(1f3iwc[:iwc]), st=:heatmap, color=:lapaz);
    p2=plot(liq[:time], liq[:height], log10.(liq[:lwc]), st=:heatmap, color=:tokyo);
    plot(p1,p2, layout=(2,1));
end

# *********************************************************************
# Function to plot the Classification data
"""
# Function to plot the Classification data

> show_classific(cnt::Dict; 

WHERE:
* cnt::Dict dictionary ouput from read_CNTfile(cloudnet_file),
* SITENANE::String (optional) string with name of site,
* maxhgt::Number (optional) indicating the maximum height in km, default=8,
* showlegend::Bool (optional) show Cloudnet legend colors, default=true
* addons::Bool to add meteo data to the plot, "wind" & "temp", default=true
Output:
* plt::Plot output plot object.
"""
function show_classific(cnt::Dict; SITENAME="", maxhgt=8, showlegend=true,
                        showatm=(:wind=true, :isoT=true), savefig=:none)

    # defining time axis ticks:
    tm_tick = cnt[:time][1]:Minute(90):cnt[:time][end];

    Xstrname = "TIME UTC from "*Dates.format(cnt[:time][2], "dd.u.yyyy")

    strtitle = (tm_tick[1], 7.5, text("CloudNet Target Classification "*SITENAME, 11, halign=:left))
    strtitle = "CloudNet Target Classification "*SITENAME
    
    cldnet = CloudNetPalette("classific")
    l = grid(2,1, heights=(0.2,0.8)) #[a{0.25h}; b];
    Y_LIM = (0, maxhgt)
    
    classplt = plot(cnt[:time], 1f-3cnt[:height], cnt[:CLASSIFY],
                    st=:heatmap, colorbar = false, framestyle = :box,
                    color=palette(cldnet, 11), ylim=Y_LIM, clim=(0,10),
                    ylabel="Height A.G.L. [km]", ytickfontsize=11, minorticks=true,
                    xlabel= Xstrname, tick_direction=:out, 
                    xticks=(tm_tick, Dates.format.(tm_tick, "H:MM")), xrot=45,
                    xtickfontsize=11, xguidefont=font(12), title=strtitle);

    # Preparing to add graphs for isotherms and wind vectors:
    ihmax = findlast(1f-3cnt[:model_height] .≤ maxhgt)

    
    Xin = haskey(cnt, :model_time) ? Vector(1:length(cnt[:model_time])) : round.(Int64, range(1, stop=length(cnt[:time]), length=25))
    
    Yin = 1f-3cnt[:model_height]
    #TLEV = [5, 0, -5, -10, -15, -20, -25, -30]
    BB = bbox(0,0,1,1)

    if showatm[:wind]
        Tin = cnt[:T][1:ihmax, Xin] .- 273.15
        Tlevels = extrema(Tin) |> x->collect.(ceil(x[1]):5:floor(x[2]))
        Attach_Isotherms(classplt, Xin, Yin[1:ihmax], Tin,
                         (1, BB), 2, TLEV = Tlevels)
    end

    if showatm[:isoT]
        Attach_Windvector(classplt, Xin[1:2:end], Yin[2:4:ihmax],
                          cnt[:UWIND][2:4:ihmax, Xin[1:2:end]],
                          cnt[:VWIND][2:4:ihmax, Xin[1:2:end]],
                          (1, BB), 3)
    end
    #
    #θv = @. cnt[:T][1:ihmax,2:2:end]*(1024f0/cnt[:Pa][1:ihmax,2:2:end])^0.286;
    
    #Attach_Profile_Cascate(classplt, Xin[2:2:end], Yin[1:ihmax], θv, (1,BB), 4)
    
    if showlegend
        classleg = ShowLegendCloudNetClassification("classific")

        pltout = plot(classleg, classplt, layout=l, size=(800,700))
    else
        pltout = plot(classplt, size=(800,600))
    end

    savefig != :none && typeof(savefig) <: String && Plots.savefig(pltout, savefig)
    
    return pltout
end
# -- OR
function show_classific(cnt_file::String; SITENAME="", maxhgt=8, showlegend=true, savefig=:none)
    cnt = CloudnetTools.readCLNFile(cnt_file)
    return show_classific(cnt, SITENAME=SITENAME, maxhgt=maxhgt, showlegend=showlegend,savefig=savefig)
end
# ----/


function CloudNetPalette(ColorType::String)
    if ColorType == "classific"
        tmp = [
            (1.00,  1.0 ,  1.00);
            (0.44,  1.0 ,  0.92);
            (0.17,  0.62,  0.95);
            (0.75,  0.6 ,  1.0);
            (0.85,  0.85,  0.87); # 0.9, 0.9, 0.92
            (0.28,  0.27,  0.72);
            (0.99,  0.65,  0.0);
            (0.78,  0.98,  0.19);
            (0.8 ,  0.73,  0.53);
            (0.89,  0.29,  0.13);
            (0.7 ,  0.21,  0.34)
        ];
    elseif ColorType == "detection"
        tmp = [
            (0.247,  0.704,  0.43);
            (0.44 ,  0.926,  0.34); 
            (0.996,  0.91 ,  0.24);
            (0.8  ,  0.96 ,  0.96); 
            (0.455,  0.51 ,  0.41); 
            (0.89 ,  0.292,  0.13);
        ];
    else
        error("$ColorType not supported!")
    end
    cldnet = map(x->RGB(x...), tmp);
    return cldnet
end

function ShowLegendCloudNetClassification(LegendType::String)
    if LegendType == "classific"
        txt_labels = [
            #X, Y, "short_name", "long_name"
            (40, 1,"Clear sky", "Clear sky"),
            (2, 1,"Liquid drops", "Cloud liquid droplets only"), 
            (20, 1,"Drizzle | Rain", "Drizzle or rain"), 
            (20, 2,"Drizzle & Cloud", "Drizzle or rain coexisting with cloud liquid droplets"),
            (2, 3,"Ice", "Ice particles"), 
            (2, 2, "SLC & Ice", "Ice coexisting with supercooled liquid droplets"),
            (20, 3,"Melting", "Melting ice particle"),
            (20, 4,"Melting & cloud", "Melting ice particles coexisting with cloud liquid droplets"), 
            (40, 2,"Aerosol", "Aerosol particles/no cloud or precipitation"), 
            (40, 3,"Insects", "Insects/no cloud or precipitation"), 
            (40, 4,"Aerosol & Insects", "Aerosol coexisting with insects/no cloud or precipitation")
        ];
    elseif LegendType == "detection"
        txt_labels = [
            (1, 1, "Ra + Li", "Good Radar & Lidar");
            (1, 2, "Ra", "Good Radar Echo"); 
            (4, 1, "Li", "Lidar Echo Only");
            (4, 2, "Ra + Att", "Radar corrected \nfor liquid attenuation");
            (7, 1, "Ra - Att", "Radar uncorrected \nfor liquid attenuation");
            (7, 2, "Ra clutter", "Radar ground clutter")
        ];
    else
        error("$LegendType not supported!")
    end 
    cldnet = CloudNetPalette(LegendType)

    scatter(map(x->(x[1], 1.5x[2]), txt_labels),
            marker=:square, grid=:false, markerstroke=:none,
            framestyle=:none, widen=true,
            markersize=8, markercolor=cldnet, legend=false, axis=false,
            clipping=:false, xlim=(-0.1, 60), ylim=(0, 8)); #1+maximum(txt_labels)[2]));
    
    tmp = map(x->(x[1]+1.5, 1.5x[2], text(x[3], 11, :left)), txt_labels) |> annotate!;

    return tmp
end
# ----/

function Attach_Isotherms(pltin, Xin, Yin, Zvar, BB, sp; maxhgt=(0,8), TLEV=10)
    #TLEV = [5, 0, -5, -10, -15, -20, -25, -30]
    
    Plots.contour!(pltin, Xin, Yin, Zvar,
                   levels = TLEV, ylim=maxhgt, xlim = extrema(Xin), ticks=:none,
                   linecolor=:black, alpha=:.5, lw=1, contour_labels = true,
                   axis=false, colorbar=:none, subplot=sp, inset=BB, 
                   background_color_subplot=:transparent);

    return pltin
end
# ----/

# **********************************************************
# Subroutine to attach the wind vector to plot
function Attach_Windvector(pltin, Xin, Yin, Uin, Vin, BB, sp; maxhgt=(0,8))

    TT, HH, UU, VV = prepare_quiver(Xin, Yin, Uin, Vin)
    
    Plots.quiver!(pltin, TT[:], HH[:], quiver=(UU[:],VV[:]),
                  ylim=maxhgt, axis=false, color=:gray, alpha=0.3,
                  ticks=:none, shape=:none, inset=BB, subplot=sp,
                  background_color_subplot=:transparent);

    return pltin
end
# ----/                                
                                
# ***********************************************************
# Subroutine to create cascate plot of profile variables
function Attach_Profile_Cascate(pltin, ts, hkm, Z, BB, sp; maxhgt=(0, 8), δts=1.5)

    Znn = normalize_2Dvar(Z, δx=δts);
    Plots.plot!(pltin, Znn .+ ts', hkm, xlim=extrema(ts), ylim=maxhgt,
                color=:lightblue, legend=false, axis=:none, ticks=:none,
                inset=BB, background_color_subplot=:transparent, subplot=sp)

    return pltin
end
# ----/

# ************************************************************
# ++++++++++ AUXILIARY FUNCTIONS FOR VISUALIZATION +++++++++++
# ************************************************************

# ************************************************************
# Prepare vectors to be used with quiver
using ..CloudnetTools

function prepare_quiver(x::Vector, y::Vector, u::Matrix, v::Matrix; factor=2)

    # 2D meshing of coordinates:
    X, Y = CloudnetTools.mymeshgrid(x, y)

    # normalizing the wind variables
    UU, VV = let WS = @. factor*sqrt(u^2 + v^2)
    
        u./WS, v./WS
    end

    return X, Y, UU, VV
    
end
# ----/

# ****************
function normalize_2Dvar(X::Matrix; δx=1.5)

    δx < 1 && error("Optional input δx needs to be ≥ 1")
    T01 = [];
    for T ∈ eachcol(X)
        let lims = extrema(T)
            (2(T .- lims[1])./(lims[2] - lims[1]) .- 1)δx |> x->push!(T01, x)
            #push!(T01, (T .- lims[1])./(lims[2] - lims[1]) |> v-> (2δx.*v) .- δx)
        end
    end
    return hcat(T01...)
end
# ----/


end  # end of Module Vis
#
# end of file
