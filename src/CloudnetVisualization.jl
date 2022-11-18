# Part of CloudnetTools.jl
# ******************************************************************
# +++++++++++++++ MODULE Vis  for visulaization ++++++++++++++++++++
# ******************************************************************

module Vis

using Plots
using Dates
using Printf

function show_LWC_IWC(LWC::Dict, IWC::Dict; SITENAME="", mxhgt=10, twoplots=false, showisoT=true, savefig=:none, cnt::Dict, extras=Dict())
    
    # defining parameters for LWC:
    CLiqLIM = (1f-2, 0.5f1) #(-2, 1);
    CLiqCOL = cgrad(:autumn1, 7, categorical=true);

    # defining parameters for IWC:
    CIceLIM = (1f-2, 0.5f1) #(-3, -1);
    CIceCOL = cgrad(:Paired_3, 7, categorical=true);

    # defining common parameters:
    tm_tick = IWC[:time][1]:Minute(90):IWC[:time][end];
    str_tick = Dates.format.(tm_tick, "H:MM");

    Xstrname = "TIME UTC from $(SITENAME) "*Dates.format(IWC[:time][2], "dd.u.yyyy")

    strtitle = "log₁₀ water content / g m⁻³";

    XLIM = extrema(IWC[:time])

    pice = plot(IWC[:time], 1f-3IWC[:height], 1f3IWC[:iwc], st=:heatmap, color=CIceCOL, clim=CIceLIM, ylim=(0, mxhgt), xlim=XLIM, xticks=(tm_tick, ""), yticks=:none, colorbar=twoplots, colorbar_scale=:log10);

    # plot both together:
    if twoplots
        mxplt = plot(LWC[:time], 1f-3LWC[:height], 1f3LWC[:lwc], st=:heatmap, color=CLiqCOL, clim=CLiqLIM, ylim=(0, mxhgt), colorbar=twoplots, colorbar_scale=:log10, xlim=XLIM, xlabel=Xstrname, xticks=(tm_tick, str_tick), xrot=45, xtickfontsize=11, xguidefont=font(12), title=strtitle); #, background_color_subplot=:transparent);

    else
        mxplt = plot(pice, colorbar=:none);
    
        plot!(mxplt, LWC[:time], 1f-3LWC[:height], 1f3LWC[:lwc], st=:heatmap, color=CLiqCOL, inset=(1,bbox(0,0,1,1)), subplot=2, colorbar=:none, colorbar_scale=:log10, background_color_subplot=:transparent, ylim=(0, mxhgt), xlim=XLIM, xlabel=Xstrname, xticks=(tm_tick, str_tick), xrot=45, xtickfontsize=11, xguidefont=font(12), ylabel="Altitude / km", title=strtitle, tick_dir=:out, framestyle=:box, gridstyle=:dash);

    end

    # Preparing to add graphs for isotherms and wind vectors:
    ihmax = findlast(1f-3cnt[:model_height] .≤ mxhgt)
    
    Xin = haskey(cnt, :model_time) ? Vector(1:length(cnt[:model_time])) : round.(Int64, range(1, stop=length(cnt[:time]), length=25))
    
    Yin = 1f-3cnt[:model_height]
    #TLEV = [5, 0, -5, -10, -15, -20, -25, -30]
    BB = bbox(0,0,1,1)

    if showisoT
        # converting to Celcius in case Temperature is in K
        Tin = let T = cnt[:T][1:ihmax, Xin]
            any(T .> 100) ? T .- 273.15 : T
        end
        
        Tlevels = extrema(filter(!isnan,Tin)) |> x->ceil.(range(ceil(x[1]), stop=floor(x[2]), length=10))
        Attach_Isotherms(mxplt, Xin, Yin[1:ihmax], Tin,
                         (1, BB), 2, TLEV = Tlevels, maxhgt=(0, mxhgt))
    end

    
    # creating colorbars for merged or gridded plots:
    if twoplots
        outplt = plot(mxplt, pice, layout=grid(2,1), ylabel="Altitude / km", tick_dir=:out, framestyle=:box, gridstyle=:dash)# , background_color_subplot=:transparent);
    else
        return mxplt
        # * liquid
        liq_val = range(CLiqLIM[1], stop=CLiqLIM[2], length=7);
        cmliq = heatmap([-1], liq_val, liq_val', xlim=(1,2), color=CLiqCOL, clim=CLiqLIM, colorbar_title="liquid", colorbar_scale=:log10, tick_dir=:out, xticks=:none, xaxis=false, box=false, title="LWC");

        # * ice
        ice_val = range(CIceLIM[1], stop=CIceLIM[2], length=7);
        cmice = heatmap!(cmliq, [-1], ice_val, ice_val', xlim=(1,2), color=CIceCOL, clim=CIceLIM, colorbar_title="ice", colorbar_scale=:log10, tick_dir=:out, xticks=:none, xaxis=false, box=false, inset=(1, bbox(0,0,0.87,1)), subplot=2) #, title="IWC", right_margin=2Plots.mm);

        # creating output layout plot:
        ll = @layout [a{0.96w} b{0.02w} c{0.02w}];
        outplt = plot!(cmliq, mxplt, inset=(1, bbox(0,0,0.8,1)), subplot=3, size=(900, 600), dpi=600, margin=15Plots.mm; extras...);
        #outplt = plot(mxplt, cmliq, cmice, layout=ll, size=(900,600), dpi=300, bottom_margin=15Plots.mm, left_margin=[15Plots.mm 0Plots.mm 0Plots.mm 0Plots.mm]);
    end
    
    return outplt
end
# ----/


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
* showatm::Dict(:wind, :isoT, :procas) to add meteo data to the plot, "wind" & "iso-temp" or "Profile cascade of variable", default=(true, true, false)
Output:
* plt::Plot output plot object.
"""
function show_classific(cnt::Dict; SITENAME="", maxhgt=8, showlegend=true,
                        atmosplot=Dict(), showatm=Dict(:wind=>true, :isoT=>true, :procas=>false), savefig=:none, extras=Dict())

    # defining time axis ticks:
    tm_tick = cnt[:time][1]:Minute(120):cnt[:time][end];
    tm_lims = (cnt[:time][1], cnt[:time][end])
    
    Xstrname = let tmp = Date.(tm_lims) |> unique
        formatstr = if length(tmp)==1
            @sprintf("Time UTC [hour] from %s", tmp[1])
        else
            @sprintf("Time UTC [hour] from %s to %s", tmp[1], tmp[2])
        end
        !isempty(SITENAME) && (formatstr *= ", "*SITENAME)
        formatstr
    end #"TIME UTC [hour] from "*Dates.format(cnt[:time][2], "dd.u.yyyy")

    strtitle = (tm_tick[1], 7.5, text("Cloudnet Target Classification "*SITENAME, 11, halign=:left))
    strtitle = !isempty(SITENAME) ? "Cloudnet Target Classification "*SITENAME : ""
    
    cldnet = CloudNetPalette("classific")
    l = grid(2,1, heights=(0.2,0.8)) #[a{0.25h}; b];
    Y_LIM = (0, maxhgt)
    
    classplt = plot(cnt[:time], 1f-3cnt[:height], cnt[:CLASSIFY],
                    st=:heatmap, colorbar = false, framestyle = :box,
                    color=palette(cldnet, 11), ylim=Y_LIM, clim=(0,10),
                    ylabel="Height A.G.L. [km]", ytickfontsize=11, minorticks=true,
                    xlabel= Xstrname, tick_direction=:out, 
                    xticks=(tm_tick, Dates.format.(tm_tick, "H")), xrot=0,
                    xtickfontsize=13, xguidefontsize=18; extras...);  # xguidefont=font(12)

    # Preparing to add graphs for isotherms and wind vectors:

    atmos = isempty(atmosplot) ? cnt : atmosplot

    ihmax = findlast(1f-3atmos[:model_height] .≤ maxhgt)
    Xin = haskey(atmos, :model_time) ? Vector(1:length(atmos[:model_time])) : unique(round.(Int64, range(1, stop=length(atmos[:time]), length=25)))
    
    Yin = 1f-3atmos[:model_height]
    #TLEV = [5, 0, -5, -10, -15, -20, -25, -30]
    BB = bbox(0,0,1,1)

    Nsubplt = 2
    if showatm[:isoT]
        # converting to Celcius in case Temperature is in K
        Tin = let T = atmos[:T][1:ihmax, Xin]
            any(T .> 100) ? T .- 273.15 : T
        end
        
        #Tlevels = extrema(Tin) |> x->collect(ceil(x[1]):5:floor(x[2]))
        Tlevels = extrema(filter(!isnan,Tin)) |> x->ceil.(range(ceil(x[1]), stop=floor(x[2]), length=7))
        classplt = Attach_Isotherms(classplt, Xin, Yin[1:ihmax], Tin,
                                    (1, BB), Nsubplt, TLEV = Tlevels)
        Nsubplt += 1
    end

    if showatm[:wind]
        classplt = Attach_Windvector(classplt, Xin[1:2:end], Yin[2:4:ihmax],
                                     atmos[:UWIND][2:4:ihmax, Xin[1:2:end]],
                                     atmos[:VWIND][2:4:ihmax, Xin[1:2:end]],
                                     (1, BB), Nsubplt)
        Nsubplt += 1
    end
    #
    #θv = @. cnt[:T][1:ihmax,2:2:end]*(1024f0/cnt[:Pa][1:ihmax,2:2:end])^0.286;
    if showatm[:procas]
        Tin = cnt[:T]
        Attach_Profile_Cascate(classplt, Xin[2:2:end], Yin[1:ihmax], Tin, (1,BB), Nsubplt)
        Nsubplt += 1
    end
    
    if showlegend
        classleg = ShowLegendCloudNetClassification("classific")

        pltout = plot(classleg, classplt, layout=l, size=(800,700), title=strtitle)
    else
        pltout = plot(classplt, size=(800,600))
    end

    savefig != :none && typeof(savefig) <: String && Plots.savefig(pltout, savefig)
    
    return pltout
end
# -- OR
function show_classific(cnt_file::String; SITENAME="", maxhgt=8, atmosplot=Dict(), showlegend=true, savefig=:none)
    cnt = CloudnetTools.readCLNFile(cnt_file)
    return show_classific(cnt, SITENAME=SITENAME, maxhgt=maxhgt, atmosplot=atmosplot, showlegend=showlegend,savefig=savefig)
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
            (41, 1,"Clear sky", "Clear sky"),
            (2, 1,"Liquid droplets", "Cloud liquid droplets only"), 
            (18, 1,"Drizzle | rain", "Drizzle or rain"), 
            (18, 2,"Drizzle & liquid droplets", "Drizzle or rain coexisting with cloud liquid droplets"),
            (2, 3,"Ice", "Ice particles"), 
            (2, 2, "Ice & SCL", "Ice coexisting with supercooled liquid droplets"),
            (18, 3,"Melting ice", "Melting ice particle"),
            (18, 4,"Melting ice & liquid droples", "Melting ice particles coexisting with cloud liquid droplets"), 
            (41, 2,"Aerosol", "Aerosol particles/no cloud or precipitation"), 
            (41, 3,"Insects", "Insects/no cloud or precipitation"), 
            (41, 4,"Aerosol & insects", "Aerosol coexisting with insects/no cloud or precipitation")
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

function Attach_Isotherms(pltin, Xin, Yin, Zvar, BB, sp; maxhgt=(0,8), TLEV=10, extras=Dict())
    #TLEV = [5, 0, -5, -10, -15, -20, -25, -30]
    
    pltin = Plots.contour!(pltin, Xin, Yin, Zvar,
                           levels = TLEV, ylim=maxhgt, xlim = extrema(Xin), ticks=:none,
                           linecolor=:black, alpha=0.5, lw=1, contour_labels=true,
                           labelsfontsize=12,
                           axis=false, colorbar=:none, subplot=sp, inset=BB, 
                           labelscolor=:black, background_color_subplot=:transparent; extras...);

    return pltin
end
# ----/

# **********************************************************
# Subroutine to attach the wind vector to plot
function Attach_Windvector(pltin, Xin, Yin, Uin, Vin, BB, sp; maxhgt=(0,8), extras=Dict())

    TT, HH, UU, VV = prepare_quiver(Xin, Yin, Uin, Vin)
    
    pltin = Plots.quiver!(pltin, TT[:], HH[:], quiver=(UU[:],VV[:]),
                          ylim=maxhgt, axis=false, color=:gray, alpha=0.3,
                          ticks=:none, shape=:none, inset=BB, subplot=sp,
                          background_color_subplot=:transparent; extras...);

    return pltin
end
# ----/                                
                                
# ***********************************************************
# Subroutine to create cascate plot of profile variables
function Attach_Profile_Cascate(pltin, ts, hkm, Z, BB, sp; maxhgt=(0, 8), δts=1.5, extras=Dict())

    Znn = normalize_2Dvar(Z, δx=δts);
    Plots.plot!(pltin, Znn .+ ts', hkm, xlim=extrema(ts), ylim=maxhgt,
                color=:lightblue, legend=false, axis=:none, ticks=:none,
                inset=BB, background_color_subplot=:transparent, subplot=sp; extras...)

    return pltin
end
# ----/

# ************************************************************
# Functions to plot the data
#
function show_measurements(cln::Dict; atmosplot=Dict(), SITENAME::String="", maxhgt=8, savefig=:none, showclassific=false, extras=Dict())
    # converting to access the Cloudnet data:
    radar = Dict(K => cln[K] for K in [:time, :height, :Ze])
    lidar = Dict(K => cln[K] for K in [:time, :height, :β])
    mwr = Dict(K => cln[K] for K in [:time, :LWP])

    rs_time = haskey(cln, :model_time) ? :model_time : :time

    rs_list = [rs_time, :model_height, :T, :UWIND, :VWIND]
    rs = isempty(atmosplot) ? Dict(K => cln[K] for K in rs_list) : atmosplot
    #atmosplot && (rs[:model_height] *= 1f-3)  # converting to km

    if showclassific
        return show_measurements(radar, lidar, mwr, atmosplot=rs, SITENAME=SITENAME, maxhgt=maxhgt, savefig=savefig, cln=cln, extras=extras)
    else
        return show_measurements(radar, lidar, mwr, atmosplot=rs, SITENAME=SITENAME, maxhgt=maxhgt, savefig=savefig, extras=extras)
    end
end
# --- OR 
function show_measurements(radar::Dict, lidar::Dict, mwr::Dict; atmosplot::Dict=Dict(),
        SITENAME::String="", maxhgt=8, savefig=false, cln=nothing, extras=Dict())

    Y_LIM = (0, maxhgt)
    #X_LIM = extrema(radar[:time])
    tm_tick = mwr[:time][1]:Minute(120):mwr[:time][end]
    tm_lims = (minimum([mwr[:time][1], lidar[:time][1], radar[:time][1]]),
               maximum([mwr[:time][end], lidar[:time][end], radar[:time][end]]))
    
    # For the Radar:
    BB = bbox(0,0,.89,1)
    radarplt = Plots.plot(radar[:time], 1f-3radar[:height], radar[:Ze],
                          st=:heatmap, color=palette(:lighttest,20), clim=(-30, 10),
                          xlim=tm_lims, ylim=Y_LIM, tick_dir=:out, ytickfontsize=11,
                          colorbar_title="Radar Reflectivity [dBz]", #titlefontsize=11,
                          ylabel="Height A.G.L. [km]", xticks=(tm_tick, ""),
                          guidefontsize=13, ticksfontsize=13,
                          bottom_margin=-1.5Plots.mm, framestyle=:box);

    # adding atmospheric information
    atmos=atmosplot
    if !isempty(atmosplot)
        # Preparing to add graphs for isotherms and wind vectors:
        Xin = haskey(atmos, :model_time) ? Vector(1:length(atmos[:model_time])) : unique(round.(Int64, range(1, stop=length(atmos[:time]), length=25)))

        #Xin = round.(Int64, range(1, stop=length(atmos[:time]), length=48)) |> unique
        #X_LIM = (0, maximum(Xin))
    
       Yin = 1f-3atmos[:model_height]
        # Yin = collect(0:0.5:maxhgt)
        K_height = haskey(atmos, :model_height) ? atmos[:model_height] : atmos[:height]
        #ihmax = [argmin(abs.(x .- 1f-3K_height)) for x in Yin]
        ihmax = findlast(1f-3atmos[:model_height] .≤ maxhgt)

        # converting to Celcius in case Temperature is in K
        Tin = let T = atmos[:T][1:ihmax, Xin]
            any(T .> 100) ? T .- 273.15 : T
        end
        
        Tlevels = extrema(filter(!isnan,Tin)) |> x->ceil.(range(ceil(x[1]), stop=floor(x[2]), length=7))
        #Tlevels = [5, 0, -5, -10, -15, -20, -25, -30]
        Attach_Isotherms(radarplt, Xin, Yin[1:ihmax], Tin,
                         (1, BB), 2, TLEV = Tlevels)

        Attach_Windvector(radarplt, Xin[1:2:end], Yin[2:4:ihmax],
                          atmos[:UWIND][2:4:ihmax, Xin[1:2:end]],
                          atmos[:VWIND][2:4:ihmax, Xin[1:2:end]],
                          (1, BB), 3)
    end
    
    # Plot for LIDAR
    BB = bbox(0,0,.89,1)
    lidarplt = Plots.plot(lidar[:time], 1f-3lidar[:height], log10.(lidar[:β]),
                          st=:heatmap, color=:roma, clim=(-7, -4),
                          xlim=tm_lims, ylim=Y_LIM, tick_dir=:out, ytickfontsize=11, colorbar_width=1,
                          colorbar_title="Lidar Attenuated\nBackscattering coefficient log10 [sr⁻¹ m⁻¹]",
                          ylabel="Height A.G.L. [km]", xticks=(tm_tick, ""),
                          guidefontsize=13, subplot=1, bottom_margin=-1.5Plots.mm,
                          framestyle=:box);
    
    # adding atmospheric information
    if !isempty(atmosplot)
        Attach_Isotherms(lidarplt, Xin, Yin[1:ihmax], Tin,
                         (1, BB), 2, TLEV = Tlevels)
    
        Attach_Windvector(lidarplt, Xin[1:2:end], Yin[2:4:ihmax],
                          atmos[:UWIND][2:4:ihmax, Xin[1:2:end]],
                          atmos[:VWIND][2:4:ihmax, Xin[1:2:end]],
                          (1, BB), 3)
    end
    
    # For Radiometer LWP
    titletext = let tmp = Date.(tm_lims) |> unique
        formatstr = if length(tmp)==1
            @sprintf("Time UTC [hour] from %s", tmp[1])
        else
            @sprintf("Time UTC [hour] from %s to %s", tmp[1], tmp[2])
        end
        !isempty(SITENAME) && (formatstr *= ", "*SITENAME)
        formatstr
        #@sprintf(formatstr, (Dates.format.(tmp, "dd-uuu-yyyy")...), SITENAME);
    end
    
    mwrplt = Plots.plot(1,1, axis=nothing, border=:none, label=false);
    
    Plots.plot!(mwrplt, mwr[:time], mwr[:LWP],
                ylim=(0, max(200, maximum(mwr[:LWP]))),
                ylabel="LWP [g m⁻²]", l=2,
                label=false,  tick_dir=:out,
                yguidefont=font(:steelblue), ytickfontsize=11,
                xticks = (tm_tick, !isnothing(cln) ? "" : Dates.format.(tm_tick, "H")),
                xlim = tm_lims, inset=(1, BB), subplot=2,
                tickfontsize=13, xguidefontsize=15, framestyle=:box, #font(15),
                xlabel = !isnothing(cln) ? "" : titletext, xrot=0);
    #ytickfontcolor=:steelblue,

    # Composite figure:

    finplt = if !isnothing(cln) && typeof(cln)<:Dict
        ll = @layout [a0{0.3h}; a1{0.3h}; b{0.1h}; a2{0.9w} b1];
        clnplt = show_classific(cln; SITENAME="", maxhgt=maxhgt, showlegend=false)
 
        Plots.plot(radarplt, lidarplt,  mwrplt, clnplt, layout=ll,  link=:y,
                        size=(1000,1100), yguidefontsize=19, ytickfontsize=12,
                        left_margin =10Plots.mm, rigth_margin=10Plots.mm,
                        bottom_margin=[-3 -3 -3 -3 50].*Plots.mm)

    elseif !isnothing(cln) && typeof(cln)<:Plots.Plot
       ll = @layout [a0{0.3h}; a1{0.3h}; b{0.1h}; a2];
 
        Plots.plot(radarplt, lidarplt,  mwrplt, cln, layout=ll,  link=:y,
                        size=(1000,1100), yguidefontsize=19, ytickfontsize=12,
                        left_margin =10Plots.mm, rigth_margin=10Plots.mm,
                        bottom_margin=[-3 -3 -3 -3 50].*Plots.mm)
 
    else
    ll = @layout [a0{0.38h}; a{0.38h}; b{0.12h}];

    Plots.plot(radarplt, lidarplt,  mwrplt, layout=ll,  link=:y,
                        size=(1000,1000), yguidefontsize=13, ytickfontsize=12,
                        left_margin =10Plots.mm, rigth_margin=13Plots.mm,
                        bottom_margin=[-3 -3 50].*Plots.mm; extras...)
    #                        title = titletext, titlefontsize=15,
    #                        legend=[false false false], 
    end
    savefig != :none && typeof(savefig)<:String && Plots.savefig(finplt, savefig)
    
    return finplt

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
        let lims = filter(!isnan, T) |> extrema
            (2(T .- lims[1])./(lims[2] - lims[1]) .- 1)δx |> x->push!(T01, x)
            #push!(T01, (T .- lims[1])./(lims[2] - lims[1]) |> v-> (2δx.*v) .- δx)
        end
    end
    return hcat(T01...)
end
# ----/

# ****************
function convert_wind_to_U_V(WSP::Matrix, WDIR::Matrix)

    # converting from weather wind direction to math angle:
    MANG = 270.0 .- WDIR
    U = @. -abs(WSP)*sind(WDIR)  #cosd(MANG)
    V = @. -abs(WSP)*cosd(WDIR)  #sind(MANG)
        
    return U, V
end
# ----/

end  # end of Module Vis
#
# end of file
