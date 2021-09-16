# CloudNet classification file plotting
function 
CLOUDNET_PATH = "/home/psgarfias/LIM/data/CloudNet";
SUPERSITE = Dict(:nsa=>"output", :mosaic=>"MOSAiC/TROPOS");
LIDAR_TYPE = "HSRL";
MODEL_TYPE = "ECMWF";
PRODUCT = Dict(:class => "classific",
               :cate=>"categorize",
               :lwc=>"lwc",
               :iwc=>"iwc",
               :drizzle=>"drizzle"
               )

site_source = :mosaic;
year = 2020;
month = 7;
day = 30;

if site_source == :nsa
    classific_file = @sprintf("%04d/%02d/nsa_%s_classific_%04d%02d%02d.nc",
    year, month, MODEL_TYPE, year, month, day);
    classific_file = joinpath(CLOUDNET_PATH, SUPERSITE[site_source], LIDAR_TYPE, classific_file);
elseif site_source == :mosaic
    classific_file = @sprintf("%04d/%02d/%04d%02d%02d_polarstern_classification.nc",
    year, month, year, month, day);
    classific_file = joinpath(CLOUDNET_PATH, SUPERSITE[site_source], classific_file);
else
    println("$site_source not a supported site!")
end

isfile(classific_file)

# Opening the classification file:
ncid = NCDataset(classific_file, "r");
classific = ncid["target_classification"][:,:];
detection = ncid["detection_status"][:,:];
time = ncid["time"];
height = ncid["height"];

using Plots

function CloudNetPalette(ColorType::String)
    if ColorType == "classific"
        tmp = [
            (1.00,  1.0 ,  1.00);
            (0.44,  1.0 ,  0.92);
            (0.17,  0.62,  0.95);
            (0.75,  0.6 ,  1.0);
            (0.9 ,  0.9 ,  0.92);
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
            (0, 1,"Clear sky"),
            (0, 2,"Cloud liquid droplets only"), 
            (0, 3,"Drizzle or rain"), 
            (0, 4,"Drizzle or rain coexisting with cloud liquid droplets"),
            (0, 5,"Ice particles"), 
            (0, 6,"Ice coexisting with supercooled liquid droplets"),
            (5, 1,"Melting ice particle"),
            (5, 2,"Melting ice particles coexisting with cloud liquid droplets"), 
            (5, 3,"Aerosol particles/no cloud or precipitation"), 
            (5, 4,"Insects/no cloud or precipitation"), 
            (5, 5,"Aerosol coexisting with insects/no cloud or precipitation")
        ];
    elseif LegendType == "detection"
        txt_labels = [
            (1, 1, "Good Radar & Lidar");
            (1, 2, "Good Radar Echo"); 
            (4, 1, "Lidar Echo Only");
            (4, 2, "Radar corrected \nfor liquid attenuation");
            (7, 1, "Radar uncorrected \nfor liquid attenuation");
            (7, 2, "Radar ground clutter")
        ];
    else
        error("$LegendType not supported!")
    end 
    cldnet = CloudNetPalette(LegendType)
    xin = [zeros(6,1); 5*ones(5,1)]
    yin = [1:6; 1:5]
    scatter(map(x->(x[1], x[2]), txt_labels),
            marker=:square, grid=:false, markerstroke=:none,
            framestyle=:none, widen=true,
            markersize=9, markercolor=cldnet, legend=false, axis=false,
            clipping=:false, xlim=(-0.1,10), ylim=(0, 1+maximum(txt_labels)[2]));
    
    tmp = map(x->(x[1]+.3, x[2], text(x[3],9, :left)), txt_labels) |> annotate!;
    #map(x->(xin[x]+.3, yin[x], text(cld_labels[x],9, :left)), 1:11) |> annotate!;
    return tmp
end


tm_tick = time[1]:Minute(90):time[end];

xstringname = "TIME UTC from "*Dates.format(Date(year, month, day), "dd.mm.yyyy")
stringtitle = (tm_tick[1], 7.5, text("CloudNet Target Classification "*SUPERSITE[site_source], 11, halign=:left))
cldnet = CloudNetPalette("classific")
l = grid(2,1, heights=(0.2,0.8)) #[a{0.25h}; b];
classleg = ShowLegendCloudNetClassification("classific");
classplt = plot(time, height*1e-3, classific, st=:heatmap,
                colorbar = false, framestyle = :box,
                color=palette(cldnet, 11), ylim=(0, 8), clim=(0,10),
                ann = stringtitle,
                ylabel="Height A.G.L. [km]", ytickfontsize=10,
                xlabel= xstringname,
                tick_direction=:out, 
                xticks=(tm_tick, Dates.format.(tm_tick, "H:MM")),
                xtickfontsize=10, xguidefont=font(11));

plot(classleg, classplt, layout=l, size=(1000,800))

# For detection plot
stringtitle = (tm_tick[1], 7.5, text("CloudNet Target Detection "*SUPERSITE[site_source], 11, halign=:left))
cldnet = CloudNetPalette("detection")
classleg = ShowLegendCloudNetClassification("detection");
detecplt = plot(time, height*1e-3, detection, st=:heatmap,
                colorbar = false, framestyle = :box,
                color=palette(cldnet, 11), ylim=(0, 8), clim=(0,10),
                ann = stringtitle,
                ylabel="Height A.G.L. [km]", ytickfontsize=10,
                xlabel= xstringname,
                tick_direction=:out, 
                xticks=(tm_tick, Dates.format.(tm_tick, "H:MM")),
                xtickfontsize=10, xguidefont=font(11));

l = grid(2,1, heights=(0.15,0.85))
plot(classleg, detecplt, layout=l, size=(1000,800))

function ShowLegendCloudNetFlags()
tmp = [
 (0.247,  0.704,  0.43);
 (0.44 ,  0.926,  0.34); 
 (0.996,  0.91 ,  0.24);
 (0.8  ,  0.96 ,  0.96); 
 (0.455,  0.51 ,  0.41); 
 (0.89 ,  0.292,  0.13);
];

qinet = map(x->RGB(x...), tmp)

qilabels = [
    (1., 1.2, text("Good Radar & Lidar",:left));
    (1., 1.1, text("Good Radar Echo",:left)); 
    (1., 1., text("Lidar Echo Only",:left));
    (2.2, 1.2, text("Radar corrected \nfor liquid attenuation",:left));
    (2.2, 1.1, text("Radar uncorrected \nfor liquid attenuation",:left));
    (2.2, 1., text("Radar ground clutter",:left))
    ];
scatter(map(x->(x[1]-0.1,x[2]), qilabels), markersize=20, marker=:square, markercolor=qinet, legend=false, xlim=(0.8, 3.1), ylim=(0.9,1.3),axis=false);

outleg = annotate!(qilabels);
return outleg
end

# plotting backscattering vs depolarization-ratio
##histogram2d(log10.(β[:]), δ[:], xlims=(-8, -3), bins=(-8:.02:-3, 0:0.01:0.7), color=palette(:jet, 10), xlabel="log10 β [m⁻¹ sr⁻¹]", ylabel="δ")
    # end of script
