var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = YAXArraysToolbox","category":"page"},{"location":"#YAXArraysToolbox","page":"Home","title":"YAXArraysToolbox","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for YAXArraysToolbox.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [YAXArraysToolbox]\n\nModules = [BasicOperations]\n\nModules = [SpaceTimeAnalyses]","category":"page"},{"location":"#YAXArraysToolbox.aggregate_time-Tuple{Any}","page":"Home","title":"YAXArraysToolbox.aggregate_time","text":"Aggregate by time\n\nArguments:\n\ncube_in YAXArray Cube.\ntime_axis: String. Name of the time axis.\nnew_resolution: String. New temporal resolution can be \"day\", \"month\", \"year\".\nnew_time_step: Int64. Time step to be computed in the new time series. e.g. new_resolution=\"day\", new_time_step=8 will compute the function each 8 days. The new time dimension will only contain the days corresponding to the 8th day.\nfun: String. Function to be applied to aggregate the time. It can be \"median\", \"mean\", \"std\", \"var\", \"sum\", \"quant\", \"min\", \"max\".\np: Float64 in the interval [0,1]. If fun=quant p is the value of the quantile. \nskipMissing: Boolean. Skip missing values when aggregating the data. If all values are missing, NaN is returned.\nskipnan: Boolean. Skip NaN values when aggregating the data. If all values are NaN, NaN is returned.\nshowprog: Boolean. Progress Bar.\nmax_cache: String. Maximum cache to read the data. It must be in MB e.g. \"100MB\" or in GB \"10GB\".\n\nExamples\n\n\nesds = open_dataset(\"https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-184x90x90-2.1.1.zarr\")\nesdc = Cube(esds)\n\n# Estimating the monthly LAI\n\nlai_month = aggregate_time(esdc[Variable = \"leaf_area_index\"]; time_axis = \"time\", new_resolution = \"month\", new_time_step=1, fun=\"mean\", p=nothing, skipMissing=true, skipnan=true, showprog=true, max_cache=\"1GB\")\n\n\n\n\n\n\n","category":"method"},{"location":"#YAXArraysToolbox.altitude_mask_results_proc-Tuple{Any}","page":"Home","title":"YAXArraysToolbox.altitude_mask_results_proc","text":"Topographical variability processor\n\nArguments:\n\ncube_in_altitude : Altitude YAXARRAY with two variables mean, and sd.\nlon_axis_name : String. Name of the longitude axis on the input cubes. By default lon_axis_name = \"lon\"\nlat_axis_name :  String. Name of the longitude axis on the input cubes. By default lon_axis_name = \"lat\"\nvariable_name : String. Name of the Variable containing the variables \"mean\", and \"sd\".\nwinsize: Edge size of the moving window on pixels. By default winsize = 5. E.g. winsize = 5 will produce a moving window with 5^2 pixels.\nshowprog: Show progress bar. By default showprog = true\n\nOutput:\n\nThe Topographical variability processor produces a YAXARRAY.cube with three Indicators:\n\nv1: v_1=dfrac1nsum_i = 1^nsigma_hi\n\nHigh values of v_1 indicate hilly terrain over the considered scale, which should be discarded from the analysis.\n\nv2: v_2=vert mu_h - dfrac1n sum_i = 1^nmu_hi vert\n\nv_2 indicates how different the mean elevation within the central pixel is from the average elevation in the local window.\n\nv3: v_3= vert sigma_h - dfrac1nsum_i = 1^nsigma_hi vert\n\nv_3: compares the central pixel's standard deviation of elevation with the standard deviation across the moving window.\n\nSee also:\n\naltitude_masking_proc TO ADD LINK!!\n\nBibliography:\n\nDuveiller, G., Hooker, J., & Cescatti, A. (2018). A dataset mapping the potential biophysical effects of vegetation cover change. Scientific Data, 5(1), Article 1. https://doi.org/10.1038/sdata.2018.14\n\n\n\n\n\n","category":"method"},{"location":"#YAXArraysToolbox.altitude_masking_proc-Tuple{Any, Any}","page":"Home","title":"YAXArraysToolbox.altitude_masking_proc","text":"Topographical masking processor\n\nArguments:\n\ncube_in_to_mask: YAXArray Cube to be masked.\ncube_in_altitude : Altitude YAXARRAY with two variables mean, and sd.\nlon_axis_name : String. Name of the longitude axis on the input cubes. By default lon_axis_name = \"lon\"\nlat_axis_name :  String. Name of the longitude axis on the input cubes. By default lon_axis_name = \"lat\"\nvariable_name : String. Name of the Variable containing the variables \"mean\", and \"sd\".\ntime_axis_name : String or NaN. It is strongly recommended to pass this parameter if the cube to be masked contains a time dimension, otherwise nothing.\nwinsize: Edge size of the moving window on pixels. By default winsize = 5. E.g. winsize = 5 will produce a moving window with 5^2 pixels.\nv1_thr : Float. Threshold to mask values using v_1 indicator. All values higer or equal to v1_thr are set to NaN. By default v1_thr = 50 \nv2_thr : Float. Threshold to mask values using v_2 indicator. All values higer or equal to v2_thr are set to NaN. By detault v2_thr = 100\nv3_thr : Float. Threshold to mask values using v_3 indicator. All values higer or equal to v3 are set to NaN. By default v3_thr = 100\nshowprog: Boolean. Show progress bar. By default showprog = true.\n\nOutput:\n\nYAXArray Datase with two variables:\ncube masked: YAXArray Cube with same dimensions as cube_in_to_mask.\nmasked_pixels: YAXArray Cube with same lat, lon, dimensions as cube_masked but with a single boolean variable indicating if the pixel was masked or not.\n\nTopographical variability indicators\n\nv1: v_1=dfrac1nsum_i = 1^nsigma_hi\n\nHigh values of v_1 indicate hilly terrain over the considered scale, which should be discarded from the analysis.\n\nv2: v_2=vert mu_h - dfrac1n sum_i = 1^nmu_hi vert\n\nv_2 indicates how different the mean elevation within the central pixel is from the average elevation in the local window.\n\nv3: v_3= vert sigma_h - dfrac1nsum_i = 1^nsigma_hi vert\n\nv_3: compares the central pixel's standard deviation of elevation with the standard deviation across the moving window.\n\nSee also\n\naltitude_mask_results_proc function TO ADD LINK!!\n\nBibliography\n\nDuveiller, G., Hooker, J., & Cescatti, A. (2018). A dataset mapping the potential biophysical effects of vegetation cover change. Scientific Data, 5(1), Article 1. https://doi.org/10.1038/sdata.2018.14\n\n\n\n\n\n","category":"method"},{"location":"#YAXArraysToolbox.masking_proc-Tuple{YAXArrays.Cubes.YAXArray}","page":"Home","title":"YAXArraysToolbox.masking_proc","text":"Masking processor\n\nArguments:\n\ncube_in_to_mask: YAXArray cube to be masked.\ncube_rsquare: Nothing, or YAXArray cube with the R^2 variable. If set to nothing no mask is applied\nrsquare_thr: Float64. R^2 threshold. All values lower than rsquare_thr are set to NaN\ncube_co_occurrence: Nothing, or YAXArray cube with the co-occurrence variable. If set to nothing no mask is applied.\nco_occurence_thr: Float64. Co-occurence threshold. All values lower than co_occurence_thr are set to NaN\ncube_delta: Nothing, or YAXArray cube with delta variable. If set to nothing no mask is applied.\nminmax_delta: Tuple. Minimum and maximum thresholds of delta variable. Values lower and higher than the thresholds are set to NaN. It is also possible to set any of the thresholds as nothing e.g. (-1, nothing) or (nothing, 1) in these cases only one threshold is applied.\ntime_dim: Nothing, or String. Name of the time dimension. This dimensions needs to be present in all the cubes. If set to nothing no time dimension considered (It can result in slower computation time!). By default time_dim = time\nshowprog: Boolean. Show progress bar. By default showprog = true\n\nOutput:\n\nYAXArray cube masked.\n\n\n\n\n\n","category":"method"},{"location":"#YAXArraysToolbox.masking_space-Tuple{Any, Any}","page":"Home","title":"YAXArraysToolbox.masking_space","text":"Masking using spatial dimension\n\nThe masked vales are set as NaN!!.\n\nArguments:\n\ncube_in YAXArray Cube to be masked.\n\n-mask YAXArray Cube without time dimension and with a single variable to be used as mask. All values equal to NaN or missing will be masked in cubein. The mask will be applied to all the variables and time steps presented in ```cubein```.\n\nlat_axis: String. Name of the latitude axis.\nlon_axis: String. Name of the longitude axis.\nval_mask: NaN or missing. Value present in mask to be used as reference to mask cube_in. Must be NaN or missing.\nshowprog: Boolean. Progress Bar.\nmax_cache: String. Maximum cache to read the data. It must be in MB e.g. \"100MB\" or in GB \"10GB\".\n\nExamples\n\nusing YAXArrays, Zarr\n\n\naxlist = [\nRangeAxis(\"time\", range(1, 20, length = 20)),\nRangeAxis(\"x\", range(1, 10, length = 10)),\nRangeAxis(\"y\", range(1, 5, length = 15)),\nCategoricalAxis(\"Variable\", [\"var1\", \"var2\"]),\n]\n\n\ndata = rand(20, 10, 15, 2)\n\n\nds = YAXArray(axlist, data, props)\n\naxlist = [\nRangeAxis(\"x\", range(1, 10, length = 10)),\nRangeAxis(\"y\", range(1, 5, length = 15)),\nCategoricalAxis(\"Variable\", [\"var1\"]),\n]\n\n\ndata = rand(10, 15, 1)\n\nds_mask = YAXArray(axlist, data)\n\nmasking_space(ds, ds_mask; lat_axis = \"x\", lon_axis = \"y\")\n\n\n\n\n\n","category":"method"},{"location":"#YAXArraysToolbox.masking_time-Tuple{Any}","page":"Home","title":"YAXArraysToolbox.masking_time","text":"Masking using time dimension.\n\nThe function implements two methods:\n\nMasking based on a threshold value for one of the variables presented in the cube. e.g., masking the values of all the variables presented in the cube where radiation is lower than X.\nMasking based on the quantile threshold, where the quantile is estimated using the time series for each one of the variables presented in the cube.\n\nThe masked vales are set as NaN.\n\nArguments:\n\ncube_in YAXArray Cube.\ntime_axis: String. Name of the time axis.\nvar_axis: String. Name of the axis containing the variables.\nvar_mask: String or nothing. Name of the variable to be used to mask the other variables. If String val must be an Int64 or Float64 number. If nothing, val must be nothing and p must be a Float64 in the interval [0,1].\nval: Float64 or nothing. The value of the threshold in var_mask to be used to mask all the variables in the cube. If var_mask = nothing then, val=nothing\np: Float64 or nothing. Quantile value used as a threshold to mask the variables.\ncomp: String. Standard comparison operation between the threshold value and each one of the elements. comp Must be one of the following: \"==\", \"!=\" \"<\", \"<=\", \">\", \">=\".\nshowprog: Boolean. Progress Bar.\nmax_cache: String. Maximum cache to read the data. It must be in MB e.g. \"100MB\" or in GB \"10GB\".\n\nExamples\n\nusing YAXArrays, Statistics, Zarr, NetCDF, YAXArraysToolbox\n\nesds = open_dataset(\n    \"https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-184x90x90-2.1.1.zarr\",\n)\nesdc = Cube(esds)\n\nesdc_small = esdc[\n    lon = (-86, -35),\n    lat = (-56, 14),\n    time = (Date(2010), Date(2014)),\n    Variable = [\"leaf_area_index\", \"sensible_heat\", \"potential_evaporation\"],\n]\n\ntest = masking_time(\n    esdc_small;\n    time_axis = \"time\",\n    var_axis = \"Variable\",\n    var_mask = \"leaf_area_index\",\n    val = 0.2,\n    comp = \"<\",\n    showprog = true,\n    max_cache = \"1GB\",\n)\n\nplot_time(esdc_small; time_axis=\"time\", var_axis=\"Variable\", var = \"leaf_area_index\", lat_axis = \"lat\", lon_axis=\"lon\", fun = \"min\")\n\nplot_time(test; time_axis=\"time\", var_axis=\"Variable\", var = \"leaf_area_index\", lat_axis = \"lat\", lon_axis=\"lon\", fun = \"min\")\n\n\n\n\n\n","category":"method"},{"location":"#YAXArraysToolbox.plot_space-Tuple{YAXArrays.Cubes.YAXArray}","page":"Home","title":"YAXArraysToolbox.plot_space","text":"Plot Space/Maps\n\nArguments\n\ncube_in: YAXArray Cube.\ntime_axis: String. Name of the time axis.\nvar_axis: String. Name of the axis containing the variables.\nvar: String or nothing. Name of the variable to be plotted. If nothing all the variables presented in the cube are plotted.\nlat_axis: String. Name of the latitude axis.\nlon_axis: String. Name of the longitute axis.\nfun: String. Name of the function used to collapse the spatial dimensions. It must be \"median\", \"mean\", \"std\", \"var\", \"sum\", \"quant\", \"min\", or \"max\".\np: Float64. in the interval [0,1]. If fun=quant p is the value of the quantile.\ncolormap: Color Map. By default: colormap = Reverse(:batlow)\nresolution: Plot resolution. By default resolution = (800, 300).\nncol: Number of plots by column. By default ncol = 1.\nnrow: Number of plots by row. By default ncol = 1.\nshowprog: Boolean. Progress Bar.\nmax_cache: String. Maximum cache to read the data. It must be in MB e.g. \"100MB\" or in GB \"10GB\".\n\nExamples\n\n\ncube_in = open_dataset(\n    \"https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-1x720x1440-2.1.1.zarr\",\n)\n\ncube_in = Cube(cube_in)\n\n\ncube_in = cube_in[\n    lon = (-9.0, 0.0),\n    lat = (35, 40),\n    time = (Date(2010), Date(2014)),\n    Variable = [\"leaf_area_index\", \"sensible_heat\"],\n]\n\nplot_space(cube_in; time_axis = \"time\", resolution = (900, 600), var_axis = \"Variable\", var =  \"leaf_area_index\", fun = \"median\")\n\n\nmetric = [\"median\", \"mean\", \"std\", \"var\", \"sum\", \"quant\", \"min\", \"max\"]\n\n\nfor i in eachindex(metric)\n    println(metric[i])\n    plot_space(\n        cube_in;\n        time_axis = \"time\",\n        var_axis = \"Variable\",\n        lon_axis = \"lon\",\n        lat_axis = \"lat\",\n        var = \"sensible_heat\",\n        fun = metric[i],\n        p = 0.2,\n        showprog = true,\n        max_cache = \"100MB\",\n    )\nend\n\n\n\nplot_space(\n    cube_in;\n    time_axis = \"time\",\n    var_axis = \"Variable\",\n    lon_axis = \"lon\",\n    lat_axis = \"lat\",\n    var = nothing,\n    fun = \"median\",\n    resolution = (1200, 600),\n    p = 0.2,\n    showprog = true,\n    max_cache = \"100MB\",\n    ncol = 2,\n)\n\n\n\n\n\n\n","category":"method"},{"location":"#YAXArraysToolbox.plot_time-Tuple{YAXArrays.Cubes.YAXArray}","page":"Home","title":"YAXArraysToolbox.plot_time","text":"Plot time\n\nThe function allow to plot the time series of a given variables in a cube or all the variables present in a cube. As is expected that cubes contain spatial dimensions the spatial dimensions are collapsed using a function e.g., estimating the mean of the variable using the pixels of a certain area for each time step.\n\nArguments:\n\ncube_in YAXArray Cube.\ntime_axis: String. Name of the time axis.\nvar_axis: String. Name of the axis containing the variables.\nvar: String or nothing. Name of the variable to be plotted. If nothing all the variables presented in the cube are plotted.\nlat_axis: String. Name of the latitude axis.\nlon_axis: String. Name of the longitute axis.\nfun: String. Name of the function used to collapse the spatial dimensions. It must be \"median\", \"mean\", \"std\", \"var\", \"sum\", \"quant\", \"min\", or \"max\".\np: Float64. in the interval [0,1]. If fun=quant p is the value of the quantile. \nresolution: Tuple. Plot resolution. By default resolution = (600, 400). \nncol: Number of plots by column. By default ncol = 1.\nnrow: Number of plots by row. By default ncol = 1.\nshowprog: Boolean. Progress Bar.\nmax_cache: String. Maximum cache to read the data. It must be in MB e.g. \"100MB\" or in GB \"10GB\".\n\nExamples\n\ncube_in = open_dataset(\n    \"https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-1x720x1440-2.1.1.zarr\",\n)\n\ncube_in = Cube(cube_in)\ncube_in.Variable\ncube_in = cube_in[\n    lon = (-9.0, 0.0),\n    lat = (35, 40),\n    time = (Date(2010), Date(2014)),\n    Variable = [\"leaf_area_index\", \"sensible_heat\"],\n]\n\nplot_time(\n    cube_in;\n    time_axis = \"time\",\n    var_axis = \"Variable\",\n    lon_axis = \"lon\",\n    lat_axis = \"lat\",\n    var = nothing,\n    fun = \"median\",\n    resolution = (900, 600),\n    p = 0.2,\n    showprog = true,\n    max_cache = \"100MB\",\n    ncol = 2\n)\n\nplot_time(\n    cube_in;\n    time_axis = \"time\",\n    var_axis = \"Variable\",\n    lon_axis = \"lon\",\n    lat_axis = \"lat\",\n    var = \"sensible_heat\",\n    fun = \"median\",\n    p = 0.2,\n    showprog = true,\n    max_cache = \"100MB\",\n)\n\n\n\n\n\n","category":"method"},{"location":"#YAXArraysToolbox.s4time-NTuple{6, Any}","page":"Home","title":"YAXArraysToolbox.s4time","text":"space4time(climatecube, pftscube, pft_list::Vector{String}, winsize = 5, minpxl = 100, minDiffPxlspercentage = 40)\n\nCompute the space for time analysis for a given climate variable.     ...\n\nArguments\n\nclimate_cube: YAXARRAY cube with dimenssions: lon, lat, time.\npfts_cube: YAXARRAY cube with dimenssions: pfts, lat,lon, time.\n...\n\nOutput\n\nThree output cubes are generated.\nout1: Summary statistics. YAXARRAY cube where summary_stat axis contains: \n* \"rsquare\": XXXX \n* \"cumulative_variance\": XXXX\n* \"predicted\": Prediction of Z for the central pixel with its real PFT combination.\nout2:\n\n#Examples\n\n\n\n\n\n","category":"method"},{"location":"#YAXArraysToolbox.space4time_proc-Tuple{Any, Any}","page":"Home","title":"YAXArraysToolbox.space4time_proc","text":"Space for time processor\n\nArguments:\n\ncube_con : YAXARRAY with the continous variable to be analyized.\n\ncube_classes: YAXARRAY with the discrete classes to be used in the space4time.\n\ntime_axis_name : String or nothing. Name of the time axis on the input cubes. By default time_axis_name = \"time\". if time_axis_name = nothing, not time dimension considered.\n\nlon_axis_name : String. Name of the longitude axis on the input cubes. By default lon_axis_name = \"lon\"\n\nlat_axis_name :  String. Name of the longitude axis on the input cubes. By default lon_axis_name = \"lat\"\n\nclasses_var_name : String. Name of the Variable containing the discrete classes. By default classes_var_name = \"classes\".\n\nwinsize: Edge size of the moving window on pixels. By default winsize = 5. E.g. winsize = 5 will produce a moving window with 5^2 pixels.\n\nminpxl : Minimum number of pixels in the moving window. By default minpxl = 25. Change accordindly to your winsize parameter.\n\nminDiffPxlspercentage: Percentage of minimum number pixels in the moving window that must have different compositions. Must be any value in the interval 30-100. By default minDiffPxlspercentage = 40\n\nclasses_vec: A string vector with the names of the classes on cube_classes to be used. e.g. from MPI-BGC internal structure classes_vec = [\"Evergreen_Needleleaf_Forests\", \"Evergreen_Broadleaf_Forests\", \"Deciduous_Needleleaf_Forests\", \"Deciduous_Broadleaf_Forests\", \"Mixed_Forests\", \"Closed_Shrublands\", \"Open_Shrublands\", \"Woody_Savannas\", \"Savannas\", \"Grasslands\", \"Permanent_Wetlands\", \"Croplands\", \"Urban_and_Built-up_Lands\", \"Cropland/Natural_Vegetation_Mosaics\", \"Permanent_Snow_and_Ice\", \"Barren\", \"Water_Bodies\"]\nmax_value: Indicates if the scale of the presence of the discrete classes if from 0 to 1 or 0 to 100 if max_value = 100 then the data is re-scaled from 0 to 1. By default max_value = 1\nshowprog: Show progress bar. By default showprog = true\nmax_cache: Size of the cache to allocate temporarily sections of the cubes. By default max_cache = 1e8\n\nOutput:\n\nThe space4time_proc produces a YAXARRAY.Dataset with three cubes:\n\nSummaryStats cube has one axis summary_stat, and three variables:\nrsquared:  \ncumulative_variance:\npredicted:\nmetrics_for_classes cube has one axis Values of Z for pure classes, and two variables:\nestimated:\nestimated_error:\nmetricsfortransitions has two axis transitions (all the transitions by pairs between the different classes), and Differences with three variables:\ndelta:\ndelta_error:\ncoocurence:\n\n\n\n\n\n\n\n","category":"method"}]
}