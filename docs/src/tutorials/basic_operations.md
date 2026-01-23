# [Basic Operations](@id basic_operations)

```@raw html
<div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 2em; border-radius: 8px; margin-bottom: 2em;">
  <h2 style="margin: 0; color: white;">📊 Getting Started with YAXArraysToolbox.jl</h2>
  <p style="margin: 0.5em 0 0 0; opacity: 0.9;">Learn the core functions for visualizing and processing Earth System data cubes</p>
</div>
```

**Author:** Daniel E. Pabon-Moreno

---

## What You'll Learn

In this tutorial, you will learn how to:

- ✅ Load data from the Earth System Data Cube (ESDC)
- ✅ Create time series plots with [`plot_time`](@ref)
- ✅ Generate spatial maps with [`plot_space`](@ref)
- ✅ Aggregate data to different temporal resolutions with [`aggregate_time`](@ref)

!!! tip "Prerequisites"
    Basic familiarity with Julia and the concept of multi-dimensional arrays is helpful but not required.

---

## Setup

First, let's load all the required packages:

```@example basic_ops
using Pkg
Pkg.instantiate()
using YAXArrays
Pkg.add(url="https://github.com/dpabon/YAXArraysToolbox.jl")
using YAXArraysToolbox
using CairoMakie
using GeoMakie
using Statistics
using Zarr
using Dates
using PythonCall
using DimensionalData
```

---

## Loading Data

We'll use the [Earth System Data Cube (ESDC)](https://www.earthsystemdatalab.net/), an analysis-ready data cube containing dozens of climate and Earth observation variables at global scale.

```@example basic_ops
esdc = open_dataset("https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-184x90x90-2.1.1.zarr")
esdc = Cube(esdc)
esdc
```

!!! info "About the ESDC"
    The Earth System Data Cube provides harmonized Earth observation and climate data at 0.25° spatial resolution and 8-day temporal resolution, spanning from 1979 to present.

---

## 1. Time Series Plotting with `plot_time`

The [`plot_time`](@ref) function creates time series plots by collapsing spatial dimensions using a specified statistic (mean, median, std, etc.).

### Selecting a Region

Let's focus on South America and select a time period:

```@example basic_ops
cube_to_plot = esdc[
    lon = -86 .. -35,
    lat = -56 .. 14,
    time = Date(2010) .. Date(2014),
    Variable = At("leaf_area_index", "sensible_heat"),
]
cube_to_plot
```

### Plotting All Variables

By default, `plot_time` plots all variables in the cube:

```@example basic_ops
plot_time(
    cube_to_plot;
    time_axis = :time,
    var_axis = :Variables,
    lon_axis = :lon,
    lat_axis = :lat,
    var = nothing,
    fun = "std",
    resolution = (900, 600),
    p = 0.2,
    showprog = true,
    max_cache = "1GB",
    ncol = 1,
    nrow = 2
)
```

### Plotting a Single Variable

You can also target a specific variable:

```@example basic_ops
plot_time(
    cube_to_plot;
    time_axis = :time,
    var_axis = :Variable,
    lon_axis = :lon,
    lat_axis = :lat,
    var = "sensible_heat",
    fun = "std",
    resolution = (900, 600),
    p = 0.2,
    showprog = true,
    max_cache = "1GB",
    ncol = 1,
    nrow = 2
)
```

### Available Statistics

The `fun` parameter supports multiple aggregation statistics:

| Statistic | Description |
|:----------|:------------|
| `"mean"` | Arithmetic mean |
| `"median"` | Median (50th percentile) |
| `"std"` | Standard deviation |
| `"var"` | Variance |
| `"sum"` | Sum of values |
| `"min"` | Minimum value |
| `"max"` | Maximum value |
| `"quant"` | Quantile (requires `p` parameter) |

```@example basic_ops
metrics = ["median", "mean", "std", "var", "sum", "quant", "min", "max"]

for metric in metrics
    fig = plot_time(
        cube_to_plot;
        time_axis = :time,
        var_axis = :Variable,
        lon_axis = :lon,
        lat_axis = :lat,
        var = "sensible_heat",
        fun = metric,
        resolution = (900, 600),
        p = 0.2,
        showprog = true,
        max_cache = "1GB",
        ncol = 1,
        nrow = 2
    )
    display(fig)
end
```

---

## 2. Spatial Mapping with `plot_space`

The [`plot_space`](@ref) function creates spatial maps by collapsing the time dimension.

### Single Variable Map

```@example basic_ops
cube_to_plot = esdc[
    lon = -86 .. -34,
    lat = -56 .. 14,
    time = Date(2010) .. Date(2014),
    Variable = At("leaf_area_index", "sensible_heat"),
]

plot_space(
    cube_to_plot;
    time_axis = :time,
    resolution = (900, 600),
    var_axis = :Variable,
    var = "leaf_area_index",
    fun = "median"
)
```

### Multiple Variables Side by Side

Set `var = nothing` to plot all variables:

```@example basic_ops
plot_space(
    cube_to_plot;
    time_axis = :time,
    resolution = (900, 600),
    var_axis = :Variables,
    var = nothing,
    ncol = 2,
    nrow = 1,
    fun = "median"
)
```

### Comparing Different Statistics

```@example basic_ops
metrics = ["median", "mean", "std", "var", "sum", "quant", "min", "max"]

for metric in metrics
    fig = plot_space(
        cube_to_plot;
        time_axis = :time,
        var_axis = :Variable,
        lon_axis = :lon,
        lat_axis = :lat,
        var = "sensible_heat",
        fun = metric,
        p = 0.2,
        showprog = true,
        max_cache = "100MB"
    )
    display(fig)
end
```

---

## 3. Temporal Aggregation with `aggregate_time`

The [`aggregate_time`](@ref) function allows you to resample data to different temporal resolutions.

### Monthly Aggregation Example

The ESDC has 8-day temporal resolution. Let's aggregate to monthly means:

```@example basic_ops
lai_month = aggregate_time(
    esdc[Variable = At("leaf_area_index")];
    time_axis = :time,
    new_resolution = "month",
    new_time_step = 1,
    fun = "mean",
    p = nothing,
    skipMissing = true,
    skipnan = true,
    showprog = true,
    max_cache = "1GB"
)
lai_month
```

### Checking the Result

Let's verify the new time axis:

```@example basic_ops
lookup(lai_month, :Ti)
```

!!! success "Result"
    The original 8-day temporal resolution has been aggregated to monthly values. The time axis now contains approximately 480 values (one per month over the full time range) instead of the original ~1800+ 8-day values.

---

## Quick Reference

| Function | Purpose | Key Parameters |
|:---------|:--------|:---------------|
| [`plot_time`](@ref) | Time series plot | `fun`, `var`, `time_axis` |
| [`plot_space`](@ref) | Spatial map | `fun`, `var`, `time_axis` |
| [`aggregate_time`](@ref) | Temporal aggregation | `new_resolution`, `fun` |

### Common Parameters

```julia
fun = "mean"           # Aggregation statistic
var = "temperature"    # Variable name (or nothing for all)
time_axis = :time      # Name of time dimension
showprog = true        # Show progress bar
max_cache = "1GB"      # Memory limit for caching
```

---

## Next Steps

Now that you've mastered the basics, explore more advanced analysis:

- 📖 **[Space-for-Time Method](@ref space4time)** — Learn how to estimate land cover change impacts
- 📚 **[API Reference](@ref api)** — Complete documentation of all functions

---

```@raw html
<div style="background: #f0f7ff; border-left: 4px solid #3b82f6; padding: 1em; margin-top: 2em;">
  <strong>💡 Tip:</strong> All examples in this tutorial use lazy evaluation. Data is only loaded when needed, making it efficient to work with large datasets.
</div>
```