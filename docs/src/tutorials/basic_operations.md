# [Basic Operations](@id basic_operations)

*Getting started with YAXArraysToolbox.jl*

**Author:** Daniel E. Pabon-Moreno

## Introduction

This tutorial demonstrates the core functionality of YAXArraysToolbox.jl for working with spatio-temporal data cubes. We'll cover:

- Loading data from the Earth System Data Cube
- Plotting time series with `plot_time`
- Creating spatial maps with `plot_space`
- Aggregating data over time with `aggregate_time`

## Setup

First, load the required packages:

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

## Loading Data

We'll use the [Earth System Data Cube (ESDC)](https://www.earthsystemdatalab.net/), a analysis-ready data cube containing many climate and Earth observation variables.

```@example basic_ops
esdc = open_dataset("https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-184x90x90-2.1.1.zarr")
esdc = Cube(esdc)
esdc
```

## Basic Functions

### Plot Time

The `plot_time` function creates time series plots by collapsing spatial dimensions using a specified statistic (mean, median, std, etc.).


:::
#### Subsetting the Data

Let's select a region (South America) and time period:

```@example basic_ops
cube_to_plot = esdc[
    lon = -86 .. -35,
    lat = -56 .. 14,
    time = Date(2010) ..  Date(2014),
    Variable = At("leaf_area_index", "sensible_heat"),
]
cube_to_plot
```

#### Plotting All Variables

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

#### Plotting a Single Variable

You can also plot a specific variable:

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

#### Available Statistics

The `fun` parameter supports multiple statistics:

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

### Plot Space

The `plot_space` function creates spatial maps by collapsing the time dimension.

```@example basic_ops
#| output: false
@doc plot_space
```

#### Single Variable Map

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

#### Multiple Variables

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

#### All Statistics for a Variable

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

### Aggregate by Time

The `aggregate_time` function allows you to aggregate data to different temporal resolutions (e.g., daily → monthly).

```@example basic_ops
#| output: false
@doc aggregate_time
```

#### Monthly Aggregation

Let's compute the monthly mean of the Leaf Area Index:

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

Check the new time axis:

```@example basic_ops
lookup(lai_month, :Ti)
```

The original 8-day temporal resolution has been aggregated to monthly values. The time axis now contains 480 values (one per month) instead of the original ~1800+ values.

## Summary

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `plot_time` | Time series plot | `fun`, `var`, `time_axis` |
| `plot_space` | Spatial map | `fun`, `var`, `time_axis` |
| `aggregate_time` | Temporal aggregation | `new_resolution`, `fun` |

## Next Steps

- Learn about the [Space4Time methodology](@ref space4time)
