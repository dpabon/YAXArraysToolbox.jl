# YAXArraysToolbox.jl

*High-performance spatio-temporal analysis for Earth System data cubes*

## Overview

**YAXArraysToolbox.jl** is a Julia package that extends [YAXArrays.jl](https://github.com/JuliaDataCubes/YAXArrays.jl) with high-level functions for analyzing spatio-temporal data cubes.

!!! note "Why YAXArraysToolbox?"
    Because laziness is not only good when reading big data — it should also apply to common analysis tasks!

## Features

| Feature | Description |
|---------|-------------|
| 📈 **Time Series Plotting** | Visualize temporal evolution with spatial aggregation |
| 🗺️ **Spatial Mapping** | Create maps with temporal aggregation |
| ⏱️ **Time Aggregation** | Resample to monthly, yearly, or custom periods |
| 🔄 **Space-for-Time Analysis** | Estimate land cover change impacts |
| 😷 **Masking** | Flexible spatial, temporal, and altitude masking |

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/dpabon/YAXArraysToolbox.jl")
```

## Quick Start

```julia
using YAXArraysToolbox
using YAXArrays
using CairoMakie
using Dates

# Load Earth System Data Cube
esdc = Cube(open_dataset(
    "https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-184x90x90-2.1.1.zarr"
))

# Select region and variables
cube = esdc[
    lon = -10 .. 0,
    lat = 35 .. 45,
    time = Date(2010) .. Date(2012),
    Variable = At("leaf_area_index")
]

# Plot time series (spatial mean)
plot_time(cube; fun="mean")

# Create spatial map (temporal median)
plot_space(cube; fun="median")

# Aggregate to monthly resolution
monthly = aggregate_time(cube; new_resolution="month", fun="mean")
```

## Main Functions

### Basic Operations

```julia
# Time series visualization
plot_time(cube; fun="mean", var="temperature")

# Spatial mapping  
plot_space(cube; fun="median", var="lai")

# Temporal aggregation
aggregate_time(cube; new_resolution="month", fun="mean")
```

### Masking

```julia
# Spatial masking
masking_space(cube, mask_cube; threshold=0.5)

# Temporal masking
masking_time(cube; start_date=Date(2010), end_date=Date(2015))
```

### Spatio-Temporal Analysis

```julia
# Space-for-time analysis
space4time_proc(climate_cube, landcover_cube, class_list)
```

## Tutorials

Get started with our tutorials:

1. **[Basic Operations](@ref basic_operations)**: Learn `plot_time`, `plot_space`, and `aggregate_time`
2. **[Space-for-Time Method](@ref space4time)**: Understand and apply the space4time methodology

## Package Structure

```
YAXArraysToolbox
├── Basic Operations
│   ├── plot_time
│   ├── plot_space
│   └── aggregate_time
├── Masking
│   ├── masking_time
│   ├── masking_space
│   └── masking_altitude
└── Spatio-Temporal Analysis
    └── space4time_proc
```

## Dependencies

YAXArraysToolbox builds on these excellent packages:

- [YAXArrays.jl](https://github.com/JuliaDataCubes/YAXArrays.jl) - Data cube handling
- [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) - Plotting
- [GeoMakie.jl](https://github.com/MakieOrg/GeoMakie.jl) - Geographic plotting

## Citation

If you use YAXArraysToolbox in your research, please cite:

[![DOI](https://zenodo.org/badge/617361484.svg)](https://zenodo.org/badge/latestdoi/617361484)

## Acknowledgements

This project was funded by:

- [Open-Earth-Monitor](https://earthmonitor.org/)
- [NFDI4Earth](https://www.nfdi4earth.de/)

## License

[MIT License](https://github.com/dpabon/YAXArraysToolbox.jl/blob/main/LICENSE)
