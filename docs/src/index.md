<h1 align="center">YAXArraysToolbox.jl</h1>

<p align="center">
  <strong>High-level spatio-temporal analysis tools for Earth System data cubes</strong>
</p>

<p align="center">
  <em>Because laziness is not only good when reading big data — it should apply to analysis too!</em>
</p>

<p align="center">
  <a href="https://dpabon.github.io/YAXArraysToolbox.jl/stable/"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Stable Documentation"></a>
  <a href="https://dpabon.github.io/YAXArraysToolbox.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev Documentation"></a>
  <a href="https://github.com/dpabon/YAXArraysToolbox.jl/actions/workflows/CI.yml?query=branch%3Amain"><img src="https://github.com/dpabon/YAXArraysToolbox.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="Build Status"></a>
  <a href="https://codecov.io/gh/dpabon/YAXArraysToolbox.jl"><img src="https://codecov.io/gh/dpabon/YAXArraysToolbox.jl/branch/main/graph/badge.svg" alt="Coverage"></a>
  <a href="https://zenodo.org/badge/latestdoi/617361484"><img src="https://zenodo.org/badge/617361484.svg" alt="DOI"></a>
</p>

---

## ✨ Features at a Glance

YAXArraysToolbox extends [YAXArrays.jl](https://github.com/JuliaDataCubes/YAXArrays.jl) with high-level functions for climate and Earth observation data analysis:

| Feature | Description |
|:--------|:------------|
| 📈 **Time Series Plotting** | Visualize temporal evolution with automatic spatial aggregation |
| 🗺️ **Spatial Mapping** | Create publication-ready maps with temporal aggregation |
| ⏱️ **Temporal Aggregation** | Resample data to monthly, yearly, or custom periods |
| 🔄 **Space-for-Time Analysis** | Estimate land cover change impacts on climate variables |
| 😷 **Flexible Masking** | Apply spatial, temporal, and altitude-based masks |

## 📦 Installation

```julia
using Pkg
Pkg.add(url="https://github.com/dpabon/YAXArraysToolbox.jl")
```

Or using the package manager:

```julia
julia> ]
pkg> add https://github.com/dpabon/YAXArraysToolbox.jl
```

## 🚀 Quick Start

```julia
using YAXArraysToolbox
using YAXArrays
using CairoMakie
using Dates

# Load Earth System Data Cube
esdc = Cube(open_dataset(
    "https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-184x90x90-2.1.1.zarr"
))

# Select a region and variable
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

## 📚 Tutorials

Get hands-on experience with our step-by-step tutorials:

```@raw html
<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1em; margin: 2em 0;">
  <div style="border: 1px solid #ddd; border-radius: 8px; padding: 1.5em; background: #f9f9f9;">
    <h3 style="margin-top: 0;">📊 Basic Operations</h3>
    <p>Learn the core functions: <code>plot_time</code>, <code>plot_space</code>, and <code>aggregate_time</code>.</p>
    <a href="tutorials/basic_operations/">Start Tutorial →</a>
  </div>
  <div style="border: 1px solid #ddd; border-radius: 8px; padding: 1.5em; background: #f9f9f9;">
    <h3 style="margin-top: 0;">🌍 Space-for-Time Method</h3>
    <p>Understand and apply the space4time methodology for land cover change analysis.</p>
    <a href="tutorials/space4time_proof_of_concept/">Start Tutorial →</a>
  </div>
</div>
```

### Tutorial Overview

1. **[Basic Operations](@ref basic_operations)** — Essential functions for data visualization and aggregation
2. **[Space-for-Time Method](@ref space4time)** — Advanced methodology for estimating land cover change impacts


## 🏗️ Package Architecture

```
YAXArraysToolbox
├── Basic Operations
│   ├── plot_time      → Time series visualization
│   ├── plot_space     → Spatial mapping
│   └── aggregate_time → Temporal resampling
├── Masking
│   ├── masking_time     → Filter by time period
│   ├── masking_space    → Spatial filtering
│   └── masking_altitude → Elevation-based filtering
└── Spatio-Temporal Analysis
    └── space4time_proc  → Land cover change impact analysis
```
## 🔧 Main Functions

### Basic Operations

```julia
# Time series visualization with spatial aggregation
plot_time(cube; fun="mean", var="temperature")

# Spatial mapping with temporal aggregation
plot_space(cube; fun="median", var="lai")

# Temporal aggregation/resampling
aggregate_time(cube; new_resolution="month", fun="mean")
```

### Masking

```julia
# Mask by spatial extent
masking_space(cube, mask_cube; threshold=0.5)

# Mask by time period
masking_time(cube; start_date=Date(2010), end_date=Date(2015))
```

### Spatio-Temporal Analysis

```julia
# Space-for-time analysis for land cover change impacts
space4time_proc(climate_cube, landcover_cube, class_list)
```

## 📖 Dependencies

YAXArraysToolbox builds on these excellent Julia packages:

| Package | Purpose |
|:--------|:--------|
| [YAXArrays.jl](https://github.com/JuliaDataCubes/YAXArrays.jl) | Data cube handling and lazy operations |
| [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) | High-quality plotting |
| [GeoMakie.jl](https://github.com/MakieOrg/GeoMakie.jl) | Geographic projections and mapping |


## 🤝 Contributing

Contributions are welcome! Here's how you can help:

1. **Report bugs** — Open an [issue](https://github.com/dpabon/YAXArraysToolbox.jl/issues) describing the problem
2. **Suggest features** — Share your ideas in the issues section
3. **Submit PRs** — Fork the repo and submit pull requests

Please check existing [issues](https://github.com/dpabon/YAXArraysToolbox.jl/issues) before creating new ones.

## 📄 Citation

If you use YAXArraysToolbox in your research, please cite:

[![DOI](https://zenodo.org/badge/617361484.svg)](https://zenodo.org/badge/latestdoi/617361484)

## 📜 License
[MIT License](https://github.com/dpabon/YAXArraysToolbox.jl/blob/main/LICENSE)


## 🙏 Acknowledgements

This project was funded by:

<p align="center">
  <a href="https://earthmonitor.org/">
    <img src="https://earthmonitor.org/wp-content/uploads/2022/04/OEM_Logo_Horizontal_Dark_Transparent_Background_205x38.png" alt="Open Earth Monitor" width="200">
  </a>
  &nbsp;&nbsp;&nbsp;&nbsp;
  <a href="https://www.nfdi4earth.de/">
    <img src="https://www.nfdi4earth.de/templates/nfdi4earth/images/NFDI4Earth_logo.png" alt="NFDI4Earth" width="200">
  </a>
</p>

This project has received funding from the [Open-Earth-Monitor Cyberinfrastructure](https://earthmonitor.org/) project that is part of European Union's Horizon Europe research and innovation programme under grant [101059548](https://cordis.europa.eu/project/id/101059548).