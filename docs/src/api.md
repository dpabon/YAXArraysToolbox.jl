# [API Reference](@id api)

This page provides comprehensive documentation for all public functions in YAXArraysToolbox.jl.

---

## Basic Operations

These are the core functions for everyday data analysis tasks.

### Time Series Plotting

Visualize how variables change over time by aggregating spatial dimensions.

```@docs
plot_time
```

**Example:**

```julia
plot_time(
    cube;
    fun = "mean",           # Aggregation function
    var = "temperature",    # Variable to plot
    time_axis = :time,      # Name of time dimension
    resolution = (900, 600) # Figure size
)
```

---

### Spatial Mapping

Create maps by aggregating the temporal dimension.

```@docs
plot_space
```

**Example:**

```julia
plot_space(
    cube;
    fun = "median",         # Aggregation function
    var = "lai",            # Variable to plot
    time_axis = :time,      # Name of time dimension
    resolution = (900, 600) # Figure size
)
```

---

### Temporal Aggregation

Resample data to different temporal resolutions (e.g., 8-day to monthly).

```@docs
aggregate_time
```

**Example:**

```julia
# Aggregate from 8-day to monthly means
monthly_cube = aggregate_time(
    cube;
    new_resolution = "month",  # Target resolution
    fun = "mean",              # Aggregation function
    skipMissing = true         # Handle missing values
)
```

**Supported resolutions:** `"day"`, `"month"`, `"year"`

**Supported functions:** `"mean"`, `"median"`, `"std"`, `"var"`, `"sum"`, `"min"`, `"max"`, `"quant"`

---

## Masking Functions

Functions for filtering and subsetting data based on various criteria.

### Temporal Masking

Filter data by time period.

```@docs
masking_time
```

**Example:**

```julia
# Keep only data from 2010-2015
masked = masking_time(
    cube;
    start_date = Date(2010, 1, 1),
    end_date = Date(2015, 12, 31)
)
```

---

### Spatial Masking

Apply spatial masks based on another data cube.

```@docs
masking_space
```

**Example:**

```julia
# Mask using a land/water mask
masked = masking_space(
    data_cube,
    land_mask_cube;
    threshold = 0.5  # Minimum land fraction
)
```

---

### General Masking

Apply combined masks with multiple criteria.

```@docs
masking_proc
```

---

## Space-for-Time Analysis

Functions for analyzing land cover change impacts using spatial variability as a proxy for temporal change.

### Main Processing Function

```@docs
space4time_proc
```

**Example:**

```julia
results = space4time_proc(
    climate_cube,           # Climate variable (e.g., LST)
    landcover_cube,         # Land cover fractions
    altitude_cube;          # Altitude data (optional)
    classes_vec = ["forest", "grassland", "cropland"],
    winsize = 5,            # Moving window size
    showprog = true
)
```

---

### Spatial Processing Variant

```@docs
space4time_proc_space
```

---

### Legacy Function

```@docs
space4time_proc_old
```

!!! warning "Deprecated"
    This function is provided for backward compatibility. Use [`space4time_proc`](@ref) for new projects.

---

## Function Index

```@index
```

---

## Type Reference

All functions in YAXArraysToolbox work with `YAXArray` objects from [YAXArrays.jl](https://github.com/JuliaDataCubes/YAXArrays.jl).

### Common Parameters

| Parameter | Type | Description |
|:----------|:-----|:------------|
| `cube` | `YAXArray` | Input data cube |
| `fun` | `String` | Aggregation function: `"mean"`, `"median"`, `"std"`, `"var"`, `"sum"`, `"min"`, `"max"`, `"quant"` |
| `time_axis` | `Symbol` | Name of the time dimension (typically `:time` or `:Ti`) |
| `var` | `String` or `Nothing` | Variable name to process, or `nothing` for all |
| `showprog` | `Bool` | Show progress bar |
| `max_cache` | `String` | Maximum memory cache (e.g., `"1GB"`) |