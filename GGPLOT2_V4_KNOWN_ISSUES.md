# ggplot2 v4.0.0 Known Issues in Seurat

## Overview
ggplot2 version 4.0.0 (released September 2025) introduced breaking changes by migrating from S3 to S7 object system. This document tracks known compatibility issues and their status.

## Fixed Issues

### 1. DarkTheme() - Deprecated `size` parameter ✅
**Status**: FIXED  
**Files**: `R/visualization.R` (lines 6476-6491)  
**Changes**: 
- Replaced `element_rect(fill = 'black', size = 0)` with `linewidth = 0`
- Replaced `element_line(colour = 'white', size = 1)` with `linewidth = 1`
- Replaced `element_line(size = 0)` with `linewidth = 0`

### 2. LabelClusters() - S7 object color retrieval ✅
**Status**: FIXED  
**Issues**: #10156, #10176  
**Files**: `R/visualization.R` (lines 6134-6227)  
**Problem**: With S7 objects, `pb$data` returns numeric coordinates instead of valid colors  
**Solution**: Extract colors directly from plot scales using `pb$plot$scales$get_scales()`
**Changes**:
- Added color scale extraction from plot scales
- Created group_colors mapping with validation
- Fixed data.medians$color assignment to use explicit group_colors mapping

## Outstanding Issues (Require ggplot2 Package Fix)

### 3. VlnPlot - S4SXP Error 🔴
**Status**: NEEDS GGPLOT2 FIX  
**Issues**: #10101, #10188, #10160  
**Files**: Not Seurat-specific code  
**Error**: 
```
Error in deparse(substitute(e2, env = caller_env(2))) :
'S4SXP': should not happen - please report
```
**Root Cause**: Deep incompatibility between ggplot2 v4 S7 objects and R's `substitute()` function calls within ggplot2 internals  
**Workaround**: Downgrade to ggplot2 3.5.2
```r
remotes::install_version("ggplot2", version = "3.5.2", repos = "https://cran.r-project.org")
```
**Note**: Some packages (e.g., clusterProfiler >= 4.12.0) require ggplot2 >= 4.0.0, creating dependency conflicts

## Potential Issues (To Monitor)

### 4. ggplot_build()$plot$data Access Pattern ⚠️
**Status**: MONITORING  
**Files**: `R/visualization.R` (lines 5948, 7795)  
**Functions**: `GGpointToPlotly()`, `GGpointToPlotlyBuild()`  
**Concern**: Direct `$plot$data` access may behave differently with S7 objects  
**Impact**: Unknown - needs testing with real-world usage

### 5. plot$layers Access Pattern ⚠️
**Status**: MONITORING  
**Files**: `R/visualization.R` (multiple locations)  
**Concern**: `plot$layers[[i]]$data` and `plot$layers[[i]]$mapping` access may need S7 updates  
**Impact**: Unknown - backward compatibility should handle this, but monitor for issues

## Migration Notes

### What Changed in ggplot2 v4
- **S7 Object System**: ggplot objects now use S7 instead of S3
- **Deprecated Parameters**: 
  - `size` → `linewidth` in `element_rect()` and `element_line()`
  - Keep `size` in `element_text()` (NOT deprecated)
- **S7 Backward Compatibility**: `$data` and `$layers` access still works, but internal data structures changed
- **New Features**: 
  - Theme improvements (ink/paper/accent, palette.*, theme_sub_*)
  - `geom_label()` new aesthetics (linewidth, border.colour, text.colour)

### Testing Recommendations
1. Test all visualization functions with ggplot2 v4
2. Verify label.box functionality in DimPlot variants
3. Check VlnPlot alternatives if S4SXP persists
4. Monitor Plotly integration functions for S7 issues

### Dependencies
- **Minimum ggplot2 version**: 3.5.2 (for compatibility)
- **Target ggplot2 version**: 4.0.0+ (once S4SXP issue resolved upstream)

## References
- ggplot2 v4 blog post: https://www.tidyverse.org/blog/2025/09/ggplot2-4-0-0/
- Issue #10101: VlnPlot S4SXP error
- Issue #10156: DimPlot label.box color error
- Issue #10176: label.box returning numeric values
- Issue #10188: Temporary S4SXP solution thread

## Last Updated
2025-11-27
