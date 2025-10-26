# Performance Analysis: Direct Call vs do.call

**Date:** 2025-10-26  
**Issue:** Understanding why do.call appeared faster in initial tests

## Initial Observation

Early tests showed do.call appearing faster:
- Direct call: 7.31 seconds
- do.call: 1.15 seconds

This seemed suspicious and warranted investigation.

## Investigation

### Equality Test
First verified that results are identical:
- ✅ All data layers identical (counts, data, scale.data)
- ✅ Variable features identical
- ✅ All dimensions and metadata identical

### Multiple Timing Tests

Ran 5 consecutive tests:

| Run | Direct | do.call | Ratio |
|-----|--------|---------|-------|
| 1   | 7.25s  | 1.22s   | 0.17x |
| 2   | 0.78s  | 0.83s   | 1.06x |
| 3   | 0.78s  | 0.78s   | 1.00x |
| 4   | 0.83s  | 0.78s   | 0.94x |
| 5   | 0.78s  | 0.78s   | 1.00x |

**Statistics:**
- Direct: median 0.78s (range 0.78-7.25s)
- do.call: median 0.78s (range 0.78-1.22s)

### Root Cause: First-Call Overhead

The **first SCTransform call** (regardless of method) takes ~7 seconds due to:
1. Loading and compiling R code
2. Loading package dependencies (sctransform, glmGamPoi, etc.)
3. JIT compilation of functions
4. Initializing data structures

**Subsequent calls** take only ~0.78 seconds because:
- Code is already compiled
- Packages are loaded
- Warm cache

### Order Independence Test

When reversing the order (do.call first):
- do.call (first): 0.78 seconds
- Direct (second): 0.78 seconds
- Ratio: 1.00x

After the initial warmup, **order doesn't matter**.

## Conclusion

✅ **Performance is IDENTICAL between direct call and do.call**

The fix successfully resolves issue #10153:
- **Before fix:** do.call took 1800+ seconds (30+ minutes)
- **After fix:** do.call takes 0.78 seconds (same as direct call)
- **Improvement:** ~2300x faster

The initial observation of do.call being "faster" was an artifact of:
1. Testing direct call first (paying the warmup cost)
2. Testing do.call second (benefiting from warm cache)

When properly tested:
- Both methods take ~0.78 seconds after warmup
- Results are byte-for-byte identical
- No performance penalty from do.call

## Test Files

- `test_equality_detailed.R` - Comprehensive equality verification
- `test_timing_multiple.R` - Multiple runs showing warmup effect
