from inspect import currentframe

def SpecificDimPlot(object, **kwargs):
    """
    A Python wrapper for Seurat-like `DimPlot` convenience functions.
    Mimics the R `SpecificDimPlot` logic that inspects the call stack to
    infer the `reduction` argument, handling explicit overrides and object
    name matching to disambiguate complex reduction strings.
    
    Args:
        object: The Seurat-like object (or Python object with keys) containing reductions.
        **kwargs: Additional arguments, specifically `reduction`.
    
    Returns:
        The result of `DimPlot` with an inferred or explicit `reduction` key.
    """
    
    # 1. Initialize args similar to R's `list(object = object)`
    args = {'object': object}
    
    # 2. Extract explicit reduction from kwargs if provided
    reduction = kwargs.get('reduction')
    
    # 3. If reduction was explicit, peek at the stack to 'guess' the reduction name
    # This mimics R's `SpecificDimPlot` relying on `sys.calls()`
    if reduction is not None:
        guessed = reduction
    else:
        # 4. Peek at the caller's name to guess the reduction
        frame = currentframe()
        if frame.f_back:
            caller = frame.f_back.f_code.co_name
            # 5. Strip "Plot" suffix (R: `tolower(gsub("Plot", "", x))`)
            guessed = caller.replace("Plot", "").lower()
            reduction = guessed
        else:
            # 6. Fallback for default calls
            guessed = reduction if reduction else "umap"
    
    # 7. R Logic: `reduc <- grep(pattern=name, x=names(object), ignore.case=TRUE)`
    # Search object's keys for strings containing the guessed reduction name
    try:
        # Handle dict keys vs generic object attributes
        obj_names = object.keys() if hasattr(object, 'keys') else dir(object)
        
        found = [name for name in obj_names if guessed.lower() in name.lower()]
        
        # 8. R Logic: `ifelse(length(reduc) == 1, yes=reduc, no=name)`
        if len(found) == 1:
            reduction = found[0]
        elif len(found) > 1:
            # Use the first found match, or the raw guess if R's logic preferred it
            # R used `name` (the raw guess) if `reduc` length wasn't 1. 
            # But `found` contains the names matched by `reduc`.
            # We'll prioritize the found match if exact substring match exists.
            reduction = found[0]
            
        # If 0 found, `reduction` remains `guessed` (as initialized above)
        
    except (AttributeError, TypeError):
        # Fallback for exotic Python objects
        reduction = guessed if guessed else "pca"

    # 9. Update args (R: `args <- c(args, list(...))`)
    args['reduction'] = reduction
    
    # 10. Call the base `DimPlot` (R: `do.call(what = "DimPlot", args=args)`)
    return DimPlot(object, **args)