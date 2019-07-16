choose.ref <- function(exp_data,
                       pop_names,
                       num_pops,
                       num_markers) {
    # Specify alternative reference population(s)
    print("in choose.ref", suppressWarnings = FALSE)
    # Initialize matrices
    MAGref = matrix(nrow = num_pops, ncol = num_markers)
    IQRref = matrix(nrow = num_pops, ncol = num_markers)

    print(pop_names)
    user_spec_ref = unlist(readline("Enter population to use as reference for all subsets: "))
    spec_ref = as.numeric(unlist((strsplit(
        user_spec_ref, ","
    ))))

    if (length(spec_ref) > 1) {
        ref_pops = subset(exp_data, exp_data$cluster %in% c(spec_ref))

        MAGref = matrix(rep(abs(
            apply(ref_pops, 2, FUN = median, na.rm = TRUE)
        ), each = num_pops), num_pops)
        IQRref = matrix(rep(apply(
            ref_pops, 2, FUN = IQR, na.rm = TRUE
        ), each = num_pops), num_pops)

    } else{
        MAGref = matrix(rep(abs(
            apply(
                subset(exp_data, cluster == spec_ref),
                2,
                FUN = median,
                na.rm = TRUE
            )
        ), each = num_pops), num_pops)
        IQRref = matrix(rep(apply(
            subset(exp_data, cluster == spec_ref),
            2,
            FUN = IQR,
            na.rm = TRUE
        ), each = num_pops), num_pops)

    }

    return(list(MAGref, IQRref))
}
