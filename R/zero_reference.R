zero.ref <- function(exp_data, num_pops, num_markers) {
    MAGref = matrix(nrow = num_pops, ncol = num_markers)
    IQRref = matrix(nrow = num_pops, ncol = num_markers)

    MAGref = matrix(0, num_pops, num_markers)
    IQRs = colIQRs(as.matrix(exp_data))
    medIQRref = median(IQRs)
    IQRref = matrix(medIQRref, num_pops, num_markers)

    return(list(MAGref, IQRref))
}
