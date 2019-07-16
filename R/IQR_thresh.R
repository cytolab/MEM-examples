IQR_thresh <- function(MAGpop,
                       MAGref,
                       IQRpop,
                       IQRref,
                       num_markers) {
    # Automatically calculate an IQR threshold
    #   Use universal IQR threshodling
    IQR_thresh_pop = matrix()
    IQR_thresh_ref = matrix()
    for (i in 1:(num_markers - 1)) {
        MAGpop_belowThresh <- MAGpop[, i] <= quantile(MAGpop[, i])[2]
        IQR_thresh_pop[i] <- min(IQRpop[, i][MAGpop_belowThresh])
        MAGref_belowThresh <- MAGref[, i] <= quantile(MAGref[, i])[2]
        IQR_thresh_ref[i] <- min(IQRref[, i][MAGref_belowThresh])
    }
    IQR.thresh = mean(c(IQR_thresh_pop, IQR_thresh_ref))

    return(IQR.thresh)
}
