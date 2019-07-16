addclusterID <- function(pop_list, names) {
    # Get position in list
    clusterID = as.numeric(names)
    cluster = rep.int(clusterID, nrow(pop_list))

    newMatrix = cbind(pop_list, cluster)
    return(newMatrix)
}
