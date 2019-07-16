format.data <- function(exp_data, file.is.clust, add.fileID) {
    # Get all files into list of matrices
    pop_list = get.files(exp_data)

    # Add cluster and file IDs according to user input
    if (file.is.clust == TRUE) {
        pop_list = mapply(
            addclusterID,
            pop_list,
            USE.NAMES = TRUE,
            names = names(pop_list),
            SIMPLIFY = FALSE
        )
    }

    if (add.fileID == TRUE) {
        pop_list = mapply(addfileID, pop_list1, names = names(pop_list))
    }

    # Put all cells into one matrix
    all_cells = do.call(rbind, pop_list)
    return(all_cells)
}
