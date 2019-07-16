rename.markers <- function(exp_data, marker_names) {
    print(colnames(exp_data[, c(1:(ncol(exp_data) - 1))]))
    user_input_names = readline(
        "Enter new marker names, in same order they appear above, separated by commas.\n No spaces allowed in name.\n"
    )
    new_marker_names = as.character(unlist(strsplit(user_input_names, ",")))

    if (length(new_marker_names) != (length(marker_names) - 1)) {
        warning(
            "Number of new marker names does not match number of markers.",
            call. = FALSE,
            immediate. = TRUE
        )
        new_marker_names = rename.markers(exp_data, marker_names)
    }

    # Add cluster column name
    new_marker_names = c(new_marker_names, "cluster")
    return(new_marker_names)
}
