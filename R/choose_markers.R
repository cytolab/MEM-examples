choose.markers <- function(exp_data) {
    print("Numbered column names, in order they appear in file: ")
    print(paste(c(1:(
        ncol(exp_data) - 1
    )), ": ", colnames(exp_data[, c(1:(ncol(exp_data) - 1))]), sep = ""))
    markers = readline("Enter column numbers to include (e.g. 1:5,6,8:10).\n")

    sep_vals = unlist(strsplit(markers, ","))
    list_vals = vector()
    for (i in 1:length(sep_vals)) {
        val = sep_vals[i]
        if (length(unlist(strsplit(val, ":"))) > 1) {
            new_val = as.numeric(unlist(strsplit(val, ":"))[1]):as.numeric(unlist(strsplit(val, ":"))[2])
        } else{
            new_val = as.numeric(sep_vals[i])
        }
        list_vals = c(list_vals, new_val)
    }

    markerList = c(list_vals, ncol(exp_data))

    return(markerList)
}
