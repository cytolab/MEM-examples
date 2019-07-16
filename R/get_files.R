get.files <- function(exp_data) {
    # Determine extension type and read in multiple files
    ext = file_ext(exp_data[[1]])

    # For FCS files
    if (ext == "fcs") {
        fileList = list.files(pattern = "*.fcs")
        allFiles = lapply(fileList, read.FCS)
        allMatrices = lapply(allFiles, exprs)
    }

    # For TXT files
    if (ext == "txt") {
        fileList = list.files(pattern = "*.txt")
        allMatrices = lapply(fileList, read.table, header = TRUE)

    }
    # For CSV Files

    if (ext == "csv") {
        fileList = list.files(pattern = "*.csv")
        allMatrices = lapply(fileList, read.table)
    }

    # Return ordered list of files for user reference
    dir.create("./output files/", showWarnings = FALSE)
    write.table(exp_data,
                paste(
                    "./output files/",
                    strftime(Sys.time(), "%Y-%m-%d_%H%M%S"),
                    " File order.txt"
                ),
                col.names = FALSE)

    names(allMatrices) = c(1:length(allMatrices))
    return(allMatrices)

}
