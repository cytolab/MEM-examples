MEM_RMSD <-
    function(MEM_vals,
             format = NULL,
             newWindow.heatmaps = FALSE,
             output.matrix = FALSE) {
        ## Use RMSD to compare MEM scores. Calculates a similarity measure by computing the percent maximum RMSD between populations by comparing their MEM scores. Returns the similarity matrix (n x n) matrix for n populations.

        if (isTRUE(format == "pop files")) {
            infiles <- dir(MEM_vals)
            sample_names = as.character(c(1:length(infiles)))
            RMSD = matrix(nrow = length(infiles), ncol = length(infiles))
            similarity = matrix(nrow = length(infiles), ncol = length(infiles))

            for (m in 1:length(infiles)) {
                for (n in 1:length(infiles)) {
                    popA1 = read.table(paste(MEM_vals, "/", infiles[[m]], sep = ""))
                    if (dim(popA1)[2] == 2) {
                        popA1 = read.table(paste(MEM_vals, "/", infiles[[m]], sep = ""),
                                           row.names = 1)
                    }
                    popB1 = read.table(paste(MEM_vals, "/", infiles[[n]], sep =
                                                 ""))
                    if (dim(popB1)[2] == 2) {
                        popB1 = read.table(paste(MEM_vals, "/", infiles[[n]], sep = ""),
                                           row.names = 1)
                    }
                    common_markers = Reduce(intersect, list(row.names(popA1), row.names(popB1)))
                    popA2 = popA1[c(common_markers), ]
                    popB2 = popB1[c(common_markers), ]
                    popA = as.vector(t(popA2))
                    popB = as.vector(t(popB2))

                    sum_squares = 0
                    for (i in 1:length(common_markers)) {
                        sum_squares = (popA[i] - popB[i]) ^ 2 + sum_squares
                        #       print(sum_squares)
                    }
                    RMSD[m, n] = sqrt(sum_squares / length(common_markers))
                    similarity[m, n] = 100 - (RMSD[m, n] / 20 * 100)
                }
            }
        } else{
            if (is(MEM_vals)[1] == "list") {
                MEM_scores = MEM_vals[[5]][[1]]
            } else{
                MEM_scores = MEM_vals
            }

            sample_names = as.character(row.names(MEM_scores))

            if (length(sample_names) == 0) {
                sample_names = as.character(1:nrow(MEM_scores))
            }

            RMSD = matrix(nrow = nrow(MEM_scores),
                          ncol = nrow(MEM_scores))
            similarity = matrix(nrow = nrow(MEM_scores),
                                ncol = nrow(MEM_scores))

            for (m in 1:nrow(MEM_scores)) {
                for (n in 1:nrow(MEM_scores)) {
                    popA = MEM_scores[m, ]
                    popB = MEM_scores[n, ]
                    common_markers = colnames(MEM_scores)

                    sum_squares = 0
                    for (i in 1:length(common_markers)) {
                        sum_squares = (popA[i] - popB[i]) ^ 2 + sum_squares
                    }
                    RMSD[m, n] = sqrt(sum_squares / length(common_markers))
                    similarity[m, n] = 100 - (RMSD[m, n] / 20 * 100)
                }
            }
        }
        row.names(similarity) = sample_names
        colnames(similarity) = sample_names

        heat_palette <-
            colorRampPalette(c("blue", "green", "yellow", "red"))
        scale_min = min(similarity)
        scale_max = max(similarity)
        scale_range = 100 - scale_min

        pairs.breaks <-
            c(
                seq(scale_min, (scale_min + scale_range / 4), by = 0.01),
                seq((scale_min + scale_range / 3.9),
                    (scale_min + scale_range / 2),
                    by = 0.01
                ),
                seq((scale_min + scale_range / 1.9),
                    (scale_max - scale_range / 3.9),
                    by = 0.01
                ),
                seq((scale_max - scale_range / 4), 100, by = 0.01)
            )

        if (newWindow.heatmaps == TRUE) {
            X11()
        }
        clustered <-
            heatmap.2(
                as.matrix(similarity),
                main = "MEM RMSD",
                dendrogram = "both",
                trace = "none",
                key = TRUE,
                col = heat_palette,
                breaks = pairs.breaks,
                margins = c(7, 15)
            )

        reorder_similarity = as.matrix(similarity[rev(clustered$rowInd), clustered$colInd])

        if (output.matrix == TRUE) {
            dir.create(file.path(getwd(), "output files"), showWarnings = FALSE)

            write.table(
                as.matrix(reorder_similarity),
                paste(
                    "./output files/",
                    strftime(Sys.time(), "%Y-%m-%d_%H%M%S"),
                    " MEM_RMSD1.txt",
                    sep = ""
                ),
                sep = "\t",
                row.names = TRUE
            )
        }

        return(similarity)
    }
