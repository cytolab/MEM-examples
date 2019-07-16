create.labels.txt <-
    function(MEM_vals_scale,
             display.thresh,
             heatmap_data) {
        posRownamesMatrix = matrix()
        negRownamesMatrix = matrix()
        for (i in 1:nrow(heatmap_data)) {
            #Get markers above threshold
            posMarkers = MEM_vals_scale[i, ][MEM_vals_scale[i, ] >= display.thresh]
            negMarkers = MEM_vals_scale[i, ][MEM_vals_scale[i, ] < 0 &
                                                 abs(MEM_vals_scale[i, ]) >= display.thresh]
            #Order markers by abs(MEM)
            posMarkersOrdered = posMarkers[order(posMarkers, decreasing = TRUE)]
            negMarkersOrdered = negMarkers[order(abs(negMarkers), decreasing = TRUE)]

            #Get + or - sign for each corresponding marker
            posVector = rep("+", length(posMarkers))
            negVector = rep("-", length(negMarkers))

            if (length(posVector) == 0) {
                posVector = rep("None", 1)
            }

            if (length(negVector) == 0) {
                negVector = rep("None", 1)
            }

            #Put together character vectors for rownames
            posSignedVec = as.matrix(paste(posVector, abs(round(
                posMarkersOrdered, 2
            )), sep = ""))
            posMarkersRanked = as.matrix(names(posMarkersOrdered))

            posRownames = t(paste(posMarkersRanked, posSignedVec, sep = ""))
            posRownamesMatrix[i] = apply(posRownames, 1, paste, collapse = " ")

            negSignedVec = as.matrix(paste(negVector, abs(round(
                negMarkersOrdered, 2
            )), sep = ""))
            negMarkersRanked = as.matrix(names(negMarkersOrdered))

            negRownames = t(paste(negMarkersRanked, negSignedVec, sep = ""))
            negRownamesMatrix[i] = apply(negRownames, 1, paste, collapse = " ")
        }
        dot = rawToChar(as.raw(149))

        new_rownames = paste(
            row.names(MEM_vals_scale),
            ":",
            'UP',
            posRownamesMatrix,
            paste(dot),
            'DN',
            negRownamesMatrix,
            sep = " "
        )
        return(new_rownames)
    }
