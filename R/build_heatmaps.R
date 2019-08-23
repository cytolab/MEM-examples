build.heatmaps <-
    function(exp_data,
             cluster.MEM = "both",
             cluster.medians = "none",
             cluster.IQRs = "none",
             display.thresh = 1,
             newWindow.heatmaps = FALSE,
             output.files = FALSE,
             labels = FALSE,
             only.MEMheatmap = FALSE) {
        dendro_var_MEM = cluster.MEM
        dendro_var_med = cluster.medians
        dendro_var_IQR = cluster.IQRs
        ##Build MEM heatmap with natural language output
        #Get MEM values and the max and min MEM value
        heatmap_data = (exp_data[[5]])[[1]]
        
        if (length(which(apply(heatmap_data, 1, function(row) any(row < 0)))) == 0){
            scale_max = 10
            heat_palette_MEM <-
                colorRampPalette(c("black", "yellow", "#FAF7C9"))
            
            pairs.breaks_MEM <-
                c(seq(0, scale_max / 6.6, by = 0.1),
                    seq(scale_max / 6.6, scale_max / 3.3, by = 0.1),
                    seq(scale_max / 3.3, scale_max, by = 0.1)
                )
            pairs.breaks_MEM <- unique(pairs.breaks_MEM)
        }else{
            scale_max = 10
            scale_min = -10
        #Initialize heatmap color palette
        heat_palette_MEM <-
            colorRampPalette(c("#A8C3F4", "#17499B", "black", "yellow", "#FAF7C9"))
        pairs.breaks_MEM <-
            c(
                seq(scale_min, scale_min / 3.3, by = 0.1),
                seq(scale_min / 3.3, scale_min / 6.6, by = 0.1),
                seq(scale_min / 6.6, 0, by = 0.1),
                seq(0, scale_max / 6.6, by = 0.1),
                seq(scale_max / 6.6, scale_max / 3.3, by = 0.1),
                seq(scale_max / 3.3, scale_max, by = 0.1)
            )
        pairs.breaks_MEM <- unique(pairs.breaks_MEM)}

        #Initialize heatmap for medians
        medians_exp_data = (exp_data[[1]])[[1]]
        scale_max = max(medians_exp_data)
        heat_palette_med <-
            colorRampPalette(c("black", "yellow", "#FAF7C9"))
        pairs.breaks_med <-
            c(seq(0, scale_max / 6.6, by = 0.1),
                seq(scale_max / 6.6, scale_max / 3.3, by = 0.1),
                seq(scale_max / 3.3, scale_max, by = 0.1)
            )
        pairs.breaks_med <- unique(pairs.breaks_med)

        #Initialize heatmap for IQR
        medians_exp_data = (exp_data[[1]])[[1]]
        if (max(medians_exp_data) <= 2.8) {
            scale_max = 2.81
        }else{
            scale_max = max(medians_exp_data)
        }
        heat_palette_IQR <-
            colorRampPalette(c("black", "yellow", "#FAF7C9"))
        pairs.breaks_IQR <-
            c(seq(0.5, (scale_max + 0.5) / 6.6, by = 0.1),
              seq((scale_max + 0.5) / 6.6, (scale_max + 0.5) / 3.3, by = 0.1),
              seq((scale_max + 0.5) / 3.3, scale_max, by = 0.1))
        pairs.breaks_IQR <- unique(pairs.breaks_IQR)

        #Round MEM enrichment vals
        MEM_vals_scale = as.matrix(round(heatmap_data, 0))
        #Initialize variables
        #Call create.labels
        new_rownames = create.labels(MEM_vals_scale, display.thresh, heatmap_data)
        new_rownames_txt = create.labels.txt(MEM_vals_scale, display.thresh, heatmap_data)
        #Specify heatmap clustering parameters
        if (dendro_var_MEM == "both") {
            Colv_var_MEM = TRUE
            Rowv_var_MEM = TRUE
        }
        if (dendro_var_MEM == "row") {
            Colv_var_MEM = FALSE
            Rowv_var_MEM = TRUE
        }
        if (dendro_var_MEM == "col") {
            Colv_var_MEM = TRUE
            Rowv_var_MEM = FALSE
        }
        if (dendro_var_MEM == "none") {
            Colv_var_MEM = FALSE
            Rowv_var_MEM = FALSE
        }
        if (dendro_var_med == "both") {
            Colv_var_med = TRUE
            Rowv_var_med = TRUE
        }
        if (dendro_var_med == "row") {
            Colv_var_med = FALSE
            Rowv_var_med = TRUE
        }
        if (dendro_var_med == "col") {
            Colv_var_med = TRUE
            Rowv_var_med = FALSE
        }
        if (dendro_var_med == "none") {
            Colv_var_med = FALSE
            Rowv_var_med = FALSE
        }
        if (dendro_var_IQR == "both") {
            Colv_var_IQR = TRUE
            Rowv_var_IQR = TRUE
        }
        if (dendro_var_IQR == "row") {
            Colv_var_IQR = FALSE
            Rowv_var_IQR = TRUE
        }
        if (dendro_var_IQR == "col") {
            Colv_var_IQR = TRUE
            Rowv_var_IQR = FALSE
        }
        if (dendro_var_IQR == "none") {
            Colv_var_IQR = FALSE
            Rowv_var_IQR = FALSE
        }


        if (newWindow.heatmaps == TRUE) {
            X11()
        }

        # Print MEM heatmap according to cluster spec

        if (labels == TRUE) {
            table <-
                heatmap.2(
                    heatmap_data,
                    main = "   MEM Heatmap",
                    dendrogram = dendro_var_MEM,
                    Rowv = Rowv_var_MEM,
                    Colv = Colv_var_MEM,
                    breaks = pairs.breaks_MEM,
                    revC = FALSE,
                    symm = FALSE,
                    symkey = FALSE,
                    symbreaks = FALSE,
                    scale = "none",
                    cexRow = 0.8,
                    cexCol = 0.8,
                    key = TRUE,
                    col = heat_palette_MEM,
                    labRow = new_rownames,
                    margins = c(5, 18),
                    trace = "none",
                    lhei = c(0.7, 1.5),
                    lwid = c(0.6, 2.2)
                )

        }else{
            table <-
                heatmap.2(
                    heatmap_data,
                    main = "MEM Heatmap",
                    dendrogram = dendro_var_MEM,
                    Rowv = Rowv_var_MEM,
                    Colv = Colv_var_MEM,
                    breaks = pairs.breaks_MEM,
                    revC = FALSE,
                    symm = FALSE,
                    symkey = FALSE,
                    symbreaks = FALSE,
                    scale = "none",
                    cexRow = 0.8,
                    cexCol = 0.8,
                    key = TRUE,
                    col = heat_palette_MEM,
                    margins = c(5, 10),
                    trace = "none"
                )
        }

        clustered_matrix = heatmap_data[rev(table$rowInd), table$colInd]
        matrix.test= as.matrix(new_rownames)
        matrix.test_txt = as.matrix(new_rownames_txt)
        
        if (length(which(apply(heatmap_data, 1, function(row) any(row < 0)))) >0){
            enrichment_score_ordered_txt = matrix.test_txt[rev(table$rowInd), ] 
        }else{
            enrichment_score_ordered_txt = matrix.test_txt[rev(table$rowInd), ]}
        
        if (length(which(apply(heatmap_data, 1, function(row) any(row < 0)))) >0){
          enrichment_score_ordered = matrix.test[rev(table$rowInd), ] 
        }else{
          enrichment_score_ordered = matrix.test[rev(table$rowInd), ]}
        

        # Print medians heatmap according to cluster spec
        reorder_medians = as.matrix(exp_data[[1]][[1]])[rev(table$rowInd), table$colInd]
        # Print IQR heatmap according to cluster spec
        reorder_IQR = as.matrix(exp_data[[3]][[1]])[rev(table$rowInd), table$colInd]

        if (only.MEMheatmap == FALSE) {
            # Print median heatmap
            if (newWindow.heatmaps == TRUE) {
                X11()
            }
            if (cluster.medians != "none") {
                medians_table <-
                    heatmap.2(
                        as.matrix(exp_data[[1]][[1]]),
                        main = "Median Heatmap",
                        dendrogram = dendro_var_med,
                        Rowv = Rowv_var_med,
                        Colv = Colv_var_med,
                        breaks = pairs.breaks_med,
                        revC = FALSE,
                        symm = FALSE,
                        symkey = FALSE,
                        symbreaks = FALSE,
                        scale = "none",
                        cexRow = 1,
                        cexCol = 0.8,
                        key = TRUE,
                        col = heat_palette_med,
                        margins = c(5, 10),
                        trace = "none"
                    )
                reorder_medians = as.matrix(exp_data[[1]][[1]])[rev(medians_table$rowInd), medians_table$colInd]
            }else{
                table2 <-
                    heatmap.2(
                        reorder_medians,
                        main = "Median Heatmap",
                        dendrogram = "none",
                        Rowv = FALSE,
                        Colv = FALSE,
                        breaks = pairs.breaks_med,
                        revC = FALSE,
                        symm = FALSE,
                        symkey = FALSE,
                        symbreaks = FALSE,
                        scale = "none",
                        cexRow = 1,
                        cexCol = 0.8,
                        key = TRUE,
                        col = heat_palette_med,
                        margins = c(5, 15),
                        trace = "none"
                    )
            }

            # Print IQR heatmap
            if (newWindow.heatmaps == TRUE) {
                X11()
            }
            if (cluster.IQRs != "none") {
                IQR_table <-
                    heatmap.2(
                        as.matrix(exp_data[[3]][[1]]),
                        main = "IQR Heatmap",
                        dendrogram = dendro_var_IQR,
                        Rowv = Rowv_var_IQR,
                        Colv = Colv_var_IQR,
                        breaks = pairs.breaks_IQR,
                        revC = FALSE,
                        symm = FALSE,
                        symkey = FALSE,
                        symbreaks = FALSE,
                        scale = "none",
                        cexRow = 1,
                        cexCol = 0.8,
                        key = TRUE,
                        col = heat_palette_med,
                        margins = c(5, 10),
                        trace = "none"
                    )
                reorder_IQR = as.matrix(exp_data[[3]][[1]])[rev(IQR_table$rowInd), IQR_table$colInd]
            }else{
                table3 <-
                    heatmap.2(
                        reorder_IQR,
                        main = "IQR Heatmap",
                        dendrogram = "none",
                        Rowv = FALSE,
                        Colv = FALSE,
                        breaks = pairs.breaks_IQR,
                        revC = FALSE,
                        symm = FALSE,
                        symkey = FALSE,
                        symbreaks = FALSE,
                        scale = "none",
                        cexRow = 1,
                        cexCol = 0.8,
                        key = TRUE,
                        col = heat_palette_IQR,
                        margins = c(5, 15),
                        trace = "none"
                    )
            }
        }
        dir.create(file.path(getwd(), "output files"), showWarnings = FALSE)
        if (output.files == TRUE) {
            # #         # Generate pdfs from MEM heatmap plot
            pdf(paste(
                "./output files/",
                strftime(Sys.time(), "%Y-%m-%d_%H%M%S"),
                " MEM heatmap.pdf"
            ),
            width = 20)
            heatmap.2(
                heatmap_data,
                main = "MEM heatmap",
                dendrogram = dendro_var_MEM,
                Rowv = Rowv_var_MEM,
                Colv = Colv_var_MEM,
                breaks = pairs.breaks_MEM,
                revC = FALSE,
                symm = FALSE,
                symkey = FALSE,
                symbreaks = FALSE,
                scale = "none",
                cexRow = 1,
                cexCol = 0.8,
                key = TRUE,
                col = heat_palette_MEM,
                trace = "none"
                ,
                lhei = c(0.7, 1.5),
                lwid = c(0.6, 2.2)
            )
            dev.off()

            # Write data to text files
            write.table(
                as.matrix(clustered_matrix),
                paste(
                    "./output files/",
                    strftime(Sys.time(), "%Y-%m-%d_%H%M%S"),
                    " MEM matrix.txt",
                    sep = ""
                ),
                sep = "\t",
                row.names = TRUE
            )
            write.table(
                as.matrix(reorder_medians),
                paste(
                    "./output files/",
                    strftime(Sys.time(), "%Y-%m-%d_%H%M%S"),
                    " Medians matrix.txt",
                    sep = ""
                ),
                sep = "\t",
                row.names = TRUE
            )
            write.table(
                as.matrix(reorder_IQR),
                paste(
                    "./output files/",
                    strftime(Sys.time(), "%Y-%m-%d_%H%M%S"),
                    " IQRs matrix.txt",
                    sep = ""
                ),
                sep = "\t",
                row.names = TRUE
            )

            if (((exp_data[[6]])[[1]]) == 0){
                write.table(
                    enrichment_score_ordered,
                    paste(
                        "./output files/",
                        strftime(Sys.time(), "%Y-%m-%d_%H%M%S"),
                        " enrichment score-rownames.txt"
                    ),
                    sep = "\t"
                )
            }else{
                filenames <- unlist(exp_data[[6]])
                matrix.filenames = as.matrix(filenames)
                filenames_ordered = matrix.filenames[rev(table$rowInd), ]
                new_rownames_filenames <-
                    cbind(filenames_ordered, enrichment_score_ordered)
                colnames(new_rownames_filenames) <- c("File", "MEM label")
                write.table(
                    new_rownames_filenames,
                    paste(
                        "./output files/",
                        strftime(Sys.time(), "%Y-%m-%d_%H%M%S"),
                        " enrichment score-rownames.txt"
                    ),
                    sep = "\t"
                )
            }
        }
        if (((exp_data[[6]])[[1]]) == 0){
            cat(enrichment_score_ordered_txt, sep = "\n")
        }else{
            filenames <- unlist(exp_data[[6]])
            matrix.filenames = as.matrix(filenames)
            filenames_ordered = matrix.filenames[rev(table$rowInd), ]
            new_rownames_filenames <-
                cbind(filenames_ordered, enrichment_score_ordered_txt)
            colnames(new_rownames_filenames) <- c("File", "MEM label")
            print(new_rownames_filenames, sep = "\n")
        }
    }
