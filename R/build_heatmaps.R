build.heatmaps <-
    function(exp_data,
             cluster.MEM = "both",
             cluster.medians = "none",
             cluster.IQRs = "none",
             display.thresh = 1,
             newWindow.heatmaps = FALSE,
             output.files = FALSE,
             labels = FALSE,
             only.MEMheatmap = FALSE,
             output.dir = "./output files/") {
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
                    seq(scale_min, scale_min / 3.3, by = 0.2),
                    seq(scale_min / 3.3, scale_min / 6.6, by = 0.2),
                    seq(scale_min / 6.6, 0, by = 0.2),
                    seq(0, scale_max / 6.6, by = 0.2),
                    seq(scale_max / 6.6, scale_max / 3.3, by = 0.2),
                    seq(scale_max / 3.3, scale_max, by = 0.2)
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
        
        
        if (newWindow.heatmaps) X11()
        
        # Print MEM heatmap according to cluster spec
        if (length(which(apply(heatmap_data, 1, function(row) any(row < 0)))) == 0){
            title_MEM = "   MEM Heatmap*"
        }else{
            title_MEM = "   MEM Heatmap"}
        
        
        args_heatmap_MEM = list(
            heatmap_data,
            main = title_MEM,
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
            col = heat_palette_MEM,
            key = TRUE,
            key.title = NA, key.xlab = NA, key.ylab = NA,
            lhei = c(0.5, 1.5),
            lwid = c(0.5, 2.2),
            trace = "none"
        )
        
        if (labels == TRUE) {
            args_heatmap_MEM = c(
                args_heatmap_MEM,
                list(
                    labRow = new_rownames,
                    margins = c(5, 15)
                    # margins = c(5, 18),
                )
            )
        }else{
            args_heatmap_MEM = c(
                args_heatmap_MEM, 
                list(
                    margins = c(5, 10)
                )
            )
        }
        table <-
            do.call(heatmap.2, args_heatmap_MEM)
        
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
            if (newWindow.heatmaps) X11()
            
            args_heatmap_MEDs = list(
                main = "Median Heatmap",
                breaks = pairs.breaks_med,
                revC = FALSE,
                symm = FALSE,
                symkey = FALSE,
                symbreaks = FALSE,
                scale = "none",
                cexRow = 0.8,  #1,
                cexCol = 0.8,
                col = heat_palette_med,
                key = TRUE,
                key.title = NA, key.xlab = NA, key.ylab = NA,
                lhei = c(0.5, 1.5),
                lwid = c(0.5, 2.2),
                trace = "none"
            )
            
            if (cluster.medians != "none") {
                args_heatmap_MEDs = c(
                    list(as.matrix(exp_data[[1]][[1]])),
                    args_heatmap_MEDs,
                    list(dendrogram = dendro_var_med,
                         Rowv = Rowv_var_med,
                         Colv = Colv_var_med,
                         margins = c(5, 10))
                )
                medians_table <- do.call(heatmap.2, args_heatmap_MEDs)
                reorder_medians = as.matrix(exp_data[[1]][[1]])[rev(medians_table$rowInd), medians_table$colInd]
            }else{
                args_heatmap_MEDs = c(
                    list(reorder_medians),
                    args_heatmap_MEDs,
                    list(dendrogram = "none",
                         Rowv = FALSE,
                         Colv = FALSE,
                         margins = c(5, 15))
                )
                table2 <- do.call(heatmap.2, args_heatmap_MEDs)
            }
            
            # Print IQR heatmap
            if (newWindow.heatmaps) X11()
            
            args_heatmap_IQRs = list(
                main = "IQR Heatmap",
                breaks = pairs.breaks_IQR,
                revC = FALSE,
                symm = FALSE,
                symkey = FALSE,
                symbreaks = FALSE,
                scale = "none",
                cexRow = 0.8,  #1,
                cexCol = 0.8,
                col = heat_palette_IQR,
                key = TRUE,
                key.title = NA, key.xlab = NA, key.ylab = NA,
                lhei = c(0.5, 1.5),
                lwid = c(0.5, 2.2),
                trace = "none"
            )
            
            if (cluster.IQRs != "none") {
                args_heatmap_IQRs = c(
                    list(as.matrix(exp_data[[3]][[1]])),
                    args_heatmap_IQRs,
                    list(dendrogram = dendro_var_IQR,
                         Rowv = Rowv_var_IQR,
                         Colv = Colv_var_IQR,
                         margins = c(5, 10))
                )
                IQR_table <- do.call(heatmap.2, args_heatmap_IQRs)
                reorder_IQR = as.matrix(exp_data[[3]][[1]])[rev(IQR_table$rowInd), IQR_table$colInd]
            }else{
                args_heatmap_IQRs = c(
                    list(reorder_IQR),
                    args_heatmap_IQRs,
                    list(dendrogram = "none",
                         Rowv = FALSE,
                         Colv = FALSE,
                         margins = c(5, 15))
                )
                table3 <- do.call(heatmap.2, args_heatmap_IQRs)
            }
        }
        
        dir.create(file.path(getwd(), output.dir), showWarnings = FALSE)
        
        if (output.files == TRUE) {
            
            # file attributes
            file_prefix = file.path(output.dir, strftime(Sys.time(), "%Y-%m-%d_%H%M%S"))
            
            # ----- pdf output ----------------------------
            
            cairo_pdf(paste0(file_prefix, " MEM heatmap.pdf"), width = 15, onefile = TRUE)
            # pdf(paste0(file_prefix, " MEM heatmap.pdf"), width = 15)
            
            # MEM heatmap plot
            do.call(heatmap.2, args_heatmap_MEM)
            
            if (only.MEMheatmap == FALSE) {
                
                # Plot median heatmap
                do.call(heatmap.2, args_heatmap_MEDs)
                
                # Plot IQR heatmap
                do.call(heatmap.2, args_heatmap_IQRs)
                
            }
            
            dev.off()
            
            # ----- end of pdf output ---------------------
            
            # Write data to text files
            write.table(
                as.matrix(clustered_matrix),
                paste0(file_prefix, " MEM matrix.txt"),
                sep = "\t",
                row.names = TRUE
            )
            write.table(
                as.matrix(reorder_medians),
                paste0(file_prefix, " Medians matrix.txt"),
                sep = "\t",
                row.names = TRUE
            )
            write.table(
                as.matrix(reorder_IQR),
                paste0(file_prefix, " IQRs matrix.txt"),
                sep = "\t",
                row.names = TRUE
            )
            
            if (((exp_data[[6]])[[1]]) == 0){
                write.table(
                    enrichment_score_ordered,
                    paste0(file_prefix, " enrichment score-rownames.txt"),
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
                    paste0(file_prefix, " enrichment score-rownames.txt"),
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
