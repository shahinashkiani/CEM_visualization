# CEM unfolding visualization
# /Users/Shaahin/Downloads/japan_with_names.csv
efficiency_barplot <- function(dataset,dmu_labels=TRUE,num_of_inputs,return_to_scale = "crs",orientation = "in"){
        
        if(dmu_labels){
                dmu_labels <-dataset[,1]
                inputs <- dataset[,2:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+2):ncol(dataset)] %>% as.matrix()
        } else {
                dmu_labels <- data.frame(dmu = 1:nrow(dataset))
                inputs <- dataset[,1:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+1):ncol(dataset)] %>% as.matrix()
        }
        
        colnames(dmu_labels) <- "DMU"
        
        eff_vec<- Benchmarking::dea(X = inputs, Y = outputs, RTS = return_to_scale , DUAL = TRUE, ORIENTATION = orientation ) 
        eff <- eff_vec$eff
        df <- data.frame(dmu_labels, crs_efficiency = eff ) %>% arrange(desc(crs_efficiency))
        p<-df %>%
                plot_ly(y= ~ reorder(DMU,crs_efficiency) , x= ~ crs_efficiency , type = "bar" , orientation = 'h') %>% 
                layout(title = paste0("Efficiency Barplot-RTS: ",return_to_scale,"-OR: ",orientation))
        
        return(p)
        
}

efficiency_comparison_plot <- function(dataset,dmu_labels=TRUE,num_of_inputs,orientation = "in", plot_type = "bar"){
        if(dmu_labels){
                dmu_labels <-dataset[,1]
                inputs <- dataset[,2:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+2):ncol(dataset)] %>% as.matrix()
        } else {
                dmu_labels <- data.frame(dmu = 1:nrow(dataset))
                inputs <- dataset[,1:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+1):ncol(dataset)] %>% as.matrix()
        }
        
        colnames(dmu_labels) <- "DMU"
        
        eff_vec_crs<- Benchmarking::dea(X = inputs, Y = outputs, RTS = "crs" , DUAL = TRUE, ORIENTATION = orientation ) 
        eff_vec_vrs<- Benchmarking::dea(X = inputs, Y = outputs, RTS = "vrs" , DUAL = TRUE, ORIENTATION = orientation ) 
        
        eff_crs <- eff_vec_crs$eff
        eff_vrs <- eff_vec_vrs$eff
        
        df <- data.frame(dmu_labels, crs = eff_crs , vrs = eff_vrs ) 
        
        if(plot_type == "bar"){
                p<-
                        df %>%
                        plot_ly(y= ~ reorder(DMU,crs) , x= ~ crs , type = "bar" , orientation = 'h' , name = "CRS") %>% 
                        add_trace(x = ~ vrs , name = "VRS") %>% 
                        layout(title = paste0("Efficiency CRS-VRS Barplot- OR: ",orientation))
     
        }else if(plot_type == "scatter"){
                p<-
                        df %>%
                        plot_ly(y= ~crs , x= ~ vrs , type = "scatter", mode = "markers", text = ~DMU ) %>% 
                        layout(title = paste0("Efficiency CRS-VRS Scatterplot- OR: ",orientation))
 
        }

        return(p)   
        
}


creff_plot <- function(dataset, num_of_inputs, dmu_labels = TRUE){
        # dataset <- jap_df 
        # number_of_inputs <- 3 
        # dmu_labels = TRUE 
        if(dmu_labels){
                dmu_labels <-dataset[,1]
                inputs <- dataset[,2:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+2):ncol(dataset)] %>% as.matrix()
                dataset <- dataset[,-1] %>% as.matrix()
        } else {
                dmu_labels <- data.frame(dmu = 1:nrow(dataset))
                inputs <- dataset[,1:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+1):ncol(dataset)] %>% as.matrix()
                dataset <- dataset %>% as.matrix()
        }
        
        colnames(dmu_labels) <- "DMU"
       
        cem <- CEM(dataset = dataset , number_of_inputs = num_of_inputs )
        avg_creff <- apply(cem, MARGIN = 2 , FUN = mean)
        
        cem_df <- data.frame(dmu_labels, cem %>% data.frame())
        colnames(cem_df) <- c("DMU", paste0("col_",dmu_labels$DMU))
        
        self_eff <- diag(cem)
        creff_df <- data.frame(dmu_labels, avg_creff , self_eff)
        
        
        
        p_barplot <-
                creff_df %>% 
                plot_ly(x = ~ avg_creff, y = ~reorder(DMU,avg_creff), type = "bar", orientation = "h") %>% 
                layout(title = "Average Cross-Efficiency")
        
        p_grouped_barplot <- 
                creff_df %>%
                plot_ly(x = ~ avg_creff, y = ~reorder(DMU,self_eff), type = "bar", orientation = "h" , name = "Cross-Eff") %>%
                add_trace(x = ~ self_eff , name = "Simple-Eff") %>% 
                layout(title = "Average Cross-Efficiency vs Simple Efficiency")
        
        p_scatterplot <- 
                creff_df %>%
                plot_ly(x = ~ self_eff, y = ~avg_creff, type = "scatter", text = ~DMU) %>%
                layout(title = "Avg-CrossEfficiency vs Simple Efficiency")
                
        # p_box <- 
        #         data.frame(creff_df , cem) %>%
        #         plot_ly(x = ~ DMU, y = ~ avg_creff, type = "scatter", text = ~DMU , type = "box") %>%
        #         layout(title = "Avg-CrossEfficiency vs Simple Efficiency")
                
        return(list(p_barplot,p_grouped_barplot,p_scatterplot))
        
}



sammon_mapping <- function(dataset , num_of_inputs, return_to_scale = "CRS", edge_treshold = 0.1, edge_transparency_factor = 4){
        set.seed(7)
        dataset = data.frame(dataset)

        Porembsky_data = scale(dataset)
        Porembsky_sammon = MASS::sammon(d = dist(Porembsky_data),k=2,niter = 10000 )

        points <- Porembsky_sammon$points

        porembski_crs <- crs_eff(dataset,num_of_inputs)
        porembski_vrs <- vrs_eff(dataset,num_of_inputs)


        Porembski_eff <- if(return_to_scale == "CRS")  porembski_crs$eff else porembski_vrs$eff
        Porembski_eff_lambda <-  if( return_to_scale == "CRS") porembski_crs$lambda else porembski_vrs$lambda


        Porembsky_df=data.frame(DMU=1:nrow(dataset),points,Efficiency = Porembski_eff)
        Porembsky_df$color = as.character(ifelse(Porembsky_df$Efficiency==1,"Efficient","Inefficient"))


        Porembski_links_df = data.frame("start.DMU"= NA,"end.DMU"=NA,"x1"=NA,"y1"=NA,"x2"=NA,"y2"=NA,"alpha"=NA)

        index = 0
        for (r in 1:nrow(dataset)) {

                for (c in 1:nrow(dataset)){
                        index = index + 1
                        Porembski_links_df[index,] = c(r,c,points[r,1],points[r,2],points[c,1],points[c,2],Porembski_eff_lambda[r,c])
                }
        }


        total_min <-  min(Porembsky_df$X1,Porembsky_df$X2)
        total_max <- max( Porembsky_df$X1,Porembsky_df$X2)




        g = ggplot2::ggplot() +
                ggplot2::geom_point(data = Porembsky_df , aes(x = X1 , y = X2,  color =color ),
                                    size = 1,
                                    alpha = 0.5) +
                ggplot2::scale_color_manual(name = paste("DMU",return_to_scale,sep = ":"),
                                            values = c("Efficient" = "gold", "Inefficient"= "skyblue")) +
                ggplot2::geom_text(data = Porembsky_df , aes(x = X1 , y = X2 , label = DMU), alpha = 0)



        for (index in 1:nrow(Porembski_links_df)) {
                if (Porembski_links_df[index,"alpha"] > edge_treshold) {
                        g = g  + ggplot2::geom_segment(data = Porembski_links_df[index,],
                                                       aes(x = x1 , y = y1 , xend = x2 , yend = y2),
                                                       alpha = round(Porembski_links_df[index,"alpha"],2)/edge_transparency_factor ,
                                                       color = "red" )
                }
        }
        g = g + ggplot2::coord_fixed(ratio = 1,  expand = TRUE)

        g<- g+ ggplot2::theme_linedraw()  + ggplot2::annotate("text", x = Inf, y = -Inf, label = "© DEA-Viz",
                                                              hjust=1.1, vjust=-1.1, col="blue", cex=6,
                                                              fontface = "bold", alpha = 0.4) +
                ggplot2::xlim(total_min , total_max) +
                ggplot2::ylim(total_min , total_max ) +
                ggplot2::coord_fixed(ratio = 1) +
                ggplot2::ggtitle("DEA Porembski Graph")


        g <- plotly::ggplotly(g)

        return(g)



}


pca_biplot <- function(dataset , number_of_inputs, return_to_scale = "CRS" , vector_size = 2 , text_size = 2 ){
        set.seed(7)
        dataset = data.frame(dataset)
        t <- scale(dataset)
        prcomp_result <- prcomp(t)

        coordinates <- prcomp_result$x[,1:2]
        crs_efficiency <- crs_eff(dataset,number_of_inputs)$eff
        vrs_efficiency <- vrs_eff(dataset,number_of_inputs)$eff

        supplement_data <- data.frame(coordinates,crs_efficiency, vrs_efficiency)

        x <- prcomp_result
        supp <- supplement_data
        data <- dataset

        data$shape <- if (return_to_scale == "CRS"){
                ifelse(supplement_data$crs_efficiency == 1 , 1, 19)
        } else{
                ifelse(supplement_data$vrs_efficiency == 1 , 1, 19)
        }

        z1 <- data.frame(DMU = 1:nrow(data), x$x[, 1:2], crs_color = NA , vrs_color = NA)
        z1$crs_color <- ifelse(test = supp$crs_efficiency ==1 , "Efficient", "Inefficient")

        z1$vrs_color <- ifelse(test = supp$vrs_efficiency ==1 , "Efficient", "Inefficient")


        z2 <- data.frame(Variables = (rownames(x$rotation)), x$rotation[, 1:2])


        total_min <-  min(z1$PC1,z1$PC2)
        total_max <- max( z1$PC1,z1$PC2)


        if (return_to_scale == "CRS") {
                g <-
                        ggplot2::ggplot(z1, aes(x = PC1,y = PC2)) +
                        ggplot2::geom_point(aes(color = z1$crs_color) ,size=1, alpha = 0.5 ) +
                        ggplot2::scale_colour_manual(name = "DMU: CRS Eff.", values =  c("Efficient"="gold" , "Inefficient"="skyblue")) +
                        ggplot2::geom_segment(data=z2, aes(PC1*vector_size, PC2*vector_size, xend=0, yend=0), col="red",alpha = 0.7 ) +
                        ggrepel::geom_text_repel(data=z2 , aes(x = PC1*vector_size,y =  PC2*vector_size, label = Variables ), col="red", alpha = 0.7, size = text_size )

        } else {
                ggplot(z1, aes(x = PC1,y = PC2)) +
                        ggplot2::geom_point(aes(color = z1$vrs_color) ,size=1, alpha = 0.5) +
                        ggplot2::scale_colour_manual(name = "DMU: VRS Eff.", values =  c("Efficient"="gold" , "Inefficient"="skyblue")) +
                        ggplot2::geom_segment(data=z2, aes(PC1*vector_size, PC2*vector_size, xend=0, yend=0), col="red",alpha = 0.7 ) +
                        ggrepel::geom_text_repel(data=z2 , aes(x = PC1*vector_size,y =  PC2*vector_size, label = Variables ), col="red", alpha = 0.7, size = text_size )

        }



        g<- g+ ggplot2::theme_linedraw()  + ggplot2::annotate("text", x = Inf, y = -Inf, label = "© DEA-Viz",
                                                              hjust=1.1, vjust=-1.1, col="blue", cex=6,
                                                              fontface = "bold", alpha = 0.4) +
                ggplot2::xlim(total_min, total_max) +
                ggplot2::ylim(total_min, total_max) +
                ggplot2::coord_fixed(ratio = 1) +
                ggplot2::ggtitle("DEA PCA Biplot")


        return(g)



}

SOM_plot <- function(dataset, number_of_inputs, horizontal_nodes = 8 , vertical_nodes = 8, return_to_scale = "CRS") {
        dataset = data.frame(dataset)
        set.seed(7)


        t<- scale(dataset)
        d<- dataset



        crs_efficiency <- crs_eff(d,number_of_inputs)$eff
        vrs_efficiency <- vrs_eff(d,number_of_inputs)$eff


        total_nodes = horizontal_nodes*vertical_nodes

        som_factor = kohonen::som(X = t,
                                  grid = kohonen::somgrid(xdim =horizontal_nodes, ydim = vertical_nodes, topo = "hexagonal" ),
                                  rlen = 10000 ,
                                  init = t[sample(x = nrow(t), size = total_nodes, replace = TRUE),],
                                  keep.data = TRUE)

        coolBlueHotRed <- function(n, alpha = 1) {grDevices::rainbow(n, end=4/6, alpha=alpha)[n:1]}
        jitter = matrix(rnorm(n= nrow(t)*2, mean = 0,sd = 0.1),ncol = 2)

        som_eff_vec = vector(length = total_nodes)


        if (return_to_scale == "CRS"){
                som_eff_vec_sapply = sapply(X = 1:total_nodes, function(x) mean(crs_efficiency[which(som_factor$unit.classif==x)] )  )
        } else {
                som_eff_vec_sapply = sapply(X = 1:total_nodes, function(x) mean(vrs_efficiency[which(som_factor$unit.classif==x)] )  )

        }


        g<- plot(som_factor, type = "property", property = som_eff_vec_sapply, palette.name = coolBlueHotRed , main = "SOM DEA" ) +
                text(x = som_factor$grid$pts[som_factor$unit.classif,1] + jitter[,1], y = som_factor$grid$pts[som_factor$unit.classif,2] + jitter[,2],labels = 1:nrow(t) , col = alpha("white",1)    )

        #g <- plotly::ggplotly(g)
        return(g)

}

SOM_properties_plot <- function(dataset, number_of_inputs, horizontal_nodes = 8 , vertical_nodes = 8, return_to_scale = "CRS" ){
        set.seed(7)


        dataset = data.frame(dataset)


        t<- scale(dataset)
        d<- dataset



        total_nodes = horizontal_nodes*vertical_nodes

        som_factor = kohonen::som(X = t,
                                  grid = kohonen::somgrid(xdim =horizontal_nodes, ydim = vertical_nodes, topo = "hexagonal" ),
                                  rlen = 10000 ,
                                  init = t[sample(x = nrow(t), size = total_nodes, replace = TRUE),],
                                  keep.data = TRUE)

        coolBlueHotRed <- function(n, alpha = 1) {grDevices::rainbow(n, end=4/6, alpha=alpha)[n:1]}
        jitter = matrix(rnorm(n= nrow(t)*2, mean = 0,sd = 0.1),ncol = 2)

        plot_titles <- colnames(d)

        num_of_plots<-ncol(t)
        plot_rows =as.integer(num_of_plots/3) + 1

        par(mfrow=c(plot_rows,3))
        for (index in 1:ncol(d)) {
                plot(som_factor, type = "property", property = som_factor$codes[[1]][,index], palette.name = coolBlueHotRed, main = plot_titles[index]  )
                #plots[[index]] <- p1
        }
}

CostaFrontier_plot <- function(dataset, number_of_inputs, model = "CRS-In", point_size = 1 , point_transparency = 0.5 ){
        data <- data.frame(dataset)


        if (model == "CRS-In") {

                # efficiency canlculations
                eff_crs_in <- crs_eff(dataset = data, num_of_inputs = number_of_inputs , orientation = "in")

                #Normalization Factor
                S_crs_inorient <- apply(X = eff_crs_in$ux,MARGIN = 1 , FUN = sum)

                # modified virtual factors
                modified_vi_inorient <- eff_crs_in$ux/ S_crs_inorient
                modified_vo_inorient <-  eff_crs_in$vy/ S_crs_inorient

                # virtual weight sums
                wsum_inputs = round(apply(X = modified_vi_inorient * data[,1:number_of_inputs], MARGIN = 1 , FUN = sum),3)
                wsum_outputs = round(apply(X = modified_vo_inorient * data[,(number_of_inputs+1):ncol(data)], MARGIN = 1 , FUN = sum),3)

                #final df
                Costa_df = data.frame(DMU = c(1:nrow(data)), I =wsum_inputs, O = wsum_outputs, efficiency_binary = ifelse(eff_crs_in$eff==1,"Efficient","Inefficient") ,efficiency =  eff_crs_in$eff)


        }
        else if (model == "CRS-Out") {

                eff_crs_out <- crs_eff(dataset = data, num_of_inputs = number_of_inputs , orientation = "out")

                S_crs_outorient <- apply(X = eff_crs_out$vy,MARGIN = 1 , FUN = sum)

                modified_vi_outorient <- eff_crs_out$ux / S_crs_outorient
                modified_vo_outorient <- eff_crs_out$vy / S_crs_outorient

                wsum_inputs = round(apply(X = modified_vi_outorient * data[,1:number_of_inputs], MARGIN = 1 , FUN = sum),3)
                wsum_outputs = round(apply(X = modified_vo_outorient * data[,(number_of_inputs+1):ncol(data)], MARGIN = 1 , FUN = sum),3)

                Costa_df = data.frame(DMU = c(1:nrow(data)), I =wsum_inputs, O = wsum_outputs, efficiency_binary = ifelse(eff_crs_out$eff==1,"Efficient","Inefficient"), efficiency =  eff_crs_out$eff)

        }

        t <- Costa_df

        total_min <-  min(t$I,t$O)
        total_max <- max(t$I,t$O)

        g = ggplot2::ggplot() +
                ggplot2::geom_point(data = t, aes(x = I , y = O, colour = efficiency_binary ),
                                    size = point_size ,
                                    alpha = point_transparency ) +
                ggplot2::scale_colour_gradient(low = "red", high = "green") +
                ggplot2::scale_colour_manual(name = "DMU Color",
                                             values = c("Efficient"="gold","Inefficient" = "skyblue"))

        g = g + ggplot2::geom_abline(xintercept = 0 , yintercept = 0 , slope = 1, color = "red")
        g = g + ggplot2::coord_fixed(ratio = 1,  expand = TRUE)



        g<- g+ ggplot2::theme_linedraw()  + ggplot2::annotate("text", x = Inf, y = -Inf, label = "© DEA-Viz",
                                                              hjust=1.1, vjust=-1.1, col="blue", cex=6,
                                                              fontface = "bold", alpha = 0.4) +
                ggplot2::xlim(total_min , total_max) +
                ggplot2::ylim(total_min , total_max) +
                ggplot2::coord_fixed(ratio = 1) +
                ggplot2::ggtitle("DEA Costa Frontier")

        g <- ggplotly(g)

        return(g)
}

# MDS color-encoding

ratio_names <- function(dataset, number_of_inputs ){

        number_of_outputs <- ncol(dataset)-number_of_inputs

        index = 0
        ratio_vector = vector(length = number_of_inputs*number_of_outputs )

        for (o in 1:number_of_outputs ){
                for (i in 1:number_of_inputs){
                        index <- index + 1

                        ratio_vector[index] <- paste0(colnames(dataset)[number_of_inputs+o],"/",colnames(dataset)[i])
                }
        }

        return(ratio_vector)
}

mds_data_ratio <- function(dataset, number_of_inputs){

        number_of_outputs <- ncol(dataset)-number_of_inputs

        index = 0
        ratio_vector <- vector(length = number_of_inputs*number_of_outputs )

        number_of_ratios <- number_of_inputs*number_of_outputs*nrow(dataset)
        ratio_mat <- matrix(data = rep(NA,number_of_ratios),nrow = nrow(dataset))


        for (o in 1:number_of_outputs ){
                for (i in 1:number_of_inputs){
                        index <- index + 1
                        #ratio_vector[index]<- paste0("O",o,"/","I",i)
                        ratio_vector[index] <- paste0(colnames(dataset)[number_of_inputs+o],"/",colnames(dataset)[i])
                        ratio_mat[,index] <- dataset[,(number_of_inputs+o)]/dataset[,i]
                }
        }

        ratio_df <- as.data.frame(ratio_mat)
        colnames(ratio_df) = ratio_vector

        return(ratio_df)
}

mds_plot <- function(dataset, number_of_inputs, variables_transformation = "none", dist_method = "euclidean", mds_type = "ratio", encoding_column = "crs.efficiency", point_size = 2 , point_transparency = 0.5){
        # variables can be "ratio" or "original"
        # dist_method can be "euclidean" or "manhattan"
        # mds_type can be "ratio", "interval", "ordinal"
        set.seed(7)
        dataset <- data.frame(dataset)
        if (variables_transformation == "none"){
                mds_data <- dataset
                coloring_choices = c(colnames(dataset),"fdh.efficiency","vrs.efficiency","crs.efficiency","drs.efficiency","irs.efficiency","frh.efficiency")
                #print(coloring_choices)
                if (!(encoding_column %in% coloring_choices )) {
                        encoding_column <- "crs.efficiency"
                        print("The encoding column is not in the feasible list of inputs,outputs,and efficiencies. It has replaced by CRS efficiency scores")
                }
        } else if (variables_transformation == "ratio"){
                mds_data <- mds_data_ratio(dataset = dataset, number_of_inputs = number_of_inputs)
                coloring_choices = c(ratio_names(dataset,number_of_inputs),"fdh.efficiency","vrs.efficiency","crs.efficiency","drs.efficiency","irs.efficiency","frh.efficiency")
                #print(coloring_choices)
                if (!(encoding_column %in% coloring_choices )) {
                        encoding_column <- "crs.efficiency"
                        print("The encoding column is not in the feasible list of output/input ratios,and efficiencies. It has replaced by CRS efficiency scores")
                }
        }

        # can be euclidean, manhattan, etc
        mds_dist <- dist(x = mds_data, method = dist_method)

        mds_model <- smacof::smacofSym(delta = mds_dist, ndim = 2 , type = mds_type  )

        org_dataset <- dataset
        core_data <- mds_data

        coordinates <- mds_model$conf

        fdh.efficiency <- Benchmarking::dea(X = org_dataset[,1:number_of_inputs], Y = org_dataset[,(number_of_inputs+1):ncol(org_dataset)], RTS = "fdh")$eff
        vrs.efficiency <- Benchmarking::dea(X = org_dataset[,1:number_of_inputs], Y = org_dataset[,(number_of_inputs+1):ncol(org_dataset)], RTS = "vrs")$eff
        drs.efficiency <- Benchmarking::dea(X = org_dataset[,1:number_of_inputs], Y = org_dataset[,(number_of_inputs+1):ncol(org_dataset)], RTS = "drs")$eff
        crs.efficiency <- Benchmarking::dea(X = org_dataset[,1:number_of_inputs], Y = org_dataset[,(number_of_inputs+1):ncol(org_dataset)], RTS = "crs")$eff
        irs.efficiency <- Benchmarking::dea(X = org_dataset[,1:number_of_inputs], Y = org_dataset[,(number_of_inputs+1):ncol(org_dataset)], RTS = "irs")$eff
        frh.efficiency <- Benchmarking::dea(X = org_dataset[,1:number_of_inputs], Y = org_dataset[,(number_of_inputs+1):ncol(org_dataset)], RTS = "add")$eff


        final_dataset <- data.frame(core_data,coordinates, fdh.efficiency,vrs.efficiency,drs.efficiency,crs.efficiency,irs.efficiency,frh.efficiency)
        all_colnames <- c(colnames(core_data),"D1","D2","fdh.efficiency","vrs.efficiency","drs.efficiency","crs.efficiency","irs.efficiency","frh.efficiency")


        index <- which(all_colnames==encoding_column)

        total_min <-  min(final_dataset$D1,final_dataset$D2)
        total_max <- max( final_dataset$D1,final_dataset$D2)


        g<- ggplot2::ggplot(data = final_dataset) +
                ggplot2::geom_point(aes(D1,D2, color = final_dataset[,index]),

                                    size = point_size,
                                    alpha = point_transparency
                ) +
                ggplot2::scale_colour_gradient(low = "blue", high = "red", name = colnames(final_dataset)[index] )



        g<- g+ ggplot2::theme_linedraw()  + ggplot2::annotate("text", x = Inf, y = -Inf, label = "© DEA-Viz",
                                                              hjust=1.1, vjust=-1.1, col="blue", cex=6,
                                                              fontface = "bold", alpha = 0.4) +


                ggplot2::xlim(total_min, total_max) +
                ggplot2::ylim(total_min,total_max) +
                ggplot2::coord_fixed(ratio = 1) +
                ggplot2::ggtitle(paste("MDS Color-Plot of",colnames(final_dataset)[index], sep = " "))


        g <- plotly::ggplotly(g)
        return(g)

}


parcoo_plot <- function(dataset, number_of_inputs,  efficiency_coordinate = "none"){
        dataset <- data.frame(dataset)
        coords_list = list()

        if (efficiency_coordinate != "none"){
                if (efficiency_coordinate == "CRS") {
                        eff_crs_in <- crs_eff(dataset = dataset, num_of_inputs = number_of_inputs , orientation = "in")$eff
                        dataset$crs_efficiency <- eff_crs_in
                }
                else if (efficiency_coordinate == "VRS") {
                        eff_vrs_in <- vrs_eff(dataset = dataset, num_of_inputs = number_of_inputs , orientation = "in")$eff
                        dataset$vrs_efficiency <- eff_vrs_in

                }
        }

        for (c in colnames(dataset)){
                l = list(label = c , values = as.formula(paste0("~",c)))
                coords_list = c(coords_list,list(l))
        }
        p <- dataset %>%
                plotly::plot_ly(type = "parcoords",
                                dimensions = coords_list,
                                hovertext = 1:nrow(dataset))

        return(p)
}


threeD_plot <- function(dataset, number_of_inputs,dmu_labels = T , dim_x , dim_y, dim_z , color = "crs_efficiency"  ){
        num_of_inputs <- number_of_inputs
        if(dmu_labels){
                dmu_labels <-dataset[,1]
                inputs <- dataset[,2:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+2):ncol(dataset)] %>% as.matrix()
                dataset <- dataset[,-1] %>% as.matrix()
        } else {
                dmu_labels <- data.frame(dmu = 1:nrow(dataset))
                inputs <- dataset[,1:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+1):ncol(dataset)] %>% as.matrix()
                dataset <- dataset %>% as.matrix()
        }
        
        colnames(dmu_labels) <- "DMU"

        # dim_x,y,z and color can be a column name or crs_efficiency, vrs_efficiency in string format
        #dataset <- data.frame(dataset)
        #dataset$label <- dmu_labels
        eff_crs_in <- crs_eff(dataset = dataset, num_of_inputs = number_of_inputs , orientation = "in")$eff
        #print(eff_crs_in)
       
        eff_vrs_in <- vrs_eff(dataset = dataset, num_of_inputs = number_of_inputs , orientation = "in")$eff
        
        dataset <- as.data.frame(dataset)
        dataset$vrs_efficiency <- eff_vrs_in
        dataset$crs_efficiency <- eff_crs_in
        dataset$DMU <- dmu_labels$DMU

        #dataset$uniform <-  "blue"

        p <-
                plotly::plot_ly(dataset, x = as.formula(paste0("~",dim_x)),
                                y = as.formula(paste0("~",dim_y)),
                                z = as.formula(paste0("~",dim_z)) ,
                                color = as.formula(paste0("~",color))) %>%
                plotly::add_markers(text = ~DMU ) %>%
                plotly::layout(scene = list(xaxis = list(title = dim_x),
                                            yaxis = list(title = dim_y),
                                            zaxis = list(title = dim_z)))

        return(p)

}


cem_viz <- function(dataset, number_of_inputs, dmu_labels = TRUE ){
        # dataset <- jap_mat
        # number_of_inputs <- 3
        num_of_inputs <- number_of_inputs
        if(dmu_labels){
                dmu_labels <-dataset[,1]
                inputs <- dataset[,2:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+2):ncol(dataset)] %>% as.matrix()
                dataset <- dataset[,-1] %>% as.matrix()
        } else {
                dmu_labels <- data.frame(dmu = 1:nrow(dataset))
                inputs <- dataset[,1:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+1):ncol(dataset)] %>% as.matrix()
                dataset <- dataset %>% as.matrix()
        }
        
        colnames(dmu_labels) <- "DMU"
        
        c <- CEM(dataset = dataset, number_of_inputs =number_of_inputs )
        #t <- cem_unfolding()
        #c <- round(c,4)
        t <- smacof::unfolding(delta = round((1-c),2),ndim = 2)

        row_df <- data.frame(t$conf.row, "DMU" =dmu_labels, Type = "Rating", Shape = 19)
        col_df <- data.frame(t$conf.col, "DMU" = dmu_labels, Type = "Rated", Shape = 24)


        unfolded_df  <- rbind(row_df,col_df)
        unfolded_df$alpha = NA
        unfolded_df$point_size = NA

        unfolded_df$alpha[unfolded_df$Type == "Rating"] <-  0.5
        unfolded_df$alpha[unfolded_df$Type == "Rated"] <-  0.5

        unfolded_df$point_size[unfolded_df$Type == "Rating"] <- 1
        unfolded_df$point_size[unfolded_df$Type == "Rated"] <- 1


        total_min <-  min(unfolded_df$D1,unfolded_df$D2)
        total_max <- max( unfolded_df$D1,unfolded_df$D2)
        
        p <- unfolded_df %>% 
                plot_ly(x = ~D1 , y = ~ D2 , text = ~DMU , symbol = ~Type)
                
        
        # g<- ggplot2::ggplot(data = unfolded_df) +
        #         ggplot2::geom_point(aes(x = D1 , y = D2 ,
        #                        color = Type,
        #                        shape = factor(unfolded_df$Shape) ,
        #                        alpha = unfolded_df$alpha ,
        #                        size = unfolded_df$point_size)
        #         )+
        #         ggplot2::scale_alpha(guide = FALSE) +
        #         ggplot2::scale_size(guide = FALSE) +
        #         ggplot2::scale_size_identity() +
        #         ggplot2::scale_alpha_identity() +
        #         ggplot2::scale_color_manual(name = "Object Type", values = c("Rating"="blue","Rated"="red"), labels = c("Rating","Rated"))+
        #         ggplot2::scale_shape_manual(name = "Object Type", values = c(19,24), labels = c("Rating","Rated"))+
        #         ggplot2::coord_fixed(ratio = 1,  expand = TRUE) +
        #         ggplot2::ggtitle("Cross-Efficiency Unfolding")
        # 
        # g <- g  + ggplot2::theme_linedraw() + ggplot2::annotate("text", x = Inf, y = -Inf, label = "© DEA-Viz",
        #                                       hjust=1.1, vjust=-1.1, col="blue", cex=6,
        #                                       fontface = "bold", alpha = 0.4) +
        #         ggplot2::xlim(total_min , total_max) +
        #         ggplot2::ylim(total_min , total_max )

        #g <- g  +
        #        ggplot2::coord_cartesian(xlim = cem_ranges$x, ylim = cem_ranges$y, expand = TRUE)


        # g <- plotly::ggplotly(g)
        return(p)

}

ben_weights_heatmap <- function(dataset , number_of_inputs ){
        weight_dataset <- ben_weights(dataset = dataset , number_of_inputs = number_of_inputs)
        std_w <- weight_standardization(weight_dataset = weight_dataset,number_of_inputs = number_of_inputs)
        input_names <- paste0("I_",colnames(std_w)[1:number_of_inputs])
        output_names <- paste0("O_",colnames(std_w)[(1+number_of_inputs):ncol(std_w)])
        colnames(std_w) <- c(input_names, output_names)
        std_w_gathered <- std_w %>% mutate(DMU = row.names(.)) %>% gather(key = "Variable" , value = "std_weight", - "DMU")
        
        std_w_gathered %>% plot_ly(y = ~ DMU , x = ~ Variable , z = ~ std_weight, type = "heatmap")
        
}

paneldata_biplot <- function(panel_data, number_of_inputs, color = "crs_efficiency"){
        # panel data must be a tidy dataframe, i.e. long format
        # panel data must have the first column as DMU labels called "Label", the second column as time period called "Period"

        panel_scaled <- scale(panel_data[,-c(1,2)])
        prcomp_result <- prcomp(panel_scaled)

        number_of_inputs = 2



        unique_periods = unique(panel_data$Period)

        crs_efficiency <- c()
        vrs_efficiency <- c()

        for (period in unique_periods){
                #dataset <- panel_data %>% filter(Year = period) %>% select(-c("Year","Bank"))
                #print(period)

                dataset <- panel_data[panel_data$Period == period , -c(1,2)]
                dataset <- data.frame(dataset)

                #print(dataset)

                t_crs_efficiency <- crs_eff(dataset,number_of_inputs)$eff
                t_vrs_efficiency <- vrs_eff(dataset,number_of_inputs)$eff
                crs_efficiency <- round(c(crs_efficiency,t_crs_efficiency),2)
                vrs_efficiency <- round(c(vrs_efficiency,t_vrs_efficiency),2)


        }

        panel_data$crs_efficiency <- floor(crs_efficiency*100)
        panel_data$vrs_efficiency <- floor(vrs_efficiency*100)

        #print("------")
        #print(panel_data)

        gg <- ggplot2::autoplot(prcomp_result, data = panel_data,
                                loadings = TRUE, loadings.colour = 'blue',
                                loadings.label = TRUE, loadings.label.size = 4, loadings.label.colour = "red")

        gg<-
                gg + ggplot2::geom_text(aes(x = PC1 , y = PC2 , label = Label , color = get(color), size = get(color))) +
                ggplot2::geom_text(aes(x = PC1 , y = PC2 , label = get(color), color = get(color), size = get(color), hjust = 0 , vjust = -0.5)) +
                ggplot2::scale_color_gradient(low = "blue", high = "red") +
                ggplot2::theme_light() +
                #ggplot2::geom_point(aes(x = PC1 , y = PC2 , color = crs_efficiency), size = 3) +
                gganimate::transition_time(Period) +  labs(title = "Period: {frame_time}")

        #gg <- animate(gg, nframes = 100, fps=2)
        return(gg)
}

correlation_subplots <- function(dataset, dmu_labels = T , num_of_inputs ){
        
        if(dmu_labels){
                dmu_labels <-dataset[,1]
                inputs <- dataset[,2:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+2):ncol(dataset)] %>% as.matrix()
        } else {
                dmu_labels <- data.frame(dmu = 1:nrow(dataset))
                inputs <- dataset[,1:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+1):ncol(dataset)] %>% as.matrix()
        }
        
        colnames(dmu_labels) <- "DMU"
        
        total_variables = ncol(dataset)
        cols <- colnames(dataset)
        p <- list()
        counter <- 0
        for (i in 1:total_variables){
                if (i<total_variables ) {
                        for (j in (i+1):(total_variables)){
                                counter <- counter +  1

                                p[[counter]]<- plotly::plot_ly(type = "scatter",
                                                               data = data.frame(dataset,DMU = dmu_labels$DMU),
                                                               x = as.formula(paste0("~",cols[i])),
                                                               y = as.formula(paste0("~",cols[j])), 
                                                               text = ~ DMU

                                )


                        }

                }


        }

        p <- plotly::subplot(nrows = 2 , p,shareX = FALSE , shareY = FALSE , titleX = TRUE , titleY = TRUE)

        return(p)

}

corrMat_plot <- function(dataset){
        PerformanceAnalytics::chart.Correlation(dataset, histogram = TRUE , pch = 15)
}

dotplot_histograms <- function(dataset){
        dataset <- data.frame(dataset)
        number_of_variables <- ncol(dataset)
        all_plots<- lapply(X = 1:number_of_variables, function(x) plotly::ggplotly(ggplot2::ggplot()+
                                   ggplot2::geom_dotplot(data = dataset, aes(x = dataset[,x]) , alpha = 0.6)+
                                   ggplot2::theme_light()+
                                   ggplot2::xlab(colnames(dataset)[x]) +
                                   #ggplot2::scale_fill_manual(name = "DMUs", values =  c("Others"="blue" , "Selected"="orange")) +
                                   ggplot2::guides(fill = FALSE )))


        plotly::subplot(all_plots, nrows = 2, shareX = FALSE , shareY = FALSE , titleX = TRUE , titleY = TRUE)

}

histograms <- function( dataset , dmu_labels = T , num_of_inputs){
        
        if(dmu_labels){
                dmu_labels <-dataset[,1]
                inputs <- dataset[,2:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+2):ncol(dataset)] %>% as.matrix()
                scaled_df <- dataset[,-1] %>% scale()
                gathered_df <- data.frame(DMU = dataset[,1] ,scaled_df) %>% gather(key = "variable" , value = "value", -1) 
                
                
        } else if(dmu_labels) {
                dmu_labels <- data.frame(dmu = 1:nrow(dataset))
                inputs <- dataset[,1:num_of_inputs] %>% as.matrix()
                outputs <- dataset[,(num_of_inputs+1):ncol(dataset)] %>% as.matrix()
                gathered_df <- dataset %>% scale() %>% gather(key = "variable" , value = "value") 
                gathered_df <- gathered_df %>% mutate(DMU = row_number())
        }
        
        
        gathered_df %>%
                ggplot() + 
                geom_histogram(aes(value, fill = factor(variable)) ) + 
                facet_grid(vars(variable))+ theme(legend.position = "none")
        
}





