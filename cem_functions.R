library(lpSolve)
library(lpSolveAPI)
library(Benchmarking)


crs_eff <- function(dataset, num_of_inputs, orientation = "in" ){
    inputs <- dataset[,1:num_of_inputs]
    outputs <- dataset[,(num_of_inputs+1):ncol(dataset)]
    Benchmarking::dea(X = inputs, Y = outputs, RTS = "crs" , DUAL = TRUE, ORIENTATION = orientation )
}



vrs_eff <- function(dataset, num_of_inputs, orientation = "in"){
    inputs <- dataset[,1:num_of_inputs]
    outputs <- dataset[,(num_of_inputs+1):ncol(dataset)]
    Benchmarking::dea(X = inputs, Y = outputs, RTS = "vrs" , DUAL = TRUE , ORIENTATION = orientation)
}


A_mat = function(unit, dataset, number_of_inputs) {
    #this function returns A matrix - matrix of constraints - for DMUunit formulation.
    # unit can be a number from 1 to nrow(input_mat)
    # This function can be developed in a way that it returns all A matrices for
    # all units in a data structure such as list. It would be faster but
    # pre-mature optimization is the root of all devils!
    
    number_of_units = nrow(dataset)
    number_of_inputs = number_of_inputs
    number_of_outputs = ncol(dataset)- number_of_inputs
    number_of_variables = number_of_inputs+number_of_outputs
    
    
    A = matrix(nrow = number_of_units, ncol= number_of_variables)
    
    
    
    #first constraint is always the numerator constraint
    A[1,]= unlist(c(rep(0,number_of_outputs),dataset[unit,1:number_of_inputs]))
    ##print(dim(A))
    
    
    input_mat_without_unit <- as.matrix(dataset[-unit,1:number_of_inputs])
    output_mat_without_unit <- as.matrix(dataset[-unit,(number_of_inputs+1):number_of_variables])
    
    
    d = number_of_units-1
    
    ##print(d)
    for (i in 1:d) {
        ##print(i)
        A[(i+1), ]=c(output_mat_without_unit[i,],-1*input_mat_without_unit[i,])
        
    }
    
    
    return(A)
    
} # A_mat


dea_4cem = function(dataset, number_of_inputs, epsilon = 0.000001) {
    # returns a list , the first element is the efficiency of the DMUs (CRS)
    # and the second element is the matrix of the optimum weights
    
    # because of using lp() without lpSolveAPI, the lower bound is ignored in some cases
    # so better to round the results in 4 digits or so, then there would be no negative value
    
    # this function is supposed to replace dea() of Benchmarking
    # This function supposed to give back the optimum [ or one optimum set of] weights
    
    # since dea() is giving bulshit results for Colombian hospital case
    # indeed, the summation of the VIs is not equal to one!
    # so it means the formulation that dea() is using is not the one that I need
    
    
    
    number_of_units = nrow(dataset)
    number_of_inputs = number_of_inputs
    number_of_outputs = ncol(dataset)- number_of_inputs
    number_of_variables = number_of_inputs+number_of_outputs
    
    
    eff_weight_mat = matrix(nrow = number_of_units, ncol = (number_of_variables+1))
    
    lprec_eff_scores <- vector(length = number_of_units)
    lprec_opt_weights <- matrix(nrow = number_of_units, ncol = number_of_variables)
    
    for (unit in 1:number_of_units){
        
        OF = unlist(c(as.vector(dataset[unit,(number_of_inputs+1):number_of_variables]),
                      rep(0,number_of_inputs)))
        
        A_upper = unlist(c(rep(0,number_of_outputs),dataset[unit,1:number_of_inputs]))
        A_middle = cbind(dataset[-unit,(number_of_inputs+1):number_of_variables],-dataset[-unit,1:number_of_inputs])
        A_bottom = OF
        A_upper = unname(A_upper)
        A_middle = unname(A_middle)
        A_bottom = unname(A_bottom)
        
        A = as.matrix(rbind(A_upper, A_middle, A_bottom))
        
        const_directions = c("==",rep("<=",(number_of_units)))
        
        RHS_values = c(1,rep(0,(number_of_units-1)),1)
        
        lp_model=lpSolve::lp(direction = "max",
                             objective.in = OF ,
                             const.mat = A ,
                             const.dir = const_directions,
                             const.rhs = RHS_values)
        
        eff_weight_mat[unit,] = c(lp_model$objval,lp_model$solution)
        
        #----- new implementation
        
        #OF
        OF <- dataset[unit,]
        
        
        lprec <- lpSolveAPI::make.lp(nrow = 0 , ncol = number_of_variables )
        lpSolveAPI::set.objfn(lprec, OF)
        lprec_const_directions <- gsub(x = const_directions, pattern = "==", replacement = "=")
        for (constraint_no in 1:nrow(A)) {
            add.constraint(lprec, A[constraint_no,], lprec_const_directions[constraint_no], RHS_values[constraint_no])
        }
        lpSolveAPI::set.bounds(lprec, lower = rep(x = epsilon, number_of_variables))
        ColNames <- c(colnames(dataset)[(number_of_inputs+1):number_of_variables],colnames(dataset)[1:number_of_inputs])
        dimnames(lprec) <- list(1:nrow(A), ColNames)
        
        #lprec
        lpSolveAPI::lp.control(lprec,sense='max')
        solve(lprec)
        lprec_eff_scores[unit] <- lpSolveAPI::get.objective(lprec)
        lprec_opt_weights[unit, ] <- lpSolveAPI::get.variables(lprec)
        
        
    }
    eff_scores <-   eff_weight_mat[,1]
    weight_mat <- eff_weight_mat[,-1]
    
    Input_weights_colnames <- paste0(colnames(dataset)[1:number_of_inputs]," Weight")
    Output_weights_colnames <-
        paste0(colnames(dataset)[(number_of_inputs+1):number_of_variables]," Weight")
    
    colnames(weight_mat) <- c(Output_weights_colnames,Input_weights_colnames)
    
    # ----- using dea()
    bench_model <- Benchmarking::dea(X = dataset[,1:number_of_inputs],Y = dataset[,(number_of_inputs+1):number_of_variables] , RTS = "crs", DUAL = TRUE )
    bench_scores <- bench_model$eff
    bench_weights <- cbind(bench_model$vy, bench_model$ux)
    
    eff_weight_list = list(eff_scores, weight_mat, lprec_eff_scores, lprec_opt_weights,bench_scores,bench_weights)
    
    return(eff_weight_list)
    
    
    
} # dea_4cm

CEM_unit = function(dataset , unit , epsilon = 0.000001, number_of_inputs ){
    #this function must return the the benevolent optimum weights for the given unit
    # because of using lp() without lpSolveAPI, the lower bound is ignored in some cases
    # so better to round the results in 4 digits or so, then there would be no negative value
    
    set.seed(7)
    
    number_of_units = nrow(dataset)
    number_of_inputs = number_of_inputs
    number_of_outputs = ncol(dataset)- number_of_inputs
    number_of_variables <- ncol(dataset)
    
    number_of_variables = ncol(dataset)
    
    
    #eff_weight_list <- dea_4cem(dataset = dataset, number_of_inputs = number_of_inputs)
    ##4th or 2nd component of dea_4cem?
    #eff_weight_mat <- eff_weight_list[[2]]
    #simple_eff <- eff_weight_list[[1]]
    
    bench_model <- Benchmarking::dea(X = dataset[,1:number_of_inputs],Y = dataset[,(number_of_inputs+1):number_of_variables] , RTS = "crs", DUAL = TRUE )
    simple_eff <- bench_model$eff
    eff_weight_mat <- cbind(bench_model$vy, bench_model$ux)
    
    
    #simple_eff = eff_weight_mat[,1]
    unit_simple_eff = simple_eff[unit]
    
    
    
    output_mat_refined<- as.matrix(dataset[-unit,(number_of_inputs+1):number_of_variables])
    #output_mat_refined = output_mat[-unit,]
    
    
    OF = apply(output_mat_refined,2,sum)
    OF = c(OF,rep(0,number_of_inputs))
    
    
    #preparing the A matrix
    A_middle = diag(number_of_variables)
    
    
    #---- needs debugging A_mat()
    A_upper=A_mat(unit=unit, dataset=dataset, number_of_inputs = number_of_inputs)
    
    A_last = unlist(c(dataset[unit,(number_of_inputs+1):number_of_variables],
                      -1*unit_simple_eff*dataset[unit,1:number_of_inputs]))
    
    
    A = rbind(A_upper,A_middle,A_last)
    
    #two new lines
    
    #preparing the constraint directions
    C_upper = c("==",rep("<=",number_of_units-1))
    C_middle = rep(">=",number_of_variables)
    C_last = "=="
    C = c(C_upper,C_middle,C_last)
    
    #print(C)
    
    #preparing the RHS values
    rhs_upper = c(1,rep(0,number_of_units-1))
    rhs_middle = rep(0,number_of_variables)
    rhs_last = 0
    rhs_total = c(rhs_upper, rhs_middle,rhs_last)
    
    #benevolent or aggressive CEM?
    #approach <- ifelse(input$cem_approach=="Aggressive","min","max")
    #print(approach)
    approach = "max"
    t = lpSolve::lp(direction = approach, objective.in = OF , const.mat = A , const.dir = C , const.rhs = rhs_total)
    
    
    if (sum(t$solution) == 0 ) {
        cem_weight <- eff_weight_mat[unit,]
    } else {
        cem_weight <- t$solution
    }
    
    #return(t$solution)
    #print(cem_weight)
    return(cem_weight)
} # CEM_unit


CEM = function(dataset, number_of_inputs) {
    #this function returns the benevolent CEM
    
    
    number_of_units = nrow(dataset)
    number_of_inputs = number_of_inputs
    number_of_outputs = ncol(dataset)- number_of_inputs
    number_of_variables = number_of_inputs+number_of_outputs
    # print("==============")
    # print("checkpoint 1")
    # print("==============")
    eff_weight_list <- dea_4cem(dataset = dataset, number_of_inputs = number_of_inputs)
    # print("==============")
    # print("checkpoint 2")
    # print("==============")
    
    eff_weight_mat <- eff_weight_list[[4]]
    # print("==============")
    # print("checkpoint 3")
    # print("==============")
    CEM_opt_weights = matrix(nrow = number_of_units, ncol = number_of_variables )
    # print("==============")
    # print("checkpoint 4")
    # print("==============")
    
    for (unit in 1:number_of_units) {
        
        CEM_opt_weights[unit,] = CEM_unit(dataset , unit , epsilon = 0.00001, number_of_inputs = number_of_inputs)
        
        
    }
    
    
    
    CEM = matrix (nrow = number_of_units, ncol = number_of_units)
    # print("==============")
    # print("checkpoint 5")
    # print("==============")
    # matrix multiplication solution
    w_outputs = CEM_opt_weights[,1:number_of_outputs]
    w_inputs = CEM_opt_weights[,(number_of_outputs+1):number_of_variables]
    inputs <- dataset[,1:number_of_inputs]
    outputs <- dataset[,(number_of_inputs+1):number_of_variables]
    # print("==============")
    # print("checkpoint 6")
    # print("==============")
    cem_output <- w_outputs %*% t(outputs)
    cem_input <- w_inputs %*% t(inputs)
    # print("==============")
    # print("checkpoint 7")
    # print("==============")
    CEM <- cem_output/ cem_input
    # print("==============")
    # print("checkpoint 8")
    # print("==============")
    return(CEM)
    
} #CEM

MI <- function(dataset, number_of_inputs, dmu_labels = TRUE, dist_method = "euclidean"){
    # dataset <- jap_df 
    # number_of_inputs <- 3 
    # dmu_labels = TRUE 
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
    
    cem <- CEM(dataset = dataset , number_of_inputs = num_of_inputs )
    avg_creff <- apply(cem, MARGIN = 2 , FUN = mean)
    
    cem_df <- data.frame(dmu_labels, cem %>% data.frame())
    colnames(cem_df) <- c("DMU", paste0("col_",dmu_labels$DMU))
    
    self_eff <- diag(cem)
    creff_df <- data.frame(dmu_labels, avg_creff , self_eff)
    
    # Doyle and Green MI
    creff_df_with_DnG <- 
        creff_df %>% 
        mutate(DnG_MI = (self_eff-avg_creff)/avg_creff)
    
    # The new MI 
    row_to_row_dist <-dist(x =  cem_df[,-1] , method = dist_method, diag = T , upper = T) %>% as.matrix()
    col_to_row_dist <- round(1-cem,4) 
    RII <- apply(row_to_row_dist , MARGIN = 2 , FUN = mean)
    CII <- apply(col_to_row_dist, MARGIN = 2 , FUN = mean)
    #unfolding_res<- smacof::unfolding(delta = round(1-cem,2),ndim = nrow(dataset)*2 , type = "ratio"  )
    creff_df_MIs <- 
        creff_df_with_DnG %>% 
        mutate(RII= RII , CII = CII , new_MI = RII * CII)
    
    return(creff_df_MIs)
}


ben_weights <- function(dataset, number_of_inputs, dmu_labels = TRUE) {
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
    #this function is supposed to return the benevolent optimum weights
    #first the outputs, then the inputs
    number_of_units = nrow(dataset)
    
    number_of_outputs = ncol(dataset)-number_of_inputs
    #number_of_factors = ncol(dataset)
    
    number_of_variables <- ncol(dataset)
    
    
    bench_model <- Benchmarking::dea(X = inputs,Y = outputs , RTS = "crs", DUAL = TRUE )
    
    eff_weight_mat <- cbind(bench_model$vy, bench_model$ux)
    
    CEM_opt_weights = matrix(nrow = number_of_units, ncol = number_of_variables )
    
    
    # max normalization for getting reasonable weights scale
    # it is possible to use other normalizations
    dataset <- apply(X = dataset, MARGIN = 2 , FUN = function(x) x/max(x))
    
    for (unit in 1:number_of_units) {
        CEM_opt_weights[unit,] = CEM_unit(dataset , unit , epsilon = 0.00001, number_of_inputs = number_of_inputs)
        
        if (sum(CEM_opt_weights[unit,]) == 0 ) {
            CEM_opt_weights[unit,] <- eff_weight_mat[unit,]
            print("warning! zero weight vector detected")
            print("replaced with:")
            print(CEM_opt_weights[unit,])
        }
    }
    
    
    factor_labels = vector(length = number_of_variables)
    factor_labels = c(colnames(dataset[,(number_of_inputs+1):number_of_variables]),
                      colnames(dataset[,1:number_of_inputs]))
    
    factor_labels <- paste(factor_labels,"Weight",sep = " ")
    
    row_labels <- dmu_labels$DMU
    
    CEM_opt_weights = as.data.frame(CEM_opt_weights)
    colnames(CEM_opt_weights) = factor_labels
    row.names(CEM_opt_weights) = row_labels
    
    
    
    
    
    return(round(CEM_opt_weights,5))
} # benevolent weights

weight_standardization <- function(weight_dataset , number_of_inputs){
    num_of_inputs <- number_of_inputs
    number_of_variables <- ncol(weight_dataset)
    number_of_outputs <- number_of_variables - number_of_inputs
    
    #output weights standardization
    if (number_of_outputs != 1 ) {
        output_standardization_factor <-
            apply(X = weight_dataset[,1:number_of_outputs] , MARGIN = 1 , FUN = sum )
        standardized_outputs <- weight_dataset[,1:number_of_outputs]/output_standardization_factor
        
        
    } else {
        
        standardized_outputs <- weight_dataset[,1:number_of_outputs]/weight_dataset[,1:number_of_outputs]
        standardized_outputs <- as.matrix(standardized_outputs)
        colnames(standardized_outputs) <- colnames(weight_dataset)[1:number_of_outputs]
    }
    
    #input weights standardization
    if (number_of_inputs != 1){
        input_standardization_factor <-
            apply(X = weight_dataset[,(number_of_outputs+1):number_of_variables]
                  , MARGIN = 1 , FUN = sum )
        
        standardized_inputs <-
            weight_dataset[,(number_of_outputs+1):number_of_variables]/input_standardization_factor
        
    } else {
        standardized_inputs <-
            weight_dataset[,(number_of_outputs+1):number_of_variables]/
            weight_dataset[,(number_of_outputs+1):number_of_variables]
        
        standardized_inputs <- as.matrix(standardized_inputs)
        colnames(standardized_inputs) <- colnames(weight_dataset)[(number_of_outputs+1):number_of_variables]
    }
    
    standardized_weights_list <- list(round(standardized_inputs,5),round(standardized_outputs,5))
    
    
    return(data.frame(round(standardized_inputs,5),round(standardized_outputs,5)))
} #weight standardization



