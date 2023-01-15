################_______________ SIMULATION STUDY _______________################
# code for the final simulation, with more or less outlying obs

rm(list = ls())

#------------------------------------------
# import Tarabelloni code
source( 'rfda/basisOp/functional_basis.R')
source( 'rfda/geometry/geometry.R')
source( 'rfda/fdata/fData.R')
source( 'rfda/stats/stats.R' )
source( 'rfda/stats/generateFD.R')
source( 'rfda/stats/expCovariance.R')
source( 'rfda/stats/robustEstimators.R')

library(Matrix) # needed for geometry.R
#------------------------------------------


#------------------------------------------
# my defined fbplot function
source("def_fbplot.R")
#------------------------------------------


#------------------------------------------
data_path = "/u/archive/laureandi/rossi/simulated_datasets_final"
metrics_path = "/u/archive/laureandi/rossi/metrics_final"
#------------------------------------------


####_____________________________ contamination ____________________________####

contamination_types = 'shape' #c('c1','c2')

add_magnitude_outliers <- function(f_data, outlier_share,
                                   contamination_type,
                                   show_plot = TRUE)
{
  # function to generate magnitude outliers
  
  n_curves = dim(f_data$coefficients)[1]
  n_outliers = n_curves*outlier_share
  
  # flag as outliers the last outlier_share (flag = 2)
  out_highlighter <- rep(c(1,2), c(n_curves-n_outliers, n_outliers))
  
  # separate two sets of curves, using only the coefficients.
  f_data_temp <- f_data$coefficients[1:(n_curves*(1-outlier_share)),]
  
  # inflate the outlying curves according to contamination type
  if (contamination_type == 'c1'){
    if (outlier_share == 0){
      # otherwise it would go from n+1:n
      mag_temp = NULL
    }
    else{
      mag_temp <- f_data$coefficients[(n_curves*(1-outlier_share)+1):n_curves,] * runif(n_outliers,2,3)
    }
  }
  if (contamination_type == 'shape'){
    if (outlier_share == 0){
      # otherwise it would go from n+1:n
      mag_temp = NULL
    }
    else{
      expCovRule.Gauss = distroRule.Gaussian( mu = rep( 0, geometry.L2$basis$nelements ),
                                              Sigma = diag( spectrum$eigen_values ) )
      mag_temp = generateFD( n_outliers, expCovRule.Gauss, geometry.L2 ) + cos( 4 * pi * geometry.L2$basis$timeStr$grid )
      mag_temp = mag_temp$coefficients
    }
  }
  else{
    if (outlier_share == 0){
      # otherwise it would go from n+1:n
      mag_temp = NULL
    }
    else{
      
      # flag with 3 half of the existing outliers
      # in the odd case, we have one mild outliers more
      out_highlighter[(n_curves-ceiling(n_outliers/2)+1):n_curves] = 3
      
      # gross outliers
      mag_temp1 <- f_data$coefficients[(n_curves*(1-outlier_share)+1):(n_curves-n_outliers/2),] * runif(n_outliers/2,2,3) 
      # mild outliers
      mag_temp2 <- f_data$coefficients[(n_curves-n_outliers/2 + 1):n_curves,] * runif(n_outliers/2,0,0.9)  
      
      mag_temp <- rbind(mag_temp1, mag_temp2)
    }
  }

  # Add them back together and create again an fData object
  f_data_mag <- as.fData(rbind(f_data_temp, mag_temp), geometry = f_data$geometry)
  
  
  # plot
  if (show_plot == TRUE){
    par(mfrow=c(1,2))
    plot(f_data, main = "Gaussian functional data",cex.main =2)
    plot(f_data_mag, col=out_highlighter, main ="Contaminated data",cex.main =2)
  }
  
  return(list(data = f_data_mag,
              label = out_highlighter))
}

#-------------------------------------------------------------------------------


library(pbapply)

num_cores = parallel::detectCores() 
cl=parallel::makeCluster(num_cores)

####____________________________ hyperparameters ___________________________####

n_curves <- 100  # number of functions generated
P <- c(200,400)
outlier_proportions <- c(0,0.05,0.1,0.15)

# Amplitude factor - overall variability around mean
alph = 0.12

# Correlation decay factor - high beta = low correlation (noisy)
beta = 0.4

# Number of principal components to retain
L = 10


####______________________________ SIMULATION ______________________________####

cov_estimators <- c('No-adjustment', 'Ledoit_Wolf','OGK_mad', 'OGK_Qn',
                    'kMRCD', 
                    'MRCD_0.5', 'MRCD_0.75',
                    'Spherical', 'Median', 'Kendall')
# the 'No-adjustment' case refers to no adjustment in the fbplot

# all metrics of interest
metrics_list = list('all_TP','all_TN','all_FP','all_FN',
                    'all_TPR', 'all_FPR', 'all_F', 'all_time')

# repetitions
B <- num_cores*2


# simulation

simulation <- function(n_curves, p, outlier_share, c, geometry.L2, spectrum, cov_estimators){
  
  expCovRule.Gauss = distroRule.Gaussian( mu = rep( 0, geometry.L2$basis$nelements ),
                                          Sigma = diag( spectrum$eigen_values ) )
  f_data = generateFD( n_curves, expCovRule.Gauss, geometry.L2 ) + sin( 4 * pi * geometry.L2$basis$timeStr$grid )
  
  # generate the contaminated dataset
  f_data_mag <- add_magnitude_outliers(f_data, outlier_share,
                                       contamination_type = c,
                                       show_plot = FALSE)
  
  # save data for the current simulation
  df = f_data_mag$data$coefficients
  #file_name = paste(k,p,c,outlier_share,'.csv',sep = '_')
  #write.csv(df, file.path(data_path, file_name), row.names=FALSE)
  
  true_out <- which(f_data_mag$label == 2)
  
  F_row <- c()
  TP_row <- c()
  TN_row <- c()
  FP_row <- c()
  FN_row <- c()
  TPR_row <- c()
  FPR_row <- c()
  time_row <- c()
  
  for(cov_estimator in cov_estimators){
    
    time_start = Sys.time()
    
    # build the fbplot and extract F* and the outliers
    if( cov_estimator == 'No-adjustment' ){
      fb_obj <- my_fbplot.fData(Data = f_data_mag$data, Cov_estimator = cov_estimator,
                                display = FALSE)
    } else{
      fb_obj <- my_fbplot.fData(Data = f_data_mag$data, Cov_estimator = cov_estimator,
                              adjust = list(TPR = 2 * stats::pnorm( 4 * stats::qnorm( 0.15 ) )),
                              display = FALSE)
    }
    
    time_end = Sys.time()
    
    id_out <- fb_obj$ID_outliers
    
    TP <- length(intersect(true_out, id_out))
    FP <- length(setdiff(id_out, true_out))
    FN <- length(setdiff(true_out,id_out))
    TN <- n_curves - (TP + FP + FN)
    
    # save results in an array
    TP_row <- c(TP_row, TP)
    TN_row <- c(TN_row, TN)
    FP_row <- c(FP_row, FP)
    FN_row <- c(FN_row, FN)
    
    TPR_row <- c(TPR_row, TP/(TP + FN))
    FPR_row <- c(FPR_row, FP/(FP + TN))
    F_row <- c(F_row, fb_obj$Fvalue)
    
    time_row <- c(time_row, difftime(time_end, time_start, units = 'secs')[[1]])
    
  }
  
  return(list(TP = TP_row, 
              TN = TN_row, 
              FP = FP_row, 
              FN = FN_row, 
              TPR = TPR_row, 
              FPR = FPR_row, 
              Fvalue = F_row, 
              comput_time = time_row,
              simulated_data = df))
}



ls_names <- list()

# for each possible dimensionality
for(p in P){

  # data generating process according to P
  timeStr = EvenTimeStr( seq( 0, 1, length.out = p ) )

  expCov = ExpCovfunction( alph, beta )

  spectrum = ExpCov_eigencouples( timeStr, L, expCov )
  eigen.basis = generateBasis( 'Custom', timeStr, N = length( spectrum$eigen_values),
                               elements = spectrum$eigen_functions )

  geometry.L2 = Geometry( eigen.basis, space = 'L2' )

  for(c in contamination_types){

    # for each possible contamination proportion
    for (outlier_share in outlier_proportions){

      print(paste('Dimensionality = ',p,', Contamination = ',c,
                  ', Outlier_share = ',outlier_share))

      #----------------- parallel run of the simulation ----------------
      #cl=parallel::makeCluster(num_cores, type = 'FORK')
      parallel::clusterExport(cl=cl, varlist=c(ls(),'.my_fbplot_fData','Lwls2D','Matrix'))
      
      res = pbreplicate(n = B, 
                        expr = simulation(n_curves, p, outlier_share, c, geometry.L2, 
                                          spectrum, cov_estimators), 
                        simplify = TRUE, 
                        cl = cl)
      #parallel::stopCluster(cl)
      #-----------------------------------------------------------------
      
      j <- 0
      for (i in metrics_list){
        j <- j + 1
        m <- t(matrix(unlist(res[j,]), nrow = length(cov_estimators)))
        colnames(m) <- cov_estimators
        assign(i, data.frame(m))
      }

      # save results for the current combination of the parameters
      file_name = paste(p,c,outlier_share,'.RData',sep = '_')
      save(all_F, all_TP, all_TN, all_FP, all_FN, all_TPR, all_FPR, all_time,
           file = file.path(metrics_path,file_name) )
      ls_names <- append(ls_names, file_name)
      
      # save datasets from each repetition, for the current combination of parameters
      for (j in 1:B){
        df = res[,j]$simulated_data
        file_name = paste(j,p,c,outlier_share,'.csv',sep = '_')
        write.csv(df, file.path(data_path, file_name), row.names=FALSE)
      }

    }

  }

}

# save the list of names of parameters' combinations
write.table(ls_names, file = 'ls_names.txt', sep = ',', row.names = FALSE)

