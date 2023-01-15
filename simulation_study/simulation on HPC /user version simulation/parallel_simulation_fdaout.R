################_______________ SIMULATION STUDY _______________################
# NB: DATA GENERATION DONE THROUGH PACKAGE fdaoutlier

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
source("user_fbplot.R")
library(fdaoutlier)
#------------------------------------------


#------------------------------------------
data_path = "/u/archive/laureandi/rossi/simulated_datasets_final"
 #"/Users/annachiararossi/github/inflation-factor-tuning-fbplot/simulation_study/simulated_datasets" 
metrics_path = "/u/archive/laureandi/rossi/metrics_final"
  #"/Users/annachiararossi/github/inflation-factor-tuning-fbplot/simulation_study/metrics"
#------------------------------------------


####_____________________________ contamination ____________________________####

contamination_types = 'c1' 

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
                    'kMRCD', 'MRCD',
                    'Spherical', 'Median', 'Kendall')
# the 'No-adjustment' case refers to no adjustment in the fbplot

# all metrics of interest
metrics_list = list('all_TP','all_TN','all_FP','all_FN',
                    'all_TPR', 'all_FPR', 'all_F', 'all_time')

# repetitions
B <- num_cores*2

# simulation

simulation <- function(n_curves, p, outlier_share, c, grid, cov_estimators){
  
  s <- fdaoutlier::simulation_model2(n = n_curves, p = p, outlier_rate = outlier_share)
  
  f_data_mag <- roahd::fData(grid, s$data)
  
  # save data for the current simulation
  #df = f_data_mag$data$values
  #file_name = paste(k,p,c,outlier_share,'.csv',sep = '_')
  #write.csv(df, file.path(data_path, file_name), row.names=FALSE)
  
  # NB this has to be changed when there are both gross and mild outliers
  true_out <- s$true_outliers
  
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
      fb_obj <- my_fbplot.fData(Data = f_data_mag, adjust = FALSE, Fvalue = 1.5, display = FALSE)
      
    } else{
      fb_obj <- my_fbplot.fData(Data = f_data_mag, adjust = list(Cov_estimator = cov_estimator, 
                                                                      TPR = 2 * stats::pnorm( 4 * stats::qnorm( 0.15 ) )),
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
  
  grid <-  seq( 0, 1, length.out = p)

  for(c in contamination_types){

    # for each possible contamination proportion
    for (outlier_share in outlier_proportions){

      print(paste('Dimensionality = ',p,', Contamination = ',c,
                  ', Outlier_share = ',outlier_share))

      #----------------- parallel run of the simulation ----------------
      #cl=parallel::makeCluster(num_cores, type = 'FORK')
      parallel::clusterExport(cl=cl, varlist=c(ls(),'.my_fbplot_fData','Matrix','sparseMatrix','Lwls2D')) 
      
      res = pbreplicate(n = B, 
                        expr = simulation(n_curves, p, outlier_share, c, grid, cov_estimators), 
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
      #for (j in 1:B){
      #  df = res[,j]$simulated_data
      #  file_name = paste(j,p,c,outlier_share,'.csv',sep = '_')
      #  write.csv(df, file.path(data_path, file_name), row.names=FALSE)
      #}

    }

  }

}

# save the list of names of parameters' combinations
write.table(ls_names, file = 'ls_names.txt', sep = ',', row.names = FALSE)

