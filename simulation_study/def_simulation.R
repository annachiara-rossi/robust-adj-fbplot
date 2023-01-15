################_______________ SIMULATION STUDY _______________################

# NB: this version is for runnning the simulation locally

setwd("~/github/inflation-factor-tuning-fbplot/simulation_study")

#------------------------------------------
rm(list = ls())

#------------------------------------------


#------------------------------------------
# import Tarabelloni code
source( 'rfda/basisOp/functional_basis.R')
source( 'rfda/geometry/geometry.R')
source( 'rfda/fdata/fData.R')
source( 'rfda/stats/stats.R' )
source( 'rfda/stats/generateFD.R')
source( 'rfda/stats/expCovariance.R')
source( 'rfda/stats/robustEstimators.R' )

library(Matrix) # needed for geometry.R
#------------------------------------------


#------------------------------------------
# my defined fbplot function
source("def_fbplot.R")
#------------------------------------------


# Simulate functional data

####____________________________ hyperparameters ___________________________####

n_curves <- 100  # number of functions generated
P <- c(200)#,400)
outlier_proportions <- c(0.05,0.1) #c(0,0.05,0.1,0.15)

# Amplitude factor - overall variability around mean
alph = 0.12

# Correlation decay factor - high beta = low correlation (noisy)
beta = 0.4

# Number of principal components to retain
L = 10


####_____________________________ contamination ____________________________####

contamination_types = 'c1' #c('c1','c2')

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
  else{
    message('Create a different contamination type')
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

data_path = "/Users/annachiararossi/github/inflation-factor-tuning-fbplot/simulation_study/simulated_datasets"
metrics_path = "/Users/annachiararossi/github/inflation-factor-tuning-fbplot/simulation_study/metrics"


####______________________________ SIMULATION ______________________________####


#cov_estimators <- c('Ledoit_Wolf','OGK_mad','OGK_Qn','MRCD_0.5',
#               'MRCD_0.75', 'Spherical', 'Median', 'Kendall' )
#,'kMRCD_alpha0.5','kMRCD_alpha0.75')

cov_estimators <- c('Ledoit_Wolf','MRCD_0.5','Spherical')

# repetitions
# TODO: check parallel implementation (replicate)
B <- 2

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

      all_F <- data.frame(matrix(ncol = length(cov_estimators), nrow = 0))
      colnames(all_F) <- cov_estimators
      all_TP <- all_F
      all_TN <- all_F
      all_FP <- all_F
      all_FN <- all_F
      all_TPR <- all_F
      all_FPR <- all_F
      all_time <- all_F
      
      for( k in 1:B ){  
        
        expCovRule.Gauss = distroRule.Gaussian( mu = rep( 0, geometry.L2$basis$nelements ),
                                                Sigma = diag( spectrum$eigen_values ) )
        f_data = generateFD( n_curves, expCovRule.Gauss, geometry.L2 ) + sin( 4 * pi * timeStr$grid )
        
        # generate the contaminated dataset
        f_data_mag <- add_magnitude_outliers(f_data, outlier_share,
                                             contamination_type = c,
                                             show_plot = FALSE)
        
        # save data for the current simulation
        df = f_data_mag$data$coefficients
        
        file_name = paste(k,p,c,outlier_share,'.csv',sep = '_')
        write.csv(df, file.path(data_path, file_name), row.names=FALSE)
        
        true_out <- which(f_data_mag$label == 2) # same, for the same outlier_share
        
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
                                      adjust = list(VERBOSE = FALSE),
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
        # at this point I have one result for each estimator: add the row to the dataframe
        
        all_TP[nrow(all_TP) + 1,] = TP_row
        all_TN[nrow(all_TN) + 1,] = TN_row
        all_FP[nrow(all_FP) + 1,] = FP_row
        all_FN[nrow(all_FN) + 1,] = FN_row
        
        all_TPR[nrow(all_TPR) + 1,] = TPR_row
        all_FPR[nrow(all_FPR) + 1,] = FPR_row
        all_F[nrow(all_F) + 1,] = F_row
        
        all_time[nrow(all_time) + 1,] = time_row
        
      }
      # here I have B results for each estimator
      
      # save results for each combination of the parameters
      file_name = paste(p,c,outlier_share,'.RData',sep = '_')
      save(all_F, all_TP, all_TN, all_FP, all_FN, all_TPR, all_FPR, all_time,
           file = file.path(metrics_path,file_name) )
      ls_names <- append(ls_names, file_name)
      
    }

  }

}

# save the list of names of parameters' combinations
write.table(ls_names, file = 'ls_names.txt', sep = ',', row.names = FALSE)




