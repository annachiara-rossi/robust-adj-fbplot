####_____________________________ VISUALIZATION ____________________________####

library(ggplot2)
library(dplyr)

####__________________________ SIMULATION RESULTS __________________________####


setwd("~/github/inflation-factor-tuning-fbplot/simulation_study/final_simulation")

# copy here the list from ls_names.txt
ls_names <- c("200_c1_0_.RData","200_c1_0.05_.RData","200_c1_0.1_.RData",
              "200_c1_0.15_.RData","400_c1_0_.RData","400_c1_0.05_.RData",
              "400_c1_0.1_.RData","400_c1_0.15_.RData")


# first combination of parameters
load(ls_names[[1]])
df_F <- all_F
df_TPR <- all_TPR
df_FPR <- all_FPR
df_time <- all_time

# second combination of parameters
load(ls_names[[2]])
df_F <- rbind(df_F, all_F)
df_TPR <- rbind(df_TPR, all_TPR)
df_FPR <- rbind(df_FPR, all_FPR)
df_time <- rbind(df_time, all_time)

# ...
load(ls_names[[3]])
df_F <- rbind(df_F, all_F)
df_TPR <- rbind(df_TPR, all_TPR)
df_FPR <- rbind(df_FPR, all_FPR)
df_time <- rbind(df_time, all_time)

load(ls_names[[4]])
df_F <- rbind(df_F, all_F)
df_TPR <- rbind(df_TPR, all_TPR)
df_FPR <- rbind(df_FPR, all_FPR)
df_time <- rbind(df_time, all_time)

load(ls_names[[5]])
df_F <- rbind(df_F, all_F)
df_TPR <- rbind(df_TPR, all_TPR)
df_FPR <- rbind(df_FPR, all_FPR)
df_time <- rbind(df_time, all_time)

load(ls_names[[6]])
df_F <- rbind(df_F, all_F)
df_TPR <- rbind(df_TPR, all_TPR)
df_FPR <- rbind(df_FPR, all_FPR)
df_time <- rbind(df_time, all_time)

load(ls_names[[7]])
df_F <- rbind(df_F, all_F)
df_TPR <- rbind(df_TPR, all_TPR)
df_FPR <- rbind(df_FPR, all_FPR)
df_time <- rbind(df_time, all_time)

load(ls_names[[8]])
df_F <- rbind(df_F, all_F)
df_TPR <- rbind(df_TPR, all_TPR)
df_FPR <- rbind(df_FPR, all_FPR)
df_time <- rbind(df_time, all_time)

#-------------------------------------------------------------------------------
# I need to reorganize data for facet grid: 

B <- 64  # number of repetitions of the simulation

build_dataframe <- function(df, B){
  n_estimators = dim(df)[2]
  df <- stack(df)
  df$dimensionality <- rep(c(replicate(B*4,'dim = 200'), replicate(B*4,'dim = 400')), 
                             times=n_estimators)
  
  df$outlier_prop <- rep(rep(c('proportion of outliers: 0','proportion of outliers: 0.05',
                               'proportion of outliers: 0.10','proportion of outliers: 0.15'), 
                             each = B, times = 2), times=n_estimators)
  return(df)
}


df_F <- build_dataframe(df_F, B)
df_TPR <- build_dataframe(df_TPR, B)
df_FPR <- build_dataframe(df_FPR, B)
df_time <- build_dataframe(df_time, B)
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# BOXPLOTS
#-------------------------------------------------------------------------------

# F*
g_F <- ggplot(df_F, aes(x = ind, y = values)) +
  xlab('') + ylab('') +
  ggtitle('Optimal F value') +
  geom_boxplot() + 
  geom_jitter(height = 0, width = 0.2, size = 0.25, colour = 'red') +
  coord_flip()
g_F + facet_grid(dimensionality ~ outlier_prop) #, scales = 'free_x'


# TPR
g_TPR <- ggplot(df_TPR[complete.cases(df_TPR), ], aes(x = ind, y = values)) +
  xlab('') + ylab('') +
  ggtitle('TPR') +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2, size = 0.25, colour = 'red') +
  coord_flip()
g_TPR + facet_grid(dimensionality ~ outlier_prop) #, scales = 'free_x'


# FPR
g_FPR <- ggplot(df_FPR, aes(x = ind, y = values)) +
  xlab('') + ylab('') +
  ggtitle('FPR') +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2, size = 0.25, colour = 'red') +
  coord_flip()
g_FPR + facet_grid(dimensionality ~ outlier_prop) #, scales = 'free_x'


# computational time
g_time <- ggplot(df_time, aes(x = ind, y = values)) +
  xlab('') + ylab('time in seconds') +
  ggtitle('F* tuning computational time') +
  geom_boxplot() +
  coord_flip()
g_time + facet_grid(dimensionality ~ outlier_prop) #, scales = 'free_x'

#------ to plot computational time without Kendall
df_time_sub <- df_time[-which(df_time$ind == 'Kendall'),]

g_time_sub <- ggplot(df_time_sub[complete.cases(df_time_sub), ], aes(x = ind, y = values)) +
  xlab('') + ylab('time in seconds') +
  ggtitle('F* tuning computational time') +
  geom_boxplot() +
  coord_flip()
g_time_sub + facet_grid(dimensionality ~ outlier_prop)

# or use table instead
summary_time <- df_time %>%
  group_by(ind, dimensionality, outlier_prop) %>%
  summarise(across(.cols = everything(), list(mean = mean, sd = sd)))




#------------------------------------------------------------------------------#
#------------------------------ VIOLIN PLOT -----------------------------------#
#------------------------------------------------------------------------------#


# F*
g_F <- ggplot(df_F, aes(x = ind, y = values)) +
  xlab('') + ylab('') +
  ggtitle('Optimal F value') +
  geom_violin( scale = 'width' ) +
  geom_jitter(height = 0, width = 0.2, size = 0.25, colour = 'bisque4') +
  coord_flip()
g_F + facet_grid(dimensionality ~ outlier_prop) #, scales = 'free_x'


# TPR
g_TPR <- ggplot(df_TPR[complete.cases(df_TPR), ], aes(x = ind, y = values)) +
  xlab('') + ylab('') +
  ggtitle('TPR') +
  geom_violin(scale = 'width') +
  geom_jitter(height = 0, width = 0.2, size = 0.25, colour = 'bisque4') +
  coord_flip()
g_TPR + facet_grid(dimensionality ~ outlier_prop) #, scales = 'free_x'


# FPR
g_FPR <- ggplot(df_FPR, aes(x = ind, y = values)) +
  xlab('') + ylab('') +
  ggtitle('FPR') +
  geom_violin(scale = 'width') +
  geom_jitter(height = 0, width = 0.2, size = 0.25, colour = 'bisque4') +
  coord_flip()
g_FPR + facet_grid(dimensionality ~ outlier_prop)


# log-time violin plot, without the 'no adj' case
df_time_sub <- df_time[-which(df_time$ind == 'No.adjustment'),]

g_time_sub <- ggplot(df_time_sub[complete.cases(df_time_sub), ], aes(x = ind, y = values)) +
  scale_y_continuous(trans='log10')+
  xlab('') + ylab('time in seconds') +
  ggtitle('F* tuning log-computational time') +
  geom_violin(scale = 'width') +
  #geom_jitter(height = 0, width = 0.2, size = 0.25, colour = 'bisque4') +
  coord_flip()
g_time_sub + facet_grid(dimensionality ~ outlier_prop)



#------------------------------------------------------------------------------#
#----------------------------- 2D DENSITY PLOT --------------------------------#
#------------------------------------------------------------------------------#

df_TPR_FPR <- df_TPR
colnames(df_TPR_FPR)[1] <- 'TPR'
colnames(df_TPR_FPR)[2] <- 'estimator'
df_TPR_FPR$FPR <- df_FPR$values  

# fill with 1s where the TPR is NaN, due to no outliers
df_TPR_FPR$TPR[is.na(df_TPR_FPR$TPR)] <- 1

# Let's try with F vs FPR
df_F_FPR <- df_F
colnames(df_F_FPR)[1] <- 'F'
colnames(df_F_FPR)[2] <- 'estimator'
df_F_FPR$FPR <- df_FPR$values

df_F_FPR %>%
  ggplot(aes(x = F, y = FPR)) +
  geom_density_2d_filled(alpha = 0.5) +
  facet_grid(dimensionality ~ outlier_prop)


# df_FPR %>%
#   filter(ind=='Ledoit_Wolf',outlier_prop=="outlier share: 0", dimensionality=='dim = 200') %>%
#   ggplot() +
#   geom_violin(aes(y=values,x=ind))


#------------------------------------------------------------------------------#
#------------------------------- TPR VS FPR -----------------------------------#
#------------------------------------------------------------------------------#

# create a dataframe containing the mean and the sd by groups 

d <- df_TPR_FPR
summary_df <- d %>%
            group_by(estimator, dimensionality, outlier_prop) %>%
            summarise(across(.cols = everything(), list(mean = mean, sd = sd)))


g_TPR_FPR <- ggplot(summary_df, aes(x = TPR_mean, y = FPR_mean, colour = estimator)) +
  xlab('TPR') + ylab('FPR') +
  ggtitle('TPR VS FPR') +
  geom_point(aes(size = FPR_sd )) +
  scale_color_brewer(palette="YlOrRd") +
  theme(legend.title = element_text(size = 7), 
        legend.text = element_text(size = 7),
        axis.text.x = element_text(size=7, angle=0))
g_TPR_FPR + facet_grid(dimensionality ~ outlier_prop)



