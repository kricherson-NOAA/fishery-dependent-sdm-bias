library(tidyverse)
library(philentropy) #K-L distance (binning method)
library(statip) #Hellinger function 
#library(ggrepel) 
library(RColorBrewer)
library(gridExtra)
library(readxl)

#Load the data
dat_hist <- readRDS("dat_hist_results_full_2_17_21.rds")

dat_hist <- dat_hist %>% 
  mutate(all_sampled = 1) #add a dummy column for ALL sampled if we want to use that 

#Load the RMSE data
rmse <- read_xlsx("performance_metrics_SDMs.xlsx", sheet = 1) %>% 
  slice(29:46) %>% 
  mutate(sampling_regime = names(dat_hist)[grepl("sampled", names(dat_hist)) &
                                             !grepl("all_", names(dat_hist))]) %>% 
  dplyr::select(16,15) %>% 
  rename(mean_rmse = `...15`) %>% 
  mutate(mean_rmse = as.numeric(mean_rmse))
  
  

#This function is taken directly from here: https://stackoverflow.com/questions/15436702/estimate-cohens-d-for-effect-size
# Modified on 5April2021 to not take absolute value of mean difference
#It looks correct to me and results match cohen.d() from effsize package, but using this means you don't have to install an extra package
cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  md  <- mean(x) - mean(y)        ## mean difference (numerator)
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation
  
  cd  <- md/csd                        ## cohen's d
  
  return(cd)
}

#Slightly altered from Steven's code. My change is just setting n = 1024 instead of 1000 because documentation says "it almost always makes sense to specify n as a power of two"
norm_vec <- function(x) sqrt(sum(x^2))

hell_dist <- function (p, q, from, to, n = 1024) {
  P <- density(p, kernel = "gaussian", from = from, to = to, n = n)
  p <- P$y
  p <- p / sum(p)
  Q <- density(q, kernel = "gaussian", from = from, to = to, n = n)
  q <- Q$y
  q <- q / sum(q)
  hd <- norm_vec(sqrt(p) - sqrt(q)) / sqrt(2)
  hd
}

#Function that returns Cohen's d, K-L distance, and Hellinger distance for any given pair of comparisons
#Predictor is the name of the column with the variable we want to compare, sampregime is the sampling regime column name
#Returns a list with the data frame of results as the first element and a ggplot as the second element


compare_dat <- function(predictor, #Can be "temp" "zoo_200", "mld" or "chl_surface"
                        sampregime1, #Can be any of the sampling regimes (e.g."random_sampled", "pref_sampled_1") or "all" for all data
                        sampregime2) 
{
  
  dat <- dat_hist
  
  #vectors of the data we are comparing (sampling regime 1 and 2)
  var1 <- dat[ ,predictor][dat[ ,sampregime1] == 1]
  var2 <- dat[ ,predictor][dat[ ,sampregime2] == 1]
  
  #Calculate Cohen's D
  cd <- cohens_d(var1,
                 var2)
  
  #Do we want a descriptor of the effect size for reference?
  #Note: Cohen suggested that d=0.2 be considered a 'small' effect size, 0.5 represents a 'medium' effect size and 0.8 a 'large' effect size
  cd_effect <- ifelse(cd <= 0.2, "small",
                      ifelse(cd > 0.2 & cd <= 0.5, "medium", 
                             "large"))
  
  #Calculate K-L distance with a binned approach using density
  
  n_bins <- 1024 #arbitrary
  
  bins <- seq(floor(min(c(var1, var2))), ceiling(max(c(var1, var2))), length = n_bins)
  
  var1_dens <- density(var1, n=n_bins, from=bins[1], to=bins[n_bins])
  var2_dens <- density(var2, n=n_bins, from=bins[1], to=bins[n_bins])
  
  # density output does not sum to 1, take sum of the density vector and scale so sums to 1
  dens1 <- tibble(x=var1_dens$x,
                  dens1=var1_dens$y) %>%
    mutate(total_dens1=sum(dens1),
           rel_dens1=dens1/total_dens1)
  
  dens2 <- tibble(x=var2_dens$x,
                  dens2=var2_dens$y) %>%
    mutate(total_dens2=sum(dens2),
           rel_dens2=dens2/total_dens2)
  
  
  dens_join <- full_join(dens1, dens2, by = "x") %>% 
    replace(., is.na(.), 0)
  
  kld_bin <- suppressMessages(philentropy::KL(rbind(dens_join$rel_dens1, dens_join$rel_dens2)) %>% as.numeric())
  
  #Hellinger distance (continuous version). Integration fails for comparing chl for all vs. dist_sampled_mpn, so I have it return an NA
  hell_dist_cont <- tryCatch(statip::hellinger(var1, var2, lower = -Inf, upper = Inf), error=function(err) NA) 
  
  #Discrete Hellinger distance, using same bins as K-L distance
  hell_dist_discr <- hell_dist(var1, var2, from=bins[1], to=bins[n_bins])
  
  #Output
  out <- data.frame("sampling_regime_1" = sampregime1,
                    "sampling_regime_2" = sampregime2,
                    "predictor" = predictor,
                    "cohens_d" = cd, 
                    "cohens_d_effect" = cd_effect, 
                    "kullback-leibler_dist" = kld_bin,
                    "hellinger_dist" = hell_dist_discr,
                    "hellinger_dist_discr" = hell_dist_discr)
  
  #plot
  df1 <- data.frame(var = var1, 
                    regime = sampregime1)
  
  df2 <- data.frame(var = var2, 
                    regime = sampregime2)
  
  df <- rbind(df1, df2)
  
  plotlabel <- paste0("Cohen's D = ", round(cd, 3), ", ", cd_effect, "\n",
                      "Kullback-Leibler distance = ", round(kld_bin, 3),"\n",
                      "Hellinger distance = ", round(hell_dist_discr, 3))
  
  histplot <- ggplot(df, aes(x = var, group = regime, fill = regime)) +
    #geom_histogram(position = "dodge", bins = 30)+
    stat_bin(aes(y = ..density..), position = 'dodge', bins = 30)+
    annotate(geom = "label", x = Inf, y = Inf, label = plotlabel,  vjust = "inward", hjust = "inward")+
    theme_bw()+
    ggtitle(paste0(predictor, ", ", sampregime1, " ", " ", "vs " , sampregime2, " ")) +
    labs(x = predictor)
  
  return(list(out,histplot))
  
}


#Apply to all combinations of random + fishery dependent sampling, predictors

fishdep_regimes <- names(dat_hist)[grepl("sampled", names(dat_hist)) &
                                     #!grepl("random_", names(dat_hist)) &
                                     !grepl("all_", names(dat_hist))]

predictors <- c("temp", "zoo_200", "mld", "chl_surface")


comps <- expand.grid(c("random_sampled", "all_sampled"), fishdep_regimes, predictors, stringsAsFactors = FALSE) %>% 
  rename("regime1" = Var1, "regime2" = Var2, "predictor" = Var3) %>% 
  filter(paste0(regime1, regime2) != "random_sampledrandom_sampled") %>% 
  arrange(regime1, regime2)

#Using a loop for now until I figure out how to use pmap right...
df_res <- list()
plot_res <- list()

for(i in 1:nrow(comps))
{
  res <- compare_dat(predictor = comps$predictor[i],
                     sampregime1 = comps$regime1[i],
                     sampregime2 = comps$regime2[i])
  
  df_res[[i]] <- res[[1]]
  plot_res[[i]] <- res[[2]]
  
}

save(plot_res, file = paste0(""))

#put results into a data frame
res_df <- bind_rows(df_res) %>%  #not sure why this changes kullback_leibler_dist to kullback.leibler_dist??
  rename(kullback_leibler_dist = kullback.leibler_dist) %>% 
  mutate(samp_grp = case_when(grepl("pref", sampling_regime_2) ~ "pref",
                              grepl("dist", sampling_regime_2) ~ "dist",
                              grepl("Closed", sampling_regime_2) ~ "closed",
                              grepl("BY", sampling_regime_2) ~ "bycatch",
                              grepl("random", sampling_regime_2) ~ "random")) %>% 
  arrange(samp_grp, sampling_regime_2) %>% 
  mutate(samp_reg = sampling_regime_2) %>% #Display names for plotting
  mutate(samp_reg = gsub("_sampled_", " ", samp_reg),
         samp_reg = ifelse(samp_reg == "BY_sampled", "Bycatch", samp_reg),
         samp_reg = ifelse(samp_reg == "random_sampled", "Random vs all", samp_reg),
         samp_reg = gsub("dist", "Dist", samp_reg),
         samp_reg = gsub("pref", "Pref", samp_reg)) %>% 
  rename(`Sampling regime` = samp_reg) %>% 
  mutate(disp_predictor = case_when(predictor == "temp" ~ "Temperature",
                                    predictor == "zoo_200" ~"Zooplankton",
                                    predictor == "mld" ~ "Mixed layer depth",
                                    predictor == "chl_surface" ~ "Chlorophyll")) %>% 
  left_join(rmse %>% rename(sampling_regime_2 = sampling_regime), by = "sampling_regime_2")


#generate colors for different groups
pref_pal <- colorRampPalette(brewer.pal(9, "PuBu"))
dist_pal <- colorRampPalette(brewer.pal(9, "YlOrRd"))
closed_pal <- colorRampPalette(brewer.pal(9, "Greens"))

#Make a darker version
pref_pal2 <- colorRampPalette(brewer.pal(9, "PuBu")[4:8])
dist_pal2 <- colorRampPalette(brewer.pal(9, "YlOrRd")[4:9])
closed_pal2 <- colorRampPalette(brewer.pal(9, "Greens")[4:8])


#Show an example of the plot output for one comparison
plot_res[[120]]


#Let's make sure the discrete and continuous versions of the Hellinger distance match -- seems to look good. 
ggplot(res_df, aes(x = hellinger_dist, y = hellinger_dist_discr)) +
  geom_point()+
  geom_abline(intercept = 0, slope = 1)

#Cohen's D vs K-L distance (random vs fishery dependent)
cd_vs_kl_plot <- ggplot(res_df %>% filter(sampling_regime_1 == "random_sampled" | sampling_regime_2 == "random_sampled"), 
                        aes(x = cohens_d, y = kullback_leibler_dist)) +
  geom_vline(xintercept = 0, color = "gray70") + #for reference
  geom_point(aes(color = `Sampling regime`), size = 2) +
  scale_color_manual(values = c("black", closed_pal2(3), dist_pal2(8), pref_pal2(5), "gray70"))+
  #scale_shape_manual(values = c(18, rep(17, 3), rep(15, 8), rep(19, 5), 4))+
  facet_wrap(~ disp_predictor) +
  theme_bw()+
  xlab("Cohen's D")+
  ylab("Kullback-Leibler distance")+
  ggtitle("Randomly sampled vs fishery dependent data")

ggsave(cd_vs_kl_plot, file = "cd_vs_kl_plot.pdf")

#Cohen's D vs Hellinger distance (random vs fishery dependetn)
cd_vs_hd_plot <- ggplot(res_df %>% filter(sampling_regime_1 == "random_sampled" | sampling_regime_2 == "random_sampled"), 
                        aes(x = cohens_d, y = hellinger_dist_discr)) +
  geom_vline(xintercept = 0, color = "gray70") + #for reference
  geom_point(aes(color = `Sampling regime`), size = 2) +
  scale_color_manual(values = c("black", closed_pal2(3), dist_pal2(8), pref_pal2(5), "gray70"))+
  #scale_shape_manual(values = c(18, rep(17, 3), rep(15, 8), rep(19, 5), 4))+
  facet_wrap(~ disp_predictor) +
  theme_bw()+
  xlab("Cohen's D")+
  ylab("Hellinger distance")+
  ggtitle("Randomly sampled vs fishery dependent data")

ggsave(cd_vs_hd_plot, file = "cd_vs_kl_plot.pdf")


#I'm not sure if it will make sense to use both color and shape when we add in size
ggsave(cd_vs_hd_rmse_plot, file = "cd_vs_hd_rmse_plot.pdf", width = 8, height = 8)

cd_vs_kl_rmse_plot <- ggplot(res_df %>% filter(sampling_regime_1 == "random_sampled" | sampling_regime_2 == "random_sampled"), 
                             aes(x = cohens_d, y = kullback_leibler_dist, size = mean_rmse)) +
  geom_vline(xintercept = 0, color = "gray70") + #for reference
  geom_point(aes(color = `Sampling regime`)) +
  scale_color_manual(values = c("black", closed_pal2(3), dist_pal2(8), pref_pal2(5), "gray70"))+
  #scale_shape_manual(values = c(18, rep(17, 3), rep(15, 8), rep(19, 5), 4))+
  facet_wrap(~ disp_predictor) +
  theme_bw()+
  xlab("Cohen's D")+
  ylab("Kullbeck-Leibler distance")+
  ggtitle("Randomly sampled vs fishery dependent data") +
  scale_size_continuous()

ggsave(cd_vs_kl_rmse_plot, file = "cd_vs_kl_rmse_plot.pdf", width = 8, height = 8)

cd_vs_hd_rmse_plot <- ggplot(res_df %>% filter(sampling_regime_1 == "random_sampled" | sampling_regime_2 == "random_sampled"), 
                             aes(x = cohens_d, y = hellinger_dist_discr, size = mean_rmse)) +
  geom_vline(xintercept = 0, color = "gray70") + #for reference
  geom_point(aes(color = `Sampling regime`)) +
  scale_color_manual(values = c("black", closed_pal2(3), dist_pal2(8), pref_pal2(5), "gray70"))+
  #scale_shape_manual(values = c(18, rep(17, 3), rep(15, 8), rep(19, 5), 4))+
  facet_wrap(~ disp_predictor) +
  theme_bw()+
  xlab("Cohen's D")+
  ylab("Hellinger distance")+
  ggtitle("Randomly sampled vs fishery dependent data") +
  scale_size_continuous()


#Save all data histograms
ggsave(
  filename = "data_distribution_plots.pdf", 
  plot = marrangeGrob(plot_res, nrow=1, ncol=1), 
  width = 8, height = 5
)

