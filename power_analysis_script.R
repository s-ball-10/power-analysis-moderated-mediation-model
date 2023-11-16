library('paramtest') 
library('pwr') 
library('ggplot2') 
library('knitr') 
library('dplyr')
library('mosaic') 
library('lavaan')
library('xtable')
library('stargazer')

## Power Analysis for Moderated Mediation Model
## Code based on: https://rdrr.io/github/jeff-hughes/paramtest/f/vignettes-build/Simulating-Power-build.Rmd

# Function to simulate data and to estimate model based on simulated data
moderated_mediation <- function(simNum, N, a1, a2, a3, a4, a5, a0=0, bb, cc) {
  
  # Simulating data
  signal <- scale(sample(0:1, N, replace= TRUE)) # signal
  x2 <- sample(0:2, N, replace= TRUE) # gender
  female <- scale(ifelse(x2 ==1, 1, 0))
  non_binary <- scale(ifelse(x2 ==2, 1, 0))
  mvar <- sqrt(1 - (a1^2) - (a2^2) - (a3^2) - (a4^2) - (a5^2)) # residual variance mediator
  competence <- rnorm(N, a0 + a1*signal + a2*female + a3*non_binary + a4*signal*female + a5*signal*non_binary, mvar)
  will_var <- sqrt(1-(cc^2)-(bb^2)) # residual variance of y variable 
  will <- rnorm(N, cc*signal + bb*competence, will_var)
  data <- data.frame(signal, female,non_binary, competence, will)
  data$signal_female <- data$signal * data$female
  data$signal_non_binary <- data$signal* data$non_binary
  
  # Estimating model based on simulated data 
  model <- '
  competence ~ a1*signal + a2*female + a3*non_binary + a4*signal_female + a5*signal_non_binary
  will ~ cc*signal
  will ~ bb*competence
  ab_male := a1*bb
  ab_female := (a1+a4)*bb
  ab_non_binary := (a1+a5)*bb
  diff_male_female := ab_male - ab_female
  diff_male_non_binary := ab_male - ab_non_binary
  diff_female_non_binary :=  ab_female - ab_non_binary
  '
  fit <- lavaan::sem(model, data=data)
  ests <- lavaan::parameterEstimates(fit)
  
  # pull output from model (main effects and interactions)
  a1_est <- ests[ests$label == 'a1', 'est']
  a1_p <- ests[ests$label == 'a1', 'pvalue']
  a1_sig <- a1_p < 0.05
  
  a2_est <- ests[ests$label == 'a2', 'est']
  a2_p <- ests[ests$label == 'a2', 'pvalue']
  a2_sig <- a2_p < 0.05
  
  a3_est <- ests[ests$label == 'a3', 'est']
  a3_p <- ests[ests$label == 'a3', 'pvalue']
  a3_sig <- a3_p < 0.05
  
  a4_est <- ests[ests$label == 'a4', 'est']
  a4_p <- ests[ests$label == 'a4', 'pvalue']
  a4_sig <- a4_p < 0.05
  
  a5_est <- ests[ests$label == 'a5', 'est']
  a5_p <- ests[ests$label == 'a5', 'pvalue']
  a5_sig <- a5_p < 0.05
  
  b_est <- ests[ests$label == 'bb', 'est']
  b_p <- ests[ests$label == 'bb', 'pvalue']
  b_sig <- b_p < 0.05
  
  c_est <- ests[ests$label == 'cc', 'est']
  c_p <- ests[ests$label == 'cc', 'pvalue']
  c_sig <- c_p < 0.05
  
  ab_male_est <- ests[ests$label == 'ab_male', 'est']
  ab_male_p <- ests[ests$label == 'ab_male', 'pvalue']
  ab_male_sig <- ab_male_p < 0.05
  
  ab_female_est <- ests[ests$label == 'ab_female', 'est']
  ab_female_p <- ests[ests$label == 'ab_female', 'pvalue']
  ab_female_sig <- ab_female_p < 0.05
  
  ab_non_binary_est <- ests[ests$label == 'ab_non_binary', 'est']
  ab_non_binary_p <- ests[ests$label == 'ab_non_binary', 'pvalue']
  ab_non_binary_sig <- ab_non_binary_p < 0.05
  
  diff_male_female_est <- ests[ests$label == 'diff_male_female', 'est']
  diff_male_female_p <- ests[ests$label == 'diff_male_female', 'pvalue']
  diff_male_female_sig <- diff_male_female_p < 0.05
  
  diff_male_non_binary_est <- ests[ests$label == 'diff_male_non_binary', 'est']
  diff_male_non_binary_p <- ests[ests$label == 'diff_male_non_binary', 'pvalue']
  diff_male_non_binary_sig <- diff_male_non_binary_p < 0.05
  
  diff_female_non_binary_est <- ests[ests$label == 'diff_female_non_binary', 'est']
  diff_female_non_binary_p <- ests[ests$label == 'diff_female_non_binary', 'pvalue']
  diff_female_non_binary_sig <- diff_female_non_binary_p < 0.05
  
  
  return(c(a1_est=a1_est, a1_p=a1_p, a1_sig=a1_sig,
           a2_est=a2_est, a2_p=a2_p, a2_sig=a2_sig,
           a3_est=a3_est, a3_p=a3_p, a3_sig=a3_sig,
           a4_est=a4_est, a4_p=a4_p, a4_sig=a4_sig,
           a5_est=a5_est, a5_p=a5_p, a5_sig=a5_sig,
           b_est=b_est,b_p=b_p, b_sig=b_sig,
           c_est=c_est, c_p=c_p, c_sig=c_sig,
           ab_male_est=ab_male_est, ab_male_p=ab_male_p, ab_male_sig= ab_male_sig,
           ab_female_est=ab_female_est, ab_female_p=ab_female_p, ab_female_sig= ab_female_sig,
           ab_non_binary_est=ab_non_binary_est, ab_non_binary_p=ab_non_binary_p, ab_non_binary_sig= ab_non_binary_sig,
           diff_male_female_est=diff_male_female_est, diff_male_female_p = diff_male_female_p, diff_male_female_sig = diff_male_female_sig,
           diff_male_non_binary_est=diff_male_non_binary_est, diff_male_non_binary_p=diff_male_non_binary_p, diff_male_non_binary_sig=diff_male_non_binary_sig,
           diff_female_non_binary_est=diff_female_non_binary_est, diff_female_non_binary_p=diff_female_non_binary_p, diff_female_non_binary_sig=diff_female_non_binary_sig))
}


## Power analysis/simulations

# Determine expected effect sizes
full <- expand.grid(c(0.5),c(-0.2),c(-0.25), c(-0.1), c(-0.15), c(0.4), c(0.25))
full

# Run simulations
correct_moderation_2 <- purrr::pmap(.l = list(a=full$Var1,b=full$Var2,c=full$Var3,d=full$Var4, e=full$Var5, f=full$Var6, g=full$Var7), 
                       .f = function(a,b,c,d,e,f,g){
                         power_med_int <- grid_search(moderated_mediation, 
                                                      params=list(N=seq(100,3500,100)),
                                                      n.iter=1000, output='data.frame', 
                                                      a1=a,
                                                      a2=b,
                                                      a3=c,
                                                      a4=d,
                                                      a5=e,
                                                      a0=0,
                                                      bb=f,
                                                      cc=g,
                                                      parallel='snow',
                                                      ncpus=4)
                         
                         results(power_med_int)  %>% 
                           group_by(N.test) %>%
                           summarise(
                             power_signal=mean(a1_sig), 
                             power_female=mean(a2_sig),
                             power_non_binary = mean(a3_sig),
                             power_sig_fem = mean(a4_sig),
                             power_sig_non = mean(a5_sig),
                             power_comp = mean(b_sig),
                             power_sig_will = mean(c_sig),
                             power_ab_male=mean(ab_male_sig),
                             power_ab_female=mean(ab_female_sig),
                             power_ab_non_binary=mean(ab_non_binary_sig),
                             power_diff_male_female=mean(diff_male_female_sig),
                             power_diff_male_non_binary=mean(diff_male_non_binary_sig),
                             power_diff_female_non_binary=mean(diff_female_non_binary_sig))
                       }) %>% bind_rows()

print(xtable(correct_moderation_2, type = "latex"), file = "power_table.tex")

#selected_columns <- correct_moderation_2[,c('N.test','power_sig_fem', 'power_diff_male_female', 'power_diff_male_non_binary', 'power_diff_female_non_binary')]
selected_columns <- correct_moderation_2[,c('N.test','power_diff_male_female', 'power_diff_male_non_binary', 'power_diff_female_non_binary')]

out_long <- selected_columns %>% tidyr::gather(key,value,contains("power"))


effect_size_1 <- 0.04
effect_size_2 <- 0.06
effect_size_3 <- 0.02



out_long <- out_long #%>%
  #mutate(lab = paste0("ab(female) - ab(non-binary) = ", effect_size_3,"; ab(male) - ab(female) = ",effect_size_1,"; ab(male) - ab(non-binary) = ",effect_size_2, "; ab(female) = ",-0.1))

gg2 <-
  ggplot(out_long,aes(x = N.test,y = value, color = key )) +
  geom_line(size = 1.3,alpha = 0.85) +
  #geom_line(aes(color = factor(key)), stat = "identity", size = 1,alpha = 0.85) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0.8,linetype = "dashed",color = "#FF3333",size = 1.2)+
  #geom_segment(aes(x=100,xend=3500,y=0.8, yend = 0.81), linetype = "dashed",color = "red")+
  #facet_wrap(~ lab,ncol = 2,scales = "fixed") +
  #scale_color_viridis_d(labels = c("Interaction signal female","Diff. ind. effect male female","Diff. ind. effect male non-binary", "Diff. ind. effect female non-binary")) + 
  #scale_color_brewer(labels = c("Diff. ind. effect female non-binary", "Diff. ind. effect male female","Diff. ind. effect male non-binary", "Interaction signal female"), palette = "RdYlBu", name = "Key: ") + 
  scale_color_manual(#labels = c("Diff. ind. effect female non-binary (ab = 0.02)", "Diff. ind. effect male female (ab = 0.04)","Diff. ind. effect male non-binary (ab = 0.06)", "Interaction signal female (b = -0.1)"),
                     labels = c("Difference ind. effects female non-binary (ab = 0.02)", "Difference ind. effects male female (ab = 0.04)","Difference ind. effects male non-binary (ab = 0.06)"), 
                     #values = c("#F4A582", "#4393C3", "#92C5DE","#993366"),name = "Key: ") +
                     #values = c("#FF9933", "#993366", "#4393C3", "#99CC00"),name = "Key: ") +
                     values = c("#FF9933", "#99CC00", "#4393C3"),name = "Key: ") +
  
  labs(title = NULL,
       subtitle = NULL,
       x = "\nNumber of observations",
       y = "Statistical power\n")+
  scale_x_continuous(breaks = seq(0,3500,200)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        plot.title = element_text(color = "black"),
        plot.subtitle = element_text(color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(),
        strip.text = element_text(size = 12), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))

gg2
ggsave("power_plot.pdf", plot=gg2, width=15, height=8) 
  
  

test_power_628 <- purrr::pmap(.l = list(a=full$Var1,b=full$Var2,c=full$Var3,d=full$Var4, e=full$Var5, f=full$Var6, g=full$Var7), 
                                    .f = function(a,b,c,d,e,f,g){
                                      power_med_int <- grid_search(moderated_mediation, 
                                                                   params=list(N=628),
                                                                   n.iter=1000, output='data.frame', 
                                                                   a1=a,
                                                                   a2=b,
                                                                   a3=c,
                                                                   a4=d,
                                                                   a5=e,
                                                                   a0=0,
                                                                   bb=f,
                                                                   cc=g,
                                                                   parallel='snow',
                                                                   ncpus=4)
                                      
                                      results(power_med_int)  %>% 
                                        group_by(N.test) %>%
                                        summarise(
                                          power_signal=mean(a1_sig), 
                                          power_female=mean(a2_sig),
                                          power_non_binary = mean(a3_sig),
                                          power_sig_fem = mean(a4_sig),
                                          power_sig_non = mean(a5_sig),
                                          power_comp = mean(b_sig),
                                          power_sig_will = mean(c_sig),
                                          power_ab_male=mean(ab_male_sig),
                                          power_ab_female=mean(ab_female_sig),
                                          power_ab_non_binary=mean(ab_non_binary_sig),
                                          power_diff_male_female=mean(diff_male_female_sig),
                                          power_diff_male_non_binary=mean(diff_male_non_binary_sig),
                                          power_diff_female_non_binary=mean(diff_female_non_binary_sig))
                                    }) %>% bind_rows()


test_power_527 <- purrr::pmap(.l = list(a=full$Var1,b=full$Var2,c=full$Var3,d=full$Var4, e=full$Var5, f=full$Var6, g=full$Var7), 
                              .f = function(a,b,c,d,e,f,g){
                                power_med_int <- grid_search(moderated_mediation, 
                                                             params=list(N=527),
                                                             n.iter=1000, output='data.frame', 
                                                             a1=a,
                                                             a2=b,
                                                             a3=c,
                                                             a4=d,
                                                             a5=e,
                                                             a0=0,
                                                             bb=f,
                                                             cc=g,
                                                             parallel='snow',
                                                             ncpus=4)
                                
                                results(power_med_int)  %>% 
                                  group_by(N.test) %>%
                                  summarise(
                                    power_signal=mean(a1_sig), 
                                    power_female=mean(a2_sig),
                                    power_non_binary = mean(a3_sig),
                                    power_sig_fem = mean(a4_sig),
                                    power_sig_non = mean(a5_sig),
                                    power_comp = mean(b_sig),
                                    power_sig_will = mean(c_sig),
                                    power_ab_male=mean(ab_male_sig),
                                    power_ab_female=mean(ab_female_sig),
                                    power_ab_non_binary=mean(ab_non_binary_sig),
                                    power_diff_male_female=mean(diff_male_female_sig),
                                    power_diff_male_non_binary=mean(diff_male_non_binary_sig),
                                    power_diff_female_non_binary=mean(diff_female_non_binary_sig))
                              }) %>% bind_rows()

