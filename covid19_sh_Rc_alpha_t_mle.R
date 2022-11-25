##################
## 2022-10-24 MLE
## Shanghai epidemic
#################
# set the work directory
setwd("/Users/baoyinyuan/Documents/withProTang/script") 
libraries = c("readxl","dplyr","plotly","tidyverse","ggplot2","tidyr","magrittr","lattice",
              "plotrix","R2jags","MASS","splines","tmvtnorm",
              "tidyquant","gridExtra", "dfoptim", "matrixStats", "Matrix", "coda")
for (lib in libraries){
  library(lib, character.only = TRUE)
}

#####
#---- 1. Load the original data of COVID-19 in Shanghai
data0 <- read_excel("/Users/baoyinyuan/Documents/withProTang/script/shanghai_covid19_social.xlsx") #.name_repair = "universal"
data0 %>% head(3)
colnames(data0)

df_sh0 <-  data0 %>% dplyr::select("date", "local_total", "controlled_total", "social_total") %>% as.data.frame %>% na.omit
df_sh0 %>% head(3)
df_sh1 <- df_sh0[nrow(df_sh0):1,]
T0_sh <- as.Date("2022-02-20") # the first epidemic day
T2_sh = as.Date("2022-06-30") # the day of the last reported case
df_sh <- df_sh1 %>% filter(date >= T0_sh & date <= T2_sh)
df_sh %<>% mutate(date = as.Date(date), 
                  day_num = as.numeric(date - T0_sh + 1),
                  local_total = as.numeric(local_total),
                  controlled_total = as.numeric(controlled_total),
                  social_total = as.numeric(social_total)) %>% na.omit
df_sh %>% head(3)
df_sh %>% tail(3)

(T2_num_sh = df_sh$day_num[which(df_sh$date == T2_sh)]) # the day number of the last reported case
df_sh %>% head(10)

(n_sh <- nrow(df_sh))  # Find number of data points in x 

#####
#---- 2. Set the parameters for the serial interval
si.para = c(2.75, 2.53) # (mean, std) # Omicron variant
g.mu = si.para[1]  # mean
g.var = si.para[2]^2 # variance
shap = g.mu^ 2/g.var # shape
rat = g.mu/g.var # rate

# the probability mass function of gamma distribution 
gam = function(t, shap, rat){
  gt =  pgamma(t, shape = shap, rate = rat) - pgamma((t-1),shape = shap, rate = rat)
  return(gt)
}

# right truncated distribution at 12 days
gam_rt = function(t, shap, rat){
  num_elem = length(t)
  gt = rep(0, num_elem)
  total_gam <- sum(gam(1:12, shap, rat))
  for(i in seq(num_elem)){
    if(t[i] > 12){
      gt[i] = 0}else{
        gt<- gam(t, shap, rat)/total_gam 
      }
  }
  return(gt)
}


########
# B-spline model
## To create the B-spline basis functions
bbase <- function(x, xl = min(x), xr = max(x), nseg = 20, deg = 3) {
  # Construct B-spline basis
  dx <- (xr - xl) / nseg
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  return(B)
}

# A function that uses the bs() function to generate the B-spline basis functions
# following Eilers and Marx 'Craft of smoothing' course. This bs_bbase() function
# is equivalent to the bbase() function available at http://statweb.lsu.edu/faculty/marx/

bs_bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3) {
  # Compute the length of the partitions
  dx <- (xr - xl) / nseg
  # Create equally spaced knots
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  # Use bs() function to generate the B-spline basis
  get_bs_matrix <- matrix(bs(x, knots = knots, degree = deg, Boundary.knots = c(knots[1], knots[length(knots)])), nrow = length(x))
  # Remove columns that contain zero only
  bs_matrix <- get_bs_matrix[, -c(1:deg, ncol(get_bs_matrix):(ncol(get_bs_matrix) - deg))]
  
  return(bs_matrix)
}

## Simulate data
set.seed(123)
inc_social_report_sh <- df_sh$social_total
inc_controlled_sh <- df_sh$controlled_total

n <- n_sh  # Find number of datapoints in x 
x <- 1:n
B <- bs_bbase(x, nseg = 8)

######
# Use MLE to estimate alpha_t and Rc
######
dim(B)
N_knots <- 11
param <- c(1.5, rep(0.2, N_knots))
Negloglik <- function(param){
  par_num <- length(param)
  Rc <- param[1]
  beta <- param[2:par_num]
  alpha_t <- rep(0, n_sh)
  inc_control_social <- rep(0, n_sh)
  inc_control_exp <- rep(0.1, n_sh)
  for(t in 2 : n_sh){
    alpha_t[t] <- exp(sum(B[t,]*beta))
    inc_control_social[t] <- inc_social_report_sh[t]*alpha_t[t]
  }
  
  sum_iter <- rep(0, 12)
  for(t in 2 : n_sh){
    for (i in (t-1) : max(t-12, 1)){
      sum_iter[t-i] <- inc_controlled_sh[i] * gam_rt(t-i, shap, rat)
      }
    inc_control_exp[t] <- max(inc_control_social[t] + Rc * sum(sum_iter[]), 0.1)  # mu_t
  } 
  
  # likelihood based two poisson distributions
  (negloglik = sum(- inc_controlled_sh * log(inc_control_exp[1:n_sh]) + inc_control_exp[1:n_sh]))
  aic = 2*par_num + 2*negloglik
  aicc = aic + (2*par_num^2 + 2*par_num)/(n_sh - par_num - 1)
  
  return(aicc) # Return the AICc
}


######
param0 <-  c(0.9, rep(0.01, N_knots))
dfopt_lps <- optim(
  par = param0,
  fn = Negloglik,
  method = "Nelder-Mead",
  hessian = TRUE,
  control = list(maxit = 5000)
)
dfopt_lps
dfopt_lps$par
dfopt_lps$value

#######
#*---Item (1): calculate confidence interval of the estimated parameters 
#*----------   by the bootstrapping method and rejection re-sampling
(esti_para <- dfopt_lps$par)
(tot_num_para <- length(param0))
hessian_matrix_sh = dfopt_lps$hessian
sigmat<-solve(hessian_matrix_sh[, 1:tot_num_para])
sigmat2<-forceSymmetric(sigmat,uplo="L")
nsample = 100
par_ini_low <- c(0.1, rep(-10, N_knots)) 
par_ini_upp <- c(5, rep(10, N_knots))
hess_sam_sh <- rtmvnorm(n = nsample, 
                     mean = as.vector(esti_para), 
                     sigma = sigmat2,
                     lower = par_ini_low,
                     upper = par_ini_upp,
                     algorithm= "rejection")
hess_sam_sh
# write.matrix(hess_sam_sh, file = "hess_sam_sh_2SI.csv")
hess_sam_sh2 = read.csv("hess_sam_sh.csv", header = F, sep='')
hess_sam_sh3 <- as.matrix(hess_sam_sh2) 


ci_value = c(0.025, 0.5, 0.975)
ci_name = c("Lower",'Median', "Upper")
nci = length(ci_value)

##*---Item (2): calculate CI for alpha_t
alpha_t.mat <- matrix(0, ncol = nsample, nrow = n_sh)
for(i in seq(nsample)){
  beta_est <- hess_sam_sh[i, 2:(N_knots+1)]
  alpha_t <- rep(0, n_sh)
  for(t in 2 : n_sh){
    alpha_t.mat[t, i] <- exp(sum(B[t,]*beta_est))
  }
}


ci_alpha = data.frame(matrix(0, nrow = n_sh, ncol = nci))
colnames(ci_alpha) = ci_name
for(i in seq(n_sh)){
  ci_alpha[i, ] = as.numeric(quantile(alpha_t.mat[i,], probs = ci_value))
}
ci_alpha

plot(x= seq(n_sh) ,ci_alpha[,1], type = "p", col = "red", ylim=c(0, 10))
lines(ci_alpha[,2], type = "p", col = "green")
lines(ci_alpha[,3], type = "p", col = "blue")

##*---Item (3) Method(i): transform ci_alpha to ci_beta, i.e., proportion of I_t^{sr} into I_t{s}
ci_beta <- data.frame(matrix(0, nrow = n_sh, ncol = nci))
for(j in seq(nci)){
  for(t in 2 : n_sh){
    ci_beta[t,j] <- 1/(1 + as.numeric(ci_alpha[t, j]))
  }
}
colnames(ci_beta) = c("Upper", "Median", "Lower")

plot(ci_beta[, 1], type = "p", col = "red")
lines(ci_beta[, 2], type = "p", col = "green")
lines(ci_beta[, 3], type ="p", col = "blue")

##*---Item (4): calculate CI for Rc as const
ci_Rc = matrix(0, nrow = 1, ncol = nci)
colnames(ci_Rc) = ci_name
for(i in seq(n_sh)){
  ci_Rc[1,] = as.numeric(quantile(hess_sam_sh[,1], probs = ci_value))
}
ci_Rc

colnames(ci_Rc) <- c("2.5%", "50%", "97.5%")
rownames(ci_Rc) <- c("Rc ")
(ci_Rc.table <- as.table(round(ci_Rc, 3)))

#####

##*---Item (5): calculate CI for Rs_t
hess_sam_sh

#*---Item (5.1) calculate alpha_t.mat based on the bootstrapped parameters
#* this step is the replication of Item(2), but not calculate CI for alpha_t
alpha_t.mat <- matrix(0, ncol = nsample, nrow = n_sh)
for(i in seq(nsample)){
  beta_est <- hess_sam_sh[i, 2:(N_knots+1)]
  alpha_t <- rep(0, n_sh)
  for(t in 2 : n_sh){
    alpha_t.mat[t, i] <- exp(sum(B[t,]*beta_est))
  }
}
alpha_t.mat
dim(alpha_t.mat)
#*---Item (5.2): calculate Rc_t based on the bootstrapped parameters
#* Each element of Rc_t here corresponds to the above-calculated alpha_t.mat by column 
Rc_t.mat <- hess_sam_sh[, 1]
Rc_t.mat

#*---Item (5.3): reconstruct the epidemic curve by each combination of Rc_t and alpha_t 
#*               which are calculated by the bootstrapping 
inc_control_exp.mat <- matrix(0, nrow = n_sh, ncol = nsample)
inc_control_social.mat <- matrix(0, nrow = n_sh, ncol = nsample)
sum_iter <- rep(0, 12)
for(j in seq(nsample)){
  # j = 1
  for(t in 2 : n_sh){
    inc_control_social.mat[t, j] <- inc_social_report_sh[t] * alpha_t.mat[t,j]
    for (i in (t-1):max(t-12, 1)){
      sum_iter[t-i] <- inc_controlled_sh[i] * gam_rt(t-i, shap, rat)
    }
    inc_control_exp.mat[t, j] <- inc_control_social.mat[t, j] + Rc_t.mat[j] * sum(sum_iter[])  # mu_t
  }
}


#*---Item (5.4): calculate CI of reconstructed epidemic curve, i.e., the expected incidence in control
ci_value = c(0.025, 0.5, 0.975)
ci_name = c("Lower",'Median', "Upper")
nci = length(ci_value)

ci_inc_control_exp_sh = data.frame(matrix(NA, nrow = n_sh, ncol = nci))
colnames(ci_inc_control_exp_sh) = ci_name
for(i in (1:n_sh)){
  ci_inc_control_exp_sh[i, ] = as.numeric(quantile(inc_control_exp.mat[i,], probs = ci_value))
}
ci_inc_control_exp_sh

plot(inc_controlled_sh, type = 's', col = "black")
lines(ci_inc_control_exp_sh[,1], type = 's', col = "red")
lines(ci_inc_control_exp_sh[,2], type = 's', col = "green")
lines(ci_inc_control_exp_sh[,3], type = 's', col = "blue")


#*---Item (5.5): use reconstructed epidemic curve to calculate the Rs_t with CI
inc_social_sh.mat <- inc_control_social.mat + matrix(replicate(nsample, inc_social_report_sh), nrow = n_sh)
Rs_t_sh.mat <- matrix(NA, nrow = n_sh, ncol = nsample)
Rs_t_movaverage_sh.mat <- matrix(NA, nrow = n_sh, ncol = nsample)
for(j in seq(nsample)){
  for(t in 2 : n_sh){
    renew <- inc_social_report_sh[(t-1) : 1] * gam_rt(1 : (t-1), shap, rat) 
    Rs_t_sh.mat[t, j] <- inc_social_sh.mat[t, j]/sum(renew)
  } 
  # moving average
  Rs_t_movaverage_sh.mat[2:n_sh,j]<- rollmean(Rs_t_sh.mat[2:n_sh, j], k = 7, fill = NA, align = "right") # moving average over 7 days
}

plot(Rs_t_movaverage_sh.mat[, 1], ylim = c(0, 15), type = "l")
abline(h = 1, col = "red") 

ci_value = c(0.025, 0.5, 0.975)
ci_name = c("Lower",'Median', "Upper")
nci = length(ci_value)

ci_Rs_t_sh = data.frame(matrix(NA, nrow = n_sh, ncol = nci))
colnames(ci_Rs_t_sh) = ci_name
for(i in (17:n_sh)){
  ci_Rs_t_sh[i, ] = as.numeric(quantile(Rs_t_movaverage_sh.mat[i,], probs = ci_value))
}
ci_Rs_t_sh

plot(ci_Rs_t_sh[,1], type = "l", col = "red", ylim = c(0, 20))
lines(ci_Rs_t_sh[,2], type = "l", col = "green")
lines(ci_Rs_t_sh[,3], type = "l", col = "blue")
abline(h = 1, col = "black")

######
#*---Item (6): visualization of all the above results

#*---Item (6.1): visualization of estimated parameters: Rc as const and beta_t in ts

ci_Rc.table

ci_beta$date <- df_sh$date
ci_beta[1:10,1:3] <- NA
ci_beta_plt <- ci_beta %>% filter(date < as.Date("2022-06-01")) %>% filter(date > as.Date("2022-03-09"))
ci_beta %>% head(15)

ci_inc_control_exp_sh$date <- df_sh$date
ci_inc_control_exp_sh_plt <- ci_inc_control_exp_sh %>% filter(date < as.Date("2022-06-01")) %>% filter(date > as.Date("2022-03-09"))
ci_inc_control_exp_sh %>% head(3)

df_sh
df_sh %>% head(3)
df_sh_plt <- df_sh %>% filter(date < as.Date("2022-06-01")) %>% filter(date > as.Date("2022-03-09")) 
df_sh_plt

color.label = c("Predicted case in control" = "red", 
                "Reported case in control" = "blue", 
                "Reported case in society" = "green")

Rc_beta_est_sh <- ggplot() +
  geom_step(data = df_sh_plt, aes(x= date, y = controlled_total, colour = "Reported case in control"), linetype="solid",size=0.6) +
  geom_step(data = df_sh_plt, aes(x= date, y = social_total, colour = "Reported case in society"), linetype="solid",size=0.6) +
  RcmdrPlugin.KMggplot2::geom_stepribbon(data = ci_inc_control_exp_sh_plt, aes(x= date, ymin = round(Lower), ymax = round(Upper)), fill = "violet") +
  geom_step(data = ci_inc_control_exp_sh_plt, aes(x= date, y = round(Median), colour = "Predicted case in control"), linetype="solid", size=0.6) + 
  geom_ribbon(data = ci_beta_plt, mapping = aes(x = date, ymin = Lower*30000 ,ymax = Upper*30000), fill = "gray60", alpha = 0.4) + 
  geom_line(data = ci_beta_plt, mapping=aes(x = date, y = Median*30000), colour = "black", linetype="solid",size=0.8) +
  scale_color_manual(values = color.label) +
  scale_x_date(date_breaks="7 day", date_labels="%Y-%m-%d") +
  scale_y_continuous(expand = c(0, 0), 
                     breaks=seq(0, 30000, by = 5000),
                     limits = c(0, 30000),
                     sec.axis = sec_axis(~./30000, name = expression(beta(t)), breaks=seq(0, 1, 0.1))) +
  xlab("Date") + ylab("Number of cases") +
  ggtitle("Shanghai") +
  guides(color=guide_legend("Case type")) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(size=10, angle = 90, colour = "black", vjust=0.5, hjust = 1),
        axis.text.y=element_text(size=10, colour = "black"),
        axis.line.y.right = element_line(color = "black"),
        axis.title=element_text(size=12),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "grey95"),
        plot.title = element_text(size = 12,face="bold"),
        legend.background = element_rect(fill="transparent"),
        legend.box.background = element_rect(colour = "black", fill = "transparent"),
        legend.title=element_text(size= 10),
        legend.text=element_text(size= 10),
        legend.position = c(0.85, 0.86),
        # legend.justification = c("right", "top"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  annotation_custom(tableGrob(ci_Rc.table), xmin = as.Date("2021-12-30"), ymin = 25000) 

Rc_beta_est_sh


ggsave(filename = "../figure/Rc_beta_est_sh_275_253.tiff",
       plot = Rc_beta_est_sh,
       width = 25,  # <=19.05cm
       height = 15, # <=22.225cm
       units= "cm",
       dpi= 330,
       compression = "lzw")


#*---Item (6.2): visualization of calculated Rs_t
ci_Rs_t_sh
ci_Rs_t_sh$date <- df_sh$date
ci_Rs_t_sh %>% head(5)
ci_Rs_t_sh_plt <- ci_Rs_t_sh %>% filter(date < as.Date("2022-06-01")) %>% filter(date > as.Date("2022-03-15")) 
ci_Rs_t_sh_plt %>% head(5)


df_sh
df_sh %>% head(3)
df_sh_plt <- df_sh %>% filter(date < as.Date("2022-06-01")) %>% filter(date > as.Date("2022-03-09")) 
df_sh_plt

inc_sh <- data.frame(matrix(nrow = n_sh, ncol = 0))
inc_sh$inc_control <- df_sh$controlled_total
inc_sh$inc_social <- df_sh$social_total
inc_sh$date <- df_sh$date
inc_stackbar_sh <- inc_sh %>% gather(case_type, case_num, -date)
inc_stackbar_sh_plt <- inc_stackbar_sh %>% filter(date < as.Date("2022-06-01")) %>% filter(date > as.Date("2022-03-09"))  


Rs_inc_stack_sh  <- ggplot() +
  geom_bar(data = inc_stackbar_sh_plt, aes(y = case_num, x = date, fill = case_type), position = "stack", stat = "identity") + 
  geom_ribbon(data = ci_Rs_t_sh_plt, mapping = aes(x = date, ymin = Lower*2500 ,ymax = Upper*2500), fill = "gray30", alpha = 0.5) + 
  geom_line(data = ci_Rs_t_sh_plt, mapping=aes(x = date, y = Median*2500), colour = "black", linetype="solid", size=0.6) +
  geom_hline(yintercept = 1*2500, linetype = "dashed", col = "red", size = 1) +
  scale_fill_manual(
    labels = c("Reported case in control", "Reported case in society"),
    values = c("gray80", "gray60"), 
    name = "Case type", 
    guide = guide_legend(reverse = TRUE)) +
  scale_x_date(date_breaks="7 day", date_labels="%Y-%m-%d") +
  scale_y_continuous(expand = c(0, 0), 
                     breaks=seq(0, 30000, by = 5000),
                     limits = c(0, 30000),
                     sec.axis = sec_axis(~./2500, name="Rs(t)", breaks=seq(0, 15, 1))) +
  xlab("Date") + ylab("Number of cases") +
  ggtitle("Shanghai") +
  guides(color=guide_legend("Incidence type")) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(size=10, angle = 90, colour = "black", vjust=0.5, hjust = 1),
        axis.text.y=element_text(size=10, colour = "black"),
        axis.line.y.right = element_line(color = "black"),
        axis.title=element_text(size=12),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "grey95"),
        plot.title = element_text(size = 12,face="bold"),
        legend.background = element_rect(fill="transparent"),
        legend.box.background = element_rect(colour = "black", fill = "transparent"),
        legend.title=element_text(size= 10),
        legend.text=element_text(size= 10),
        legend.position = c(0.85, 0.86),
        # legend.justification = c("right", "top"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
Rs_inc_stack_sh

ggsave(filename = "../figure/Rs_inc_stack_sh_275_253.tiff",
       plot = Rs_inc_stack_sh,
       width = 25,  # <=19.05cm
       height = 15, # <=22.225cm
       units= "cm",
       dpi= 330,
       compression = "lzw")



