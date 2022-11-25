##################
## 2022-10-24 MLE
## Estimate effective reproduciton number following Scenario (a)
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
# T2_sh = as.Date("2022-05-17") # the day of the last reported case
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

#####
#---- 2. Set the parameters for the serial interval
si.para = c(3.42, 7.93)
g.mu = si.para[1]  # mean
g.var = si.para[2]^2 # variance
shap = g.mu^ 2/g.var # shape
rat = g.mu/g.var # rate

# the probability mass function of gamma distribution 
gam = function(t, shap, rat){
  gt =  pgamma(t, shape = shap, rate = rat) - pgamma((t-1),shape = shap, rate = rat)
  return(gt)
}

#####
#---- 3. Estimate Rt by using only social cases
inc_social_report_sh <- df_sh$social_total
inc_controlled_sh <- df_sh$controlled_total
inc_total <- inc_social_report_sh + inc_controlled_sh
(n_sh <- nrow(df_sh))  # Find number of data points in x 

#####
#*---- 3.1 Scenario(1): no internal infection in control; all the infections occur in society;
#*--- $I_t^{sc}$ are not infectious, i.e., all the cases are caused
#*--- by the reported case in society $I_t^{sr}$.
Rt_sh = rep(NA, n_sh)
Rt_movaverage_sh = rep(NA, n_sh)
for(t in 2 : n_sh){
  renew <- inc_social_report_sh[(t-1) : 1] * gam(1 : (t-1), shap, rat) 
  Rt_sh[t] <- inc_total[t]/sum(renew)
}
# moving average
Rt_movaverage_sh[2:n_sh]<- rollmean(Rt_sh[2:n_sh], k = 7, fill = NA, align = "right") # moving average over 7 days

plot(Rt_movaverage_sh, type = "l", ylim = c(0,60))
abline(h= 1, col="red")

Rt_sh_df <- data.frame(Rt = Rt_movaverage_sh, date= df_sh$date)
Rt_sh_df %>% head(5)
Rt_sh_plt <- Rt_sh_df %>% filter(date < as.Date("2022-06-01")) %>% filter(date > as.Date("2022-03-15")) 
Rt_sh_plt %>% head(5)


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

Rt_inc_stack_sh  <- ggplot() +
  geom_bar(data = inc_stackbar_sh_plt, aes(y = case_num, x = date, fill = case_type), position = "stack", stat = "identity") + 
  geom_line(data = Rt_sh_plt, mapping=aes(x = date, y = Rt*500), colour = "black", linetype="solid", size=0.6) +
  geom_hline(yintercept = 1*500, linetype = "dashed", col = "red", size = 1) +
  scale_fill_manual(
    labels = c("Reported case in control", "Reported case in society"),
    values = c("gray80", "gray60"), 
    name = "Case type", 
    guide = guide_legend(reverse = TRUE)) +
  scale_x_date(date_breaks="7 day", date_labels="%Y-%m-%d") +
  scale_y_continuous(expand = c(0, 0), 
                     breaks=seq(0, 30000, by = 5000),
                     limits = c(0, 30000),
                     sec.axis = sec_axis(~./500, name="R(t) when no internal infections in control",
                                         breaks=seq(0, 60, 5))) +
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
        legend.position = c(0.86, 0.88),
        # legend.justification = c("right", "top"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
Rt_inc_stack_sh

ggsave(filename = "../figure/Rt_inc_stack_sh.tiff",
       plot = Rt_inc_stack_sh,
       width = 25,  # <=19.05cm
       height = 15, # <=22.225cm
       units= "cm",
       dpi= 330,
       compression = "lzw")

#
################

#---- 1. Load the original data of COVID-19 in Xi'an
data0_xian <- read_excel("/Users/baoyinyuan/Documents/withProTang/script/Xianguankong.xlsx") #.name_repair = "universal"
data0_xian %>% head(3)
colnames(data0_xian)

df_xian <-  data0_xian %>% dplyr::select("date", "local_total", "controlled_total", "social_total") %>% as.data.frame %>% na.omit
df_xian %>% head(3)
df_xian %>% tail(3)

(day_ini <- as.Date(df_xian$date[1]))
df_xian %<>% mutate(date = as.Date(date), 
                    day_num = as.numeric(date - day_ini + 1),
                    local_total = as.numeric(local_total),
                    controlled_total = as.numeric(controlled_total),
                    social_total = as.numeric(social_total)) %>% na.omit
df_xian %>% head(3)
df_xian %>% tail(3)

T2_xian = as.Date(df_xian$date[nrow(df_xian)]) # the day of the last reported case
(T2_num_xian = df_xian$day_num[which(df_xian$date == T2_xian)]) # the day number of the last reported case
df_xian %>% head(10)
(n_xian <- nrow(df_xian))  # Find number of data points in x 

#####
#---- 2. Set the parameters for the serial interval
si.para = c(3.00, 2.48) # Delta variant
g.mu = si.para[1]  # mean
g.var = si.para[2]^2 # variance
shap = g.mu^ 2/g.var # shape
rat = g.mu/g.var # rate

# the probability mass function of gamma distribution 
gam = function(t, shap, rat){
  gt =  pgamma(t, shape = shap, rate = rat) - pgamma((t-1),shape = shap, rate = rat)
  return(gt)
}


#####
#---- 3. Estimate Rt by using only social cases
inc_social_report_xian <- df_xian$social_total
inc_controlled_xian <- df_xian$controlled_total
inc_total <- inc_social_report_xian + inc_controlled_xian
(n_xian <- nrow(df_xian))  # Find number of data points in x 

#####
#*---- 3.1 Scenario(1): no internal infection in control; all the infections occur in society;
#*--- $I_t^{sc}$ are not infectious, i.e., all the cases are caused
#*--- by the reported case in society $I_t^{sr}$.
Rt_xian = rep(NA, n_xian)
Rt_movaverage_xian = rep(NA, n_xian)
for(t in 2 : n_xian){
  renew <- inc_social_report_xian[(t-1) : 1] * gam(1 : (t-1), shap, rat) 
  Rt_xian[t] <- inc_total[t]/sum(renew)
}
# moving average
Rt_movaverage_xian[2:n_xian]<- rollmean(Rt_xian[2:n_xian], k = 5, fill = NA, align = "right") # moving average over 7 days

plot(Rt_movaverage_xian, type = "l", ylim = c(0,12))
abline(h=1, col="red")

Rt_xian_df <- data.frame(Rt = Rt_movaverage_xian, date= df_xian$date)
Rt_xian_df %>% head(10)
Rt_xian_plt <- Rt_xian_df %>% filter(date < as.Date("2022-01-17")) %>% filter(date > as.Date("2021-12-09")) 
Rt_xian_plt %>% head(5)


df_xian
df_xian %>% head(3)
df_xian_plt <- df_xian %>% filter(date < as.Date("2022-01-17")) %>% filter(date > as.Date("2021-12-09")) 
df_xian_plt

inc_xian <- data.frame(matrix(nrow = n_xian, ncol = 0))
inc_xian$inc_control <- df_xian$controlled_total
inc_xian$inc_social <- df_xian$social_total
inc_xian$date <- df_xian$date
inc_stackbar_xian <- inc_xian %>% gather(case_type, case_num, -date)
inc_stackbar_xian_plt <- inc_stackbar_xian %>% filter(date < as.Date("2022-01-17")) %>% filter(date > as.Date("2021-12-09"))  


Rt_inc_stack_xian  <- ggplot() +
  geom_bar(data = inc_stackbar_xian_plt, aes(y = case_num, x = date, fill = case_type), position = "stack", stat = "identity") + 
  geom_line(data = Rt_xian_plt, mapping=aes(x = date, y = Rt*10), colour = "black", linetype="solid", size=0.6) +
  geom_hline(yintercept = 1*10, linetype = "dashed", col = "red", size = 1) +
  scale_fill_manual(
    labels = c("Reported case in control", "Reported case in society"),
    values = c("gray80", "gray60"), 
    name = "Case type", 
    guide = guide_legend(reverse = TRUE)) +
  scale_x_date(date_breaks="7 day", date_labels="%Y-%m-%d") +
  scale_y_continuous(expand = c(0, 0), 
                     breaks=seq(0, 200, by = 20),
                     limits = c(0, 200),
                     sec.axis = sec_axis(~./10, name="R(t) when no internal infections in control",
                                         breaks=seq(0, 15, 1))) +
  xlab("Date") + ylab("Number of cases") +
  ggtitle("Xi'an") +
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
        legend.position = c(0.86, 0.88),
        # legend.justification = c("right", "top"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
Rt_inc_stack_xian

ggsave(filename = "../figure/Rt_inc_stack_xian.tiff",
       plot = Rt_inc_stack_xian,
       width = 25,  # <=19.05cm
       height = 15, # <=22.225cm
       units= "cm",
       dpi= 330,
       compression = "lzw")


#














