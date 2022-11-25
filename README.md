# EpidemicZoning


# List of data in excel:
1. Xi'an_Epidemic.xlsx <br />
   In Xi'an epidemic, 
   the column of date is the day when cases were reported;
   the column of social_total is the number of cases who were reported in social zone;
   the column of controlled_total is the number of cases who were reported in control zone;
   the column of local_total is the total number of cases who were reported in the city.
2. Shanghai_Epidemic.xlsx <br />
   In Shanghai epidemic,
   the column of date is the day when cases were reported;
   the column of social_total is the number of cases who were reported in social zone;
   the column of controlled_total is the number of cases who were reported in control zone;
   the column of local_total is the total number of cases who were reported in the city.


# List of R code:
1. covid19_Rt_xian_shanghai.R <br />
    Assume that there are no internal infections within the control zone. 
    The effective reproduction numbers over time for both Xi'an epidemic and Shanghai epidemic were calculated and visualized as in Figure 2 in the main text.  <br />
2. covid19_xian_Rc_alpha_t_mle.R <br />
    Assume that there are internal infections within the control zone in Xi'an epidemic.
    (1) The control reproduction number as constant and the time-varying $beta_t$ was estimated and visualized as in the top panel in Figure 3 in the main text.
    (2) The effective reproduction number in social zone was calculated and visualized as in the top panel in Figure 4 in the main text.<br />
3. covid19_sh_Rc_alpha_t_mle.R <br />
    Assume that there are internal infections within the control zone in Shanghai epidemic.
    (1) The control reproduction number as constant and the time-varying $beta_t$ was estimated and visualized as in the top panel in Figure 3 in the main text.
    (2) The effective reproduction number in social zone was calculated and visualized as in the top panel in Figure 4 in the main text.<br />

