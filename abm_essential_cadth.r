# Agent Based Model
# - Built based on 10 states microsimulation model: Susceptible, Ho
# - Quarantine Policy added 

# last update Sep 14 by Jasper Zhang


#load required packages
if (!require('pacman')) install.packages('pacman'); library(pacman) # use this package to conveniently install other packages
# load (install if required) packages from CRAN
p_load("diagram") 
#install_github("DARTH-git/darthtools", force = TRUE) #Uncomment if there is a newer version
#p_load_gh("DARTH-git/darthtools")

library(darthtools)
require(dplyr)
library(here)
library(readr)
library(lubridate)
library(tidyr)
library(mvtnorm)
library(readxl)
library(tictoc)
source(here("abm_functions.R"))
#global variables




#set seed
#g_seed = 1
df_list = list()
time_list = list()
for (g_seed in c(1:2)){

set.seed(g_seed)  # set the seed 

#starting and end date of simulation
date_start                  <- as.Date("2020-11-03")
date_cali                   <- as.Date("2020-11-19")
date_end                    <- as.Date("2020-12-17")

# quarantine period
q_lock                      <- 14
# incubation period
latentperiod                <- 4 # time from receiving virus to being infectious
incubation                  <- 3 # symptomatic: time from being infectious to develop symptoms

# Model structure
#states
v_n                         <- c("susceptible", 
                                 "infected_asymptomatic",
                                 "infected_presymptomatic",
                                 "infected_symptomatic",
                                 "susceptible_isolated", 
                                 "infected_asymptomatic_isolated",
                                 "infected_presymptomatic_isolated",
                                 "infected_symptomatic_isolated",
                                 "hospitalized",
                                 "recovered")  # vector with state names

n_states                    <- length(v_n)   # number of states
n_i                         <- 1000000       # number of individuals
seed_inf                    <- 2000          # number of infected people at day 0 


### function definitions
#not in
`%!in%`                      <- Negate(`%in%`)


#define groups of states

ending_states                <- c("hospitalized","recovered")

infect_states                <- c("infected_asymptomatic",
                                 "infected_presymptomatic",
                                 "infected_symptomatic")

infect_states_iso            <- c("infected_asymptomatic_isolated",
                                 "infected_presymptomatic_isolated",
                                 "infected_symptomatic_isolated")


all_infect_states            <- c(infect_states,infect_states_iso)

all_presympt_states          <- c("infected_presymptomatic_isolated",
                                 "infected_presymptomatic")

iso_states                   <- c("susceptible_isolated",
                                 "infected_asymptomatic_isolated",
                                 "infected_presymptomatic_isolated",
                                 "infected_symptomatic_isolated")

sus_states                   <- c("susceptible","susceptible_isolated")



#for (n_i in seq(10000, 100000,10000)){
tic() #start time record
  
#calendar date set up
length_of_study <- as.numeric(date_end  - date_start) + 1
cali_start      <- as.numeric(date_cali - date_start) + 1
cali_len        <- as.numeric(date_end  - date_cali)  + 1
#PP: ARE YOU SURE THAT CALI_END IS CORRECT ABOVE?
#JZ: in david's recent result, he ran the model till end of Dec.


n_t             <- length_of_study
all_dates       <- as.POSIXlt(seq(date_start, date_start + n_t - 1, by="day")) # convert to appropriate format
holidays        <- read_excel(here::here("holidays","holiday.xlsx"))           # load dates of holidays
holidays_date   <- as.POSIXlt(holidays$Date)                                   # convert to appropriate format
not_holiday     <- all_dates %!in% holidays_date                               # identify days not a holiday
all_wday        <- all_dates$wday                                              # day of the week number (0-6)
all_mon         <- all_dates$mon                                               # month number


#PP: Why not 5 to 8? please provide appropriate documenation throughout the code so i can debug
#JZ: in posixlt mon starts from 0 School year is from Sep code: 9 to June code 5
schoolday_ind   <- ((all_mon >= 8 | all_mon <=5) & (all_wday > 0 & all_wday < 6) & (not_holiday))# test if today is a school day mon 
work_dates      <- all_dates[all_wday  > 0 & all_wday  < 6]       # only select 5 working days PP: this is in agreement to david's model? # JZ: Yes
part_time_date  <- work_dates[seq(1, length(work_dates), by = 2)] # part time implies only working half the dates
fulltime_ind    <- all_dates %in% work_dates     & (not_holiday)  #identify days working for full time workers
parttime_ind    <- all_dates %in% part_time_date & (not_holiday)  # identify dats working for part time workers
designated_ind  <- schoolday_ind & parttime_ind                   # not sure what this is


# create a data frame with each individual's 
# ID number, treatment effect modifier, age and initial time in infected state 
df_X <- create_abmpop_real(n_i)

#sample from prob conditional prob of detection
#are you detected or not
#infected are you 

head(df_X) # print the first rows of the dataframe

MicroSim <- function(n_i, df_X, seed = 1) {
  # Arguments:  
  # n_i: number of individuals
  # df_X: data frame with individual data 
  # seed: seed for the random number generator, default is 1
  # Returns:
  # results: data frame with states tracing
  set.seed(seed) # set a seed to be able to reproduce the same results
  
  # create three matrices called m_M, 
  # number of rows is equal to the n_i, the number of columns is equal to n_t 
  # (the initial state and all the n_t cycles)
  # m_M is used to store the health state information over time for every individual

  m_M    <-  matrix(nrow = n_i, ncol = n_t + 1, 
                               dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                               paste("cycle", 0:n_t, sep = " ")))  
  
  m_vts    <-  matrix(nrow = n_i, ncol = n_t + 1, 
                    dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                    paste("cycle", 0:n_t, sep = " ")))  
  
  
  # PP: how is that different # JZ: this is designed to keep track of exposures everyday m_M is for states record
  m_expl <-  matrix(nrow = n_i, ncol = n_t, 
                 dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                 paste("cycle", 1:n_t, sep = " ")))  
 
  
  m_M[, 1]         <- df_X$state     # initial health state
  m_vts[,1]        <- df_X$v_Ts  
  v_Ts             <- df_X$v_Ts      # initialize time since illness onset
  v_T_presympt     <- rep(0, n_i)    # initialize time since being in presymptomatic state
  v_inf            <- rep(0, n_i)    # initialize number of exposures in a family


  v_total_inf_daily <- rep(0,n_t)    # initialize the number of total infections per cycle
  v_new_inf_daily   <- rep(0,n_t)    # initialize the number of total new infections per day
  v_detected_daily  <- rep(0,n_t)    # initialize the number of total infections per cycle
  
  # open a loop for time running cycles 1 to n_t 
  for (t in 1:n_t) {
 # for (t in 1:9) {
    print(paste("running day: ",as.character(t),sep = ""))
    # check if current running cycle is 
    
    is_school_day   <-  schoolday_ind[t] # if students will go to school today
    is_designated   <- designated_ind[t] # if students from designated school will go to school
    is_parttime_day <-   parttime_ind[t] # if parttime worker will go to work
    is_fulltime_day <-   fulltime_ind[t] # if fulltime worker
    
    pop <- df_X
    #find people who are not at home
    # univ students, on work day:  healthy essential workers ,school day students and teachers, sick essential workers without pay sick leave
    
    ft_student <- (is_school_day & pop$inclass & pop$designated == FALSE) & (pop$state %!in% iso_states)
    pt_student <- (is_designated & pop$inclass & pop$designated == TRUE)  & (pop$state %!in% iso_states)
    
    ft_worker  <- ((pop$essential == 1) & is_fulltime_day & pop$ftpt2 == 2) & (pop$state %!in% iso_states)
    pt_worker <- ((pop$essential == 1)  & is_parttime_day & pop$ftpt2 == 1) & (pop$state %!in% iso_states)
    
    # PP: please breakdown the statement below to simpler steps so that we are sure nothing is wrong # JZ: i'll breakdown
    left_home  <- ft_student |  #full time student #(pop$inuniv) |  # university students wont stay at home PP: you or david made this assumption?
                  pt_student |  #partime student                    # also i dont see how you sleect the uni students here. 
        !(pop$state %in% iso_states) & # people in quarantine are all home
         (pop$state %in% sus_states &( #healthy
        ((pop$essential == 1) & is_parttime_day & pop$ftpt2 == 1) |   # part time essential worker
        ((pop$essential == 1) & is_fulltime_day & pop$ftpt2 == 2))) | # full time essential worker PP: FT is 1 or 2? 
         (pop$state %in% all_infect_states &                          # infected essential worker without pay sick leave
          pop$essential == 1  & pop$getspsl == 0)
    
    #stay_home <- !left_home
    
    stay_home <- rep(TRUE,nrow(pop))#!(pop$inuniv)  # JZ: everyone will have at home interaction with others PP: not sure why you assume that? 
    
    
    
    #check exposure
    df_X$expl            <- 0                                 # reset expl level
    df_X$expl_home       <- 0                                 # reset expl level
    df_X$contact_type    <- 0                                 # reset contact type
    df_X <- family_exposure(df_X, stay_home)
    df_X <- school_exposure(df_X, ft_student, pt_student)
    df_X <- work_exposure(  df_X, ft_worker, pt_worker)
    
    
    m_expl[,t] <- df_X$expl
    # calculate the transition probabilities for the cycle based on health state t
    m_P <- Probs(m_M[, t], df_X, v_Ts)
    # check if transition probabilities are between 0 and 1
    check_transition_probability(m_P, verbose = TRUE)
    # check if each of the rows of the transition probabilities matrix sum to one
    check_sum_of_transition_array(m_P, n_states = n_i, n_cycles = n_t, verbose = TRUE)
    # sample the next health state and store that state in matrix m_M
    prob_new_state = samplev(m_P) #generate new states based on the prob
    df_X$state = prob_new_state
    
    #if (t > 1){
    df_X = detection(df_X,v_Ts)
    #}
    prob_new_state = df_X$state
    
    v_detected_daily[t] = df_X$new_detection[1] # sum((m_M[, t] %in% infect_states  & prob_new_state %in% infect_states_iso ))
    v_total_inf_daily[t] = sum(prob_new_state %in% all_infect_states)
    v_new_inf_daily[t] = sum((m_M[, t] %in% sus_states   & prob_new_state %in% all_infect_states ))
 
    
    
    #update Quarantine Status
    df_X = Quarantine_Check(df_X,m_M[, t]) #m_M[,t] is the previous day's state
    qc_state = df_X$state # state_after quarantine check
    

    # update time since illness onset for t + 1 
    v_Ts <- if_else(qc_state %in% all_infect_states, v_Ts + 1, 0)
    m_vts[,t + 1] =  v_Ts 
    
    df_X$inf <- if_else(v_Ts >= latentperiod, TRUE, FALSE)# infectious people only people who are in infected states more than latent days are infectious
    
    #check_presymptomatic state time
    v_T_presympt <-if_else(qc_state %in% all_presympt_states, v_T_presympt + 1, 0)
    v_go_to_sympt <- if_else(v_T_presympt >=3, TRUE, FALSE)
    df_X$state[v_go_to_sympt] =  paste(substr(df_X$state[v_go_to_sympt],start = 1, stop = 9) ,substring(df_X$state[v_go_to_sympt],13),sep = "") #this will remove "pre"
    
    v_go_to_recover <- if_else(v_Ts >= 17, TRUE, FALSE)
    df_X$state[v_go_to_recover] = "recovered"
    
    m_M[, t + 1]  <- df_X$state  # new day's state
    
    
    # Display simulation progress
    if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n_t * 100, "% done", sep = " "))
    }
   
    
  } # close the loop for the time points 
  

  # store the results from the simulation in a list
  results <- list(m_M = m_M ,
                  pop = df_X, 
                  v_new_inf_daily = v_new_inf_daily, 
                  v_total_inf_daily =   v_total_inf_daily ,
                  v_detected_daily =  v_detected_daily)   
  
  return(results)  # return the results
  
} 
# end of the `MicroSim` function  


# By specifying all the arguments in the `MicroSim()` the simulation can be started

# Run the simulation model

outcomes <- MicroSim(n_i = n_i, df_X = df_X, seed = 1)

toc(log = TRUE, quiet = TRUE)
# Show results
#}





# running time recorded by tictoc
log.txt <- tic.log(format = TRUE)
log.nf <- tic.log(format = FALSE)
tic.clearlog()

# time in txt format
writeLines(unlist(log.txt))
# time in numerical format
timings <- unlist(lapply(log.nf, function(x) x$toc - x$tic))



result <-
  data.frame(
    dates = as.Date(all_dates) ,
    v_new_inf_daily = outcomes$v_new_inf_daily, v_detected_daily =  outcomes$v_detected_daily)#,
    #v_total_inf_daily =   outcomes$v_total_inf_daily)


df_list[[g_seed]] = result
time_list[[g_seed]] = timings 
}    
saveRDS(df_list, file = "100iters.rds")



library("reshape2")
library("ggplot2")
result$v_new_inf_daily = round(result$v_new_inf_daily* 14.5)
#result$v_total_inf_daily = result$v_total_inf_daily * 14.5
result$v_detected_daily = round(result$v_detected_daily * 14.5)

running = result %>% slice(2:16)
cali = result %>% slice(17:n())

run_data_long <- melt(running, id.vars="dates")
run_data_long$days = rep(c(1:15),2)# convert to long format
cali_data_long <- melt(cali, id.vars="dates")  
cali_data_long$days = rep(c(1:29),2)






ggplot(data=run_data_long,
       aes(x=days, y=value, colour=variable)) +
  geom_line() + ggtitle("R ABM model running period1104-1118") + 
  labs(color='Variables')  +
  scale_color_manual(labels = c("Daily New Infections", "Daily Detected Cases")
                     ,values = c("blue", "red"))+ ylim(0, 15000) 

ggplot(data=cali_data_long,
       aes(x=days, y=value, colour=variable)) +
  geom_line() + ggtitle("R ABM model calibration period1119-1217") + labs(color='Variables')  +
  scale_color_manual(labels = c("Daily New Infections", "Daily Detected Cases")
                     ,values = c("blue", "red")) + ylim(0, 15000)



