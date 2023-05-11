
# Functions definitions for Agent Based Model 
# last update Sep 14 by Jasper Zhang

#function that helps change column name of a dataframe
#args:
# data: input dataset
# name: the name to be changed
# change_to: the objective column name
change_col_name <- function(data,name,change_to){
  names(data)[names(data) == name] = change_to
  return(data)
}

# convert char to num
convert_num_char <- function(x) {
  return(as.numeric(as.character(x)))
}

#PP:  can you please provide a tiny example that proves these functions work?

#prob function of symptomatic
fSympt <- function(age, sex){
  # Arguments:
  # age: age of individual
  # sex: gender of individual male 0 female 1
  # return: the probability of being symptomatic
  
  male                      = as.integer(1-sex)
  b0                        = -0.711118121955952
  b1019                     = -0.472050382359539
  b2029                     = -0.764169562107295
  b3039                     = -0.868515165981507
  b4049                     = -0.926066409326952
  b5059                     = -1.02056140823593
  b6069                     = -0.778262795462909
  b7079                     = -0.609852069159838
  b80plus                   = -0.102291803009561
  bmale                     = 0.0674273429920238
  
  v_b1 = c(0,b1019 , 
           b2029 ,
           b3039 ,
           b4049 ,
           b5059 ,
           b6069 ,
           b7079 ,
           b80plus, 
           b80plus,
           b80plus,
           b80plus,
           b80plus,
           b80plus)
  
  bsum                      = b0 + v_b1[floor(age / 10) + 1] + (male * bmale) 
  odds                      = exp(bsum)
  prob                      = odds / (1 + odds)
  return (1 - prob)
}




#create abm population from data files
create_abmpop_real <- function(npop){
  # Arguments:
  # npop: population size
  # Returns: 
  # initialized new population dataframe
  
  # Step 1 : read data tables  
  tbl.individual            <- read_tsv(here::here("individual_data","tblIndivData.txt"))
  tbl.classrooms            <- read_tsv(here::here("individual_data","tblClassrooms.txt"))
  tbl.custudents            <- read_tsv(here::here("individual_data","tblCUstudents.txt"))
  tbl.hhdata                <- read_tsv(here::here("individual_data","tblHHdata.txt"))
  tbl.patches               <- read_tsv(here::here("individual_data","tblPatches.txt"))
  tbl.regions               <- read_tsv(here::here("individual_data","tblRegions.txt"))
  tbl.schools               <- read_tsv(here::here("individual_data","tblSchools.txt"))
  tbl.workers               <- read_tsv(here::here("individual_data","tblWorkers.txt"))
  #vaccination status (future)
  #tbl.vaccschedage         <- read_tsv(here::here("individual_data","tblVaccSchedAge.txt"))
  #tbl.vaccschedincome      <- read_tsv(here::here("individual_data","tblVaccSchedIncome.txt"))  
  
  # Step 2: select columns needed
  #PP what is 2 and what is 10? # JZ: updated
  classinfo                 <- read_tsv(here::here("individual_data","tblClassrooms.txt"))
  sample_ppl                <- tbl.individual[1:npop,]
  sample_ppl                <- sample_ppl[,names(sample_ppl)[1:20]]
  sample_ppl                <- change_col_name(sample_ppl,"famid","fid")
  sample_ppl                <- change_col_name(sample_ppl,"workplaceid","wid")

  
  cols_select               <- c("classrmid","collunivID","fid","wid","regionid","patchid","essential","ftpt2")
  sample_ppl[cols_select]   <- sapply(sample_ppl[cols_select], convert_num_char)  # convert character to numeric
  sample_ppl$atwork         <- as.logical(sample_ppl$atwork)    # indicator: at work
  sample_ppl$inclass        <- sample_ppl$classrmid  != 0 #indicator: if individual is  inclass or not
  sample_ppl$inuniv         <- sample_ppl$collunivID != 0 #indicator: if individual is in university or not
  sample_ppl$designated     <- FALSE # designated school
  sample_ppl$qtime          <- 0     # quarantine time clock
  sample_ppl$time           <- 0     # total time put in quarantine
  sample_ppl$v_Ts           <- 0     # time of being infected
  sample_ppl$inf            <- FALSE # infectious indicator
  sample_ppl$age_group      <- ifelse(sample_ppl$age <= 10,1,2)                      # young child/oldchild
  sample_ppl$age_group      <- ifelse(sample_ppl$age <= 17, sample_ppl$age_group, 3) # adult
  
  sample_ppl$new_detection  <- 0
  
  # Step 3: initialize states and other variables
  sample_ppl$state          <- "susceptible"
  #previously infected people by age group: 0-9, 10-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+
  #PP: where is this information coming and is it per million? in other words you can't use that with a different npop value
  #JZ: Problem Solved
  pinfect_prop            <- npop/1000000 #below num in different age group is 1million population based
  vec_pinfect               <- round(pinfect_prop * c(377,1920,4450,3322,2964,3168,2056,1155,2093))
  age_break                 <- c(0,10,20,30,40,50,60,70,80,200)
  seed_inf                    <- round(2000 * pinfect_prop )     #update seed_inf based on the propotion of npop/1million
  
  
  
  #PP please break such statements to  smaller concepts
  for (i in c(1:length(vec_pinfect))){ # set them to ending state "recovered"
    pi_num                  <- vec_pinfect[i]
    within_age   <- (sample_ppl$age >= age_break[i]) & (sample_ppl$age < age_break[i + 1])  # find those in th age range 
    age_id       <- sample_ppl$ptid[within_age]                                             # identify their IDs
    sampl_id_age <- sample(age_id , pi_num)                                                 # sample those previously infected
    
    sample_ppl$state[sample_ppl$ptid %in% sampl_id_age] = "recovered"                       # set them as recovered. 
  }
  
  
  sample_ppl$expl           <- 0 # exposure level: number of infected people contacted within a day
  sample_ppl$expl_home      <- 0 # exposure level: number of infected people contacted at home within a day
  sample_ppl$contact_num    <- 0 # number of contacts a day
  sample_ppl$contact_max    <- 0 # max number of contact
  sample_ppl$contact_type   <- 0 # contact type
  # JZ: assign probability of being symptomatic based on individual's age and gender
  sample_ppl$prob_sympt     <- fSympt(sample_ppl$age, sample_ppl$sex)   # PP EXPLANATION WHAT IS THIS? # probability of being symptomatic?
  m_Asymp_Symp              <- as.matrix(data.frame(Asympt = 1-sample_ppl$prob_sympt ,
                                                     Sympt =   sample_ppl$prob_sympt)   )
  
  # PP Can't you do that with a rbinom?
  
  sample_ppl$sympt          <- samplev(m_Asymp_Symp) # asymptomatic and symptomatic directions are pre-defined.
  
  # infected seed pop seed_inf 
  #PP break the statements to simpler steps
  susc_ids                          <- sample_ppl[sample_ppl$state == "susceptible",]$ptid       #the IDs of those who are susceptible
  seed_ids                          <- sample(susc_ids, seed_inf) # sample seed infection pt ids
  select_seed_inf                   <- sample_ppl$ptid %in% seed_ids #boolean
  seed_inf_states                   <- ifelse(sample_ppl$sympt[select_seed_inf] == "Sympt",  # if symptomatic and infected
                                              "infected_presymptomatic",                     # set as presymptomatic
                                              "infected_asymptomatic")                       # otherwise set as asymptomatic
  sample_ppl$state[select_seed_inf] <- seed_inf_states                                       # assign to those seed infected a state
  
  
  seedpop                           <- sample_ppl$state %in% all_infect_states              #logical variable indicating if someone is in the seed population 
  prepop                            <- sample_ppl$state == "infected_presymptomatic"        # logical variable indicating if someone is infected.
  
  
  sample_ppl$v_Ts[seedpop] <- round(runif(sum(seed_inf),1,17))                             #duration of infection by model start
  
  sample_ppl$inf <- if_else(sample_ppl$v_Ts >=latentperiod, TRUE, FALSE)# infectious people only people who are in infected states more than latentperiod days are infectious
  
  sample_ppl$state[prepop & sample_ppl$v_Ts >=(latentperiod + incubation)] = "infected_symptomatic"# infectious people only people who are in infected states more than latentperiod days are infectious
  
  
  #probability of being detected
  p_detect_sympt_overall  <- 0.635  # symptomatic PP: source?
  p_detect_asympt_overall <- 0.435  # asymptomatic PP: source?
  
  sample_ppl$detected = 0
  sample_ppl$detect_time = 9999
  
  n_sympt = sum(sample_ppl$sympt == "Sympt")
  
  n_asympt = npop - n_sympt
  
  n_detect_sympt = round(p_detect_sympt_overall*n_sympt)
  n_detect_asympt = round(p_detect_asympt_RT*n_asympt)
  
  id_detect_sympt = sample((sample_ppl[sample_ppl$sympt == "Sympt" ,]$ptid),n_detect_sympt )
  id_detect_asympt= sample((sample_ppl[sample_ppl$sympt == "Asympt" ,]$ptid),n_detect_asympt )
  
  sample_ppl$detected[id_detect_sympt] = 1
  sample_ppl$detected[id_detect_asympt] = 1
  
  sample_ppl$detect_time[id_detect_sympt] = floor(runif(n_detect_sympt, min=1, max=5))
  sample_ppl$detect_time[id_detect_asympt] = floor(runif(n_detect_asympt, min=1, max=12))# 12 is infectious period -1
  
  
  #designated school, students study on part-time base, go to school every the other day # JZ: updated 
  designated                <- select(tbl.schools,"index","designated" )  #PP not sure hat designated school means
  names(designated)[1]      <- "schoolid"
  class                     <- merge(tbl.classrooms,designated,by = c("schoolid"))
  designated_class          <- class[class$designated == 1,]
  designated_id             <- designated_class[,5:28] # PP what is 5 nd 28?
  #change dataframe structure to vector
  # PP there is a "as.vector" command no?
  m                         <- as.matrix(designated_id)
  q                         <- matrix(m, 1)
  
  designated_id_vec         <- sort(unique(as.numeric(as.character(q))),decreasing = FALSE) # sort school ids?
  designated_id_vec         <- designated_id_vec[-1]                                        # remove first 0 
  designated_id_sample      <- sample_ppl$ptid %in% designated_id_vec                       # identify the pt ids who are in the designated vector 
  sample_ppl$designated[designated_id_sample] = TRUE
  
  # JZ: count teacher as extra member of the class
  tbl.classrooms$nclass     = tbl.classrooms$nclass + 1                  #   #update teachers PP : why +1? 
  #teachers' working schedules are the same as their students
  teacher_id                = sample_ppl$ptid %in% tbl.classrooms$teacherid
  teachers                  = sample_ppl[teacher_id,]
  sample_ppl[teacher_id,]$atwork = FALSE
  
  
  #initialize number of contacts from multivariate normal
  #work, school, patch, region
  # PP where is this being used?
  # JZ 
  v_mean                    <- c(2.282382386,1.609437912,0.405465108,0.587786665)
  v_sd                      <- c(0.456476477,0.321887582,0.081093022,0.117557333)
  corr_matrix               <- matrix(data = c(1,0.9,0.9,0.9,
                                              0.9,1,0.9,0.9,
                                              0.9,0.9,1,0.9,
                                              0.9,0.9,0.9,1),
                                     nrow = 4,byrow = TRUE)
  cov_matrix                <- diag(v_sd) %*% corr_matrix %*% diag(v_sd)
  
  #round(exp(rmvnorm(n = 1000, mean = v_mean, sigma = cov_matrix,method = "chol")),0)[,1] #work
  #round(exp(rmvnorm(n = 1000, mean = v_mean, sigma = cov_matrix,method = "chol")),0)[,2] #school
  
  
  #set number of contacts in class
  #tbl.classrooms$nclass = tbl.classrooms$nclass + 1 # add teacher to the class size count
  #PP check if this is done right and add more documentation
  class_and_size            <- tbl.classrooms[,4:28] %>% gather(key = "id_order" ,value = "ptid", -nclass) %>% select(ptid,nclass)
  sid_in_sample             <- which(sample_ppl$ptid %in% class_and_size$ptid)
  
  
  ##double booked students PP: isnt that the reverse of the above? and whey do you have double booked students?
  student_class_match       = which(class_and_size$ptid %in% sid_in_sample)
  appear_count              = table(class_and_size$ptid[student_class_match])
  double_booked             = which(table(class_and_size$ptid[student_class_match]) == 2)
  double_id                 = as.numeric(as.character(names(appear_count[double_booked])))
  double_teacher            = which(tbl.classrooms$teacherid %in% double_id)
  tbl.classrooms[double_teacher, ]$teacherid = 0
  tbl.classrooms[double_teacher, ]$nclass = tbl.classrooms[double_teacher, ]$nclass - 1
  
  #PP are you doing this to remove the dublicates? 
  class_and_size            = tbl.classrooms[,4:28] %>% gather(key = "id_order" ,value = "ptid", -nclass) %>% select(ptid,nclass)
  sid_in_sample             = which(sample_ppl$ptid %in% class_and_size$ptid)
  n_students                = length(sid_in_sample) #total number of students and teachers
  #update contact num
  
  # PP ?You sample teh mean contacts per student here right? These are daily contacts? you need to measeure the number of contacts on david's data to see if you are having fewer contacts
  sample_ppl$contact_num[sid_in_sample]  <- round(exp(rmvnorm(n = n_students, mean = v_mean, sigma = cov_matrix,method = "chol")),0)[,2] #school
  sample_ppl$contact_max[sid_in_sample]  <- class_and_size$nclass[which(class_and_size$ptid %in% sid_in_sample)] -1
  sample_ppl$contact_type[sid_in_sample] <- 1 #in class
  
  
  #set number of contacts in university/college
  s_univ                     <- sample_ppl$ptid[which(sample_ppl$inuniv == TRUE)]
  s_class                    <- sample_ppl$ptid[which(sample_ppl$inclass == TRUE)]
  double_identity            <- intersect(s_univ,s_class) #people can be a teacher at highschool and univ students PP: unlikely and probably not true
  sample_ppl[double_identity,]$inuniv     = FALSE
  sample_ppl[double_identity,]$collunivID = 0             # ok so you remove univ and teacher
  
  univ_in_sample             <- sample_ppl$ptid[which(sample_ppl$inuniv == TRUE)]    # identify who is a uni student
  univ_size                  <- table(sample_ppl[univ_in_sample,]$collunivID)        # university size by uni ID
  univ_members               <- sample_ppl[univ_in_sample,]                          # all characteristics of uni students
  univ_info                  <- data.frame(collunivID = as.numeric(as.character(names(univ_size))),usize =  as.numeric(as.matrix(univ_size)))
  univ_members_size          <- merge(univ_members,univ_info, by = "collunivID")
  univ_members_size          <- univ_members_size[order(univ_members_size$ptid),]
  
  sample_ppl$contact_num[univ_in_sample]  <- 15                            # PP: what does the 15 mean here? why not sample for university? and is that assumption the saem as davids?
  sample_ppl$contact_max[univ_in_sample]  <- univ_members_size$usize - 1
  sample_ppl$contact_type[univ_in_sample] <- 2                             # in university/college
  
  #set number of contacts atwork
  work_in_sample             <- sample_ppl$ptid[which(sample_ppl$atwork == TRUE)]
  work_size                  <- table(sample_ppl[work_in_sample,]$wid) #workplacesize
  work_ppl                   <- sample_ppl[work_in_sample,]
  work_info                  <- data.frame(wid = as.numeric(as.character(names(work_size))),wsize =  as.numeric(as.matrix(work_size)))
  work_place_size            <- merge(work_ppl, work_info, by = "wid")
  work_place_size            <- work_place_size[order(work_place_size$ptid),]
  #PP whay not sampling once?
  sample_ppl$contact_num[work_in_sample]  <- round(exp(rmvnorm(n = length(work_in_sample), mean = v_mean, sigma = cov_matrix,method = "chol")),0)[,1]
  sample_ppl$contact_max[work_in_sample]  <- work_place_size$wsize -1
  sample_ppl$contact_type[work_in_sample] <- 3 #atworks
  
  
  #update_contact_num
  sample_ppl$contact_num     =  apply(X=data.frame(sample_ppl$contact_num,sample_ppl$contact_max), MARGIN=1, FUN=min) # PP: is such restriction sthgn that you saw on dvid's model too? 
  
  return(sample_ppl)
}



#check
detection <- function(pop,v_Ts){
  detection_period = v_Ts - latentperiod
  detection_period[pop$sympt == "Sympt"] = detection_period[pop$sympt == "Sympt"] - incubation# individual won't get detection before develop symptoms.
  detected_individuals = (detection_period == pop$detect_time) & pop$state %!in% ending_states
  pop$new_detection = sum(detected_individuals)
  new_iso = detected_individuals & pop$state %!in% iso_states
  pop$state[new_iso ] = paste(pop$state[new_iso ],"isolated",sep = "_")
  return(pop)
  }




#Quarantine_Check can update the quarantine clock at current time cycle
Quarantine_Check <- function(pop, pre_state){
  # Arguments:
  # pop: population dataframe of current time cycle t
  # pre_state: states at previously cycle t-1
  # Returns: 
  # new population dataframe with updated states
  
  # Step1: find people in different quarantine conditions
  
  # People in quarantine:
  # condition 1: reset quarantine clock for those who become infected or have close contact with infected people during quarantine
  new_infected_in_Q          = which(pop$state %in% infect_states_iso & pre_state == "susceptible_isolated") #got infected in quarantine
  
  reset_family               = unique(pop$fid[new_infected_in_Q])
  reset_people               = pop$fid %in% reset_family #find all the people who had expose in family quarantine
  
  
  # condition 2: find people are going to be released from quarantine
  # check if quarantine clock hit q_lock
  full_q                     = pop$qtime == q_lock & !reset_people & pop$state %!in% ending_states 
  
  #find susceptible people and infected people who meet quarantine time
  sus_full_q                 = full_q & pop$state == "susceptible_isolated"
  inf_full_q                 = full_q & pop$state %in% infect_states_iso 
  
  # condition 3: find people who stay in quarantine and the quarantine clock has to be incremented.
  else_q                     = (pop$qtime < q_lock & pop$qtime >0) & !reset_people & pop$state %!in% ending_states  
  
  # Other people:
  # condition 4: find newly infected people entering quarantine NOTE: people enter quarantine by detection
  new_infected               = which(pop$state %in% iso_states & pop$qtime == 0) #new infected generated from infected
  
  # condition q_lock: find newly infected people's dependencies in families and workpalces
  new_infected_family        = unique(pop$fid[new_infected]) #families have newly infected in quarantine cases
  
  #new_infected_work= unique(pop$wid[new_infected]) #workplaces have newly infected in quarantine cases #comment out if same workplace contact don't need to be put in quarantine
  
  new_qurantine              = (((pop$fid %in% new_infected_family & pop$qtime == 0) | #(pop$wid %in% new_infected_work & pop$qtime == 0) |
                                   (pop$fid %in% reset_family)) & (pop$state %!in% ending_states) & (!(pop$essential == 1 & pop$getspsl == 0))  ) #people whose quarantine clocks have to be set to 1 including
  new_qurantine_add = new_qurantine & pop$state %!in% iso_states #find people need to be put in quarantine state
  
  
  # Step2: Update quarantine clock
  # set new_quarantine clock to 1 
  pop$qtime[new_qurantine]    = 1
  pop$state[new_qurantine_add] = paste(pop$state[new_qurantine_add],"_isolated",sep = "") #update quarantine states
  pop$time[new_qurantine]     = pop$time[new_qurantine]+1 #increment total quarantine time
  
  
  # release people with meet q_lock-days quarantine requirement
  # - isolated susceptible people will go to susceptible
  # - infected people will go to recovered
  pop$qtime[full_q]           = 0
  pop$state[sus_full_q]       = "susceptible" #when quarantine ends, state goes back to susceptible
  pop$state[inf_full_q]       = "recovered" #when quarantine ends, state goes back to recovered
  
  #increment quarantine clock and total quarantine time for everyone else in quarantine
  pop$qtime[else_q]           = pop$qtime[else_q] + 1
  pop$time[else_q]            = pop$time[else_q] + 1
  return (pop)
}




#Exposure_Check family
family_exposure <- function(pop,stay_home){
  # Arguments:
  # pop: population dataframe of current time cycle t
  # stay_home: boolean values that TRUE for people who stay at home in a day/ or at the end of the day
  # Return:
  # population dataframe with updated family exposure level
  
  homeid                      <- pop$ptid[stay_home]        # ID of individuals at home
  infectious                  <- pop$inf[stay_home]         # who is infectious
  fid                         <- pop$fid[stay_home]         # family ID of individuals at home
  
  current_state               <- pop$state[stay_home]       # 
  expl_list                   <- pop$expl_home[stay_home]   # PP What does expl_list mean?
  
  df_exp                      <- data.frame(fid,expl_list) # PP: Why is this needed?
  
  # exposure check at home
  infectious_family           <- data.frame(table(fid[infectious])) # find number of infectious people in each family
  
  # if an infection has been identified in at least one family
  if(nrow(infectious_family) > 0){
    names(infectious_family)    <- c("fid","Freq_f")
    infectious_family$fid       <- convert_num_char(infectious_family$fid)
    infectious_family$Freq_f    <- convert_num_char(infectious_family$Freq_f)
    df_exp                      <- left_join(df_exp,infectious_family, by = "fid")   #why left join? to accomodate for the zero infections n families??
    df_exp$Freq_f[is.na(df_exp$Freq_f)] = 0
    df_exp$expl_list            <- df_exp$expl_list + df_exp$Freq_f
  }
  pop$expl_home[stay_home]      <- df_exp$expl_list   # PP: is that updating exposure?
  #pop$contact_type[stay_home]  <- 0
  return(pop)
  
}

#Exposure_Check School
school_exposure <- function(pop,ft_student,pt_student){
  allpop                      <- pop
  if(sum(ft_student) > 0){
    select_pop                <- allpop[ft_student,]
    expl_list                 <- select_pop$expl
    infectious                <- select_pop$inf
    cid_list                  <- select_pop$classrmid
    out_of_iso                <- select_pop$state %!in% iso_states
    ptid_list                 <- select_pop$ptid[out_of_iso & select_pop$state == "susceptible"]
    inclass                   <- select_pop$inclass
    classtable                <- data.frame(table(cid_list[inclass]))
    names(classtable)         <- c("classrmid","classsize")
    
    inf_student               <- select_pop[select_pop$inf,]
    inf_class                 <- inf_student$classrmid
    
    if (nrow(inf_student) > 0){
      contact                   <- c()
      for( i in c(1:nrow(inf_student))){
        select_student          = inf_student[i,]
        classptid               = select_pop[select_pop$classrmid == select_student$classrmid,]$ptid
        to_contact              = classptid[classptid %!in% select_student$ptid]                      # remove own student
        if (length(to_contact) > 1){
          sample_contacts <- sample(to_contact, min(length(to_contact),select_student$contact_num)) # PP why do you need teh min statement here and not just lenght(to_contact?)
          contact                 = c(contact,sample_contacts)
        }else if (length(to_contact) == 1){
          new_contact = to_contact
          contact                 = c(contact,new_contact)
        }
        
      }
      
      if (length(contact) != 0){
        contact_table             <- data.frame(table(contact))
        names(contact_table)      = c("ptid","expl")
        
        contact_table$ptid        <- as.numeric(as.character(contact_table$ptid))
        contact_table$expl        <- as.numeric(as.character(contact_table$expl))
        expl_list[select_pop$ptid %in% contact_table$ptid] = expl_list[select_pop$ptid %in% contact_table$ptid] + contact_table$expl
        allpop[(allpop$ptid %in% contact_table$ptid ),]$contact_type = 1
      }
      
      
    }
    allpop[ft_student,]$expl  <- expl_list
    
    
    
  }
  
  if(sum(pt_student) > 0){
    select_pop                = allpop[pt_student,]
    expl_list                 = select_pop$expl
    infectious                = select_pop$inf
    cid_list                  = select_pop$classrmid
    out_of_iso                = select_pop$state %!in% iso_states
    ptid_list                 = select_pop$ptid[out_of_iso & select_pop$state == "susceptible"]
    inclass                   = select_pop$inclass
    classtable                = data.frame(table(cid_list[inclass]))
    names(classtable)         = c("classrmid","classsize")
    inf_student               = select_pop[select_pop$inf,]
    inf_class                 = inf_student$classrmid
    
    if (nrow(inf_student) > 0){
      contact <- c()
      for( i in c(1:nrow(inf_student))){
        select_student        = inf_student[i,]
        classptid             = select_pop[select_pop$classrmid == select_student$classrmid,]$ptid
        to_contact            = classptid[classptid %!in% select_student$ptid]
        if (length(to_contact) > 1){
          contact               = c(contact,sample(to_contact,
                                                   min(length(to_contact),select_student$contact_num)))
        }else if (length(to_contact) == 1){
          new_contact = to_contact
          contact                 = c(contact,new_contact)
        }}
      
      if (length(contact) != 0){
        contact_table           <- data.frame(table(contact))
        names(contact_table)    = c("ptid","expl")
        
        contact_table$ptid      <- as.numeric(as.character(contact_table$ptid))
        contact_table$expl      <- as.numeric(as.character(contact_table$expl))
        expl_list[select_pop$ptid %in% contact_table$ptid] = expl_list[select_pop$ptid %in% contact_table$ptid] + contact_table$expl
        allpop[(allpop$ptid %in% contact_table$ptid ),]$contact_type = 1
      }}
    
    allpop[pt_student,]$expl  <- expl_list
    
  }
  
  return(allpop) 
}



#Exposure_Check Work
work_exposure <- function(pop,ft_worker,pt_worker){
  allpop                      = pop
  if(sum(ft_worker) > 0){
    select_pop                = allpop[ft_worker,]
    expl_list                 = select_pop$expl
    infectious                = select_pop$inf
    wid_list                  = select_pop$wid
    out_of_iso                = select_pop$state %!in% iso_states
    ptid_list                 = select_pop$ptid[out_of_iso & select_pop$state == "susceptible"]
    atwork                    = select_pop$atwork
    
    workplacetable                = data.frame(table(wid_list[atwork]))
    
    names(workplacetable)         = c("wid","wsize")
    
    inf_worker               = select_pop[select_pop$inf,]
    inf_workplace                = inf_worker$wid
    
    if (nrow(inf_worker) > 0){
      contact                   <- c()
      
      for( i in c(1:nrow(inf_worker))){
        select_worker           = inf_worker[i,]
        workplaceptid           = select_pop[select_pop$wid == select_worker$wid,]$ptid
        to_contact              = workplaceptid[workplaceptid %!in% select_worker$ptid]
        if (length(to_contact) > 1){
          new_contact = sample(to_contact,
                               min(length(to_contact),select_worker$contact_num))
          contact                 = c(contact,new_contact)
        }else if (length(to_contact) == 1){
          new_contact = to_contact
          contact                 = c(contact,new_contact)
        }
      }
      
      if (length(contact) != 0){
        contact_table             <- data.frame(table(contact))
        names(contact_table)      = c("ptid","expl")
        
        contact_table$ptid        <- as.numeric(as.character(contact_table$ptid))
        contact_table$expl        <- as.numeric(as.character(contact_table$expl))
        
        expl_list[select_pop$ptid %in% contact_table$ptid] = expl_list[select_pop$ptid %in% contact_table$ptid] + contact_table$expl
        allpop[(allpop$ptid %in% contact_table$ptid ),]$contact_type = 3
      }}
    
    allpop[ft_worker,]$expl  <- expl_list
    
  }
  
  
  if(sum(pt_worker) > 0){
    select_pop                = allpop[pt_worker,]
    expl_list                 = select_pop$expl
    infectious                = select_pop$inf
    wid_list                  = select_pop$wid
    
    out_of_iso                = select_pop$state %!in% iso_states
    ptid_list                 = select_pop$ptid[out_of_iso & select_pop$state == "susceptible"]
    atwork                    = select_pop$atwork
    
    workplacetable                = data.frame(table(wid_list[atwork]))
    names(workplacetable)         = c("wid","wsize")
    
    inf_worker               = select_pop[select_pop$inf,]
    inf_workplace                = inf_worker$wid
    
    if (nrow(inf_worker) > 0){
      contact                   <- c()
      for( i in c(1:nrow(inf_worker))){
        select_worker           = inf_worker[i,]
        workplaceptid           = select_pop[select_pop$wid == select_worker$wid,]$ptid
        to_contact              = workplaceptid[workplaceptid %!in% select_worker$ptid]
        
        
        if (length(to_contact) > 1){
          contact                 = c(contact,sample(to_contact,
                                                     min(length(to_contact),select_worker$contact_num)))
        }else if (length(to_contact) == 1){
          new_contact = to_contact
          contact                 = c(contact,new_contact)
        }}
      
      if (length(contact) != 0){
        contact_table             <- data.frame(table(contact))
        names(contact_table)      = c("ptid","expl")
        
        contact_table$ptid        <- as.numeric(as.character(contact_table$ptid))
        contact_table$expl        <- as.numeric(as.character(contact_table$expl))
        expl_list[select_pop$ptid %in% contact_table$ptid] = expl_list[select_pop$ptid %in% contact_table$ptid] + contact_table$expl
        allpop[(allpop$ptid %in% contact_table$ptid ),]$contact_type = 3
      }}
    allpop[pt_worker,]$expl  <- expl_list
    
  }
  
  return(allpop) 
}



#### Deterministic analysis ####
#Transition probabilities
#(all terminal probabilities are conditional on not being admitted/recovered)

#probability of being detected
p_detect_sympt_overall  <- 0.635  # symptomatic PP: source?
p_detect_asympt_overall <- 0.435  # asymptomatic PP: source?

p_detect_sympt <-  1 - (1 - p_detect_sympt_overall) ^  (1 / 10)  # PP: what is 10 ?
p_detect_asympt <- 1 - (1 - p_detect_asympt_overall) ^ (1 / 13)  # PP: what is 13 ?

rrCT <- 0.752835329151619  #PP: source? explanation?
p_detect_asympt_CT    <- p_detect_asympt_overall * rrCT
p_detect_asympt_CTadj <- min(1, p_detect_asympt_CT / p_detect_sympt)
p_detect_asympt_RT    <-   p_detect_asympt_overall * (1 - rrCT)


#probability of being admitted PP: how  is that needed? this is probbility of admission
p_Admit_base <- 0.0524597414138787# 1 - (1 - 0.0524597414138787) ^  (1 / 10)

#infection probabilities vary by location
#household
p_transhh_young_child <- 0.0160176427409341  #  <= 10 PP: is that from young to child? different child age groups
p_transhh_old_child   <- 0.077154345991161   #  <=17  PP: is that from old   to child?
p_transhh_adult       <- 0.076831649685427   #  > 18

#vector of household transmission probabilites
p_transhh = c(p_transhh_young_child, p_transhh_old_child,  p_transhh_adult)

#atwork \ university
p_TransmitOther        <- 0.0220595686685414
#school
orSchoolMeasures       <-  0.15 #PP: what is this? odds ratio
p_TransmitOther_school <- odds_to_prob(prob_to_odds(p_TransmitOther)*orSchoolMeasures) #from treeage
#colegeUnivs #p_TransmitOther 



infect_prob <- function(prob,exposure) { 
  # Arguments:
  # prob:     the probability of infection in one time contact
  # exposure: the time of being exposed to infectious individuals
  # Returns: 
  # the probability of being infected based on certain exposure time
  return(1 - (1 - prob) ^ (exposure))   # 1 - prob is the prob of being safe, 
  # 1 minus the prob of always being safe is the prob of getting infection
}

Probs <- function(M_t, df_X, v_Ts) { 
  # Arguments:
  # M_t:  health state occupied at cycle t (character variable)
  # df_X: data frame with individual characteristics data 
  # v_Ts: vector with the duration of being infected
  # Returns: 
  # transition probabilities for that cycle
  
  # create matrix of state transition probabilities
  m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  
  # give the state names to the rows
  rownames(m_p_t) <-  v_n                               
  
  
  #Calculate the probability of being infected.
  sus <- df_X[df_X$state == "susceptible",] #identify those susceptible
  
  #PP: plese provide documentation
  sus_prob_home  <- rep(0, nrow(sus)) 
  one_prob       <- rep(1, nrow(sus))
  contact_active <- sus[sus$expl > 0 | sus$expl_home > 0, ]
  at_home        <- contact_active$expl_home > 0
  
  if(sum(at_home) > 0){
    
    home_infect = infect_prob(p_transhh[contact_active[at_home,]$age_group],contact_active[at_home,]$expl_home)
    
    sus_prob_home[sus$expl_home > 0] = home_infect
    
  }
  
  
  sus_prob = rep(0,nrow(sus))
  at_school = contact_active$contact_type == 1
  if(sum(at_school) > 0){
    school_infect = infect_prob(p_TransmitOther_school,contact_active[at_school,]$expl)
    sus_prob[sus$expl > 0 & sus$contact_type == 1] = school_infect
    
  }
  
  # at_univ = contact_active$contact_type == 2
  # if(sum(at_univ) > 0){
  #   univ_infect = infect_prob(p_TransmitOther,contact_active[at_univ,]$expl)
  #   sus_prob[sus$expl > 0 & sus$contact_type == 2] = univ_infect
  #   
  # }
  
  at_work = contact_active$contact_type == 3
  if(sum(at_work) > 0){
    work_infect = infect_prob(p_TransmitOther,contact_active[at_work,]$expl)
    sus_prob[sus$expl > 0 & sus$contact_type == 3] = work_infect
  }
  
  asy_prob = rep(0,length(sus$sympt ))
  asy_prob[sus$sympt == "Asympt"] = 1
  p_Sus_Asy <- asy_prob * (1- ((1- sus_prob_home) * (1 - sus_prob))) # 1 minus the prob of being safe all day
  
  pre_prob = rep(0,length(sus$sympt ))
  pre_prob[sus$sympt == "Sympt"] = 1
  p_Sus_Pre <- pre_prob * (1- ((1- sus_prob_home) * (1 - sus_prob)))
  
  
  
  sus_iso = df_X[df_X$state == "susceptible_isolated",]
  sus_iso_prob_home = rep(0,nrow(sus_iso))
  one_prob = rep(1,nrow(sus_iso))
  contact_active <- sus_iso[sus_iso$expl > 0 | sus_iso$expl_home > 0 ,]
  at_home = contact_active$expl_home > 0
  
  if(sum(at_home) > 0){
    home_infect = infect_prob(p_transhh[contact_active[at_home,]$age_group],contact_active[at_home,]$expl_home)
    sus_iso_prob_home[sus_iso$expl_home > 0] = home_infect
  }
  
  #other contact locations
  sus_iso_prob = rep(0,nrow(sus_iso))
  at_school = contact_active$contact_type == 1
  if(sum(at_school) > 0){
    school_infect = infect_prob(p_TransmitOther_school,contact_active[at_school,]$expl)
    sus_iso_prob[sus_iso$expl > 0 & sus_iso$contact_type == 1] = school_infect
  }
  
  
  at_work = contact_active$contact_type == 3
  if(sum(at_work) > 0){
    work_infect = infect_prob(p_TransmitOther,contact_active[at_work,]$expl)
    sus_iso_prob[sus_iso$expl > 0 & sus_iso$contact_type == 3] = work_infect
  }
  
  asy_prob_iso = rep(0,length(sus_iso$sympt ))
  asy_prob_iso[sus_iso$sympt == "Asympt"] = 1
  p_Sus_iso_Asy <- asy_prob_iso* (1- (1 - sus_iso_prob) * (1 - sus_iso_prob_home))
  
  pre_prob_iso = rep(0,length(sus_iso$sympt ))
  pre_prob_iso[sus_iso$sympt == "Sympt"] = 1
  p_Sus_iso_Pre <- pre_prob_iso * (1- (1 - sus_iso_prob) * (1 - sus_iso_prob_home))
  
  
  
  #to susceptible 
  m_p_t["susceptible"                         , M_t == "susceptible"]                   <- 1 - p_Sus_Asy - p_Sus_Pre
  m_p_t["infected_asymptomatic"               , M_t == "susceptible"]                   <- p_Sus_Asy
  m_p_t["infected_presymptomatic"             , M_t == "susceptible"]                   <- p_Sus_Pre 
  m_p_t["infected_symptomatic"                 , M_t == "susceptible"]                  <- 0
  m_p_t["susceptible_isolated"                , M_t == "susceptible"]                   <- 0
  m_p_t["infected_asymptomatic_isolated"      , M_t == "susceptible"]                   <- 0
  m_p_t["infected_presymptomatic_isolated"    , M_t == "susceptible"]                   <- 0
  m_p_t["infected_symptomatic_isolated"        , M_t == "susceptible"]                  <- 0
  m_p_t["hospitalized"                        , M_t == "susceptible"]                   <- 0 
  m_p_t["recovered"                           , M_t == "susceptible"]                   <- 0
  
  
  #to infected asy  
  m_p_t["susceptible"                       , M_t == "infected_asymptomatic"   ]        <- 0
  m_p_t["infected_asymptomatic"             , M_t == "infected_asymptomatic"   ]        <- 1 #(1 -  p_detect_asympt)
  m_p_t["infected_presymptomatic"           , M_t == "infected_asymptomatic"   ]        <- 0
  m_p_t["infected_symptomatic"               , M_t == "infected_asymptomatic"   ]       <- 0
  m_p_t["susceptible_isolated"              , M_t == "infected_asymptomatic"   ]        <- 0
  m_p_t["infected_asymptomatic_isolated"    , M_t == "infected_asymptomatic"   ]        <- 0 #p_detect_asympt # detection
  m_p_t["infected_presymptomatic_isolated"  , M_t == "infected_asymptomatic"   ]        <- 0
  m_p_t["infected_symptomatic_isolated"      , M_t == "infected_asymptomatic"   ]       <- 0
  m_p_t["hospitalized"                      , M_t == "infected_asymptomatic"   ]        <- 0
  m_p_t["recovered"                         , M_t == "infected_asymptomatic"   ]        <- 0
  
  #infected_pre
  m_p_t["susceptible"                       , M_t == "infected_presymptomatic"  ]       <- 0
  m_p_t["infected_asymptomatic"             , M_t == "infected_presymptomatic"  ]       <- 0
  m_p_t["infected_presymptomatic"           , M_t == "infected_presymptomatic"  ]       <- 1 #(1 - p_Admit_base)# * (1 -  p_detect_asympt)
  m_p_t["infected_symptomatic"               , M_t == "infected_presymptomatic"  ]      <- 0 
  m_p_t["susceptible_isolated"              , M_t == "infected_presymptomatic"  ]       <- 0
  m_p_t["infected_asymptomatic_isolated"    , M_t == "infected_presymptomatic"  ]       <- 0
  m_p_t["infected_presymptomatic_isolated"  , M_t == "infected_presymptomatic"  ]       <- 0 #(1 - p_Admit_base) # * p_detect_asympt
  m_p_t["infected_symptomatic_isolated"      , M_t == "infected_presymptomatic"  ]      <- 0
  m_p_t["hospitalized"                      , M_t == "infected_presymptomatic"  ]       <- 0 # p_Admit_base 
  m_p_t["recovered"                         , M_t == "infected_presymptomatic"  ]       <- 0
  
  
  #infected_sym
  m_p_t["susceptible"                       , M_t == "infected_symptomatic"  ]           <- 0
  m_p_t["infected_asymptomatic"             , M_t == "infected_symptomatic"  ]           <- 0
  m_p_t["infected_presymptomatic"           , M_t == "infected_symptomatic"  ]           <- 0
  m_p_t["infected_symptomatic"               , M_t == "infected_symptomatic"  ]          <- (1 - p_Admit_base) #* (1 -  p_detect_sympt)
  m_p_t["susceptible_isolated"              , M_t == "infected_symptomatic"  ]           <- 0
  m_p_t["infected_asymptomatic_isolated"    , M_t == "infected_symptomatic"  ]           <- 0
  m_p_t["infected_presymptomatic_isolated"  , M_t == "infected_symptomatic"  ]           <- 0
  m_p_t["infected_symptomatic_isolated"      , M_t == "infected_symptomatic"  ]          <- 0 #(1 - p_Admit_base) #* p_detect_sympt
  m_p_t["hospitalized"                      , M_t == "infected_symptomatic"  ]           <- p_Admit_base 
  m_p_t["recovered"                         , M_t == "infected_symptomatic"  ]           <- 0
  
  #sus_iso
  m_p_t["susceptible"                       , M_t == "susceptible_isolated"   ]         <- 0
  m_p_t["infected_asymptomatic"             , M_t == "susceptible_isolated"   ]         <- 0
  m_p_t["infected_presymptomatic"           , M_t == "susceptible_isolated"   ]         <- 0
  m_p_t["infected_symptomatic"               , M_t == "susceptible_isolated"   ]         <- 0
  m_p_t["susceptible_isolated"              , M_t == "susceptible_isolated"   ]         <- 1-p_Sus_iso_Asy -p_Sus_iso_Pre
  m_p_t["infected_asymptomatic_isolated"    , M_t == "susceptible_isolated"   ]         <- p_Sus_iso_Asy 
  m_p_t["infected_presymptomatic_isolated"  , M_t == "susceptible_isolated"   ]         <- p_Sus_iso_Pre
  m_p_t["infected_symptomatic_isolated"      , M_t == "susceptible_isolated"   ]         <- 0
  m_p_t["hospitalized"                      , M_t == "susceptible_isolated"   ]         <- 0
  m_p_t["recovered"                         , M_t == "susceptible_isolated"   ]         <- 0
  
  #infected_asy iso
  m_p_t["susceptible"                       , M_t == "infected_asymptomatic_isolated" ]       <- 0
  m_p_t["infected_asymptomatic"             , M_t == "infected_asymptomatic_isolated" ]       <- 0
  m_p_t["infected_presymptomatic"           , M_t == "infected_asymptomatic_isolated" ]       <- 0
  m_p_t["infected_symptomatic"               , M_t == "infected_asymptomatic_isolated" ]       <- 0
  m_p_t["susceptible_isolated"              , M_t == "infected_asymptomatic_isolated" ]       <- 0
  m_p_t["infected_asymptomatic_isolated"    , M_t == "infected_asymptomatic_isolated" ]       <- 1
  m_p_t["infected_presymptomatic_isolated"  , M_t == "infected_asymptomatic_isolated" ]       <- 0
  m_p_t["infected_symptomatic_isolated"      , M_t == "infected_asymptomatic_isolated" ]       <- 0
  m_p_t["hospitalized"                      , M_t == "infected_asymptomatic_isolated" ]       <- 0
  m_p_t["recovered"                         , M_t == "infected_asymptomatic_isolated" ]       <- 0
  
  #infected presym iso
  m_p_t["susceptible"                       , M_t == "infected_presymptomatic_isolated" ]     <- 0
  m_p_t["infected_asymptomatic"             , M_t == "infected_presymptomatic_isolated" ]     <- 0
  m_p_t["infected_presymptomatic"           , M_t == "infected_presymptomatic_isolated" ]     <- 0
  m_p_t["infected_symptomatic"               , M_t == "infected_presymptomatic_isolated" ]     <- 0
  m_p_t["susceptible_isolated"              , M_t == "infected_presymptomatic_isolated" ]     <- 0
  m_p_t["infected_asymptomatic_isolated"    , M_t == "infected_presymptomatic_isolated" ]     <- 0
  m_p_t["infected_presymptomatic_isolated"  , M_t == "infected_presymptomatic_isolated" ]     <- 1 - p_Admit_base
  m_p_t["infected_symptomatic_isolated"      , M_t == "infected_presymptomatic_isolated" ]     <- 0
  m_p_t["hospitalized"                      , M_t == "infected_presymptomatic_isolated" ]     <- p_Admit_base
  m_p_t["recovered"                         , M_t == "infected_presymptomatic_isolated" ]     <- 0
  
  #infected sym iso
  m_p_t["susceptible"                       , M_t == "infected_symptomatic_isolated" ]    <- 0
  m_p_t["infected_asymptomatic"             , M_t == "infected_symptomatic_isolated" ]    <- 0
  m_p_t["infected_presymptomatic"           , M_t == "infected_symptomatic_isolated" ]    <- 0
  m_p_t["infected_symptomatic"               , M_t == "infected_symptomatic_isolated" ]    <- 0
  m_p_t["susceptible_isolated"              , M_t == "infected_symptomatic_isolated" ]    <- 0
  m_p_t["infected_asymptomatic_isolated"    , M_t == "infected_symptomatic_isolated" ]    <- 0
  m_p_t["infected_presymptomatic_isolated"  , M_t == "infected_symptomatic_isolated" ]    <- 0
  m_p_t["infected_symptomatic_isolated"      , M_t == "infected_symptomatic_isolated" ]    <- 1 - p_Admit_base
  m_p_t["hospitalized"                      , M_t == "infected_symptomatic_isolated" ]    <- p_Admit_base
  m_p_t["recovered"                         , M_t == "infected_symptomatic_isolated" ]    <- 0
  
  #hospitalized
  m_p_t["susceptible"                       , M_t == "hospitalized"]                     <- 0
  m_p_t["infected_asymptomatic"             , M_t == "hospitalized"]                     <- 0
  m_p_t["infected_presymptomatic"           , M_t == "hospitalized"]                     <- 0
  m_p_t["infected_symptomatic"               , M_t == "hospitalized"]                     <- 0
  m_p_t["susceptible_isolated"              , M_t == "hospitalized"]                     <- 0
  m_p_t["infected_asymptomatic_isolated"    , M_t == "hospitalized"]                     <- 0
  m_p_t["infected_presymptomatic_isolated"  , M_t == "hospitalized"]                     <- 0
  m_p_t["infected_symptomatic_isolated"      , M_t == "hospitalized"]                     <- 0
  m_p_t["hospitalized"                      , M_t == "hospitalized"]                     <- 1
  m_p_t["recovered"                         , M_t == "hospitalized"]                     <- 0
  
  
  #recovered
  m_p_t["susceptible"                       , M_t == "recovered"]                        <- 0
  m_p_t["infected_asymptomatic"             , M_t == "recovered"]                        <- 0
  m_p_t["infected_presymptomatic"           , M_t == "recovered"]                        <- 0
  m_p_t["infected_symptomatic"               , M_t == "recovered"]                        <- 0
  m_p_t["susceptible_isolated"              , M_t == "recovered"]                        <- 0
  m_p_t["infected_asymptomatic_isolated"    , M_t == "recovered"]                        <- 0
  m_p_t["infected_presymptomatic_isolated"  , M_t == "recovered"]                        <- 0
  m_p_t["infected_symptomatic_isolated"      , M_t == "recovered"]                        <- 0
  m_p_t["hospitalized"                      , M_t == "recovered"]                        <- 0
  m_p_t["recovered"                         , M_t == "recovered"]                        <- 1
  
  
  return(t(m_p_t))
}       