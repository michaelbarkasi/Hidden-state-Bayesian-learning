
#### Preamble ##########################################################################################################
#### 
#### Title: Simulation script (in R) for Hidden-state Bayesian simulation of online motor 
####    responses to sonification feedback.
#### Description: This script contains the code for generating simulated reach trials. 
####    It ends by producing a data frame with the results of the simulation.
####    There will be a second script written which will then analyse the error, learning, 
####    and jerk in the simulation, using the LMM approach from the original RS4 data analysis. 
####    
#### Update (Jan 12, 2024): I have not made that second script yet and probably won't,
####    as I've set this project aside. Also, even this script is not finished. It runs
####    as a demo, but the code for all intended options is not yet written.
####    
#### Author: Michael Barkasi

# Want to clear R session and start fresh?
if (TRUE) rm(list = ls())

# For reproducability
set.seed(1234)

run_parallel <- TRUE # Set to TRUE if running on RIS cluster or Mac, FALSE if running on Windows
if (run_parallel) {
  if(!require(doMC)) {
    install.packages("doMC")
    library(doMC)
  }
  core_num <- detectCores() - 2
  registerDoMC(cores=core_num) 
}

#### Setup #############################################################################################################

# Begin R script by installing (if needed) and loading required packages. 
if(!require(plyr)) {
  install.packages("plyr")
  library(plyr)
}
if(!require(foreach)) {
  install.packages("foreach")
  library(foreach)
}
if(!require(Rcpp)) {
  install.packages("Rcpp")
  library(Rcpp)
}
if(!require(dtw)) {
  install.packages("dtw")
  library(dtw)
}

# Next, must source all C++ files. 
sourceCpp("RS4-Bayes-script-Cpp.cpp")

# Next, import data from reach study 4 which will serve as the basis for our simulation. 
if(!require(readr)) {
  install.packages("readr")
  library(readr)
}
motionSU <- read_csv(
  "reach-study-4-motionSU-final-p15only.csv", 
  col_types = cols(.default = col_double())
  )
pTable <- read_csv(
  "reach-study-4-pTable-final.csv", 
  col_types = cols(Condition = col_character(), 
  sex = col_character(),
  .default = col_double())
  )

#### Define simulation object classes ##################################################################################
#### Description: These are objects which will store the component parts of our simulations.

# Simulation Parameters
hypothesis_space_size <- 25 # Number of reaches to consider as possible model reaches
  # ... EVENTUALLY WILL SET TO 50 (25 FROM RS4 DATA, 25 SIMULATED), BUT LEAVE AT 10 FOR MWE AND DEBUGGING
num_of_trials <- 30 # How many reaches should each participant perform during a simulated session?
  # ... EVENTUALLY CHANGE THIS TO 75 TO INCLUDE SWITCH BLOCK, BUT LEAVE AT 50 FOR MWE AND DEBUGGING
switch_trial <- 51 # At which trial should the feedback be switched from online to terminal or vice versa? 

# Object Class #1: Hypothesis Space
slot_names_r <- paste0("r", 1:hypothesis_space_size) # Slots for reaches in the hypothesis space
slot_names_p <- paste0("P", 1:hypothesis_space_size) # Slots to hold current (prior) probabilities for each reach
setClass(
  "HypothesisSpace",
  slots = setNames(
    c(rep("matrix", hypothesis_space_size), rep("numeric", hypothesis_space_size)), 
    c(slot_names_r, slot_names_p)
    )
  )

# Example: if we have a hypothesis space H, we can access a reach in it by calling H@r1, H@r2, etc.
# Demo of how to access and modify slots in our new class (useful in loops): 
#   new_H <- new("HypothesisSpace")
#   slot(new_H, paste0("r", 1)) <- matrix(1:4, nrow = 2, ncol = 2)
#   slot_value <- slot(new_H, paste0("r", 1))

# Object Class #2: Participant
setClass(
  "Participant",
    slots = c(
      number = "numeric",
      motor_noise_level = "numeric",  # Amount of jitter to add to each reach in the hypothesis space, representing motor noise. 
      forgetfulness = "numeric",      # How quickly does the participant forget low probability reaches? 
      learning_ability = "numeric",   # How good is the participant at inferring useful new reaches to test?
      other_parameters = "character", # To be filled in as we go. 
      curiousity = "numeric",         # How much does the participant want to explore new reaches?
      hypothesis_space = "HypothesisSpace"
    )
  )

# Object Class #3: Session
# ... Objects in this class will be simulated experimental/learning sessions, with a specified simulated participant and number of trials (reaches). 
slot_names_hs <- paste0("hs", 1:num_of_trials) # Slots for saving hypothesis space at start of each trial in session
slot_names_t <- paste0("t", 1:num_of_trials)   # Slots for trials (sensor quaternions) in the session
slot_names_st <- paste0("st", 1:num_of_trials) # Slots for trials (spatial coordinates of wrist) in the session
slot_names_s <- paste0("s", 1:num_of_trials)   # Slots for sonification feedback in the session
setClass(
  "Session",
  slots = setNames(
    c("Participant", 
      "numeric", 
      "numeric", 
      "matrix", 
      "list", 
      "matrix", 
      rep("HypothesisSpace", num_of_trials), 
      rep("matrix", num_of_trials), 
      rep("matrix", num_of_trials), 
      rep("matrix", num_of_trials),
      "matrix",
      "matrix"
      ), 
    c("participant", 
      "current_trial", 
      "hyp_num_tried", 
      "model_son", 
      "spatial_model", 
      "model", 
      slot_names_hs, 
      slot_names_t, 
      slot_names_st, 
      slot_names_s,
      "hspace_error",
      "hspace_posteriors"
      )
    )
  )

#### Core (Helper) Functions ###########################################################################################

# We need to define several core functions which will be used to simulate trials and process the data. 

# Function #0A: Function for printing session status
print_session_status <- function(simulated_session, return_errors = TRUE) {
  
  cat("\n\nParticipant Number:", simulated_session@participant@number)
  cat("\nCurrent Trial Number:", simulated_session@current_trial)
  cat("\nHypothesis Number tried this Trial:", simulated_session@hyp_num_tried)
  cat("\nModel Reach Length:", nrow(simulated_session@model))
  cat("\nModel Sonification Length:", nrow(simulated_session@model_son))
  
  m <- simulated_session@model
  m_final1 <- m[nrow(m),1:4]
  m_final2 <- m[nrow(m),5:8]
  R <- slot(simulated_session, paste0("t", simulated_session@current_trial))
  final_position1 <- R[nrow(R),1:4]
  final_position2 <- R[nrow(R),5:8]
  d <- sqrt(sum((m_final1 - final_position1)^2)) + sqrt(sum((m_final2 - final_position2)^2))
  cat("\nTarget distance of current trial:", d)
  
  distances <- c()
  for ( i in 1:hypothesis_space_size ) {
    Ri <- slot(simulated_session@participant@hypothesis_space, paste0("r", i))
    final_position1 <- Ri[nrow(Ri),1:4]
    final_position2 <- Ri[nrow(Ri),5:8]
    d <- sqrt(sum((m_final1 - final_position1)^2)) + sqrt(sum((m_final2 - final_position2)^2))
    distances <- c(distances, d)
  }
  target_errors <- round(distances, 3)
  probabilities <- sapply(1:hypothesis_space_size, function(i) {slot(simulated_session@participant@hypothesis_space, paste0("P", i))})
  probabilities_rounded <- round(probabilities, 3)
  probs_and_errors <- array(NA, dim=c(2,hypothesis_space_size))
  probs_and_errors[1,] <- probabilities_rounded
  probs_and_errors[2,] <- target_errors
  colnames(probs_and_errors) <- paste0("H", 1:ncol(probs_and_errors))  # Assign column names
  rownames(probs_and_errors) <- c("Prob", "Err")  # Assign row names
  order_index <- order(probs_and_errors["Err", ])
  probs_and_errors <- probs_and_errors[, order_index]
  cat("\nHypothesis Probabilities and Target Error:\n")
  print(probs_and_errors)
  
  plot(
    probabilities, 
    ylim=c(0,1), xlab="Hypothesis Number", ylab="posterior probability of hypothesis", 
    main=paste0("Hypothesis Probabilities, Trial ", simulated_session@current_trial)
    )
  
  #cat("\nSpatial Model (Current status):")
  #print(simulated_session@spatial_model)
  
  if (return_errors) return(distances)
  
}

# Function #0B: Function for simulating sonification feedback (sources C++),
# ... for use only with RS4 data directly; for testing and debugging. 
sonification_feedback <- function(p,r) {
  
  # Inputs:
  #  p: integer, the number of the participant from the RS4 data.
  #  r: integer, the number of the reach from p from the RS4 data to sonify.
  # Output: 
  #  ARM: A matrix with number of rows = length of r (one row per motion sample), with three columns: 
  #   1. amplitude level of the simulated sonification,
  #   2. pitch level of the simulated sonification,
  #   3. matched model number from real-time time-warping algorithm, used to compute amplitude level. 
  
  # Grab participant's model
  m <- pTable$model[which(pTable$number==p)]
  model_q1_ <- as.matrix(motionSU[motionSU$participant==p & motionSU$trials==m, c("qax", "qay", "qaz", "qar")])
  model_q2_ <- as.matrix(motionSU[motionSU$participant==p & motionSU$trials==m, c("qbx", "qby", "qbz", "qbr")])

  # Grab reach
  Sensor1_Quat_ <- as.matrix(motionSU[motionSU$participant==p & motionSU$trials==r, c("qax", "qay", "qaz", "qar")])
  Sensor2_Quat_ <- as.matrix(motionSU[motionSU$participant==p & motionSU$trials==r, c("qbx", "qby", "qbz", "qbr")])

  # Compute sonification
  #   Note: process_sample is sourced from C++
  APM <- simulate_sonification(
    Sensor1_Quat = Sensor1_Quat_, 
    Sensor2_Quat = Sensor2_Quat_,
    model_q1 = model_q1_,
    model_q2 = model_q2_
    )

  return(APM)
  
}

# Function #1: Function for simulating sonification feedback (sources C++)
terminal_sonification_feedback_of_current_trial <- function(simulated_session) {
  
  # Inputs:
  #  simulated_session: an object of class Session, 
  #     which will contain the current trial number, the reach made in that trial, and the participant's model reach.
  # Output: 
  #  ARM: A matrix with number of rows = length of r (one row per motion sample), with three columns: 
  #   1. amplitude level of the simulated sonification,
  #   2. pitch level of the simulated sonification,
  #   3. matched model number from real-time time-warping algorithm, used to compute amplitude level. 
  
  # Grab participant's model
  model <- simulated_session@model
  model_q1_ <- model[,1:4]
  model_q2_ <- model[,5:8]

  # Grab reach from the current trial
  reach <- slot(simulated_session, paste0("t", simulated_session@current_trial))
  Sensor1_Quat_ <- reach[,1:4]
  Sensor2_Quat_ <- reach[,5:8]
  
  # Compute sonification
  #   Note: process_sample is sourced from C++
  APM <- simulate_sonification(
    Sensor1_Quat = Sensor1_Quat_, 
    Sensor2_Quat = Sensor2_Quat_,
    model_q1 = model_q1_,
    model_q2 = model_q2_
    )
  
  return(APM)
  
}

# Function #2: Function for initiating a hypothesis space 
initiate_hypothesis_space <- function(real_participant_num) {
  
  # Inputs: 
  #   real_participant_num: integer, the number of the participant from the RS4 data to use 
  #     as the basis for the hypothesis space. NOTE: These are real participants from the old data,
  #     not simulated participants!
  # Output: 
  #   An object of class HypothesisSpace, with each slot filled with a reach based on data 
  #     from participant number participant_num.
  
  # Initiate an empty hypothesis space
  new_H <- new("HypothesisSpace")
  
  # Now fill with reaches based on data from participant number participant_num,
  #   and set the initial (prior) probability of each reach being the model to 1/hypothesis_space_size.
  reaches <- unique(na.omit(motionSU$reach[which(motionSU$participant==real_participant_num)]))
  for ( i in 1:hypothesis_space_size ) {
    # NOTE! WE'RE JUST GRABBING THE FIRST hypothesis_sample_size AVAILABLE REACHES FROM THE PARTICIPANT;
    #  IN THE FINISHED SIMULATION, WILL GRAB FIRST 25 FROM RS4 DATA (UP TO REACH 25 ONLY), 
    #  THEN SIMULATE 25 REACHES TO FILL OUT THE HYPOTHESIS SPACE. NEED REACH BIOMECHANICS SIMULATION CODE!
    slot(new_H, paste0("r", i)) <- as.matrix(
      motionSU[motionSU$participant==real_participant_num & motionSU$trials==reaches[i], 
      c("qax", "qay", "qaz", "qar", "qbx", "qby", "qbz", "qbr")])
    slot(new_H, paste0("P", i)) <- 1/hypothesis_space_size
  }

  # Return our initiated hypothesis space
  return(new_H)
  
}

# Function #3: Function for initiating a simulated participant
initiate_simulated_participant <- function(real_participant_num, p_num = 0) {
  
  # Inputs: 
  #   real_participant_num: integer, the number of the participant from the RS4 data to use 
  #     as the basis for the simulated participant. NOTE: These are real participants from the old data,
  #     not simulated participants!
  # Output: 
  #   An object of class Participant
  
  # Initiate an empty hypothesis space
  new_P <- new("Participant")
  
  # If simulated participant number not explicitly set, use real participant number
  if (p_num==0) {
    p_num <- real_participant_num
  }
  
  # Initiate the simulated participant's hypothesis space
  new_P@hypothesis_space <- initiate_hypothesis_space(real_participant_num)
  
  # Initiate other parameters for simulated participant
  new_P@number <- p_num
  new_P@motor_noise_level <- 0.5 # replace, make a parameter that's randomly generated from a distribution we can set in simulation 
  new_P@forgetfulness <- 5       # replace, make a parameter that's randomly generated from a distribution we can set in simulation  (integer number)
  new_P@learning_ability <- 0.1  # replace, make a parameter that's randomly generated from a distribution we can set in simulation 
  new_P@curiousity <- 3          # replace, make a parameter that's randomly generated from a distribution we can set in simulation
  # ... Initiate other parameters 
  
  # Return our initiated simulated participant
  return(new_P)

}

# Function #4: Function to compute difference between time-warped reach paths,
#   used to find likelihoods in Bayesian update of hypothesis space.
#   This function (TWD) essentially quantifies the "match" between the sonified reach and sonified model. 
quat_step_distance <- function(index, matrix1) {
  row1 <- matrix1[index,]
  row2 <- c(0,0,0,1)
  if (index>1) row2 <- matrix1[index-1,]
  return(sqrt(sum((row1 - row2)^2)))
}
TWD <- function(r, m, r_son, m_son) {
  
  q1_r <- r[,1:4]
  q2_r <- r[,5:8]
  q1_m <- m[,1:4]
  q2_m <- m[,5:8]

  q1_rv <- aaply(1:nrow(q1_r),1,quat_step_distance,matrix1=q1_r,.parallel=run_parallel)
  q2_rv <- aaply(1:nrow(q2_r),1,quat_step_distance,matrix1=q2_r,.parallel=run_parallel)
  q1_mv <- aaply(1:nrow(q1_m),1,quat_step_distance,matrix1=q1_m,.parallel=run_parallel)
  q2_mv <- aaply(1:nrow(q2_m),1,quat_step_distance,matrix1=q2_m,.parallel=run_parallel)
  
  rv <- (q1_rv + q2_rv) * 1000 # put into qu/s
  mv <- (q1_mv + q2_mv) * 1000 # put into qu/s
  
  rv[which(is.na(rv))] <- 0
  mv[which(is.na(mv))] <- 0
  
  d <- dtw(mv,rv)
  mi <- d$index1
  ri <- d$index2
  
  indices <- matrix(1:length(mi), ncol = 1)
  amp_diff <- apply(indices, 1, function(i) {abs(m_son[mi[i],1] - r_son[ri[i],1])/max(1,r_son[ri[i],1],na.rm=TRUE)})
  pitch_diff <- apply(indices, 1, function(i) {abs(m_son[mi[i],2] - r_son[ri[i],2])/max(1,r_son[ri[i],2],na.rm=TRUE)})
  
  amp_diff <- mean(amp_diff, na.rm = TRUE)
  pitch_diff <- mean(pitch_diff, na.rm = TRUE)
  
  return(amp_diff + pitch_diff)
  
}

# Function #4: Function for updating a hypothesis space H between trials,
#   based on the simulated participant's motor noise, forgetfulness, and learning ability. 
base_update_hypothesis_space <- function(simulated_participant,last_hypothesis_number_tried) {
  
  # "base" = "MNFL" = Motor Noise, Forget and Learning"
  
  # Inputs: 
  #   simulated_participant: an object of class Participant, 
  #     which will contain the hypothesis space to be updated and other parameters, 
  #     such as motor noise, forgetfulness, and learning ability needed for the MNFL update.
  # Output: 
  #   A copy of that participant, with the hypothesis space updated. 
  
  # (1) First, "jitter" each reach in the hypothesis space by motor_noise, representing motor noise and forgetfulness.
  for ( i in 1:hypothesis_space_size ) {
    h_q1 <- slot(simulated_participant@hypothesis_space, paste0("r", i))[,1:4]
    h_q2 <- slot(simulated_participant@hypothesis_space, paste0("r", i))[,5:8]
    j1 <- c(
      runif(n=1,max=simulated_participant@motor_noise_level),
      runif(n=1,min=-1,max=1),
      runif(n=1,min=-1,max=1),
      runif(n=1,min=-1,max=1)
      )
    j2 <- c(
      runif(n=1,max=simulated_participant@motor_noise_level),
      runif(n=1,min=-1,max=1),
      runif(n=1,min=-1,max=1),
      runif(n=1,min=-1,max=1)
      )
    slot(simulated_participant@hypothesis_space, paste0("r", i)) <- jitter_quat_trajectory( h_q1, h_q2, j1, j2 )
    # just using jitter_quat_trajectory function as a dummy place holder, must write real function, with motor_noise as a parameter
  }
  
  # (2) Second, forget the lowest probability reach in the hypothesis space, 
  #   and replace with a new simulated reach, created by jittering the last performed reach. 

  # Step 1: Find lowest nonzero probability
  probabilities <- sapply(1:hypothesis_space_size, function(i) {slot(simulated_participant@hypothesis_space, paste0("P", i))})
  non_zero_indices <- which(probabilities != 0)
  num_of_lowest_probability_hypothesis <- non_zero_indices[which.min(probabilities[non_zero_indices])]
  # Step 2: Is this probability less than half that of the highest probability reach?  
  num_of_highest_probability_hypothesis <- order(probabilities)[hypothesis_space_size]
  if ( probabilities[num_of_lowest_probability_hypothesis] < 0.5*probabilities[num_of_highest_probability_hypothesis] && length(non_zero_indices) > 1 ) {
    # If so, forget the lowest probability reach
    slot(simulated_participant@hypothesis_space, paste0("P", num_of_lowest_probability_hypothesis)) <- 0.0
  } else if ( probabilities[num_of_highest_probability_hypothesis] < 0.3 ) { # Step 3: If not, replace it with a new reach to try. 
    # Grab the last tried reach and quantify its uncertainty. 
    last_hyp_num_tried <- slot(simulated_participant@hypothesis_space, paste0("r", last_hypothesis_number_tried))
    last_hyp_num_tried_q1 <- last_hyp_num_tried[,1:4]
    last_hyp_num_tried_q2 <- last_hyp_num_tried[,5:8]
    uncertainty <- 1 - slot(simulated_participant@hypothesis_space, paste0("P", last_hypothesis_number_tried))
    # Step 4: Replace the lowest probability reach with a new simulated reach, 
    #   gotten by jittering the last tried reach as a function of the uncertainty and scaled motor noise. Also, 
    #   set its probability to be the same as that of the last tried reach.
    exploration_factor <- uncertainty*simulated_participant@curiousity*simulated_participant@motor_noise_level
    j1e <- c(
      runif(n=1,min=-exploration_factor,max=exploration_factor),
      runif(n=1,min=-1,max=1),
      runif(n=1,min=-1,max=1),
      runif(n=1,min=-1,max=1)
      )
    j2e <- c(
      runif(n=1,min=-exploration_factor,max=exploration_factor),
      runif(n=1,min=-1,max=1),
      runif(n=1,min=-1,max=1),
      runif(n=1,min=-1,max=1)
      )
    slot(simulated_participant@hypothesis_space, paste0("r", num_of_lowest_probability_hypothesis)) <- jitter_quat_trajectory( last_hyp_num_tried_q1, last_hyp_num_tried_q2, j1e, j2e )
    slot(simulated_participant@hypothesis_space, paste0("P", num_of_lowest_probability_hypothesis)) <- slot(simulated_participant@hypothesis_space, paste0("P", last_hypothesis_number_tried))
  }
  # Renormalize probabilities
  new_probabilities <- sapply(1:hypothesis_space_size, function(i) {slot(simulated_participant@hypothesis_space, paste0("P", i))})
  for ( i in 1:hypothesis_space_size ) {
    slot(simulated_participant@hypothesis_space, paste0("P", i)) <- slot(simulated_participant@hypothesis_space, paste0("P", i))/sum(new_probabilities)
  }

  # Return simulated participant with updated hypothesis space
  return(simulated_participant)
  
}

# Function #5: Function for updating the hypothesis space using Baye's rule with 
#   sonification feedback as the evidence
bayesian_update_hypothesis_space <- function(simulated_session) {
  
  # Inputs: 
  #   simulated_session: an object of class Session, 
  #     which will contain the participant's hypothesis space to be updated, 
  #     as well as the sonification feedback from the current trial 
  #     and the session model reach, both needed for the Bayesian update.
  # Output: 
  #   A copy of that session, with the hypothesis space updated. 
  
  # CODE HERE FOR BAYESIAN UPDATE OF HYPOTHESIS SPACE BASED ON SONIFICATION FEEDBACK
  
  # Preliminary stuff needed to compute likelihoods of the evidence, for each hypothesis in the hypothesis space
  r_actual <- slot(simulated_session, paste0("t", simulated_session@current_trial))
  m_actual <- simulated_session@model
  r_actual_son_actual <- slot(simulated_session, paste0("s", simulated_session@current_trial))
  m_actual_son_actual <- simulated_session@model_son
  # Quantify evidence, i.e., match between heard sonification and heard model playback.
  TWD_actual <- TWD(r_actual, m_actual, r_actual_son_actual, m_actual_son_actual)
  if (!is.numeric(TWD_actual)) cat("\nWARNING: TWD_actual is not numeric!\n")
  
  sigmoid <- function(x,a=5,b=0.5) {1/(1+exp(-a*(x-b)))}

  # Compute unnomralized posterior probability for each reach in the hypothesis space
  prior_probabilities_H <- sapply(1:hypothesis_space_size, function(i) {slot(simulated_session@participant@hypothesis_space, paste0("P", i))})
  likelihoods_E <- c() 
  count_of_discriminable_hypotheses <- 0
  for ( i in 1:hypothesis_space_size ) {
    
    m_hypothesis <- slot(simulated_session@participant@hypothesis_space, paste0("r", i))
    m_hypothesis_son_hypothesis <- simulate_sonification(
      Sensor1_Quat = m_hypothesis[,1:4], 
      Sensor2_Quat = m_hypothesis[,5:8],
      model_q1 = m_hypothesis[,1:4],
      model_q2 = m_hypothesis[,5:8]
      )
    r_actual_son_hypothesis <- simulate_sonification(
      Sensor1_Quat = r_actual[,1:4], 
      Sensor2_Quat = r_actual[,5:8],
      model_q1 = m_hypothesis[,1:4],
      model_q2 = m_hypothesis[,5:8]
      )
    TWD_h <- TWD(r_actual, m_hypothesis, r_actual_son_hypothesis, m_hypothesis_son_hypothesis)
    if (!is.numeric(TWD_h)) cat("\nWARNING: TWD_h is not numeric! Hypothesis number: ", i, "\n")
    # Basic idea: Likelihood of a "good match" (TWD_actual=1) is 1, on the hypothesis that the actual reach made is the model. 
    #   As the match gets worse (TWD_actual decreases), the likelihood of that poor match, give the actual reach made is the model, 
    #   goes down. 
    if ( sigmoid ( abs(TWD_actual - TWD_h)/max(c(TWD_actual,TWD_h)) ) > 0.5 ) count_of_discriminable_hypotheses <- count_of_discriminable_hypotheses + 1
    likelihoods_E_h <- 1 - sigmoid( abs(TWD_actual - TWD_h)/max(c(TWD_actual,TWD_h)) )
    likelihoods_E <- c( likelihoods_E, likelihoods_E_h ) 
    
  }
  unnormalized_posterior_probabilities_H <- likelihoods_E * prior_probabilities_H
  
  # Normalize to find the posterior probabilities
  posterior_probabilities_H <- unnormalized_posterior_probabilities_H / sum(unnormalized_posterior_probabilities_H)
  
  # Save posteriors 
  simulated_session@hspace_posteriors <- rbind(simulated_session@hspace_posteriors, posterior_probabilities_H)
  
  # Update the hypothesis space with the new posterior probabilities
  for ( i in 1:hypothesis_space_size ) {
    slot(simulated_session@participant@hypothesis_space, paste0("P", i)) <- posterior_probabilities_H[i]
  }
  
  # Return simulated session with updated hypothesis space
  return(simulated_session)
  
  # Need to think about whether the same function will work for both terminal and online sessions. 
  
}

# Function #6: Function for simulating a reach and recording sonification feedback
simulate_reach <- function(simulated_session, online_feedback) {
  
  # Inputs: 
  #   simulated_session: an object of class Session, 
  #     which will contain the participant's hypothesis space from which to draw a preplanned reach,
  #     as well as the participant's model reach, needed for the simulation.
  #   online_feedback: logical, TRUE if simulating making the reach with online feedback, 
  #     FALSE if simulating making the reach with terminal terminal feedback. 
  # Output: 
  #   The updated session with the reach made saved and the heard sonification feedback recorded.
  
  # Select a reach from the hypothesis space to make for the current trial, baed on probability distribution of the hypothesis space
  prior_probabilities <- sapply(1:hypothesis_space_size, function(i) {slot(simulated_session@participant@hypothesis_space, paste0("P", i))})
  # Bias towards reaches with higher prior probabilities
  prior_probabilities <- prior_probabilities^4
  prior_probabilities <- prior_probabilities / sum(prior_probabilities) # Renormalize
  reach_number <- sample(1:hypothesis_space_size, size = 1, prob = prior_probabilities) 
  simulated_session@hyp_num_tried <- reach_number
  preplanned_reach <- slot(simulated_session@participant@hypothesis_space, paste0("r", reach_number))
  
  if (online_feedback) {
    
    # CODE HERE FOR SIMULATING ONLINE FEEDBACK
    
  } else {
   
    # Make preplanned reach without any online adjustments and save it as the reach made for the current trial
    slot(simulated_session, paste0("t", simulated_session@current_trial)) <- preplanned_reach
    
    # Simulate and record sonification feedback, i.e. "hear sonification".
    slot(simulated_session, paste0("s", simulated_session@current_trial)) <- terminal_sonification_feedback_of_current_trial(simulated_session) 
    
  }
  
  return(simulated_session)
  
}

# Function #7: Function for computing spatial trajectory from two-pivot biomechanical model
compute_spatial_reach_two_pivot <- function(simulated_session, trial_num, Hyp_space = FALSE) {
  
  # Update spatial model and save wrist position
  #   NOTE: For now, all of this is a placeholder using the simple two-pivot model, 
  #   just to get a MWE, until we have the full biomechanics simulation up and going.
  simulated_session@spatial_model <- Initiate_two_pivot_system() # This is a function from the biomechanics simulation, for now a placeholder
  #    using the basic two-pivot model, until OpenSim is up and going. 
  
  # To compute the spatial trajectory for the stored model, set trial_num to zero. 
  if ( trial_num == 0 ) {
    t_temp <- matrix(0, nrow = nrow(simulated_session@model), ncol = 3)
    for ( s in 1:nrow(t_temp) ) {
      quat1_now <- simulated_session@model[s,1:4]
      quat2_now <- simulated_session@model[s,5:8]
      quat1_prev <- c(0,0,0,1)
      quat2_prev <- c(0,0,0,1)
      if ( s>1 ) {
        quat1_prev <- simulated_session@model[s-1,1:4]
        quat2_prev <- simulated_session@model[s-1,5:8]
      }
      quat1_diff <- quat_diff(quat1_prev, quat1_now)
      quat2_diff <- quat_diff(quat2_prev, quat2_now)
      simulated_session@spatial_model <- compute_two_pivot_next_step(simulated_session@spatial_model, quat1_diff, quat2_diff) 
      t_temp[s,] <- simulated_session@spatial_model$t
    }
  } else if (!Hyp_space) {
    t_temp <- matrix(0, nrow = nrow(slot(simulated_session, paste0("t", trial_num))), ncol = 3)
    for ( s in 1:nrow(t_temp) ) {
      quat1_now <- slot(simulated_session, paste0("t", trial_num))[s,1:4]
      quat2_now <- slot(simulated_session, paste0("t", trial_num))[s,5:8]
      quat1_prev <- c(0,0,0,1)
      quat2_prev <- c(0,0,0,1)
      if ( s>1 ) {
        quat1_prev <- slot(simulated_session, paste0("t", trial_num))[s-1,1:4]
        quat2_prev <- slot(simulated_session, paste0("t", trial_num))[s-1,5:8]
      }
      quat1_diff <- quat_diff(quat1_prev, quat1_now)
      quat2_diff <- quat_diff(quat2_prev, quat2_now)
      simulated_session@spatial_model <- compute_two_pivot_next_step(simulated_session@spatial_model, quat1_diff, quat2_diff) 
      t_temp[s,] <- simulated_session@spatial_model$t
    }
  } else { 
    t_temp <- matrix(0, nrow = nrow(slot(simulated_session@hs1, paste0("r", trial_num))), ncol = 3)
    for ( s in 1:nrow(t_temp) ) {
      quat1_now <- slot(simulated_session@hs1, paste0("r", trial_num))[s,1:4]
      quat2_now <- slot(simulated_session@hs1, paste0("r", trial_num))[s,5:8]
      quat1_prev <- c(0,0,0,1)
      quat2_prev <- c(0,0,0,1)
      if ( s>1 ) {
        quat1_prev <- slot(simulated_session@hs1, paste0("r", trial_num))[s-1,1:4]
        quat2_prev <- slot(simulated_session@hs1, paste0("r", trial_num))[s-1,5:8]
      }
      quat1_diff <- quat_diff(quat1_prev, quat1_now)
      quat2_diff <- quat_diff(quat2_prev, quat2_now)
      simulated_session@spatial_model <- compute_two_pivot_next_step(simulated_session@spatial_model, quat1_diff, quat2_diff) 
      t_temp[s,] <- simulated_session@spatial_model$t
    }
  }
  
  # slot(simulated_session, paste0("st", trial_num)) <- t_temp
  
  return(t_temp)
  
}

#### Simulation ########################################################################################################

# NOTE: This is so far set up for terminal feedback online. 
#   Once we have that running, we'll need to add several parameters (the ones of real interest) 
#   as function variables, parameters which will govern when and how twitches are made in response to online feedback. 
#   At that point, might also need to revisit the terminal feedback, so that it perhaps effects these factors, to 
#   perhaps reproduce the switch effects for terminal learning feedback. 
simulate_session <- function(
    real_participant_num,           # Simulation based around this participant's random, no-feedback reaches. 
    learning_feedback = "terminal", # Must be "terminal" or "online", specifies feedback given in first 50 trials before switch. 
    print_progress = TRUE,          # If TRUE, will print out (at least) the trial number as it goes. 
    print_details = TRUE,           # If TRUE, will print out details from each trial. 
    save_hypothesis_space = FALSE   # Likely to take up a massive amount of space, but useful for debugging. 
    ) {
  
  cat("\nRunning Simulation of Participant", real_participant_num, "with", learning_feedback, "learning feedback.")
  if ( print_progress & !print_details ) cat("\nTrial: ")
  
  new_S <- new("Session")
  new_S@hspace_error <- matrix(0, nrow = 0, ncol = hypothesis_space_size) 
  new_S@hspace_posteriors <- matrix(0, nrow = 0, ncol = hypothesis_space_size)
  
  new_S@participant <- initiate_simulated_participant(real_participant_num)
  
  # The model reach (matrix) will be randomly selected from the hypothesis space
  model_number <- sample(1:hypothesis_space_size, 1)
  new_S@model <- slot(new_S@participant@hypothesis_space, paste0("r", model_number))
  new_S@current_trial <- model_number
  new_S@model_son <- simulate_sonification(
    Sensor1_Quat = new_S@model[,1:4], 
    Sensor2_Quat = new_S@model[,5:8],
    model_q1 = new_S@model[,1:4],
    model_q2 = new_S@model[,5:8]
    )
  
  for ( i in 1:num_of_trials ) {
    
    if ( print_progress & !print_details ) {
      if ( i == num_of_trials ) {
        cat(i)
      } else {
        cat(i, ", ")
      }
    } 
    
    # Keep track of current trial
    new_S@current_trial <- i
    
    # Save hypothesis space at the start of the trial
    #   Have to see: This may take up a massive amount of space, e.g., if we have 50 reaches
    #   in the hypothesis space and 75 trials, that will be 3,750 + 75 = 3,825 matrices saved, each with
    #   approximately 1000 rows and at least 8 columns. 
    if (save_hypothesis_space) {
      slot(new_S, paste0("hs", i)) <- new_S@participant@hypothesis_space
    } else if (i==1) {
      # We definitely want to always at least save the initial hs, as this is needed to compute stuff later, 
      #   such as normalized error and normalized jerk. 
      slot(new_S, paste0("hs", i)) <- new_S@participant@hypothesis_space
    }
    
    # Simulate and record a reach, i.e. "make reach". 
    #   Note: This will also simulate and record the sonification feedback, either as online or terminal. 
    if ( learning_feedback == "terminal" ) {
      if ( i < switch_trial ) {
        online_feedback <- FALSE
      } else if ( i >= switch_trial ) {
        online_feedback <- TRUE
      }
    } else if ( learning_feedback == "online" ) {
      if ( i < switch_trial ) {
        online_feedback <- TRUE
      } else if ( i >= switch_trial ) {
        online_feedback <- FALSE
      }
    } else {
      stop("learning_feedback must be either 'terminal' or 'online'.")
    }
    new_S <- simulate_reach(new_S, online_feedback = online_feedback)
    
    if (print_details) {
      errors <- print_session_status(new_S) # ... and return errors
      new_S@hspace_error <- rbind(new_S@hspace_error, errors)
    }
    
    # Update spatial model and save wrist position
    #   NOTE: For now, all of this is a placeholder using the simple two-pivot model, 
    #   just to get a MWE, until we have the full biomechanics simulation up and going.
    slot(new_S, paste0("st", i)) <- compute_spatial_reach_two_pivot(new_S, i)
    
    # Perform a Bayesian update of the hypothesis space based on the sonification feedback
    new_S <- bayesian_update_hypothesis_space(new_S)
    
    # Perform a base update of the hypothesis space to be applied for the next trial
    new_S@participant <- base_update_hypothesis_space(new_S@participant, new_S@hyp_num_tried)
    
  }
  
  if (print_details) {
    cat("\n\nSimulation Complete.\n")
    cat("\nFinal Session Status:\n")
    print_session_status(new_S)
  }
  
  # Note: The final hypothesis space can be extracted from the participant object within the session object, 
  #   at the end of the simulation.
  
  return(new_S)
  
}

start_time <- Sys.time()
test_simulation <- simulate_session(real_participant_num = 15)
end_time <- Sys.time() 

run_duration <- end_time - start_time
cat("\n\nSimulation took", run_duration, " ", units(run_duration), "to run.")

#### Print results #####################################################################################################

analyze_final_hypothesis_space <- function(simulated_session) {
  
  m <- simulated_session@model
  m_final1 <- m[nrow(m),1:4]
  m_final2 <- m[nrow(m),5:8]
  
  posteriors <- c()
  distances <- c()
  for ( i in 1:hypothesis_space_size ) {
    Pi <- slot(simulated_session@participant@hypothesis_space, paste0("P", i))
    Ri <- slot(simulated_session@participant@hypothesis_space, paste0("r", i))
    posteriors <- c(posteriors, Pi)
    final_position1 <- Ri[nrow(Ri),1:4]
    final_position2 <- Ri[nrow(Ri),5:8]
    d <- sqrt(sum((m_final1 - final_position1)^2)) + sqrt(sum((m_final2 - final_position2)^2))
    distances <- c(distances, d)
  }
  
  plot(x=distances, y=posteriors, xlab="hypothesis distance from model", ylab="posterior probability", main="Posterior Probability vs. Hypothesis Distance from Model")
  plot(x=1:length(distances), y=distances, xlab="Hypothesis Number", ylab="hypothesis distance from model", main="Hypothesis Distance from Model")
  plot(x=1:length(posteriors), y=posteriors, xlab="Hypothesis Number", ylab="posterior probability of hypothesis", main="Hypothesis Posterior Probability")
  
  num_trials <- nrow(simulated_session@hspace_posteriors)
  highest_hypothesis_confidences <- rep(NA, num_trials)
  total_error_weighted_by_confidence <- rep(NA, num_trials)
  for ( t in 1:num_trials ) {
    highest_hypothesis_confidences[t] <- simulated_session@hspace_posteriors[t,which.max(simulated_session@hspace_posteriors[t,])]
    total_error_weighted_by_confidence[t] <- sum(simulated_session@hspace_error[t,]*simulated_session@hspace_posteriors[t,])
  }
  
  par(mar = c(5, 4, 4, 5) + 0.1)
  plot(
    x=1:num_trials,
    y=highest_hypothesis_confidences,
    xlab="simulated trial number",
    ylab="highest posterior probability in hypothesis space",
    main="Highest Posteriors (blue) and Mean Weighted Error (red) by Trial",
    col="blue",
    type="l",
    ylim=c(0,1)
    )
  par(new = TRUE)
  plot(
    x=1:num_trials, 
    y=total_error_weighted_by_confidence,
    type = "l", col = "red", axes = FALSE, xlab = "", ylab = ""
    )
  axis(4)
  mtext("mean weighted error of hypothesis space", side = 4, line = 3)  # Label for the second y-axis
  
}
analyze_final_hypothesis_space(test_simulation)

if(!require(rgl)) {
  install.packages("rgl")
  library(rgl)
}
spatial_plot <- function(simulated_session, learning_feedback = "terminal") {
  
  C1_color <- "blue3"        # feedback color (reaches 26-75) for condition 1
  C1_color_light <- "azure"
  C1_color_dark <- "blue4"
  C2_color <- "red3"         # feedback color (reaches 26-75) for condition 2
  C2_color_light <- "lightpink"
  C2_color_dark <- "red4"
  
  if ( learning_feedback == "online" ) {
    colorgradFunc <- colorRampPalette(c(C1_color_light,C1_color_dark))
    print_colors <- C1_color
  } else if ( learning_feedback == "terminal" ) {
    colorgradFunc <- colorRampPalette(c(C2_color_light,C2_color_dark))
    print_colors <- C2_color
  }
  
  # Plot Model
  s_model <- compute_spatial_reach_two_pivot(simulated_session, 0)
  plot3d(
    s_model[,1],
    s_model[,2],
    s_model[,3],
    col = "orange",xlab="x (forward)",ylab="y (vertical)",zlab="z (left-right)"
  )
  
  # Plot reaches in the initial hypothesis space ("no feedback, random")
  for ( i in 1:hypothesis_space_size ) {
    s_reach <- compute_spatial_reach_two_pivot(simulated_session, i, Hyp_space=TRUE)
    plot3d(
      s_reach[,1],
      s_reach[,2],
      s_reach[,3],
      col = "grey", add=TRUE
    )
  }
  
  # Plot reaches made in the simulation (stimulated learning trials)
  colorgrad <- colorgradFunc(num_of_trials)
  for ( i in 1:num_of_trials ) {
    plot3d(
      slot(simulated_session, paste0("st", i))[,1],
      slot(simulated_session, paste0("st", i))[,2],
      slot(simulated_session, paste0("st", i))[,3],
      col = colorgrad[i], add=TRUE
    )
  }
  
  # Plot the points stored in the initial-position spatial model
  x <- unlist(lapply(Initiate_two_pivot_system(), function(x) x[1]))
  y <- unlist(lapply(Initiate_two_pivot_system(), function(x) x[2]))
  z <- unlist(lapply(Initiate_two_pivot_system(), function(x) x[3]))
  plot3d(x, y, z, col = "black", add=TRUE)
  
  # Correct aspect ratio
  aspect3d(1, 1, 1)
  
}
spatial_plot(test_simulation)

