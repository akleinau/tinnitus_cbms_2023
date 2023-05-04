library(tidyverse)
library(feasts)
library(tsibble)
library(lubridate)
library(fable)

setup_diary <- function(data_version = '22_11_28') {
  data <- read_csv(paste("data/diary_", data_version, ".csv", sep=""))
  data <- data %>% filter(user_id <= 42101) %>% filter(user_id >= 2101) # filter out test users
  data$ts <- as.Date(data$created_at)
  data <- data %>% distinct(user_id, ts, .keep_all = TRUE)
  
  #add user names
  csv_mapping <- read_csv("data/ema_to_rct.csv")
  data <- data %>% left_join(csv_mapping, by=c('user_id'))
  data <- data %>% mutate(name=ifelse(is.na(name), '???', name))

  return(data)
}

setup_clinical <- function(data_version = '22_11_11') {
  clin_data <- read_csv(paste("data/longitudinal_", data_version, ".csv", sep=""))
  clin_data <- clin_data %>% mutate(treatment=str_sub(treatment_code, start=-1) %>% str_to_upper, name=external_id)
  #clin_data$name <- str_replace(clin_data$external_id, "-[a-zA-Z0-9]*-", "-")
  # sort visits
  clin_data$visit_type <- factor(clin_data$visit_type, 
                                 c('screening', 'baseline', 'interim_visit', 
                                   'final_visit', 'followup_1', 'followup_2'))
  return(clin_data)
}

setup_treatment_codes <- function() {
  treatment_codes <- read_csv("data/TreatmentCodes.csv")
  return(treatment_codes)
}

sample_data <- function(local_data, num_samples) {
  sample_name <- local_data %>% select(name) %>% distinct() %>% slice_sample(n=num_samples)
  sample <- local_data %>% filter (name %in% sample_name$name)
  return (sample)
}

split_inter_patient <- function(data, test_prop) {
  sample_name <- data %>% dplyr::select(user_id) %>% distinct() %>% slice_sample(prop=test_prop)
  sample_test <- data %>% filter (user_id %in% sample_name$user_id)
  sample_train <- data %>% filter (! user_id %in% sample_name$user_id)
  return (list(sample_train, sample_test))
}

split_intra_patient <- function(data, test_prop) {
  data <- data %>% rowid_to_column(var="rowid")
  sample_test <- data %>% group_by(user_id)  %>% slice_tail(prop=test_prop) %>% ungroup()
  sample_train <- data %>% filter (! rowid %in% sample_test$rowid) %>% select(!c("rowid"))
  sample_test <- sample_test %>% select(!c("rowid"))
  return (list(sample_train, sample_test))
}

standardize_cumberness <- function(train_data, test_data) {
  data_summary <- train_data %>% group_by(user_id) %>% summarize("cumberness_mean"=mean(cumberness), "cumberness_stddev"=sd(cumberness))
  train_data <- train_data %>% left_join(data_summary, by="user_id") %>% mutate(std_cumberness = (cumberness-cumberness_mean)/cumberness_stddev)
  test_data <- test_data %>% left_join(data_summary, by="user_id") %>% mutate(std_cumberness = (cumberness-cumberness_mean)/cumberness_stddev)
  return(list(train_data, test_data))
}

get_correlations <- function(data, y, attributes) {
  correlations <- data %>% group_by(user_id) %>% select(all_of(append(attributes, y))) %>% 
    summarize(correlate(x=cur_data(), diagonal=1, quiet=TRUE) %>% filter(term==y)) %>% 
    select(!all_of(c("term")))
  correlations <- correlations %>% pivot_longer(cols=all_of(append(attributes, y)), names_to="term", values_to="correlation") %>% 
    mutate(absCorrelation = abs(correlation)) %>% filter(term != y)
  return(correlations)
}

treatment_similarity <- function(treatment_codes, xn, yn) {
  x <- treatment_codes %>% filter(user_id == xn)
  y <- treatment_codes %>% filter(user_id == yn)
  score <- 0
  if (x$HA & y$HA) score = score + 1
  if (x$CBT & y$CBT) score = score + 1
  if (x$SC & y$SC) score = score + 1
  if (x$ST & y$ST) score = score + 1
  return(score)
}

standardize_attributes <- function(train_data, test_data) {
  #"loudness"
  data_summary <- train_data %>% group_by(user_id) %>% summarize("mean_"=mean(loudness), "stddev"=sd(loudness))
  train_data <- train_data %>% mutate(raw_loudness = loudness) %>%
    left_join(data_summary, by="user_id") %>% 
    mutate(loudness = (loudness-mean_)/stddev) %>%
    mutate(loudness = replace_na(loudness,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  test_data <- test_data %>% mutate(raw_loudness = loudness) %>%
    left_join(data_summary, by="user_id") %>% 
    mutate(loudness = (loudness-mean_)/stddev) %>%
    mutate(loudness = replace_na(loudness,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  
  #"jawbone"
  data_summary <- train_data %>% group_by(user_id) %>% summarize("mean_"=mean(jawbone), "stddev"=sd(jawbone))
  train_data <- train_data %>% mutate(raw_jawbone = jawbone) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(jawbone = (jawbone-mean_)/stddev) %>%
    mutate(jawbone = replace_na(jawbone,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  test_data <- test_data %>% mutate(raw_jawbone = jawbone) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(jawbone = (jawbone-mean_)/stddev) %>%
    mutate(jawbone = replace_na(jawbone,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  
  #"neck"
  data_summary <- train_data %>% group_by(user_id) %>% summarize("mean_"=mean(neck), "stddev"=sd(neck))
  train_data <- train_data %>% mutate(raw_neck = neck) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(neck = (neck-mean_)/stddev) %>%
    mutate(neck = replace_na(neck,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  test_data <- test_data %>% mutate(raw_neck = neck) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(neck = (neck-mean_)/stddev) %>%
    mutate(neck = replace_na(neck,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  
  #"tin_day"
  data_summary <- train_data %>% group_by(user_id) %>% summarize("mean_"=mean(tin_day), "stddev"=sd(tin_day))
  train_data <- train_data %>% mutate(raw_tin_day = tin_day) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(tin_day = (tin_day-mean_)/stddev) %>%
    mutate(tin_day = replace_na(tin_day,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  test_data <- test_data %>% mutate(raw_tin_day = tin_day) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(tin_day = (tin_day-mean_)/stddev) %>%
    mutate(tin_day = replace_na(tin_day,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  
  #"tin_cumber"
  data_summary <- train_data %>% group_by(user_id) %>% summarize("mean_"=mean(tin_cumber), "stddev"=sd(tin_cumber))
  train_data <- train_data %>% mutate(raw_tin_cumber = tin_cumber) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(tin_cumber = (tin_cumber-mean_)/stddev) %>%
    mutate(tin_cumber = replace_na(tin_cumber,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  test_data <- test_data %>% mutate(raw_tin_cumber = tin_cumber) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(tin_cumber = (tin_cumber-mean_)/stddev) %>%
    mutate(tin_cumber = replace_na(tin_cumber,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  
  #"tin_max"
  data_summary <- train_data %>% group_by(user_id) %>% summarize("mean_"=mean(tin_max), "stddev"=sd(tin_max))
  train_data <- train_data %>% mutate(raw_tin_max = tin_max) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(tin_max = (tin_max-mean_)/stddev) %>%
    mutate(tin_max = replace_na(tin_max,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  test_data <- test_data %>% mutate(raw_tin_max = tin_max) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(tin_max = (tin_max-mean_)/stddev) %>%
    mutate(tin_max = replace_na(tin_max,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  
  #"movement"
  data_summary <- train_data %>% group_by(user_id) %>% summarize("mean_"=mean(movement), "stddev"=sd(movement))
  train_data <- train_data %>% mutate(raw_movement = movement) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(movement = (movement-mean_)/stddev) %>%
    mutate(movement = replace_na(movement,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  test_data <- test_data %>% mutate(raw_movement = movement) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(movement = (movement-mean_)/stddev) %>%
    mutate(movement = replace_na(movement,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  
  #"stress"
  data_summary <- train_data %>% group_by(user_id) %>% summarize("mean_"=mean(stress), "stddev"=sd(stress))
  train_data <- train_data %>% mutate(raw_stress = stress) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(stress = (stress-mean_)/stddev) %>%
    mutate(stress = replace_na(stress,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  test_data <- test_data %>% mutate(raw_stress = stress) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(stress = (stress-mean_)/stddev) %>%
    mutate(stress = replace_na(stress,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  
  #"emotion"
  data_summary <- train_data %>% group_by(user_id) %>% summarize("mean_"=mean(emotion), "stddev"=sd(emotion))
  train_data <- train_data %>% mutate(raw_emotion = emotion) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(emotion = (emotion-mean_)/stddev) %>%
    mutate(emotion = replace_na(emotion,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  test_data <- test_data %>% mutate(raw_emotion = emotion) %>% 
    left_join(data_summary, by="user_id") %>% 
    mutate(emotion = (emotion-mean_)/stddev) %>%
    mutate(emotion = replace_na(emotion,0)) %>% #fix NA values when there is only 0 in the data
    select(!c("mean_", "stddev"))
  
  return(list(train_data, test_data))
}