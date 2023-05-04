get_similar_data <- function(patient_data, compare_data, distance, amount = NULL, threshold = NULL) {
  if (!is.null(amount) && amount == 1) return(patient_data %>% mutate(cur_user_id = user_id))
  if (is.null(amount) && is.null(threshold)) return(patient_data %>% mutate(cur_user_id = user_id))
  cur_id <- patient_data %>% slice_head() %>% pull(user_id)
  compare_data <- compare_data %>% mutate(cur_user_id = cur_id)

  similarity_per_id <- compare_data %>% group_by(user_id) %>%
    summarize("similarity" = distance(patient_data, cur_data_all())) %>%
    arrange(similarity) 
  if (!is.null(amount)) {
    most_similar <- similarity_per_id %>% slice_head(n=amount)
  }
  else {
    most_similar <- similarity_per_id %>% filter(similarity <= threshold)
  }
  
  similar_data <- compare_data %>% inner_join(most_similar, by=c("user_id"))
  
  return(similar_data)
}

select_features <- function(data, usedGoal, usedValues) {
  if (usedValues == "date") return(data %>% select(local_date)) 
  if (usedValues == "select") {
    cur_user_data <- data %>% filter(user_id == cur_user_id)
    subset <- NULL
    if (usedGoal == "std_cumberness") {
      subset <- summary(regsubsets(std_cumberness ~ loudness + jawbone + neck + tin_day + tin_cumber + tin_max + movement + stress + emotion,
                                   data=cur_user_data, method="exhaustive"), matrix.logical=TRUE) 
    }
    else if (usedGoal == "std_residuals" ) {
      subset <- summary(regsubsets(std_residuals ~ loudness + jawbone + neck + tin_day + tin_cumber + tin_max + movement + stress + emotion,
                                   data=cur_user_data, method="exhaustive"), matrix.logical=TRUE) 
    }
    
    selected_features <- data.frame("include" =subset$which[which.max(subset$adjr2),]) %>% #subset feature selection with highest adjr2 score
      rownames_to_column(var="features") %>% 
      filter(features != "(Intercept)") %>% 
      filter(include) %>% 
      pull(features)
    
    return(data %>% select(selected_features))
  }
  
  return(NULL)
}

regression <- function(train_data, test_data, attributes, lm_model, usedGoal, usedValues, distance, amount = NULL, threshold=NULL) {
  
  
  correlations <- get_correlations(train_data, usedGoal, attributes)
  
  #add all correlations
  correlations <- correlations %>% select(user_id, term, correlation) %>% pivot_wider(names_from=term, names_prefix="corr_", values_from="correlation")
  train_data <- train_data %>% left_join(correlations, by=c("user_id"))
  
  #calculate intercept, slope
  models <- train_data %>% group_by(user_id) %>%
    summarize(get_similar_data(cur_data_all(), train_data, distance, amount = amount, threshold=threshold) %>%
                       summarize(model = list(fit_xy(lm_model, x= select_features(cur_data_all(), usedGoal, usedValues), y= cur_data_all() %>%
                                                pull(usedGoal)))))
  
  train_data <- train_data %>% mutate("split"="train") %>% right_join(models, by=c("user_id"))
  test_data <- test_data %>% mutate("split"="test") %>% right_join(models, by=c("user_id"))
  
  pred <- function(data) {
    return(predict(data %>% pull(model) %>% first(), data) %>% pull(.pred))
  }
  
  #calculate residuals
  train_data <- train_data %>% group_by(user_id) %>% mutate(prediction=pred(cur_data())) %>% ungroup() %>%
    mutate(std_residuals=std_cumberness-prediction) %>%
    mutate(residuals=std_residuals*cumberness_stddev)
  test_data <-  test_data  %>% group_by(user_id) %>% mutate(prediction=pred(cur_data())) %>% ungroup() %>%
    mutate(std_residuals=std_cumberness-prediction) %>%
    mutate(residuals=std_residuals*cumberness_stddev)
  
  train_data <- train_data %>% select(!c("model")) %>% select(!starts_with("corr_"))
  test_data <- test_data %>% select(!c("model"))
  
  return(list(train_data, test_data))
}

regression_all <- function(train_data, test_data, attributes, lm_model, usedGoal, usedValues) {
  
  
  #calculate intercept, slope
  model <- fit_xy(lm_model, x= select_features_all(train_data, usedGoal, usedValues), 
                  y= train_data %>% pull(usedGoal))
  
  train_data <- train_data %>% mutate("split"="train")
  test_data <- test_data %>% mutate("split"="test")

  pred <- function(data, model) {
    return(predict(model, data) %>% pull(.pred))
  }
  
  #calculate residuals
  train_data <- train_data %>% group_by(user_id) %>% mutate(prediction=pred(cur_data(), model)) %>% ungroup() %>%
    mutate(std_residuals=std_cumberness-prediction) %>%
    mutate(residuals=std_residuals*cumberness_stddev)
  test_data <-  test_data  %>% group_by(user_id) %>% mutate(prediction=pred(cur_data(), model)) %>% ungroup() %>%
    mutate(std_residuals=std_cumberness-prediction) %>%
    mutate(residuals=std_residuals*cumberness_stddev)
  
  return(list(train_data, test_data))
}

select_features_all <- function(data, usedGoal, usedValues) {
  if (usedValues == "date") return(data %>% select(local_date)) 
  if (usedValues == "select") {
    cur_user_data <- data
    subset <- NULL
    if (usedGoal == "std_cumberness") {
      subset <- summary(regsubsets(std_cumberness ~ loudness + jawbone + neck + tin_day + tin_cumber + tin_max + movement + stress + emotion,
                                   data=cur_user_data, method="exhaustive"), matrix.logical=TRUE) 
    }
    else if (usedGoal == "std_residuals" ) {
      subset <- summary(regsubsets(std_residuals ~ loudness + jawbone + neck + tin_day + tin_cumber + tin_max + movement + stress + emotion,
                                   data=cur_user_data, method="exhaustive"), matrix.logical=TRUE) 
    }
    
    selected_features <- data.frame("include" =subset$which[which.max(subset$adjr2),]) %>% #subset feature selection with highest adjr2 score
      rownames_to_column(var="features") %>% 
      filter(features != "(Intercept)") %>% 
      filter(include) %>% 
      pull(features)
    
    return(data %>% select(selected_features))
  }
  
  return(NULL)
}