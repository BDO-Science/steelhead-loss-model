library(tidyverse)
library(CDECRetrieve)
library(busdater)
library(zoo)
library(MASS)
library(MuMIn)
library(DHARMa)
start_date <- '2007-01-01'
end_date <- '2024-12-31'
stations <- list('MRC', 'TLG', 'SNS')

flow_list <- lapply(stations, function(station) {
  cdec_query(station, '66', 'M', start_date, end_date)
})


flows <- bind_rows(flow_list) %>%
  mutate(wy = get_fy(as.Date(datetime), opt_fy_start = '10-01')) %>%
  mutate(parameter_value = na.approx(parameter_value))

summary <- flows %>%
  group_by(location_id, wy) %>%
  summarize(af = sum(parameter_value)) %>%
  filter(wy < 2025)

summary2 <- summary %>%
  pivot_wider(names_from = 'location_id', values_from = 'af') %>%
  mutate(stan_prev = lag(SNS, 1),
         stan_prev2 = lag(SNS,2),
         tul_prev = lag(TLG, 1),
         tul_prev2 = lag(TLG, 2),
         mer_prev = lag(MRC, 1),
         mer_prev2 = lag(MRC, 2)) %>%
  filter(wy > 2008) %>%
  dplyr::select(-2:-4)

data <- readRDS('model_data.rds') %>%
  filter(wy >2008) %>%
  left_join(summary2, by = 'wy')

model_list <- list(
  #all covariates on the filtered dataframe WY > 2008
  glm.nb(loss ~ facility + OMR_mean + exports + 
           sac + sjr+
           stan_prev + stan_prev2 +
           tul_prev + tul_prev2 + mer_prev + mer_prev2, data = data,
         na.action = na.fail),
  glm.nb(loss ~ facility + OMR_mean + exports + 
           sac +
           stan_prev + stan_prev2 +
           tul_prev + tul_prev2 + mer_prev + mer_prev2, data = data,
         na.action = na.fail),
  glm.nb(loss ~ facility + OMR_mean + exports + 
           sac +
           stan_prev + stan_prev2 +
           tul_prev + tul_prev2 + mer_prev, data = data,
         na.action = na.fail),
  glm.nb(loss ~ facility + OMR_mean + exports + 
           sac +
           stan_prev + stan_prev2 +
           tul_prev2 + mer_prev, data = data,
         na.action = na.fail),
  glm.nb(loss ~ facility + OMR_mean + exports + 
           stan_prev + stan_prev2 +
           tul_prev2 + mer_prev, data = data,
         na.action = na.fail)
)


names(model_list) <- paste0("model", seq_along(model_list))
vifs <- lapply(model_list, car::vif)

full_model_no_interact <- glm.nb(loss ~ facility + OMR_mean + exports + 
                                   stan_prev + stan_prev2 +
                                   tul_prev2 + mer_prev, data = data,
                                 na.action = na.fail)

model_set_no_interact <- dredge(full_model_no_interact) #modeling all combinations of covariates

full_model_interact <- glm.nb(loss ~ facility + OMR_mean + exports + 
                                stan_prev + stan_prev2 +
                                tul_prev2 + mer_prev +
                                stan_prev:stan_prev2 + tul_prev2:mer_prev, data = data,
                              na.action = na.fail)

model_set_interact <- dredge(full_model_interact) #modeling all combinations of covariates

top_models <- get.models(model_set_interact, subset = 1:5)


#####k-fold method
cv_list <- lapply(top_models, function(m) {
  boot::cv.glm(data, m, K = 5)
})
# Extract delta values (raw CV error and adjusted CV error)
delta_df <- do.call(rbind, lapply(cv_list, function(x) x$delta))
colnames(delta_df) <- c("raw_cv_error", "adjusted_cv_error")

delta_df <- as.data.frame(delta_df) # Convert to data frame
delta_df$model <- names(top_models)  # optional: use model names if available

#####leave one out method
n <- nrow(data) # Number of observations

cv_rmse_loocv <- function(model_formula, data) { #Function to perform LOOCV
  rmse_values <- numeric(n) #return mean RMSE for a given model formula
  
  for(i in 1:n){
    train_data <- data[-i, ]
    test_data <- data[i, , drop = FALSE]
    model_cv <- glm.nb(model_formula, data = train_data) #Refit model on training data
    pred <- predict(model_cv, newdata = test_data, type = "response") #Predict on the left-out observation
    rmse_values[i] <- (test_data$loss - pred)^2 #Calculate RMSE for the single prediction
  }
  sqrt(mean(rmse_values))
}

rmse_results_loocv <- sapply(top_models, function(m) { #Apply to each model in the list
  cv_rmse_loocv(formula(m), data)
})

rmse_df <- data.frame( #put loocv results in dataframe
  model = as.character(names(rmse_results_loocv)),
  rmse = as.numeric(rmse_results_loocv)
)

#####cross validation results for top 5 models
cv_results <- delta_df %>% left_join(rmse_df, by = 'model') %>% #join k-fold and loocv dataframes
  dplyr::select(4,1,2,3) %>%
  arrange(rmse)

top_model_name <- cv_results %>% slice(1) %>%
  pull(model) #pull name of model with lowest rmse

#################################################
#model testing on top model
#################################################
top_model <- top_models[['56']] #extract top model from model list

#####drop1 to look at variance explained by covariates in top models
drop_1 <- drop1(top_model, test = "Chi") %>%
  arrange(desc(Deviance))

#####simulating and plotting residuals
res <- simulateResiduals(top_model)
plot(res)

library(here)
project <- here()
source(here(project,'functions.R'))
data_predictions <- add_predictions_with_ci(top_model, data)

top_model_predict <- data_predictions %>%
  ggplot(aes(x = wy, group = facility)) + ##graphing the predictions
  geom_point(aes(y = loss), size = 4) +
  geom_crossbar(aes(y = pred, ymin = lower_ci, ymax = upper_ci, fill = facility), alpha = 0.5) +
  labs(x = 'Water Year', y = 'Predicted Loss') +
  facet_wrap(~facility, nrow = 2, scales = 'free_y') +
  theme_bw() +
  theme(plot.margin = margin(0.5,0.5,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
top_model_predict

top_model_predict2 <- data_predictions %>%
  ggplot(aes(x = loss, y = pred)) + ##graphing the predictions
  geom_point(size = 3) +
  geom_smooth(method = 'lm') +
  labs(x = 'Observed Loss', y = 'Predicted Loss') +
  facet_wrap(~facility, nrow = 2, scales = 'free') +
  theme_bw() +
  theme(plot.margin = margin(0.5,0.5,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
top_model_predict2