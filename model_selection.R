library(tidyverse)
library(MASS)
library(MuMIn)
library(DHARMa)
library(corrplot)

project <- here()
source(here(project,'functions.R'))

data_all <- readRDS(here(project,'model_data.rds')) %>%
  mutate(era = if_else(wy <= 2008, 'Early', 'Current'))

data_filtered <- data_all %>%
  dplyr::filter(wy > 2008)

################
#quick corr plot
################
num_data <- data_filtered[, sapply(data_filtered, is.numeric)]
cor_matrix <- cor(num_data, use = "complete.obs")
corrplot <- corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45,
         addCoef.col = "darkgrey", number.cex = 0.8, number.digits = 1,
         number.font = 2,
         order = "hclust",
         addrect = 3)

######################################################################
#recursively removing covariates based on vif to reduce collinearity
#removing one covariate at a time with highest vif
#stop removing covariate when all vifs < 10
######################################################################
model_list <- list(
  #all covariates on the filtered dataframe WY > 2008
  glm.nb(loss ~ facility + OMR_mean + exports + 
                         sac + sjr + sj_index_prev + sj_index_prev2 +
                         sac_index_prev + sac_index_prev2, data = data_filtered,
         na.action = na.fail),
  #sjr had highest vif so removed on next model run
  glm.nb(loss ~ facility + OMR_mean + exports +
           sac + sj_index_prev + sj_index_prev2 +
           sac_index_prev + sac_index_prev2, data = data_filtered,
         na.action = na.fail),
  #sac_index_prev2 had highest vif so removed on next model run
  glm.nb(loss ~ facility + OMR_mean + exports +
           sac + sj_index_prev + sj_index_prev2 +
           sac_index_prev, data = data_filtered,
         na.action = na.fail)
  #all vifs are <10 after this model
)

names(model_list) <- paste0("model", seq_along(model_list))
vifs <- lapply(model_list, car::vif)


######################################################################
#using dredge on final model after removing collinearity
#included with and without interaction between sj_index_prev and sj_index_prev2
######################################################################

full_model_no_interact <- glm.nb(loss ~ facility + OMR_mean + exports + sac + 
                                   sac_index_prev + sj_index_prev + 
                                   sj_index_prev2, 
                                 data = data_filtered,
                                 na.action = na.fail)

model_set_no_interact <- dredge(full_model_no_interact) #modeling all combinations of covariates

full_model_interact <- glm.nb(loss ~ facility + OMR_mean + exports + sac + 
                                sj_index_prev + sj_index_prev2 + 
                                sac_index_prev + sj_index_prev:sj_index_prev2, 
                              data = data_filtered, na.action = na.fail)

model_set_interact <- dredge(full_model_interact) #modeling all combinations of covariates


top_models <- get.models(model_set_interact, subset = 1:5) #pulling top 5 models based on AICc from interaction model

##########################################
#cross-validation of all top models
#using both k-fold of 5 and leave one out
##########################################

#####k-fold method
cv_list <- lapply(top_models, function(m) {
  boot::cv.glm(data_filtered, m, K = 5)
})
# Extract delta values (raw CV error and adjusted CV error)
delta_df <- do.call(rbind, lapply(cv_list, function(x) x$delta))
colnames(delta_df) <- c("raw_cv_error", "adjusted_cv_error")

delta_df <- as.data.frame(delta_df) # Convert to data frame
delta_df$model <- names(top_models)  # optional: use model names if available

#####leave one out method
n <- nrow(data_filtered) # Number of observations

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
  cv_rmse_loocv(formula(m), data_filtered)
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
top_model <- top_models[[top_model_name]] #extract top model from model list

#####drop1 to look at variance explained by covariates in top models
drop_1 <- drop1(top_model, test = "Chi") %>%
  arrange(desc(Deviance))

#####simulating and plotting residuals
res <- simulateResiduals(top_model)
plot(res)

##########################################################################
#predicting top model on dataset and comparing actuals to predicted values
##########################################################################

data_predictions <- add_predictions_with_ci(top_model, data_filtered)

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
ggsave(top_model_predict, file = 'viz/prediction.png', height = 6, width = 8)

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
ggsave(top_model_predict2, file = 'viz/prediction_linear.png', height = 6, width = 8)

############################################
#model fit on full data frame of covariates
############################################
fit <- expand.grid(sj_index_prev = seq(0.5,7,0.5),
                           sj_index_prev2 = seq(0.5,7,0.5),
                         facility = c('CVP', 'SWP'),
                         exports = seq(2000,6000,500),
                         OMR_mean = seq(-4000,2000,500),
                         sac = seq(8000,50000,1000))

fit <- add_predictions_with_ci(top_model, fit)
