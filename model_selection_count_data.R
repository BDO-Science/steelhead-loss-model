library(tidyverse)
library(MASS)
library(MuMIn)
library(DHARMa)
library(corrplot)
library(here)
library(patchwork)

project <- here()
source(here(project,'functions.R'))

data_all <- readRDS(here(project,'model_data.rds')) %>%
  mutate(era = if_else(wy <= 2008, 'Early', 'Current'))

data_count <- data_all %>% 
  group_by(wy) %>%
  summarize(count = sum(count),
            sj_index_prev2 = min(sj_index_prev2),
            sj_index_prev = min(sj_index_prev),
            sac_index_prev2 = min(sac_index_prev2),
            sac_index_prev = min(sac_index_prev),
            OMR_mean = min(OMR_mean),
            sjr = min(sjr),
            sac = min(sac),
            exports = min(exports))
#####################
#all years
#####################
data_count_filtered <- data_count %>%
  filter(wy >= 2009)

model_list_filtered <- list(
  #all covariates on the filtered dataframe WY > 2008
  glm.nb(count ~ OMR_mean + exports + 
           sac + sjr + sj_index_prev + sj_index_prev2 +
           sac_index_prev + sac_index_prev2, data = data_count_filtered,
         na.action = na.fail),
  #sjr had highest vif so removed on next model run
  glm.nb(count ~ OMR_mean + exports +
           sac + sj_index_prev + sj_index_prev2 +
           sac_index_prev + sac_index_prev2, data = data_count_filtered,
         na.action = na.fail),
  #sac_index_prev2 had highest vif so removed on next model run
  glm.nb(count ~ OMR_mean + exports +
           sac + sj_index_prev + sj_index_prev2 +
           sac_index_prev, data = data_count_filtered,
         na.action = na.fail),
  #all vifs are <10 after this model
)

names(model_list_filtered) <- paste0("model", seq_along(model_list_filtered))
vifs_filtered <- lapply(model_list_filtered, car::vif)

full_model_no_interact_filtered <- glm.nb(count ~ OMR_mean + exports +
                                            sj_index_prev + sj_index_prev2 +
                                            sac_index_prev, data = data_count_filtered,
                                          na.action = na.fail)

model_set_no_interact_filtered <- dredge(full_model_no_interact_filtered)

full_model_interact_filtered <- glm.nb(count ~ OMR_mean + exports +
                                         sj_index_prev + sj_index_prev2 +
                                         sac_index_prev + sj_index_prev:sj_index_prev2, data = data_count_filtered,
                                       na.action = na.fail)

model_set_interact_filtered <- dredge(full_model_interact_filtered)

top_models_filtered <- get.models(model_set_interact_filtered, subset = 1:10)

#####leave one out method
n <- nrow(data_count_filtered) # Number of observations

cv_rmse_loocv <- function(model_formula, data) { #Function to perform LOOCV
  rmse_values <- numeric(n) #return mean RMSE for a given model formula
  
  for(i in 1:n){
    train_data <- data[-i, ]
    test_data <- data[i, , drop = FALSE]
    model_cv <- glm.nb(model_formula, data = train_data) #Refit model on training data
    pred <- predict(model_cv, newdata = test_data, type = "response") #Predict on the left-out observation
    rmse_values[i] <- (test_data$count - pred)^2 #Calculate RMSE for the single prediction
  }
  sqrt(mean(rmse_values))
}

rmse_results_loocv <- sapply(top_models_filtered, function(m) { #Apply to each model in the list
  cv_rmse_loocv(formula(m), data_count_filtered)
})

rmse_df <- data.frame( #put loocv results in dataframe
  model = as.character(names(rmse_results_loocv)),
  rmse = as.numeric(rmse_results_loocv)
)

top_model_name <- rmse_df %>% 
  arrange(rmse) %>%
  slice(1) %>%
  pull(model)

top_model <- top_models_filtered[[top_model_name]] #extract top model from model list

drop_1 <- drop1(top_model, test = "Chi") %>%
  arrange(desc(Deviance))

res <- simulateResiduals(top_model)
plot(res)

data_predictions <- add_predictions_with_ci(top_model, data_count_filtered)

top_model_predict <- data_predictions %>%
  ggplot(aes(x = factor(wy))) + ##graphing the predictions
  geom_crossbar(aes(y = pred, ymin = lower_ci, ymax = upper_ci), fill = 'lightgrey',alpha = 0.5) +
  geom_point(aes(y = count), size = 4) +
  labs(x = 'Water Year', y = 'Predicted Salvage') +
  #facet_wrap(~facility, nrow = 2, scales = 'free_y') +
  theme_bw() +
  theme(plot.margin = margin(0.5,0.5,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 13))
top_model_predict

top_model_predict2 <- data_predictions %>%
  ggplot(aes(x = count, y = pred)) + ##graphing the predictions
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(x = 'Observed Salvage', y = 'Predicted Loss') +
  #facet_wrap(~facility, nrow = 2, scales = 'free') +
  theme_bw() +
  theme(plot.margin = margin(0.5,0.5,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
top_model_predict2

final <- top_model_predict/top_model_predict2
final

ggsave(filename = here(project, 'viz/count_all_years.png'), final, height = 9, width = 10)


