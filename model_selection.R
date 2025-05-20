library(tidyverse)
library(MASS)
library(MuMIn)
library(DHARMa)
library(gridExtra)
source('functions.R')

data_all <- readRDS('model_data.rds') %>%
  mutate(era = if_else(wy <= 2008, 'Early', 'Current'))

data_filtered <- data_all %>%
  dplyr::filter(wy > 2008)

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
  #excluded: sjr
  glm.nb(loss ~ facility + OMR_mean + exports + 
           sac + sj_index_prev + sj_index_prev2 +
           sac_index_prev + sac_index_prev2, data = data_filtered,
         na.action = na.fail),
  #excluded: sjr, sac_index_prev2
  glm.nb(loss ~ facility + OMR_mean + exports + 
           sac + sj_index_prev + sj_index_prev2 +
           sac_index_prev, data = data_filtered,
         na.action = na.fail)
)

names(model_list) <- paste0("model", seq_along(model_list))
vifs <- lapply(model_list, car::vif)


######################################################################
#using dredge on final model after removing collinearity
#included interaction between sj_index_prev and sj_index_prev2
######################################################################

full_model <- glm.nb(loss ~ facility + OMR_mean + exports + 
                       sac + sj_index_prev + sj_index_prev2 +
                       sac_index_prev + sj_index_prev:sj_index_prev2, data = data_filtered,
                     na.action = na.fail)

model_set <- dredge(full_model) #modeling all combinations of covariates
top_models <- get.models(model_set, subset = 1:5) #pulling top 5 models based on AICc


#cross-validation of all top models
cv_list <- lapply(top_models, function(m) {
  boot::cv.glm(data_filtered, m, K = 5)
})

#simulating and plotting residuals for best fitting model 4 and model 8
res_list <- lapply(top_models, simulateResiduals)
lapply(res_list, plot)

#predicting on all of the top models
pred_list <- lapply(top_models, function(m) {
  add_predictions_with_ci(m, data_filtered)
})


graph_list <- lapply(pred_list, function(df){
  ggplot(df, aes(x = wy, group = facility)) + ##graphing the predictions
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
})

all_predictions <- grid.arrange(grobs = graph_list, nrow = 2)
ggsave(all_predictions, file = 'model_predictions.png', height = 8, width = 15)

