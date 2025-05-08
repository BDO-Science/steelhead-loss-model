library(tidyverse)
library(janitor)
library(busdater)
library(ggrepel)
library(ggpmisc)
library(pscl)

###import and clean loss data
sh_import <- read_csv('https://www.cbr.washington.edu/sacramento/data/php/rpt/juv_loss_detail.php?sc=1&outputFormat=csv&year=all&species=2%3Aall&dnaOnly=no&age=no') %>%
  clean_names()

sh <- sh_import %>%
  mutate(wy = get_fy(as.Date(sample_time), opt_fy_start = '10-01')) %>%
  filter(wy < 2025 & !is.na(wy),
         !is.na(adipose_clip))

###import and clean hatchery data
old_hatchery <- read_csv('CV_Steelhead_Hatchery_Release_Database.csv') %>% ###adding a few older years that aren't in hatchery_releases.csv
  clean_names() %>%
  group_by(water_year_wy) %>%
  summarize(number_released = sum(total_number_released)) %>%
  filter(water_year_wy %in% c(1994,1995,1996,1997,1998)) %>%
  rename('wy' = 'water_year_wy')
hatchery <- bind_rows(read_csv("hatchery_releases.csv"), old_hatchery) #binding all hatchery releases together

###import and clean water year type data
wytypes <- read_csv("https://www.cbr.washington.edu/sacramento/data/php/rpt/hci.php?sc=1&outputFormat=csv&classification=Official") %>%
  clean_names() %>%
  mutate(wy = as.numeric(wy))

sj_wy_types_prev_year <- wytypes %>%
  filter(basin == 'San Joaquin Valley') %>%
  select(wy = 1, sj_index_prev = 3) %>%
  bind_rows(data.frame(wy = c(1991,1992,1993,1994), sj_index_prev = c(1.96,1.56 ,4.20, 2.05))) %>%
  mutate(wy = wy + 1)  #keying in wy_type - 1 to line up with current years salvage

all_wy_types <- wytypes %>%
  filter(basin == 'San Joaquin Valley') %>%
  select(wy = 1, sj_index_prev2 = 3) %>%
  bind_rows(data.frame(wy = c(1991,1992,1993,1994), sj_index_prev2 = c(1.96,1.56 ,4.20, 2.05))) %>%
  mutate(wy = wy + 2) %>% #keying in wy_type - 1 to line up with current years salvage
  left_join(sj_wy_types_prev_year, by = 'wy') #joining to wy_type - 1

###summarizing loss data by WY and joining water year types
sh_year <- sh %>%
  group_by(wy, adipose_clip) %>%
  summarize(loss = sum(loss, na.rm = TRUE)) %>%
  pivot_wider(names_from = 'adipose_clip', values_from = 'loss') %>%
  left_join(hatchery, by = 'wy') %>%
  mutate(hatch_prop = Clipped/number_released) %>%
  left_join(all_wy_types, by = 'wy') %>%
  filter(wy > 2008)

###exploratory data analysis showing correlation between hatchery and natural steelhead loss 
hatch_v_natural <- ggplot(sh_year, aes(x = Unclipped, y = hatch_prop)) +
  geom_smooth(method = 'lm', se = FALSE, size = 1.5, color = 'steelblue3', linetype = 'dashed') +
  geom_point(size = 2) + 
  geom_text_repel(aes(label = wy), min.segment.length = 0.1) +
  stat_poly_eq(label.x = 0.1, label.y = 0.9, size = 5) +
  labs(x = 'Natural Steelhead Total Loss', y = 'Proportion of Hatchery Loss') +
  theme_bw() +
  theme(plot.margin = margin(0.5,0.5,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
hatch_v_natural

###linear model correlating annual steelhead loss to previous years water type
model1 <- glm(Unclipped ~ sj_index_prev * sj_index_prev2, #model both years with interaction
                        family = gaussian(link = "log"), 
                        data = sh_year)
summary(model1)

model2 <- glm(Unclipped ~ sj_index_prev + sj_index_prev2, #model both years without interaction
              family = gaussian(link = "log"), 
              data = sh_year)
summary(model2)

model3 <- glm(Unclipped ~ sj_index_prev, #model only previous years wy_type
              family = gaussian(link = "log"), 
              data = sh_year)
summary(model3)

model4 <- glm(Unclipped ~ sj_index_prev2, #model only wy_tpe from two years ago
              family = gaussian(link = "log"), 
              data = sh_year)
summary(model4)

aics <- data.frame(model = c('model1', 'model2', 'model3', 'model4'), #comparing AICs and pseudoR2
                   covariate = c('both years interact', 'both years noninteract', 'wytype - 1', 'wytype - 2'),
                   aic = c(model1[["aic"]], model2[["aic"]], model3[["aic"]], model4[["aic"]]),
                   pseudo_r2 = round(c(pR2(model1)["McFadden"],pR2(model2)["McFadden"],
                                 pR2(model3)["McFadden"],pR2(model4)["McFadden"]),2)) %>%
  arrange(aic)

###making predictions on existing data using the best model
preds <- predict(model1, newdata = sh_year, type = "response", se.fit = TRUE)

sh_year$pred <- preds$fit
sh_year$se <- preds$se.fit
sh_year$lower_ci <- preds$fit - 1.96 * preds$se.fit
sh_year$upper_ci <- preds$fit + 1.96 * preds$se.fit

predict_graph <- ggplot(sh_year, aes(x = wy)) + ##graphing the predictions
  geom_point(aes(y = Unclipped), size = 4) +
  geom_crossbar(aes(y = pred, ymin = lower_ci, ymax = upper_ci), fill = 'steelblue3', alpha = 0.5) +
  labs(x = 'Water Year', y = 'Predicted Loss') +
  theme_bw() +
  theme(plot.margin = margin(0.5,0.5,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
predict_graph

###graphing overall model
final_fit <- expand.grid(sj_index_prev = seq(0.5,7,0.1),
                         sj_index_prev2 = seq(1,7,1))
preds_final <- predict(model1, newdata = final_fit, type = "response", se.fit = TRUE)

final_fit$pred <- preds_final$fit
final_fit$se <- preds_final$se.fit
final_fit$lower_ci <- preds_final$fit - 1.96 * preds_final$se.fit
final_fit$upper_ci <- preds_final$fit + 1.96 * preds_final$se.fit

final_fit_graph <- ggplot(final_fit, aes(sj_index_prev, pred, 
                                         color = factor(sj_index_prev2),
                                         fill = factor(sj_index_prev2))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.5) +
  facet_wrap(~sj_index_prev2) +
  labs(x = 'San Joaquin Previous Water Year Index', y = 'Predicted Loss') +
  theme_bw() +
  theme(legend.position = 'none',
        plot.margin = margin(0.5,0.5,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
final_fit_graph  
