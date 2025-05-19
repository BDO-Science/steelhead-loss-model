

sj_wy_types_prev_year <- wytypes %>%
  filter(basin == 'San Joaquin Valley') %>%
  dplyr::select(wy = 1, sj_index_prev = 3) %>%
  bind_rows(data.frame(wy = c(1991,1992,1993,1994), sj_index_prev = c(1.96,1.56 ,4.20, 2.05))) %>%
  mutate(wy = wy + 1)  #keying in wy_type - 1 to line up with current years salvage

all_sj_wy_types <- wytypes %>%
  filter(basin == 'San Joaquin Valley') %>%
  dplyr::select(wy = 1, sj_index_prev2 = 3) %>%
  bind_rows(data.frame(wy = c(1991,1992,1993,1994), sj_index_prev2 = c(1.96,1.56 ,4.20, 2.05))) %>%
  mutate(wy = wy + 2) %>% #keying in wy_type - 1 to line up with current years salvage
  left_join(sj_wy_types_prev_year, by = 'wy') #joining to wy_type - 1

sac_wy_types_prev_year <- wytypes %>%
  filter(basin == 'Sacramento Valley') %>%
  dplyr::select(wy = 1, sac_index_prev = 3) %>%
  bind_rows(data.frame(wy = c(1991,1992,1993,1994), sac_index_prev = c(1.96,1.56 ,4.20, 2.05))) %>%
  mutate(wy = wy + 1)  #keying in wy_type - 1 to line up with current years salvage

all_sac_wy_types <- wytypes %>%
  filter(basin == 'Sacramento Valley') %>%
  dplyr::select(wy = 1, sac_index_prev2 = 3) %>%
  bind_rows(data.frame(wy = c(1991,1992,1993,1994), sac_index_prev2 = c(1.96,1.56 ,4.20, 2.05))) %>%
  mutate(wy = wy + 2) %>% #keying in wy_type - 1 to line up with current years salvage
  left_join(sac_wy_types_prev_year, by = 'wy') #joining to wy_type - 1

all_wy_type <- all_sj_wy_types %>%
  left_join(all_sac_wy_types, by = 'wy')

sh_year_all <- sh %>%
  group_by(wy, facility, adipose_clip) %>%
  summarize(loss = sum(loss, na.rm = TRUE)) %>%
  pivot_wider(names_from = 'adipose_clip', values_from = 'loss') %>%
  left_join(hatchery, by = 'wy') %>%
  mutate(hatch_prop = Clipped/number_released) %>%
  left_join(all_wy_type, by = 'wy') 

sh_year <- sh_year_all %>%
  filter(wy > 2008) %>%
  mutate(Unclipped = as.integer(Unclipped)) %>%
  mutate(wyindex = sj_index_prev/sac_index_prev,
         wyindex2 = sj_index_prev2/sac_index_prev2)

model.nm <- glm.nb(Unclipped ~ sj_index_prev * sj_index_prev2 + facility,
                   data = sh_year)
summary(model.nm)

model.nm1 <- glm.nb(Unclipped ~ wyindex * wyindex2 + facility,
                    data = sh_year)
summary(model.nm1)

preds <- predict(model.nm, newdata = sh_year, type = "response", se.fit = TRUE)

sh_year$pred <- preds$fit
sh_year$se <- preds$se.fit
sh_year$lower_ci <- preds$fit - 1.96 * preds$se.fit
sh_year$upper_ci <- preds$fit + 1.96 * preds$se.fit

predict_graph <- ggplot(sh_year, aes(x = wy, group = facility)) + ##graphing the predictions
  geom_point(aes(y = Unclipped), size = 4) +
  geom_crossbar(aes(y = pred, ymin = lower_ci, ymax = upper_ci, fill = facility), alpha = 0.5) +
  labs(x = 'Water Year', y = 'Predicted Loss') +
  facet_wrap(~facility, nrow = 2, scales = 'free_y') +
  theme_bw() +
  theme(plot.margin = margin(0.5,0.5,0.2,0.2, unit = 'cm'),
        axis.title.y = element_text(margin=margin(r=15), size = 15),
        axis.title.x = element_text(margin=margin(t=15), size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
predict_graph
