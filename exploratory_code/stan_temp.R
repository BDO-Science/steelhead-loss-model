library(tidyverse)
library(CDECRetrieve)
library(janitor)
library(busdater)

years <- factor(seq(2006,2024,1))
templist <- list()
for(i in years){
  url <- paste0('https://www.cbr.washington.edu/sacramento/data/php/rpt/mg.php?map=1&mgconfig=river&tempUnit=F&avgyear=0&consolidate=1&grid=1&y1min=&y1max=&y2min=&y2max=&size=large&outputFormat=csvSingle&data[]=WaterTemperature&loc[]=OBB&year[]=',i)
  temp <- read_csv(url) %>% clean_names()
  templist[[i]] <- temp
}

all <- bind_rows(templist) %>%
  mutate(Date = as.Date(paste0(year,'-',mm_dd))) %>%
  mutate(wy = get_fy(Date, opt_fy_start = '10-01')) %>%
  mutate(exceed = if_else(value > 63, 1, 0)) %>%
  filter(!is.na(Date))

temp_summary <- all %>%
  group_by(wy) %>%
  summarize(days = sum(exceed, na.rm = TRUE),
            temp = mean(value, na.rm = TRUE)) %>%
  mutate(temp_1year = lag(temp, n = 1),
         temp_2year = lag(temp, n = 2)) %>%
  filter(wy > 2008)

