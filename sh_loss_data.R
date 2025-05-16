setwd("~/GitHub/steelhead-loss-model")

source("functions.R")

library(tidyverse)
library(janitor)
library(busdater)
library(ggrepel)
library(ggpmisc)
library(pscl)
library(AICcmodavg)
library(dataRetrieval)
library(rvest)

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

sj_wy_types_prev_year_3 <- wytypes %>%
  filter(basin == 'San Joaquin Valley') %>%
  select(wy = 1, sj_index_prev3 = 3) %>%
  bind_rows(data.frame(wy = c(1991,1992,1993,1994), sj_index_prev3 = c(1.96,1.56 ,4.20, 2.05))) %>%
  mutate(wy = wy + 3)  #keying in wy_type - 1 to line up with current years salvage

all_wy_types <- wytypes %>%
  filter(basin == 'San Joaquin Valley') %>%
  select(wy = 1, sj_index_prev2 = 3) %>%
  bind_rows(data.frame(wy = c(1991,1992,1993,1994), sj_index_prev2 = c(1.96,1.56 ,4.20, 2.05))) %>%
  mutate(wy = wy + 2) %>% #keying in wy_type - 1 to line up with current years salvage
  left_join(sj_wy_types_prev_year, by = 'wy') %>% #joining to wy_type - 1
  left_join(sj_wy_types_prev_year_3, by = 'wy') #joining to wy_type - 1

###summarizing loss data by WY and joining water year types
sh_year <- sh %>%
  group_by(wy, adipose_clip) %>%
  summarize(loss = sum(loss, na.rm = TRUE)) %>%
  pivot_wider(names_from = 'adipose_clip', values_from = 'loss') %>%
  left_join(hatchery, by = 'wy') %>%
  mutate(hatch_prop = Clipped/number_released) %>%
  left_join(all_wy_types, by = 'wy')

###############
#Pull OMR Data

data_OMR <- calc_OMR(dateStart = "1993-01-01",dateEnd = "2024-08-01", timing = "daily", extrap = T, proof = T)
#Regression of daily middle vs old river flow - adj R = 0.96
#data from 1993-01-01 to 2024-08-01

data_OMR_sum <- data_OMR %>% filter(month %in% c(12,1:6)) %>% rename(wy=waterYear) %>% group_by(wy) %>%
  summarise(OMR_mean = mean(omrFlowExtrap,na.rm=T))

###############
#Pull Vernalis, Freeport, and Export Data


data_df<-get_dayflow()

data_flow <- data_df %>% mutate(month = month(date),wy=ifelse(month(date)>9,year(date)+1,year(date))) %>%
  group_by(wy) %>% filter(month %in% c(12,1:6)) %>% summarise(sjr=mean(sjr),sac=mean(sac),exports=mean(exports))

data_vernalis <- readNWISuv(siteNumbers = "11303500",
                            parameterCd = "00060",
                            startDate = "2023-01-01",
                            endDate = "2024-08-01")
data_vernalis_sum <- data_vernalis %>% mutate(date=as.Date(dateTime)) %>% mutate(month = month(date),wy=ifelse(month(date)>9,year(date)+1,year(date))) %>%
  group_by(wy) %>% rename(sjr=X_00060_00000 )%>% filter(month %in% c(12,1:6)) %>% summarise(sjr=mean(sjr,na.rm=T)) %>%
  filter(wy==2024)

data_freeport <- readNWISuv(siteNumbers = "11447650",
                            parameterCd = "00060",
                            startDate = "2023-01-01",
                            endDate = "2024-08-01")
data_freeport_sum <- data_freeport %>% mutate(date=as.Date(dateTime)) %>% mutate(month = month(date),wy=ifelse(month(date)>9,year(date)+1,year(date))) %>%
  group_by(wy) %>% rename(sac=X_00060_00000 )%>% filter(month %in% c(12,1:6)) %>% summarise(sac=mean(sac,na.rm=T)) %>%
  filter(wy==2024)

data_comb_sum<- left_join(data_vernalis_sum,data_freeport_sum)

data_flow <- bind_rows(data_flow,data_comb_sum)

### 
#Rejoin data
sh_year <- sh_year %>% left_join(data_OMR_sum) %>% left_join(data_flow)


