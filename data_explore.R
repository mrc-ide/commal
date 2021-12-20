dhs_data_raw <- readRDS("ignore/dhs/processed_data/processed_dhs.RDS")
library(UpSetR)

cty <- dhs_data_raw %>%
  filter(iso %in% c("KEN", "UGA", "TZA"),
         rdt %in% c("negative", "positive"),
         microscopy %in% c("negative", "positive"),
         anemia_level %in% c("non_severe", "severe"),
         !is.na(prevalence)) %>%
  mutate(pfpr_g = cut_interval(prevalence, n = 10))

table(cty$anemia_level, cty$microscopy)


cty2 <- cty %>%
  filter(microscopy == "positive")

table(cty2$anemia_level, cty2$Impaired_breathing)
table(cty$anemia_level, cty$Impaired_consciousness)
table(cty$Impaired_breathing, cty$Impaired_consciousness)


li <- list(
  sa = which(cty2$anemia_level == "severe"),
  Impaired_breathing = which(cty2$Impaired_breathing == "yes"),
  Impaired_consciousness = which(cty2$Impaired_consciousness == "yes")
)

upset(fromList(li), order.by = "freq")



pd <- readRDS("ignore/prob_hosp/paton_inferred.rds")

ggplot(pd, aes(x = pfpr, y = 1000 * (sma / py), col = distance)) +
  geom_point()


dis <- dhs_data_raw %>%
  group_by(distance) %>%
  summarise(
    n = n(),
    sma = sum(microscopy == "positive" & anemia_level == "severe", na.rm = TRUE) / n) %>%
  ungroup() %>%
  arrange(distance) %>%
  mutate(cs = cumsum(sma) / sum(sma))

ggplot(dis, aes(x = distance, y = cs)) +
  geom_line() +
  xlim(0, 25)
  