### Figure 1 - masa fit wrt PfPr ################################################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

source("R/model_functions.R")

ndraw <- 100

dhs_masa <- readRDS("ignore/prob_hosp/dhs_masa.rds") %>%
  mutate(country = ifelse(country == "Congo Democratic Republic", "DRC", country))

# Load fit
parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") %>%
  select(sample, country, global_capacity, country_capacity, shift, pfpr_beta) %>%
  mutate(country = ifelse(country == "Congo Democratic Republic", "DRC", country))

# Global fit
fit_median <- parameters %>%
  select(-country, -country_capacity) %>%
  summarise_all(median) %>%
  mutate(country_capacity = 0) %>%
  left_join(data.frame(
    pfpr = seq(0, 0.7, 0.01)),
    by = character()) %>%
  mutate(masapr = pmap_dbl(select(., -sample), gompertz))

fit_draws <- parameters %>%
  select(-country, -country_capacity) %>%
  slice_sample(n = ndraw) %>%
  mutate(country_capacity = 0) %>%
  left_join(data.frame(
    pfpr = seq(0, 0.7, 0.01)),
    by = character()) %>%
  mutate(masapr = pmap_dbl(select(., -sample), gompertz))

# Country fits
country_data <- dhs_masa %>%
  group_by(country) %>%
  mutate(pfprg = cut_number(pfpr, 4)) %>%
  group_by(country, pfprg) %>%
  summarise(
    pfprl = quantile(pfpr, 0.025),
    pfpru = quantile(pfpr, 0.975),
    pfpr = mean(pfpr),
    masal = binom::binom.exact(sum(masa), n())$lower,
    masau = binom::binom.exact(sum(masa), n())$upper,
    masa = mean(masa)
  )

country_pfpr_max_bound <- country_data %>%
  select(country, pfpru) %>%
  group_by(country) %>%
  slice_max(order_by = pfpru, n = 1, with_ties = FALSE) %>%
  mutate(pfprq = pmax(0.2, pfpru)) %>%
  select(-pfpru)

country_draws <- parameters %>%
  group_by(country) %>%
  slice_sample(n = ndraw) %>%
  ungroup() %>%
  left_join(data.frame(
    pfpr = seq(0, 0.9, 0.01)),
    by = character()
  ) %>%
  left_join(country_pfpr_max_bound, by = "country") %>%
  filter(pfpr < pfprq) %>%
  select(-pfprq) %>%
  mutate(masapr = pmap_dbl(select(., -country, -sample), gompertz))

country_median <- parameters %>%
  group_by(country) %>%
  summarise_all(median) %>%
  ungroup() %>%
  left_join(data.frame(
    pfpr = seq(0, 0.9, 0.01)),
    by = character()
  ) %>%
  left_join(country_pfpr_max_bound, by = "country") %>%
  filter(pfpr < pfprq) %>%
  select(-pfprq) %>%
  mutate(masapr = pmap_dbl(select(., -country, -sample), gompertz))

# To order country plots by mean masa prevalence
country_order <- dhs_masa %>%
  group_by(country) %>%
  summarise(masa = mean(masa)) %>%
  ungroup() %>%
  arrange(-masa)

country_draws$country <- factor(country_draws$country, levels = country_order$country)
country_median$country <- factor(country_median$country, levels = country_order$country)

### Figure 1a model fitted trend ####
fig1a <- ggplot() +
  geom_line(data = fit_draws, aes(x = pfpr, y = masapr, group = sample), alpha = 0.2, col = "#00798c") +
  geom_line(data = fit_median, aes(x = pfpr, y = masapr), col = "#edae49", size = 1) +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab(expression(MSA[0.5-5])) + 
  theme_bw() +
  coord_cartesian(ylim = c(0, 0.005))

### Figure 1b, fit to country data ###
breaks <- function(){
  function(x){
    up <- max(x)
    if(up <= 0.2){
      b <- seq(0, 0.2, 0.1)
    } 
    if(up > 0.2 & up <= 0.4){
      b <- seq(0, up, 0.1)
    }
    if(up > 0.4){
      b <- seq(0, max(x), 0.2)
    }
    return(b)
  }
}

country_plot <- function(country_median, country_draws, country_data, nc = 3){
  ggplot() +
    geom_line(data = country_draws, aes(x = pfpr, y = masapr, group = sample), alpha = 0.1, col = "#00798c") +
    geom_line(data = country_median, aes(x = pfpr, y = masapr), col = "#edae49", size = 1) +
    geom_linerange(data = country_data, aes(y = masa, xmin = pfprl, xmax = pfpru), col = "#2e4057") +
    geom_linerange(data = country_data, aes(x = pfpr, ymin = masal, ymax = masau), col = "#2e4057") +
    geom_point(data = country_data, aes(x = pfpr, y = masa), col = "#2e4057", size = 0.75) +
    xlab(expression(~italic(Pf)~Pr[2-10])) +
    ylab(expression(MSA[0.5-5])) + 
    scale_x_continuous(breaks = breaks()) +
    theme_bw() +
    facet_wrap(~ country, ncol = nc, scales = "free") +
    theme(strip.background = element_rect(fill = NA),
          strip.text = element_text(size = 6),
          axis.text = element_text(size = 6))
}

countries <- unique(country_median$country)
c1 <- countries[1:6]
c2 <- countries[7:21]

cp1 <- country_plot(filter(country_median, country %in% c1),
                    filter(country_draws, country %in% c1),
                    filter(country_data, country %in% c1))
cp2 <- country_plot(filter(country_median, country %in% c2),
                    filter(country_draws, country %in% c2),
                    filter(country_data, country %in% c2),
                    nc = 5)


fig1b <- (fig1a + theme(axis.title.x = element_blank()) | cp1 + theme(axis.title.x = element_blank())) +
  plot_layout(widths = c(2, 3.63))
fig1 <- (fig1b / cp2) + plot_layout(heights = c(2, 3.5))

ggsave("ignore/figures_tables/dhs_fit.png", fig1, width = 9, height = 7, scale = 0.75)
