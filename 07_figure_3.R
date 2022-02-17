### Figure 3 - cascade plot ####################################

library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(ggplot2)
library(patchwork)

n <- 100000
prev <- 0.001403431
# Load fit
parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") %>%
  filter(country %in% c("Uganda", "Tanzania", "Kenya")) %>%
  select(sample, prob_symptomatic, prob_recognise, hosp) %>%
  mutate(
    a = n * prev,
    b = a * prob_symptomatic,
    c = b * prob_recognise,
    d = c * hosp
  )

p1 <- parameters %>% 
  slice_sample(n = 100) %>%
  select(c(sample,a:d)) %>%
  pivot_longer(-sample, names_to = "step", values_to = "n") %>%
  mutate(step = factor(step,
                          levels = c("a", "b", "c", "d"),
                          labels = c("LM+ve &\nhb < 5g/dL",
                                     "Symptomatic\nSMA",
                                     "Symptoms\nrecognised",
                                     "Access\nhospital")))

fig3 <- ggplot(p1, aes(x = step, y = n, col = step)) +
  geom_jitter(size = 0.8) +
  geom_boxplot(fill = NA, col = "black", coef = 100) +
  ylab("Number of children (per 100,000)") +
  xlab("") + 
  theme_bw() + 
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, 4.5)) +
  geom_segment(x = 4.65, xend = 4.65, y = median(parameters$b), yend = median(parameters$d), col = "black") +
  geom_segment(x = 4.55, xend = 4.65, y = median(parameters$b), yend = median(parameters$b), col = "black") +
  geom_segment(x = 4.55, xend = 4.65, y = median(parameters$d), yend = median(parameters$d), col = "black") +
  geom_text(aes(label = "Care\ngap", x = 4.9 , y = (median(parameters$b) + median(parameters$d)) / 2), col = "black")

ggsave("figures/fig3.png", fig3, height = 6, width = 6)
