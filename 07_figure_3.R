### Figure 3 - cascade plot ####################################

library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(ggplot2)
library(patchwork)

n <- 100000
prev <- 0.0006291983
# Load fit
parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") %>%
  filter(country %in% c("Uganda", "Tanzania", "Kenya")) %>%
  select(sample, prob_recognise, hosp) %>%
  mutate(
    a = n * prev,
    b = a * prob_recognise,
    c = b * hosp,
    a = a+ c
  )

p1 <- parameters %>% 
  slice_sample(n = 100) %>%
  select(c(sample,a:c)) %>%
  pivot_longer(-sample, names_to = "step", values_to = "n") %>%
  mutate(step = factor(step,
                          levels = c("a", "b", "c"),
                          labels = c("Symptomatic\nSMA",
                                     "Symptoms\nrecognised",
                                     "Access\nhospital")))

fig3 <- ggplot(p1, aes(x = step, y = n, col = step)) +
  geom_jitter(size = 0.8) +
  geom_boxplot(fill = NA, col = "black", coef = 100) +
  ylab("Number of children (per 100,000)") +
  xlab("") + 
  theme_bw() + 
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, 3.5)) +
  geom_segment(x = 3.65, xend = 3.65, y = median(parameters$a), yend = median(parameters$c), col = "black") +
  geom_segment(x = 3.55, xend = 3.65, y = median(parameters$a), yend = median(parameters$a), col = "black") +
  geom_segment(x = 3.55, xend = 3.65, y = median(parameters$c), yend = median(parameters$c), col = "black") +
  geom_text(aes(label = "Care\ngap", x = 3.9 , y = (median(parameters$b) + median(parameters$c)) / 2), col = "black")

ggsave("figures/fig3.png", fig3, height = 6, width = 6)
