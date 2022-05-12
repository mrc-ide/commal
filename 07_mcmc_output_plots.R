### MCMC parameter plots  ######################################################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(drjacoby)

# Load fit
parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") 
mcmc <- readRDS("ignore/prob_hosp/mcmc_fits/mcmc.rds")

### Random effects #############################################################
re_plot <- ggplot(parameters, aes(x = country, fill = country, y = country_capacity)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(coef = 100) +
  xlab("") + 
  ylab("") +
  scale_fill_discrete(guide = "none") +
  coord_flip() +
  theme_bw()

ggsave("ignore/figures_tables/figS_random_effects.png", re_plot, height = 4, width = 6)
################################################################################

### Parameter plots ############################################################
p <- c(names(select(parameters, -sample, -hosp, -country, -country_capacity, -chronic)),
       "hosp_Kenya", "hosp_Tanzania", "hosp_Uganda",
       "chronic_Kenya", "chronic_Tanzania", "chronic_Uganda")

# Global parameters
pp <- plot_par(mcmc, p, display = FALSE)
trace_plots <- patchwork::wrap_plots(lapply(pp, function(x) x$trace + theme(legend.position = "none")), ncol = length(p))
hist_plots <- patchwork::wrap_plots(lapply(pp, function(x) x$hist), ncol = length(p))
acf_plots <- patchwork::wrap_plots(lapply(pp, function(x) x$acf), ncol = length(p))
parameter_plots <- trace_plots / hist_plots / acf_plots

ggsave("ignore/figures_tables/figS_parameter_plots.png", parameter_plots, height = 10, width = 24)

# Correlations between global pars
cor <- apply(combn(p, 2), 2, function(x){
  plot_cor(mcmc, x[1], x[2])
})
correlation_plots <- patchwork::wrap_plots(cor) + plot_layout(guides = "collect")

ggsave("ignore/figures_tables/figS_correlation_plots.png", correlation_plots, height = 10, width = 24)
################################################################################


