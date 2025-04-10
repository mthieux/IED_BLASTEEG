# ---
# title: "eeg_stat_analysis"
# author: "Marine Thieux"
# ---

# ----------------------------------------------------------------------------------------------------
# This R script performs generalized additive modeling (GAM) on EEG data to analyze the effect of 
# interictal epileptiform discharges (IED) over time.
# ----------------------------------------------------------------------------------------------------


# Libraries
# ----------------------------------------------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(ggeffects)
library(emmeans)
library(mgcv)
library(itsadug)
library(ggplot2)
library(purrr)

# ----------------------------------------------------------------------------------------------------

# Import data
# ----------------------------------------------------------------------------------------------------
# EEG Long-Format Data Specification
# Expected structure of the processed EEG dataset: 
# - Data must be in **long format**: one row per `time × channel × epoch`.
# - Should contain **only EEG channels**.
# Data Frame Structure : 
  # 
  # | Column     | Type         | Description |
  # |------------|--------------|-------------|
  # | `time`     | `numeric`    | Time within the epoch 
  # | `epoch`    | `factor`     | Identifier for each epoch
  # | `channel`  | `factor`     | EEG channel name 
  # | `ampl`     | `numeric`    | Amplitude 
  # | `session`  | `factor`     | BLAST Session number
  # | `IED`      | `factor`     | interictal epileptiform discharges: `"with"` or `"without"`.
  # | `patname`  | `factor`     | Patient ID
  # 
# ----------------------------------------------------------------------------------------------------

# Grand average
# ----------------------------------------------------------------------------------------------------
datas %>%
  group_by(IED, channel, time) %>%
  summarise(ampl_avg = mean(ampl)) %>%
  ggplot(aes(x = time, y = ampl_avg, color = channel)) +
  geom_hline(yintercept = 0) +
  geom_line(size = 0.3) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = -0.7, linetype = "dotdash") +
  theme(legend.position = "none") +
  facet_grid(row = vars(IED))
# ----------------------------------------------------------------------------------------------------


# GAM - Scale T distribution
# ----------------------------------------------------------------------------------------------------
datas %>%
  filter(time > -0.2 & time < 0.8) %>%
  subset(channel == c("Pz", "P3", "P4")) %>%
  droplevels() -> datas

model.gamm.randpat.autocorr.scalet <- mgcv::bam(
  ampl ~ s(
    time,
    by = IED,
    k = 20,
    bs = "cr",
    m = 1
  ) +
    IED +
    s(
      time,
      patname,
      k = 20,
      bs = "fs",
      m = 1
    ),
  data = datas,
  rho = 0.98,
  AR.start = datas$start.event,
  family = "scat",
  nthreads =  parallel::detectCores(),
  discrete = TRUE,
  method = "fREML",
  select = TRUE
)

appraise(model.gamm.randpat.autocorr.scalet)
# ----------------------------------------------------------------------------------------------------

# GAM to compare
# ----------------------------------------------------------------------------------------------------
model.gamm.randpat.autocorr.scalet.01 <- mgcv::bam(
  ampl ~ s(time, k = 20, bs = "cr", m = 1) +
    IED +
    s(time, patname, k = 20, bs = "fs", m = 1),
  data = datas,
  rho = 0.98,  
  AR.start = datas$start.event,
  family = "scat",
  nthreads =  parallel::detectCores(),
  discrete = TRUE,  
  method = "fREML", 
  select = TRUE    
)


model.gamm.randpat.autocorr.scalet.02 <- mgcv::bam(
  ampl ~ s(time,  k = 20, bs = "cr", m = 1) +
    s(time, patname, k = 20, bs = "fs", m = 1),
  data = datas,
  rho = 0.98,  
  AR.start = datas$start.event,
  family = "scat",
  nthreads =  parallel::detectCores(),
  discrete = TRUE,  
  method = "fREML", 
  select = TRUE    
)


model.gamm.randpat.autocorr.scalet.03 <- mgcv::bam(
  ampl ~ s(time, by = IED, k = 20, bs = "cr", m = 1) +
    s(time, patname, k = 20, bs = "fs", m = 1),
  data = datas,
  rho = 0.98,  
  AR.start = datas$start.event,
  family = "scat",
  nthreads =  parallel::detectCores(),
  discrete = TRUE, 
  method = "fREML", 
  select = TRUE    
)
# ----------------------------------------------------------------------------------------------------

# Comparisons
# ----------------------------------------------------------------------------------------------------
compareML(model.gamm.randpat.autocorr.scalet.01,
          model.gamm.randpat.autocorr.scalet.02,
          suggest.report = T,
          print.output = T,
          signif.stars = T)$table

compareML(model.gamm.randpat.autocorr.scalet.02,
          model.gamm.randpat.autocorr.scalet.03,
          suggest.report = T,
          print.output = T,
          signif.stars = T)$table

compareML(model.gamm.randpat.autocorr.scalet.01,
          model.gamm.randpat.autocorr.scalet.03,
          suggest.report = T,
          print.output = T,
          signif.stars = T)$table

compareML(model.gamm.randpat.autocorr.scalet.01,
          model.gamm.randpat.autocorr.scalet,
          suggest.report = T,
          print.output = T,
          signif.stars = T)$table


compareML(model.gamm.randpat.autocorr.scalet.02,
          model.gamm.randpat.autocorr.scalet,
          suggest.report = T,
          print.output = T,
          signif.stars = T)$table

compareML(model.gamm.randpat.autocorr.scalet.03,
          model.gamm.randpat.autocorr.scalet,
          suggest.report = T,
          print.output = T,
          signif.stars = T)$table
# ----------------------------------------------------------------------------------------------------

# Best model parameters 
# ----------------------------------------------------------------------------------------------------
summary(model.gamm.randpat.autocorr.scalet)

par(mfrow=c(1,2), cex=1.1)

plot_smooth(model.gamm.randpat.autocorr.scalet, 
            view = "time", rug = F, plot_all = "IED",
            rm.ranef = F,
            n.grid = 600,
            main = "with random")

plot_diff(model.gamm.randpat.autocorr.scalet, 
          view = "time", 
          comp = list(IED = c("with", "without")),
          rm.ranef = F,
          n.grid = 600,
          sim.ci = F)
# ----------------------------------------------------------------------------------------------------

# Post-Hoc
# ----------------------------------------------------------------------------------------------------
emm <- emmeans(model.gamm.randpat.autocorr.scalet,
                  ~ IED | time, 
                  at = list(time = seq(from = 0.2, to = 0.8, by = 0.05)))

pairs(emm) %>%
  as.data.frame() %>%
  knitr::kable()

des_cont <- c(list("without-with" = c(-1, 1)))
test(contrast(regrid(emm), des_cont, adjust = "FDR"))
# ----------------------------------------------------------------------------------------------------

# Bootstrap with permutation
# ----------------------------------------------------------------------------------------------------
fit_gam <- function(data) {
  mgcv::bam(
    ampl ~ s(time, by = IED, k = 20, bs = "cr", m = 1) +
      IED +
      s(time, patname, k = 20, bs = "fs", m = 1),
    data = data,
    rho = 0.98,
    AR.start = data$start.event,
    family = "scat",
    nthreads = parallel::detectCores(),
    discrete = TRUE,
    method = "fREML",
    select = TRUE
  )
}

time_grid <- seq(min(datas$time), max(datas$time), length.out = 600) # Prediction grid
pred_data <- expand.grid(
  time = time_grid,
  IED = c("with", "without"),
  patname = unique(datas$patname)[1]  # arbitrary subject
)

observed_model <- fit_gam(datas) 
observed_preds <- predict(observed_model, newdata = pred_data)
observed_df <- cbind(pred_data, fit = observed_preds) %>%
  pivot_wider(names_from = IED, values_from = fit) %>%
  mutate(diff = with - without)

n_iter <- 100 # Run permutations
perm_diffs <- replicate(n_iter, {
  perm_data <- datas
  perm_data$IED <- sample(datas$IED)  # shuffle labels
  perm_model <- fit_gam(perm_data)
  perm_preds <- predict(perm_model, newdata = pred_data)
  perm_df <- cbind(pred_data, fit = perm_preds) %>%
    pivot_wider(names_from = IED, values_from = fit) %>%
    mutate(diff = with - without)
  perm_df$diff
}, simplify = "matrix")

perm_CI <- apply(perm_diffs, 1, quantile, probs = c(0.025, 0.975))  # null distribution


observed_df <- observed_df %>% 
  mutate(
    lower_CI = perm_CI[1, ],
    upper_CI = perm_CI[2, ],
    sig = ifelse(diff < lower_CI | diff > upper_CI, TRUE, FALSE)
  )
# ----------------------------------------------------------------------------------------------------

# Figure
# ----------------------------------------------------------------------------------------------------
plotme <- ggemmeans(model.gamm.randpat.autocorr.scalet,
                    terms = c("time","IED"))

plot(plotme) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("with" = "#993356", "without" = "#8da0cb")) +
  scale_fill_manual(values = c("with" = "#993356", "without" = "#8da0cb")) +
  labs(
    title = "",
    x = "Time (s)",       
    y = "Predicted amplitude (uV)" 
  ) +
  theme_minimal(base_size = 16) + 
  scale_x_continuous(breaks = seq(-0.2, 0.8, by = 0.2)) +  
  scale_y_continuous(breaks = seq(-5, 5, by = 1.5)) + 
  theme(legend.position = c(0.9, 0.85),
        axis.text = element_text(size = 14),  
        axis.title = element_text(size = 16),  
        plot.title = element_text(size = 18)) -> FigA


# Get model-based difference with CI from plot_diff()
model_ci_df <- plot_diff(
  model.gamm.randpat.autocorr.scalet,
  view = "time",
  comp = list(IED = c("with", "without")),
  rm.ranef = FALSE,
  n.grid = 600,
  return.data = TRUE  
)

model_ci_df <- model_ci_df %>%
  mutate(
    model_fit = est,
    model_CI_low = est - CI,
    model_CI_up  = est + CI, 
    model_sig = model_CI_low > 0 | model_CI_up < 0 
  )

# Join model-based CI with permutation-based plot data
combined_df <- observed_df %>%
  left_join(model_ci_df, by = "time")


ggplot(combined_df, aes(x = time)) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "#8da0cb", alpha = 0.3) + 
  geom_ribbon(aes(ymin = model_CI_low, ymax = model_CI_up), fill = "grey", alpha = 0.4) +  
  geom_line(aes(y = diff), color = "grey30", size = 3) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, color = "black") + 
  geom_point(data = subset(combined_df, sig),
             aes(x = time, y = diff),
             color = "tomato3", size = 1, shape = 16, alpha = 0.8) +  
  geom_tile(data = filter(model_ci_df, model_sig),
            aes(y = -2.5), height = 0.02, width = 0.01,
            fill = "tomato2", alpha = 0.8) + 
  labs(x = "Time (s)", y = "Amplitude Difference (µV)") +
  theme_minimal(base_size = 16) + 
  scale_x_continuous(breaks = seq(-0.2, 0.8, by = 0.2)) +  
  scale_y_continuous(breaks = seq(0, -2.5, by = -0.5)) + 
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 16),  
        plot.title = element_text(size = 18))  -> FigB


(FigA / FigB) + 
  plot_annotation(
    tag_levels = 'A',             
    tag_suffix = ')',
    theme = theme(plot.tag = element_text(size = 16, face = "bold", hjust = -0.1)) 
  ) -> Figure

# ----------------------------------------------------------------------------------------------------
