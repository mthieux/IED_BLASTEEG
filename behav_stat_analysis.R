# ---
# title: "behav_stat_analysis"
# author: "Marine Thieux"
# ---

# ----------------------------------------------------------------------------------------------------
# This R script performs statistical analyses on behavioral data to analyze the effect of 
# interictal epileptiform discharges (IED) on attention performance.
# ----------------------------------------------------------------------------------------------------


# Libraries
# ----------------------------------------------------------------------------------------------------
library(knitr)
library(rmdformats)
library(tidyverse)
library(afex)
library(emmeans)
library(jtools)
library(ggstance)
library(ggsignif)
library(ggfortify)
library(gtsummary)
library(psych)
library(lme4)
library(glmmTMB)
library(pscl)
library(performance)
library(corrplot)
library(Hmisc)
library(RColorBrewer)
library(gridExtra)
library(showtext)
library(dagitty)
library(Rmisc) 
library(PupillometryR) 
library(see) 
# ----------------------------------------------------------------------------------------------------

# Import data
# ----------------------------------------------------------------------------------------------------
# Expected structure of the processed dataset: 
# - Data must be in **long format**: one row per `trial × session × patname`.
# Data Frame Structure : 
# 
# | Column      | Type         | Description |
# |------------ |--------------|-------------|
# | `patname`   | `factor`     | Patient ID
# | `group`     | `factor`     | Diagnostic group: "Epilepsy" or "Epilepsy + ADHD" 
# | `age`       | `numeric`    | Age at testing in months
# | `sex`       | `factor`     | Sex: "male" or "female"
# | `epi_type`  | `factor`     | Epilepsy type: "Generalized" or "Focal" 
# | `epi_synd`  | `factor`     | Epilepsy syndrome
# | `epi_onset` | `numeric`    | Age at first seizure in months 
# | `epi_sfreq` | `factor`     | Seizure frequency per year: None, 1-12, 13-52, 1/day
# | `epi_asm`   | `factor`     | Antiseizure medication: None, 1, 2+  
# | `adhd_inat` | `numeric`    | DSM inattention criteria
# | `adhd_hyp`  | `numeric`    | DSM hyperactivity criteria 
# | `session`   | `factor`     | BLAST Session number: 1, 2 or 3
# | `feedback`  | `factor`     | Accuracy per trial
# | `blast_RT`  | `numeric`    | Reaction time in ms per trial
# | `blast_mRT` | `numeric`    | Mean reaction time in ms per session
# | `blast_err` | `numeric`    | Error rate (%) per session
# | `blast_stab`| `numeric`    | Stability per session
# | `blast_int` | `numeric`    | Intensity per session
# | `ID_spiky`  | `factor`     | Presence of at least 1 IED: yes, no
# | `IED`       | `factor`     | Interictal epileptiform discharges: `"with"` or `"without"`
# | `IED_n`     | `numeric`    | Number of IED per block
# | `IED_mdur`  | `numeric`    | Mean duration of IED per block in seconds
# | `IED_dur`   | `numeric`    | Duration of IED per block in seconds (mean IED duration*number of IED)
# ----------------------------------------------------------------------------------------------------

# Descriptive analysis - Tables
# ----------------------------------------------------------------------------------------------------
datas %>% # Demographics and ADHD characteristics 
  group_by(patname) %>%
  slice(1) %>%
  ungroup() %>%
  filter(ID_spiky == "yes") %>% # subgroup with at least 1 IED
  dplyr::select(c(age, sex, adhd_inat, adhd_hyp, group)) %>%
  tbl_summary(by = "group",
              missing_text = "Missing values",
              digits = all_continuous() ~ 2,
              label = c(
                age ~ "Age (months)",
                sex ~ "Sex",
                adhd_inat ~ "Inattention criteria (DSM-V)",
                adhd_hyp ~ "Hyperactivity criteria (DSM-V)"
              ),
              statistic = list(
                all_continuous() ~ "{median} [{min}-{max}]",
                all_categorical() ~ "{p}% ({n} / {N})"
              )
  ) %>%
  add_overall() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Group**") %>%
  modify_caption("**Table 1. Demographics (children with >= 1 IED)**") %>%
  bold_labels() 

datas %>% # Epilepsy-related characteristics
  group_by(patname) %>%
  slice(1) %>%
  ungroup() %>%
  droplevels() %>%
  filter(ID_spiky == "yes") %>%
  dplyr::select(c(epi_onset, epi_type, epi_synd, epi_sfreq, epi_asm, group)) %>%
  tbl_summary(by = "group", 
              missing_text = "Missing values",
              label = c(
                epi_onset ~ "Age 1st seizure (months)", 
                epi_type ~ "Epilepsy type",
                epi_synd ~ "Epilepsy syndrome",
                epi_sfreq ~ "Seizure frequency / year",
                epi_asm ~ "Antiseizure medication"), 
              statistic = list(all_continuous() ~ "{median} [{min}-{max}]",
                               all_categorical() ~ "{p}% ({n} / {N})")) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Group**") %>%
  modify_caption("**Table 2. Epilepsy-related characteristics (children with >= 1 IED)**") %>%
  bold_labels()

datas %>% # IED-related characteristics 
  filter(ID_spiky == "yes") %>%
  droplevels() %>%
  dplyr::select(c(IED, IED_n, IED_mdur, IED_dur, group)) %>%
  tbl_summary(by = "group", 
              missing_text = "Missing values",
              digits = all_continuous() ~ 2,
              label = c(
                IED ~ "IED type", 
                IED_n ~ "IED / session",
                IED_mdur ~ "Mean IED duration / session (s)",
                IED_dur ~ "IED duration / session (s)"), 
              statistic = list(all_continuous() ~ "{mean} [{min}-{max}]",
                               all_categorical() ~ "{p}% ({n} / {N})")) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Group**") %>%
  modify_caption("**Table 3. IED characteristics (children with >= 1 IED)**") %>%
  bold_labels()
# ----------------------------------------------------------------------------------------------------

# Cumulative effect of IED - (n)
# ----------------------------------------------------------------------------------------------------
datas %>%
  filter(ID_spiky == "yes") %>%
  filter(feedback == "71") %>%
  dplyr::select(patname, session, group, age, blast_mRT, blast_err, blast_stab, blast_int, epi_asm, IED_n, epi_type, epi_sfreq, epi_onset) %>%
  group_by(patname, session) %>%
  dplyr::summarise(
    blast_mRT = unique(blast_mRT),
    blast_err = unique(blast_err) + 0.00001,
    blast_stab = unique(blast_stab) + 0.00001,
    blast_int = unique(blast_int) + 0.00001,
    IED_n = unique(IED_n),
    group = unique(group),
    ASM = unique(epi_asm),
    age = unique(age),
    epi_type = unique(epi_type),
    epi_sfreq = unique(epi_sfreq),
    epi_onset = unique(epi_onset)
  ) %>%
  ungroup() -> datas.to.model.n

## Reaction time 
mod.epiG.n.bloc.RT <- glm(blast_mRT ~ group + # Multicollinearity check
                            ASM +
                            poly(IED_n, 2) +
                            session +
                            poly(age, 2) +
                            epi_type +
                            epi_sfreq +
                            epi_onset,
                          data = datas.to.model.n,
                          family = Gamma())

performance::check_collinearity(mod.epiG.n.bloc.RT)
plot(performance::check_collinearity(mod.epiG.n.bloc.RT))

mod.epiG.n.bloc.RT <- glm(blast_mRT ~ group + # final model
                            session +
                            ASM +
                            poly(IED_n, 2) +
                            poly(age, 2) +
                            epi_type +
                            epi_sfreq +
                            epi_onset +
                            group:poly(age, 2),
                          data = datas.to.model.n,
                          family = Gamma())

car::Anova(mod.epiG.n.bloc.RT)
performance(mod.epiG.n.bloc.RT)

sjPlot::plot_model(mod.epiG.n.bloc.RT, # pos-hoc plot for IED
                   type = "pred",
                   terms = c("IED_n [all]"),
                   show.data = TRUE)

sjPlot::plot_model(mod.epiG.n.bloc.RT, # post-hoc for ASM
                   type = "pred",
                   terms = c("ASM"),
                   show.data = FALSE, 
                   ci.lvl = 0.95)

PH.emmeans <- emmeans(mod.epiG.n.bloc.RT, "ASM")
pairs(PH.emmeans, adjust = "FDR", side = "both", reverse = TRUE)

des_cont <- list("Mono - None" = c(-1, 1, 0),
                 "Poly - None" = c(-1, 0, 1),
                 "Poly - Mono" = c(0, -1, 1))

contrast(regrid(PH.emmeans), des_cont, adjust = "FDR", side = "both")

rm(PH.emmeans, des_cont)

sjPlot::plot_model(mod.epiG.n.bloc.RT, # post-hoc for age of onset
                   type = "pred",
                   terms = c("epi_onset"),
                   show.data = TRUE)


## Error percentage 
mod.epiG.n.bloc.err <- glm(PCT_error ~ group + 
                             ASM +
                             poly(IED_n, 2) +
                             session +
                             poly(age, 2) +
                             epi_type +
                             # Seizure frequency excluded due to multicollinearity
                             epi_onset,
                           data = datas.to.model.n,
                           family = Gamma())

performance::check_collinearity(mod.epiG.n.bloc.err)
plot(performance::check_collinearity(mod.epiG.n.bloc.err))

mod.epiG.n.bloc.err <- glm(PCT_error ~ # final model
                             group + 
                             session +
                             ASM +
                             poly(IED_n, 2) +
                             poly(age, 2) +
                             epi_type +
                             # Seizure frequency excluded
                             epi_onset + 
                             # session:poly(age, 2) removed
                             group:poly(age, 2),
                           data = datas.to.model.n,
                           family = Gamma())

car::Anova(mod.epiG.n.bloc.err)
performance(mod.epiG.n.bloc.err)

sjPlot::plot_model(mod.epiG.n.bloc.err, # post-hoc ASM
                   type = "pred",
                   terms = c("ASM"),
                   show.data = FALSE)

PH.emmeans <- emmeans(mod.epiG.n.bloc.err, "ASM")
pairs(PH.emmeans, adjust = "FDR", side = "both", reverse = TRUE)

des_cont <- list("Mono - None" = c(-1, 1, 0),
                 "Poly - None" = c(-1, 0, 1),
                 "Poly - Mono" = c(0, -1, 1))

contrast(regrid(PH.emmeans), 
         des_cont, 
         adjust = "FDR", 
         side = "both")

rm(PH.emmeans, des_cont)


sjPlot::plot_model(mod.epiG.n.bloc.err, # post-hoc IED
                   type = "pred",
                   terms = c("IED_n [all]"),
                   show.data = TRUE)

## Stability 
mod.epiG.n.bloc.stab <- glm(blast_stab ~ 
                              group +
                              epi_asm +
                              poly(IED_n, 2) +
                              session +
                              poly(age, 2) +
                              epi_type +
                              epi_sfreq +
                              poly(epi_onset, 2),
                            data = datas.to.model.n,
                            family = Gamma())

performance::check_collinearity(mod.epiG.n.bloc.stab)
plot(performance::check_collinearity(mod.epiG.n.bloc.stab))

mod.epiG.n.bloc.stab <- glm(blast_stab ~ 
                              group + 
                              session +
                              epi_asm +
                              poly(IED_n, 2) +
                              poly(age, 2) +
                              epi_type +
                              epi_sfreq +
                              poly(epi_onset, 2), 
                            data = datas.to.model.n,
                            family = Gamma())

car::Anova(mod.epiG.n.bloc.stab)
performance(mod.epiG.n.bloc.stab)

sjPlot::plot_model(mod.epiG.n.bloc.stab, # post-hoc for ASM
                   type = "pred",
                   terms = c("epi_asm"),
                   show.data = FALSE, 
                   ci.lvl = 0.95)

PH.emmeans <- emmeans(mod.epiG.n.bloc.stab, "epi_asm")
pairs(PH.emmeans, adjust = "FDR", side = "both", reverse = TRUE)

des_cont <- list("None - 1" = c(-1, 1, 0),
                 "None - 2+" = c(-1, 0, 1),
                 "1 - 2+" = c(0, -1, 1))

contrast(regrid(PH.emmeans), 
         des_cont, 
         adjust = "FDR", 
         side = "both")

rm(PH.emmeans, des_cont)

sjPlot::plot_model(mod.epiG.n.bloc.stab, # post-hoc for IED
                   type = "pred",
                   terms = c("IED_n [all]"),
                   show.data = TRUE)

sjPlot::plot_model(mod.epiG.n.bloc.stab, # post-hoc for age of onset
                   type = "pred",
                   terms = c("epi_onset"),
                   show.data = TRUE)

## Intensity 
mod.epiG.n.bloc.int <- glm(blast_int ~ 
                             group + 
                             epi_asm + 
                             poly(IED_n, 2) + 
                             session + 
                             poly(age, 2) + 
                             epi_type + 
                             epi_sfreq + 
                             poly(epi_onset, 2), 
                           data = datas.to.model.n,
                           family = Gamma())

performance::check_collinearity(mod.epiG.n.bloc.int)
plot(performance::check_collinearity(mod.epiG.n.bloc.int))

mod.epiG.n.bloc.int <- glm(blast_int ~ # final model
                             group + 
                             session + 
                             epi_asm + 
                             poly(IED_n, 2) + 
                             poly(age, 2) + 
                             epi_type + 
                             epi_sfreq + 
                             poly(epi_onset, 2), 
                           data = datas.to.model.n,
                           family = Gamma())

car::Anova(mod.epiG.n.bloc.int)
performance(mod.epiG.n.bloc.int)

sjPlot::plot_model(mod.epiG.n.bloc.int, # post-hoc for ASM
                   type = "pred",
                   terms = c("epi_asm"),
                   show.data = FALSE, 
                   ci.lvl = 0.95)

PH.emmeans <- emmeans(mod.epiG.n.bloc.int, "epi_asm")
pairs(PH.emmeans, adjust = "FDR", side = "both", reverse = TRUE)

des_cont <- list("None - 1" = c(-1, 1, 0),
                 "None - 2+" = c(-1, 0, 1),
                 "1 - 2+" = c(0, -1, 1))

contrast(regrid(PH.emmeans), 
         des_cont, 
         adjust = "FDR", 
         side = "both")

rm(PH.emmeans, des_cont)

sjPlot::plot_model(mod.epiG.n.bloc.int, # post-hoc for IED
                   type = "pred",
                   terms = c("IED_n [all]"),
                   show.data = TRUE)

sjPlot::plot_model(mod.epiG.n.bloc.int, # post-hoc for age of onset
                   type = "pred",
                   terms = c("epi_onset"),
                   show.data = TRUE)
# ----------------------------------------------------------------------------------------------------

# Local effect of IED 
# ----------------------------------------------------------------------------------------------------
datas %>%
  filter(ID_spiky == "yes") %>%
  filter(feedback == 71) %>%
  dplyr::select(patname, session, group, age,
                blast_RT, IED) %>%
  mutate(IED = as.factor(IED)) -> datas.to.model

mod.trial.RT <- glm(blast_RT ~ session + 
                      group + 
                      IED +
                      poly(age, 2),
                    data = datas.to.model,
                    family = Gamma())

performance::check_collinearity(mod.trial.RT)
plot(performance::check_collinearity(mod.trial.RT))

mod.trial.RT <- mixed(blast_RT ~ session + # final model
                        group + 
                        IED +
                        poly(age, 2) +
                        (1|patname),
                      data = datas.to.model,
                      family = Gamma(),
                      method = "LRT",
                      expand_re = TRUE,
                      control = glmerControl(optimizer = "bobyqa", 
                                             optCtrl = list(maxfun = 2e5)))

anova(mod.trial.RT)
performance(mod.trial.RT)

sjPlot::plot_model(mod.trial.RT$full_model, # post-hoc for IED
                   type = "emm",
                   terms = c("Event.perturb"),
                   show.data = F)

# ----------------------------------------------------------------------------------------------------


# Bootstrapping analysis 
# ----------------------------------------------------------------------------------------------------
n_boot <- 10000
boot_results <- data.frame(intercept = numeric(n_boot), slope = numeric(n_boot))
set.seed(123)  

for (i in 1:n_boot) {
  boot_sample <- data.frame()
  for (participant_id in unique(datas$patname)) {
    participant_data <- datas %>% filter(patname == participant_id)
    
    with_IED <- participant_data %>%    # Filter conditions: "with IED"
      filter(ID_spiky == "yes",
             feedback == 71,
             IED == "with")
    
    without_IED <- participant_data %>% # Filter conditions: "without IED"
      filter(ID_spiky == "yes",
             feedback == 71,
             IED == "without")
    
    sampled_without <- without_IED %>%
      sample_n(nrow(with_IED), replace = TRUE)

    participant_sample <- bind_rows(with_IED, sampled_without)
    boot_sample <- bind_rows(boot_sample, participant_sample)
  }
  
  model <- lm(blast_RT ~ IED, data = boot_sample)
  boot_results$intercept[i] <- coef(model)[1]  # RT for "without"
  boot_results$slope[i] <- coef(model)[2]      # Difference: "with" - "without"
}

mean_intercept <- mean(boot_results$intercept)
mean_slope <- mean(boot_results$slope)

conf_int_intercept <- quantile(boot_results$intercept, c(0.025, 0.975))
conf_int_slope <- quantile(boot_results$slope, c(0.025, 0.975))

cat("Bootstrap Results for blast_RT ~ IED (\"without\" as reference):\n")
cat("Mean Intercept (RT for 'without'): ", round(mean_intercept, 2), "ms\n")
cat("95% CI for Intercept: [", round(conf_int_intercept[1], 2), ", ", round(conf_int_intercept[2], 2), "]\n")
cat("Mean Slope (Difference in RT 'with' - 'without'): ", round(mean_slope, 2), "ms\n")
cat("95% CI for Slope: [", round(conf_int_slope[1], 2), ", ", round(conf_int_slope[2], 2), "]\n")
# ----------------------------------------------------------------------------------------------------

# Figures
# ----------------------------------------------------------------------------------------------------
datas %>%
  filter(ID_spiky == "yes") %>%
  filter(feedback == "71") %>%
  dplyr::select(patname, session, group, age, blast_mRT, blast_err, blast_stab, blast_int, epi_asm, IED_n, epi_type, epi_sfreq, epi_onset) %>%
  group_by(patname, session) %>%
  dplyr::summarise(
    blast_mRT = unique(blast_mRT),
    blast_err = unique(blast_err) + 0.00001,
    blast_stab = unique(blast_stab) + 0.00001,
    blast_int = unique(blast_int) + 0.00001,
    IED_n = unique(IED_n),
    group = unique(group),
    ASM = unique(epi_asm),
    age = unique(age),
    epi_type = unique(epi_type),
    epi_sfreq = unique(epi_sfreq),
    epi_onset = unique(epi_onset)
  ) %>%
  ungroup() -> graph

# Figure 2
graph %>%
  ggplot(aes(x = IED_n, y = blast_mRT)) +
  geom_smooth(method = "lm", formula =  y ~ poly(x, 2), colour = "#993356", fill = "#993356") +
  geom_point(colour = "#8da0cb",  width = .1, size = 2, alpha = .8) +
  theme_light() + 
  theme(legend.position = "none", 
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) + 
  labs(x = "IED quantity", y = "RT (mean / session)", tag = 'A)') -> RT

graph %>%
  ggplot(aes(x = IED_n, y = blast_err)) +
  geom_smooth(method = "lm", formula =  y ~ poly(x, 2), colour = "#993356", fill = "#993356") +
  geom_point(colour = "#8da0cb",  width = .1, size = 2, alpha = .8) +
  theme_light() + 
  theme(legend.position = "none", 
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) + 
  labs(x = "IED quantity", y = "Error % (/ session)", tag = 'B)') -> ERR

graph %>%
  ggplot(aes(x = IED_n, y = blast_stab)) +
  geom_smooth(method = "lm", formula =  y ~ poly(x, 2), colour = "#993356", fill = "#993356") +
  geom_point(colour = "#8da0cb",  width = .1, size = 2, alpha = .8) +
  theme_light() + 
  theme(legend.position = "none", 
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) + 
  labs(x = "IED quantity", y = "BLAST-Stability (/ session)", tag = 'C)') -> STAB

graph %>%
  ggplot(aes(x = IED_n, y = blast_int)) +
  geom_smooth(method = "lm", formula =  y ~ poly(x, 2), colour = "#993356", fill = "#993356") +
  geom_point(colour = "#8da0cb",  width = .1, size = 2, alpha = .8) +
  theme_light() + 
  theme(legend.position = "none", 
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) + 
  labs(x = "IED quantity", y = "BLAST-Intensity (/ session)", tag = 'D)') -> INT

grid.arrange(RT, ERR, STAB, INT, ncol = 2, nrow = 2)

# Figure 3
datas.to.model %>%
  ggplot(aes(x = IED, y = blast_RT, fill = IED)) +
  geom_flat_violin(aes(fill = IED), position = position_nudge(x = .25, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) + 
  geom_jitter(aes(x = IED, y = blast_RT, colour = IED), width = .2, size = 2, alpha = .4) + 
  geom_boxplot(aes(x = IED, y = blast_RT, fill = IED), outlier.shape = NA, notch = TRUE, alpha = .4, width = .4, colour = "#636363") +
  scale_fill_manual(values = c("#8da0cb", "#993356")) +
  scale_colour_manual(values = c("#8da0cb", "#993356")) +
  labs(y = "Reaction time (ms)", x = "IED") +
  scale_x_discrete(labels = c("Without", "With")) +
  theme_light() +
  theme(legend.position="none",
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  geom_signif(
    y_position = 2.5, xmin = c(1.3), xmax = c(1.7),
    annotations = c("***"), tip_length = 0.01, textsize = 10, color="black") 

# Figure S1
graph %>%
  ggplot(aes(x = epi_onset, y = blast_mRT)) +
  labs(x = "Age at epilepsy onset (months)", y = "RT (mean / session)", tag = "A)") +
  geom_smooth(method = "lm", formula =  y ~ poly(x, 2), colour = "#7fc97f", fill = "#7fc97f") + 
  geom_jitter(colour = "#8da0cb",  width = .1, size = 2, alpha = .8) +
  theme_light() + 
  scale_x_continuous(,sec.axis = sec_axis(transform = ~ . /12, breaks = c(0,3,6,9,12,15,18), name = "Age at epilepsy onset (years)")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) -> RT

graph %>%
  ggplot(aes(x = epi_onset, y = blast_stab)) +
  labs(x = "Age at epilepsy onset (months)", y = "BLAST-Stability (/ session)", tag = "B)") +
  geom_smooth(method = "lm", formula =  y ~ poly(x, 2), colour = "#7fc97f", fill = "#7fc97f") + 
  geom_jitter(colour = "#8da0cb",  width = .1, size = 2, alpha = .8) +
  theme_light() + 
  scale_x_continuous(,sec.axis = sec_axis(transform = ~ . /12, breaks = c(0,3,6,9,12,15,18), name = "Age at epilepsy onset (years)")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) -> STAB

graph %>%
  ggplot(aes(x = epi_onset, y = blast_int)) +
  labs(x = "Age at epilepsy onset (months)", y = "BLAST-Intensity (/ session)", tag = "C)") +
  geom_smooth(method = "lm", formula =  y ~ poly(x, 2), colour = "#7fc97f", fill = "#7fc97f") + 
  geom_jitter(colour = "#8da0cb",  width = .1, size = 2, alpha = .8) +
  theme_light() + 
  scale_x_continuous(,sec.axis = sec_axis(transform = ~ . /12, breaks = c(0,3,6,9,12,15,18), name = "Age at epilepsy onset (years)")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) -> INT

grid.arrange(RT, STAB, INT, ncol = 3, nrow = 1)

### Figure S2
ggplot(graph, aes(x = ASM, y = blast_mRT, fill = ASM)) +
  geom_flat_violin(aes(fill = ASM), position = position_nudge(x = .25, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) + 
  geom_jitter(aes(x = ASM, y = blast_mRT, colour = ASM), width = .1, size = 2, alpha = .4) + 
  geom_boxplot(aes(x = ASM, y = blast_mRT, fill = ASM), outlier.shape = NA, notch = TRUE, alpha = .5, width = .4, colour = "#636363") +
  scale_fill_manual(values = c("#993356", "#7fc97f", "#8da0cb")) +
  scale_colour_manual(values = c("#993356", "#7fc97f", "#8da0cb")) +
  labs(y = "RT (mean / session)", x = "ASM", tag = "A)") +
  theme_light() +
  theme(legend.position="none",
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  geom_signif(
    y_position = 2, xmin = c(2.3), xmax = c(2.8),
    annotations = c("**"), tip_length = 0.01, textsize = 6, color="black") +
  geom_signif(
    y_position = 2.4, xmin = c(1.3), xmax = c(2.8),
    annotations = c("***"), tip_length = 0.01, textsize = 6, color="black") -> RT

ggplot(graph, aes(x = ASM, y = blast_err, fill = ASM)) +
  geom_flat_violin(aes(fill = ASM), position = position_nudge(x = .25, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) + 
  geom_jitter(aes(x = ASM, y = blast_err, colour = ASM), width = .1, size = 2, alpha = .4) + 
  geom_boxplot(aes(x = ASM, y = blast_err, fill = ASM), outlier.shape = NA, notch = TRUE, alpha = .4, width = .4, colour = "#636363") +
  scale_fill_manual(values = c("#993356", "#7fc97f", "#8da0cb")) +
  scale_colour_manual(values = c("#993356", "#7fc97f", "#8da0cb")) +
  labs(y = "Error % (/ session)", x = "ASM", tag = "B)") +
  theme_light() +
  theme(legend.position="none",
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  geom_signif(
    y_position = 40, xmin = c(1.3), xmax = c(1.8),
    annotations = c("*"), tip_length = 0.01, textsize = 6, color="black") +
  geom_signif(
    y_position = 40, xmin = c(2.3), xmax = c(2.8),
    annotations = c("**"), tip_length = 0.01, textsize = 6, color="black") -> ERR

ggplot(graph, aes(x = ASM, y = blast_stab, fill = ASM)) +
  geom_flat_violin(aes(fill = ASM), position = position_nudge(x = .25, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) + 
  geom_jitter(aes(x = ASM, y = blast_stab, colour = ASM), width = .1, size = 2, alpha = .4) + 
  geom_boxplot(aes(x = ASM, y = blast_stab, fill = ASM), outlier.shape = NA, notch = TRUE, alpha = .4, width = .4, colour = "#636363") +
  scale_fill_manual(values = c("#993356", "#7fc97f", "#8da0cb")) +
  scale_colour_manual(values = c("#993356", "#7fc97f", "#8da0cb")) +
  labs(y = "BLAST-Stability (/ session)", x = "ASM", tag = "C)") +
  theme_light() +
  theme(legend.position="none",
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  geom_signif(
    y_position = 75, xmin = c(1.3), xmax = c(2.9),
    annotations = c("**"), tip_length = 0.01, textsize = 6, color="black") +
  geom_signif(
    y_position = 60, xmin = c(1.4), xmax = c(1.8),
    annotations = c("*"), tip_length = 0.01, textsize = 6, color="black") +
  geom_signif(
    y_position = 60, xmin = c(2.5), xmax = c(2.9),
    annotations = c("***"), tip_length = 0.01, textsize = 6, color="black") -> STAB

ggplot(graph, aes(x = ASM, y = blast_int, fill = ASM)) +
  geom_flat_violin(aes(fill = ASM), position = position_nudge(x = .25, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) + 
  geom_jitter(aes(x = ASM, y = blast_int, colour = ASM), width = .1, size = 2, alpha = .4) + 
  geom_boxplot(aes(x = ASM, y = blast_int, fill = ASM), outlier.shape = NA, notch = TRUE, alpha = .5, width = .4, colour = "#636363") +
  scale_fill_manual(values = c("#993356", "#7fc97f", "#8da0cb")) +
  scale_colour_manual(values = c("#993356", "#7fc97f", "#8da0cb")) +
  labs(y = "BLAST-Intensity (/ session)", x = "ASM", tag = "D)") +
  theme_light() +
  theme(legend.position="none",
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  geom_signif(
    y_position = 80, xmin = c(2.4), xmax = c(2.8),
    annotations = c("***"), tip_length = 0.01, textsize = 6, color="black") +
  geom_signif(
    y_position = 110, xmin = c(1.4), xmax = c(2.8),
    annotations = c("**"), tip_length = 0.01, textsize = 6, color="black") -> INT

grid.arrange(RT, ERR, STAB, INT, ncol = 2, nrow = 2)

# Figure S3
hist(boot_results$slope,
     main = "Bootstrap Distribution of Slope",
     xlab = "Slope (RT difference: with - without IED)",
     breaks = 30)
abline(v = mean_slope, col = "#993356", lwd = 2)
abline(v = conf_int_slope, col = "#8da0cb", lty = 2, lwd = 2)

hist(boot_results$intercept,
     main = "Bootstrap Distribution of Intercept",
     xlab = "Intercept (RT for 'without' IED)",
     breaks = 30)
abline(v = mean_intercept, col = "#993356", lwd = 2)
abline(v = conf_int_intercept, col = "#8da0cb", lty = 2, lwd = 2)
# ----------------------------------------------------------------------------------------------------
