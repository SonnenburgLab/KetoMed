
library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)   
library(emmeans)    

#Read and prepare data 
path <- "ketomed_gi_questionaire.xlsx"
raw  <- read_excel(path, sheet = "data")

bristol_col <- names(raw)[grepl("choose the type", names(raw), ignore.case = TRUE)]
stopifnot(length(bristol_col) == 1)

d <- raw %>%
  transmute(
    id        = factor(record_id),
    event     = redcap_event_name,
    phase_raw = tp,
    order     = factor(`first diet`),     # "Keto first" / "Med first"
    diagnosis = factor(diagnosis),         # "T2 Diabetes" / "Pre-diabetes"
    bristol   = suppressWarnings(as.numeric(.data[[bristol_col]]))  # 1-7; blanks -> NA
  )

# Participant 8016 has repeat visit rows (5b/6b_arm_1); keep the standard visits
d <- d %>% filter(!grepl("b_arm_1$", event))
stopifnot(all(table(d$id, d$phase_raw) <= 1))   # one row per id x phase

# Baseline is the reference level
keep_phases <- c("baseline", "keto del", "keto own", "med del", "med own")
d <- d %>%
  filter(phase_raw %in% keep_phases) %>%
  mutate(phase = factor(phase_raw, levels = keep_phases))

#Descriptive
desc <- d %>%
  group_by(phase) %>%
  summarise(n      = sum(!is.na(bristol)),
            median = median(bristol, na.rm = TRUE),
            mean   = round(mean(bristol, na.rm = TRUE), 2),
            sd     = round(sd(bristol,   na.rm = TRUE), 2),
            .groups = "drop")
print(desc)

#Paired Wilcoxon
np_vs_base <- lapply(keep_phases[-1], function(ph) {
  sub <- d %>% filter(phase %in% c("baseline", ph)) %>%
    select(id, phase, bristol) %>%
    pivot_wider(names_from = phase, values_from = bristol) %>% drop_na()
  wt <- wilcox.test(sub[["baseline"]], sub[[ph]], paired = TRUE, exact = FALSE)
  data.frame(comparison = paste(ph, "vs baseline"),
             n_pairs = nrow(sub), p_raw = wt$p.value)
}) %>% bind_rows() %>%
  mutate(p_BH = p.adjust(p_raw, method = "BH"))
print(np_vs_base)

#LMM-- each diet phase vs baseline
m_main   <- lmer(bristol ~ phase + (1 | id), data = d)
emm_main <- emmeans(m_main, ~ phase)
lmm_vs_base <- summary(contrast(emm_main, method = "trt.vs.ctrl", ref = 1),
                       infer = TRUE)
print(lmm_vs_base)

#Does the phase effect depend on diet order or diagnosis?
# Read the phase:order and phase:diagnosis rows (Pr(>F)).
m_order <- lmer(bristol ~ phase * order + (1 | id), data = d)
print(anova(m_order))

m_dx <- lmer(bristol ~ phase * diagnosis + (1 | id), data = d)
print(anova(m_dx))
