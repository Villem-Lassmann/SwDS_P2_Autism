require(readxl)
require(dplyr)
require(tidyr)
require(mice)
require(stringr)
require(ggplot2)
require(gridExtra)
require(ggpubr)
require(naniar)


tfences <- function(x, k = 1.5, one.sided = F, printValues = F){
  qtls <- quantile(x, names = F, na.rm = T)
  iqr <- IQR(x, na.rm = T)
  if (printValues){
    print(qtls[2]-k*iqr)
    print(qtls[4]+k*iqr)
  }
  if (one.sided){
    ext <- (x < qtls[2]-k*iqr)
  } else {
    ext <- (x < qtls[2]-k*iqr) | (x > qtls[4]+k*iqr)
  }
  return(ext)
}

defw <- 6.78
defh <- 5.73

##############

base <- read_xlsx("Data/Raw Data.xlsx", na = "888") %>%
  mutate(brief_raw_working_memory = as.numeric(brief_raw_working_memory),
         pvt_mean_rt = as.double(gsub("11.40.8333333", "1140.8333333", pvt_mean_rt)), # Assuming this was an entry mistake
         SCQ = na_if(SCQ, 0),
         wasi_sum_rawscores = na_if(wasi_sum_rawscores, 0),
         gender = factor(gender),
         diagnosis = factor(diagnosis),
         where_english = ifelse(where_english == 4 & age_acquisition < 5, NA, where_english),
         where_english = factor(where_english)) %>%
  mutate_at(vars(starts_with("brief")), ~ replace(., . == 0, NA)) %>%
  rename(brief_raw_self_monitor = `brief_raw_self-monitor`)

# nor <- base %>% filter(diagnosis == 0)
# aut <- base %>% filter(diagnosis == 1)


##############


mdv <- base %>%
  mutate(emo = factor(ifelse(bilec_home_output > bilec_english_output, "Secondary", "Primary"))) %>%
  select(-part_no, -SCQ, -tomi_compmean, -et_falsebelief_testtrial_preference_score, -bilec_home_input, -bilec_english_input,
         -bilec_home_output, -bilec_english_output)%>%
  mutate_at(c(4:16,19:35), ~ replace(., tfences(., 2.25), NA))
names(mdv) <- c("g", "a", "d", "bpvs", "vpst", "wasi", "te", "tb", "ta", "tt", "fti", "ftni",
             "fbc", "fbi", "lti", "lto", "aa", "we", "bi", "bsm", "bs", "bec", "bini", "bwm", "bpo",
             "btm", "bom", "fpc", "fpi", "fmc", "fmi", "pr", "pl", "plr", "pfs", "emo")

miss1 <- gg_miss_upset(mdv, nsets = 15)

miss2 <- gg_miss_fct(mdv, d) + labs(x="Diagnosis", y = "Variable") + scale_x_discrete(labels = c("Non-autistic", "Autistic"))

#ggsave("Images/Miss1.png", miss1, width = defw, height = defh)
ggsave("Images/Miss2.png", miss2, width = defw, height = defh)


##############


actual <- base %>%
  mutate(emo = factor(ifelse(bilec_home_output > bilec_english_output, "Secondary", "Primary"))) %>% #Input was always Primary
  select(-part_no, -SCQ, -tomi_compmean, -et_falsebelief_testtrial_preference_score, -bilec_home_input, -bilec_english_input,
         -bilec_home_output, -bilec_english_output) %>%
  mutate_at(c(4:16,19:35), ~ replace(., tfences(., 2.25), NA)) %>%
  mutate(tom_tb_totalscore = as.integer(tom_tb_totalscore),
         wasi_sum_rawscores = as.integer(wasi_sum_rawscores),
         bpvs_raw = as.integer(bpvs_raw),
         age_m = as.integer(age_m),
         age_acquisition = as.integer(age_acquisition),
         pvt_number_of_lapses = as.integer(pvt_number_of_lapses),
         pvt_count_falsestarts = as.integer(pvt_count_falsestarts)) %>%
  mutate_at(vars(starts_with("brief")), as.integer)

names(actual) <- c("g", "a", "d", "bpvs", "vpst", "wasi", "te", "tb", "ta", "tt", "fti", "ftni",
                   "fbc", "fbi", "lti", "lto", "aa", "we", "bi", "bsm", "bs", "bec", "bini", "bwm", "bpo",
                   "btm", "bom", "fpc", "fpi", "fmc", "fmi", "pr", "pl", "plr", "pfs", "emo")

act_predmat <- make.predictorMatrix(actual)
act_predmat[,3] <- 0

ini <- mice(actual, maxit = 0, predictorMatrix = act_predmat, remove.collinear = F)

meth <- ini$method
meth[c("vpst", "wasi", "te", "tb", "ta", "tt", "fti", "ftni", "fbc", "fbi", "fmc", "fmi")] <- "norm"

# seed <- round(runif(1,1,100000)) # So that the small and large would be from the same seed
seed <- 475621
act_imp_small <- mice(actual, m=5, maxit=100, predictorMatrix = act_predmat, remove.collinear = F, seed = seed)
act_imp_large <- mice(actual, m=60, maxit=100, predictorMatrix = act_predmat, remove.collinear = F, seed = seed)
act_comp <- complete(act_imp_small)

sp <- stripplot(act_imp_small, te + fti + bs + bpo + fmc + pl ~ .imp)
mix <- plot(act_imp_small, ta + bs + pl ~ .it | .ms) 
dp <- densityplot(act_imp_small, ~ te + fti + bec + fbc + fmi + pl)

#ggsave("Images/Stripplot.png", sp, width = defw, height = defh)
#ggsave("Images/Mixing.png", mix, width = defw, height = defh)

form <- "d ~ g + a + bpvs + vpst + wasi + te + tb + ta + tt + fti + ftni + 
                                       fbc + fbi + lti + lto + aa + we + bi + bsm + bs + bec + bini + 
                                       bwm + bpo + btm + bom + fpc + fpi + fmc + fmi + pr + pl + 
                                       plr + pfs"

# g + a + bpvs + vpst + wasi

fit_large1 <- with(act_imp_large, glm(d ~ a + bpvs, family = binomial("logit")), maxit = 1000000)
fit_large_sum1 <- summary(pool(fit_large1))
fit_large_sum1
pool(fit_large1)

# g + a + te + tb + ta + tt

fit_large2 <- with(act_imp_large, glm(d ~ a + ta, family = binomial("logit")), maxit = 1000000)
fit_large_sum2 <- summary(pool(fit_large2))
fit_large_sum2
pool(fit_large2)

# g + a + fti + ftni + fbc + fbi + lti + lto + aa +we

fit_large3 <- with(act_imp_large, glm(d ~ a + fbc + we, family = binomial("logit")), maxit = 1000000)
fit_large_sum3 <- summary(pool(fit_large3))
fit_large_sum3
pool(fit_large3)

# g + a + bi + bsm + bs + bec + bini + bwm + bpo + btm + bom

fit_large4 <- with(act_imp_large, glm(d ~ bs + bini + bwm, family = binomial("logit")), maxit = 1000000)
fit_large_sum4 <- summary(pool(fit_large4))
fit_large_sum4
pool(fit_large4)

# g + a + fpc + fpi + fmc + fmi + pr + pl + plr + pfs

fit_large5 <- with(act_imp_large, glm(d ~ a + fpi, family = binomial("logit")), maxit = 1000000)
fit_large_sum5 <- summary(pool(fit_large5))
fit_large_sum5
pool(fit_large5)

# vpst + wasi + ta + fti + bs + bini + bwm + fpi
## bs + bini + bwm + bom

fit_large2 <- with(act_imp_large, glm(d ~ bs + bini + 
                                        bwm + bom
                                        , family = binomial("logit")), maxit = 1000000)
summary(pool(fit_large2))


##############

diagn <- c("Non-autistic", "Autistic")
names(diagn) <- c("0","1")

age <- ggplot(actual, aes(x = a, y = ..density.., color = g)) +
  geom_histogram(alpha = 0.3, position = "identity", binwidth = 4, aes(fill = g)) +
  stat_density(adjust = 0.5, fill = NA, position = "identity", size = 1) +
  stat_density(aes(x=a, y=..density..), adjust = 0.5, color = "black", fill = NA, size = 1) +
  geom_vline(data = actual %>% select(g, d, a) %>% group_by(g, d) %>% summarise(ma = median(a)), 
             mapping = aes(xintercept = ma, color = g), linetype = "dashed", size = 1) +
  geom_vline(data = actual %>% select(d, a) %>% group_by(d) %>% summarise(ma = median(a)), 
             mapping = aes(xintercept = ma), linetype = "dashed", size = 1, color = "black") +
  facet_grid(d ~ ., labeller = labeller(d = diagn)) +
  scale_fill_discrete(name = "Gender", labels = c("Male", "Female")) + 
  scale_color_discrete(name = "Gender", labels = c("Male", "Female")) + 
  labs(x = "Age in months", y = "Density") +
  theme(legend.position = "bottom")

ggsave("Images/Age.png", age, width = defw, height = defh)


wasi <- ggplot(actual %>% filter(!is.na(wasi)), aes(y = wasi, fill = d)) +
  geom_boxplot() +
  scale_fill_discrete(name = "Diagnosis", labels = diagn) +
  labs(y = "WASI scores") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_flip()

#ggsave("Images/Wasi.png", wasi, width = defw, height = defh)

bpvs <- ggplot(actual %>% filter(!is.na(bpvs)), aes(y = bpvs, fill = d)) +
  geom_boxplot() +
  scale_fill_discrete(name = "Diagnosis", labels = diagn, guide = F) +
  labs(y = "BPVS scores") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_flip()

vpst <- ggplot(actual %>% filter(!is.na(vpst)), aes(y = vpst, fill = d)) +
  geom_boxplot() +
  scale_fill_discrete(name = "Diagnosis", labels = diagn) +
  labs(y = "Vocabulary processing speed (ms)") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "bottom") +
  coord_flip()

lang <- ggarrange(wasi, bpvs, vpst, nrow = 3, common.legend = T, legend = "bottom")

ggsave("Images/lang.png", lang, width = defw, height = defh*1.25)


tomi_data <- actual %>%
  select(d, te, tb, ta, tt) %>%
  rename(`TOMI Early` = te, `TOMI Basic` = tb, `TOMI Advanced` = ta, `TOMI Total` = tt) %>%
  pivot_longer(-d, names_to = "Test", values_to = "Score") %>%
  mutate(Test = factor(Test, levels = c("TOMI Early", "TOMI Basic", "TOMI Advanced", "TOMI Total"))) %>%
  filter(!is.na(Score))

tomi <- ggplot(tomi_data, aes(x = Test, y = Score, fill = d)) +
  geom_boxplot() +
  labs(x = "") +
  scale_fill_discrete(name = "Diagnosis", labels = diagn) +
  theme(legend.position = "bottom")

ggsave("Images/Tomi.png", tomi, width = defw, height = defh)


acq_data <- actual %>%
  select(d, aa, we) %>%
  filter(!is.na(aa), !is.na(we)) %>%
  mutate(we = factor(we, levels = 1:4, labels = c("Home", "Nursery", "Playgroup", "School")))

slang <- ggplot(acq_data, aes(x = aa, fill = d)) +
  geom_bar(position = "dodge", stat = "count") +
  facet_wrap(~ we) +
  scale_x_discrete(name = "Age of acquiring second language", limits = c(0:5)) +
  scale_fill_discrete(name = "Diagnosis", labels = diagn) +
  theme(legend.position = "bottom") +
  labs(y = "Number of children")

ggsave("Images/Seclang.png", slang, width = defw, height = defh)


bilec_data <- actual %>%
  select(d, emo, "Input" = lti, "Output" = lto) %>%
  pivot_longer(c("Input", "Output"), names_to = "Channel", values_to = "Score") %>%
  mutate(emo = factor(ifelse(Channel == "Input", 1, emo), labels =c("Primary", "Secondary"))) %>%
  filter(complete.cases(.))

bilec <- ggplot(bilec_data, aes(x = Score, color = d)) +
  geom_histogram(aes(y = ..density.., fill = d), binwidth = 5, position = "identity", alpha = 0.3) +
  stat_density(adjust = 0.5, fill = NA, position = "identity", size = 1) +
  stat_density(aes(x=Score, y=..density..), adjust = 0.5, color = "black", fill = NA, size = 1) +
  geom_vline(data = bilec_data %>% group_by(d, Channel, emo) %>% summarise(ma = median(Score)), 
             mapping = aes(xintercept = ma, color = d), linetype = "dashed", size = 1) +
  geom_vline(data = bilec_data  %>% group_by(Channel, emo) %>% summarise(ma = median(Score)), 
             mapping = aes(xintercept = ma), linetype = "dashed", size = 1, color = "black") +
  facet_grid(Channel + emo ~ .) +
  scale_fill_discrete(name = "Diagnosis", labels = diagn) + 
  scale_color_discrete(name = "Diagnosis", labels = diagn) + 
  labs(x = "Score", y = "Density") +
  theme(legend.position = "bottom")

ggsave("Images/Bilec.png", bilec, width = defw, height = defh)


etf_data <- actual %>%
  select(d, "At interacting" = fti, "At non-interacting" = ftni, "To correct" = fbc, 
         "To incorrect" = fbi) %>%
  pivot_longer(-d, names_to = "Test", values_to = "Time") %>%
  filter(complete.cases(.))

etf <- ggplot(etf_data, aes(x = Test, y = Time, fill = d)) + 
  geom_boxplot() +
  labs(x = "", y = "Time (ms)") +
  theme(legend.position = "bottom") +
  scale_fill_discrete(name = "Diagnosis", labels = diagn)

ggsave("Images/Etf.png", etf, width = defw, height = defh)


brief_data <- actual %>%
  select(d, "Inhibit" = bi, "Self monitor" = bsm, "Shift" = bs, "Emotional control" = bec, 
         "Initiate" = bini, "Working memory" = bwm, "Plan/organise" = bpo, "Task monitor" =  btm, "Organisation" = bom) %>%
  pivot_longer(-d, names_to = "Brief scale", values_to = "Score") %>%
  filter(complete.cases(.))

brief <- ggplot(brief_data, aes(x = `Brief scale`, y = Score, fill =d)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle=45, vjust=0.6),
        legend.position = "bottom") +
  scale_fill_discrete(name = "Diagnosis", labels = diagn)

ggsave("Images/Brief.png", brief, width = defw, height = defh)


flanker_data1 <- actual %>%
  select(d, "Congruent" = fpc, "Incongruent" = fpi) %>%
  pivot_longer(-d, names_to = "Task", values_to = "Percentage") %>%
  filter(complete.cases(.)) %>%
  mutate(Percentage = Percentage/100)

flanker1 <- ggplot(flanker_data1, aes(x = d, y = Percentage, fill = d)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.035) +
  #coord_flip() +
  facet_grid(Task ~ .) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom") +
  scale_fill_discrete(name = "Diagnosis", labels = diagn) +
  scale_y_continuous(labels = scales::percent)

flanker_data2 <- actual %>%
  select(d, "Congruent" = fmc, "Incongruent" = fmi) %>%
  pivot_longer(-d, names_to = "Task", values_to = "Time") %>%
  filter(complete.cases(.))

flanker2 <- ggplot(flanker_data2, aes(x = Time, color = d)) +
  geom_histogram(aes(y = ..density.., fill = d), binwidth = 100, position = "identity", alpha = 0.3) +
  stat_density(adjust = 0.5, fill = NA, position = "identity", size = 1) +
  stat_density(aes(x=Time, y=..density..), adjust = 0.5, color = "black", fill = NA, size = 1) +
  geom_vline(data = flanker_data2 %>% group_by(d, Task) %>% summarise(ma = median(Time)), 
             mapping = aes(xintercept = ma, color = d), linetype = "dashed", size = 1) +
  geom_vline(data = flanker_data2  %>% group_by(Task) %>% summarise(ma = median(Time)), 
             mapping = aes(xintercept = ma), linetype = "dashed", size = 1, color = "black") +
  facet_grid(Task ~ .) +
  scale_fill_discrete(name = "Diagnosis", labels = diagn) + 
  scale_color_discrete(name = "Diagnosis", labels = diagn) + 
  labs(x = "Time (ms)", y = "Density") +
  theme(legend.position = "bottom")

flanker <- ggarrange(flanker2, flanker1, nrow = 2, common.legend = T, legend = "bottom")

ggsave("Images/Flanker.png", flanker, width = defw, height = defh*1.5)


pvt_data1 <- actual %>%
  select(d, "Valid response" = pr, "Lapse response" = plr) %>%
  pivot_longer(-d, names_to = "Cat", values_to = "Time") %>%
  filter(complete.cases(.))

pvt1 <- ggplot(pvt_data1, aes(x = Time, color = d)) +
  geom_histogram(aes(y = ..density.., fill = d), binwidth = 100, position = "identity", alpha = 0.3) +
  stat_density(adjust = 0.5, fill = NA, position = "identity", size = 1) +
  stat_density(aes(x=Time, y=..density..), adjust = 0.5, color = "black", fill = NA, size = 1) +
  geom_vline(data = pvt_data1 %>% group_by(d, Cat) %>% summarise(ma = median(Time)), 
             mapping = aes(xintercept = ma, color = d), linetype = "dashed", size = 1) +
  geom_vline(data = pvt_data1  %>% group_by(Cat) %>% summarise(ma = median(Time)), 
             mapping = aes(xintercept = ma), linetype = "dashed", size = 1, color = "black") +
  facet_grid(Cat ~ .) +
  scale_fill_discrete(name = "Diagnosis", labels = diagn) + 
  scale_color_discrete(name = "Diagnosis", labels = diagn) + 
  labs(x = "Time (ms)", y = "Density") +
  theme(legend.position = "bottom")

pvt_data2 <- actual %>%
  select(d, "Lapses" = pl, "False starts" = pfs) %>%
  pivot_longer(-d, names_to = "Cat", values_to = "Count") %>%
  filter(complete.cases(.))

pvt2 <- ggplot(pvt_data2, aes(x = d, y = Count, fill = d)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2) +
  #coord_flip() +
  facet_grid(Cat ~ .) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom") +
  scale_fill_discrete(name = "Diagnosis", labels = diagn)

pvt <- ggarrange(pvt1, pvt2, nrow = 2, common.legend = T, legend = "bottom")

ggsave("Images/Pvt.png", pvt, width = defw, height = defh*1.5)
