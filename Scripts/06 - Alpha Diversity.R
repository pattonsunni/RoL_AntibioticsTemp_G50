# Author: Sunni Patton
# Last edited: 10/18/24
# Title: Alpha diversity 

## Set seed ====
set.seed(123)

## Load libraries ====
library(picante)
library(phyloseq)
library(car)
library(stats)
library(ggpubr)
library(DescTools)
library(rstatix)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(pbkrtest)
library(microViz)
library(fitdistrplus)
library(glmmTMB)
library(performance)
library(dharma)


## Load data ====
readRDS(here::here("Output Files/05 - Relative Abundance - Output/ps.rare.rds")) -> ps.rare

## Calculate alpha diversity ====
# Agglomerate to genus level
ps.rare.genus <- tax_glom(ps.rare, taxrank = "Genus", NArm=FALSE)

# Calculate
estimate_richness(ps.rare.genus, measures = c("Observed", "InvSimpson", "Shannon", "Chao1", "Simpson")) -> alphaDiv.genus

### Add Samples column
alphaDiv.genus$Samples <- rownames(alphaDiv.genus)

### Combine alpha metrics with sam_data
merge(as.matrix(ps.rare@sam_data), alphaDiv.genus, by = "Samples") -> alphaDiv.genus

write.csv(alphaDiv.genus, here::here("alphaDiv.genus.csv"))

## Look into data distributions shannon ====
# See what distibution is closest to 
descdist(alphaDiv.genus$Shannon, boot=500)   ## Shannon seems to have beta distribution (but isn't because values aren't from 0-1)

# Create plots to visualize different distribution fits
shan_wb <- fitdist(alphaDiv.genus$Shannon, "weibull")
shan_gam <- fitdist(alphaDiv.genus$Shannon, "gamma")
shan_exp <- fitdist(alphaDiv.genus$Shannon, "exp")
shan_lognorm <- fitdist(alphaDiv.genus$Shannon, "lnorm")
shan_norm <- fitdist(alphaDiv.genus$Shannon, "norm")


par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "gamma", "expo", "lognorm", "norm")
denscomp(list(shan_wb, shan_gam, shan_exp, shan_lognorm, shan_norm), legendtext = plot.legend)
qqcomp(list(shan_wb, shan_gam, shan_exp, shan_lognorm, shan_norm), legendtext = plot.legend)
cdfcomp(list(shan_wb, shan_gam, shan_exp, shan_lognorm, shan_norm), legendtext = plot.legend)
ppcomp(list(shan_wb, shan_gam, shan_exp, shan_lognorm, shan_norm), legendtext = plot.legend)

## Goodness of fit tests to compare models
gofstat(list(shan_wb, shan_gam, shan_lognorm, shan_norm))

gof = gofTest(alphaDiv.genus$Shannon,distribution = "lnorm", test = "ks")
print(gof) ## lognormal distribution
gof = gofTest(alphaDiv.genus$Shannon,distribution = "weibull", test = "ks")
print(gof) ## gamma distribution and weibull
gof = gofTest(alphaDiv.genus$Shannon,distribution = "gamma", test = "chisq")
print(gof) 

## Model testing (Gamma distribution)
### Shannon as response variable, treatment (pretreatment) and time as predictor variables, tankID as random effect, 
r1 <- glmer(Shannon ~ Pretreatment + Time + Pretreatment*Time + (1|TankID), data = alphaDiv.genus, family = Gamma(link = "inverse"))
# boundary is singular (i.e., not great -- really don't have enough data to fit this model)
summary(r1)
plot(residuals(r1))
qqnorm(residuals(r1))
plot(simulate_residuals(r1))

Anova(r1)
emmeans(r1, pairwise ~ Time | Pretreatment, adjust = "tukey")


## No interaction
r2 <- glmer(Shannon ~ Pretreatment + Time + (1|TankID), data = alphaDiv.genus, family = Gamma(link = "inverse"))
anova(r2)
emmeans(r1, pairwise ~ Pretreatment | Time, adjust = "tukey")

## No random effect
r3 <- glm(Shannon ~ Pretreatment + Time + Pretreatment*Time, data = alphaDiv.genus, family = Gamma(link = "inverse"))
Anova(r3)
emmeans(r3, pairwise ~ Pretreatment | Time, adjust = "tukey")

## seems like tankID is overparametizing the model so it gives the is singular warning; when I remove
## it, it goes away; tank ID is perfectly aligned with treatment (no two tanks share the same ID across treatments/
## tankID only exists within one treatment group)


### Shannon as response variable, treatment (pretreatment) and time as predictor variables, no tankID
r2 <- glm(Shannon ~ Pretreatment + Time + Pretreatment*Time, data = alphaDiv.genus, family = Gamma(link = "inverse"))
## no is singular warning
summary(r2) # AIC: -109
plot(residuals(r2))
qqnorm(residuals(r2))


emmeans(r2, pairwise ~ Pretreatment | Time, adjust = "tukey")
### Some differences in significance, some pretty major ones maybe?

## Look into distribution simpson ====
## probably beta bc it's bound by 0 and 1
# See what distibution is closest to 
descdist(alphaDiv.genus$Observed, boot=500)   ## Shannon seems to have beta distribution (but isn't because values aren't from 0-1)

# Create plots to visualize different distribution fits
obs_wb <- fitdist(alphaDiv.genus$Observed, "weibull")
obs_gam <- fitdist(alphaDiv.genus$Observed, "gamma")
obs_lognorm <- fitdist(alphaDiv.genus$Observed, "lnorm")
obs_norm <- fitdist(alphaDiv.genus$Observed, "norm")

simp_beta <- fitdist(alphaDiv.genus$Observed, "beta")


par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "gamma", "lognorm", "norm")
denscomp(list(obs_wb, obs_gam, obs_lognorm, obs_norm), legendtext = plot.legend)
qqcomp(list(obs_wb, obs_gam, obs_lognorm, obs_norm), legendtext = plot.legend)
cdfcomp(list(obs_wb, obs_gam, obs_lognorm, obs_norm), legendtext = plot.legend)
ppcomp(list(obs_wb, obs_gam, obs_lognorm, obs_norm), legendtext = plot.legend)

bm1 <- glmmTMB(Simpson ~ Pretreatment + Time + Pretreatment*Time + (1|TankID), data=alphaDiv.genus, family='beta_family')
plot(residuals(bm1))
qqnorm(residuals(bm1))



## Function for calculating summary data (good for line plot) ====
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Calculate mean and standard error ====
summarySE(alphaDiv.genus, measurevar="Shannon", groupvars=c("Treatment_Long","Time_TotalDays")) -> shannon.genus.stat
summarySE(alphaDiv.genus, measurevar="InvSimpson", groupvars=c("Treatment_Long","Time_TotalDays")) -> invSimp.genus.stat
summarySE(alphaDiv.genus, measurevar="Observed", groupvars=c("Treatment_Long","Time_TotalDays")) -> obs.genus.stat
summarySE(alphaDiv.genus, measurevar="Chao1", groupvars=c("Treatment_Long","Time_TotalDays")) -> chao1.genus.stat


## Plot line plot ====
pd <- position_dodge(0.1) # move them .05 to the left and right

## Set order
shannon.genus.stat$Treatment_Long <- factor(shannon.genus.stat$Treatment_Long, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
obs.genus.stat$Treatment_Long <- factor(obs.genus.stat$Treatment_Long, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
invSimp.genus.stat$Treatment_Long <- factor(invSimp.genus.stat$Treatment_Long, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
chao1.genus.stat$Treatment_Long <- factor(chao1.genus.stat$Treatment_Long, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)

### Shannon
ggplot(shannon.genus.stat, aes(x=Time_TotalDays, y=Shannon, colour=Treatment_Long, group=Treatment_Long)) + 
  scale_color_manual(values = c("darkgrey", "#e68fac", "#EF8737", "darkred")) +
  geom_errorbar(aes(ymin=Shannon-ci, ymax=Shannon+ci), width=0.5, size = 1, alpha = 0.5, position=pd) +
  geom_line(position=pd, linewidth = 2) +
  geom_point(position=pd, size=4) +
  xlab("Time (Days)") +
  ggtitle("Shannon Diversity - Genus") +
  theme_bw() +  theme(legend.position = c(0.80, 0.1)) + labs(color = "Treatment") +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 12))


## Plot boxplot ====

## facet by time

read.csv(here::here("alphaDiv.genus.csv")) -> alphaDiv.genus

alphaDiv.genus$Time_TotalDays <- factor(alphaDiv.genus$Time_TotalDays, c("0", "2", "4", "9", "14", "19", "24", "29", "34"), ordered = TRUE)

alphaDiv.genus$Pretreatment[alphaDiv.genus$Pretreatment == "Antibiotics + Temperature"] <- "Antibiotics + Temp"
alphaDiv.genus$Pretreatment <- factor(alphaDiv.genus$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temp"), ordered = TRUE)

plot_shan <- ggplot(alphaDiv.genus, aes(x = Pretreatment, y = Shannon, color = Pretreatment)) +
  geom_boxplot(lwd = 1.1, outlier.color = "NA")
plot_shan <- plot_shan + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Pretreatment, color = Pretreatment), size = 3,
                                      position = position_dodge(width = 0.9))
plot_shan <- plot_shan + geom_point(aes(color = Pretreatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 80, hjust = 1), axis.title = element_text(face = "bold", size = 14), title = element_text(face = "bold"))
plot_shan <- plot_shan + xlab("Time (Days)") + ylab("Shannon Diversity") 
plot_shan <- plot_shan + guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "none") +
  facet_wrap(~Time_TotalDays) + theme(strip.text = element_text(face = "bold", size = 14)) + scale_color_manual(values = c("#666666", "#61864c", "#c48f10", "#be5067"))


plot_simp <- ggplot(alphaDiv.genus, aes(x = Pretreatment, y = InvSimpson, color = Pretreatment)) +
  geom_boxplot(lwd = 1.1, outlier.color = "NA")
plot_simp <- plot_simp + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Pretreatment, color = Pretreatment), size = 3,
                                      position = position_dodge(width = 0.9))
plot_simp <- plot_simp + geom_point(aes(color = Pretreatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 80, hjust = 1), axis.title = element_text(face = "bold", size = 14), title = element_text(face = "bold"))
plot_simp <- plot_simp + xlab("Time (Days)") + ylab("Inverse Simpson Diversity") 
plot_simp <- plot_simp + guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "none") +
  facet_wrap(~Time_TotalDays) + theme(strip.text = element_text(face = "bold", size = 14)) + scale_color_manual(values = c("#666666", "#61864c", "#c48f10", "#be5067"))

plot_obs <- ggplot(alphaDiv.genus, aes(x = Pretreatment, y = Observed, color = Pretreatment)) +
  geom_boxplot(lwd = 1.1, outlier.color = "NA")
plot_obs <- plot_obs + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Pretreatment, color = Pretreatment), size = 3,
                                      position = position_dodge(width = 0.9))
plot_obs <- plot_obs + geom_point(aes(color = Pretreatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 80, hjust = 1), axis.title = element_text(face = "bold", size = 14), title = element_text(face = "bold"))
plot_obs <- plot_obs + xlab("Time (Days)") + ylab("Observed Richness") 
plot_obs <- plot_obs + guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "none") +
  facet_wrap(~Time_TotalDays) + theme(strip.text = element_text(face = "bold", size = 14)) + scale_color_manual(values = c("#666666", "#61864c", "#c48f10", "#be5067"))


ggarrange(plot_shan, plot_obs, nrow = 1) -> alphaplot


ggplot2::ggsave(here::here("Output Files/06 - Alpha Diversity - Output/01 - alphaplot_genus_nostats.png"), alphaplot,
                height = 450, width = 700, units = "mm",
                scale = 0.5, dpi = 1000)


## Plot alpha diversity when aquarickettsia removed
readRDS(here::here("Output Files/05 - Relative Abundance - Output/ps.noASV1.rare.rds")) -> ps.noASV1.rare

ps.noASV1.rare.genus <- tax_glom(ps.noASV1.rare, taxrank = "Genus", NArm=FALSE)

# Calculate
estimate_richness(ps.noASV1.rare.genus, measures = c("Observed", "InvSimpson", "Shannon", "Chao1", "Simpson")) -> alphaDiv.genus.noASV1

### Add Samples column
alphaDiv.genus.noASV1$Samples <- rownames(alphaDiv.genus.noASV1)

### Combine alpha metrics with sam_data
merge(as.matrix(ps.noASV1.rare.genus@sam_data), alphaDiv.genus.noASV1, by = "Samples") -> alphaDiv.genus.noASV1

write.csv(alphaDiv.genus.noASV1, here::here("Output Files/06 - Alpha Diversity - Output/alphaDiv.genus.noASV1.csv"))



alphaDiv.genus.noASV1$Pretreatment[alphaDiv.genus.noASV1$Pretreatment == "Antibiotics + Temperature"] <- "Antibiotics + Temp"
alphaDiv.genus.noASV1$Pretreatment <- factor(alphaDiv.genus.noASV1$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temp"), ordered = TRUE)

plot_shan <- ggplot(alphaDiv.genus.noASV1, aes(x = Pretreatment, y = Shannon, color = Pretreatment)) +
  geom_boxplot(lwd = 1.1, outlier.color = "NA")
plot_shan <- plot_shan + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Pretreatment, color = Pretreatment), size = 3,
                                      position = position_dodge(width = 0.9))
plot_shan <- plot_shan + geom_point(aes(color = Pretreatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 80, hjust = 1), axis.title = element_text(face = "bold", size = 14), title = element_text(face = "bold"))
plot_shan <- plot_shan + xlab("Time (Days)") + ylab("Shannon Diversity") 
plot_shan <- plot_shan + guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "none") +
  facet_wrap(~Time_TotalDays) + theme(strip.text = element_text(face = "bold", size = 14)) + scale_color_manual(values = c("#666666", "#61864c", "#c48f10", "#be5067"))


plot_simp <- ggplot(alphaDiv.genus.noASV1, aes(x = Pretreatment, y = Simpson, color = Pretreatment)) +
  geom_boxplot(lwd = 1.1, outlier.color = "NA")
plot_simp <- plot_simp + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Pretreatment, color = Pretreatment), size = 3,
                                      position = position_dodge(width = 0.9))
plot_simp <- plot_simp + geom_point(aes(color = Pretreatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 80, hjust = 1), axis.title = element_text(face = "bold", size = 14), title = element_text(face = "bold"))
plot_simp <- plot_simp + xlab("Time (Days)") + ylab("Simpson Diversity (1-D)") 
plot_simp <- plot_simp + guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "none") +
  facet_wrap(~Time_TotalDays) + theme(strip.text = element_text(face = "bold", size = 14)) + scale_color_manual(values = c("#666666", "#61864c", "#c48f10", "#be5067"))


plot_obs <- ggplot(alphaDiv.genus.noASV1, aes(x = Pretreatment, y = Observed, color = Pretreatment)) +
  geom_boxplot(lwd = 1.1, outlier.color = "NA")
plot_obs <- plot_obs + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Pretreatment, color = Pretreatment), size = 3,
                                    position = position_dodge(width = 0.9))
plot_obs <- plot_obs + geom_point(aes(color = Pretreatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 80, hjust = 1), axis.title = element_text(face = "bold", size = 14), title = element_text(face = "bold"))
plot_obs <- plot_obs + xlab("Time (Days)") + ylab("Observed Richness") 
plot_obs <- plot_obs + guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "none") +
  facet_wrap(~Time_TotalDays) + theme(strip.text = element_text(face = "bold", size = 14)) + scale_color_manual(values = c("#666666", "#61864c", "#c48f10", "#be5067"))

ggarrange(plot_shan, plot_obs, nrow = 1, labels = "AUTO") -> alphaplot

ggplot2::ggsave(here::here("Output Files/06 - Alpha Diversity - Output/01 - alphaplot_genus_nostats_noASV1.png"), alphaplot,
                height = 450, width = 700, units = "mm",
                scale = 0.5, dpi = 1000)



## MicrobiomeStat ====

## Alpha diversity is significantly driven by time in the control group (less than ideal)
## Going to use microbiomestat package which allows us to incorporate variables into the alpha diversity calculations?

devtools::install_github("cafferychen777/MicrobiomeStat")
library(MicrobiomeStat)

## convert phyloseq object to microbiomestat data object
data.obj <- mStat_convert_phyloseq_to_data_obj(ps.rare)
## refactor time
data.obj$meta.dat$Time_TotalDays <- factor(as.factor(data.obj$meta.dat$Time_TotalDays), levels = c("0", "2", "4", "9", "14", "19", "24", "29", "34"))
## calculate alpha diversity; this is based on multiple linear models??
generate_alpha_test_single(
  data.obj = data.obj,
  alpha.obj = NULL, 
  alpha.name = c("shannon", "observed_species"),
  depth = NULL,
  group.var = "Treatment_Long",
  adj.vars = "TankID")

## above didn't work; but trying to get visual on data
summary(mStat_summarize_data_obj(data.obj, "Pretreatment", "Time_TotalDays"))



## Plot Shannon ====
### Genus
plot_shan_genus <- ggplot(alphaDiv.genus, aes(x = as.factor(Time), y = Shannon, color = Treatment_Long)) + 
  geom_boxplot(lwd = 1.1, outlier.color = "NA") + stat_summary(fun = mean, geom = "line", mapping = aes(group = Treatment_Long, color = Treatment_Long),
                                                               linewidth = 1.25, position = position_dodge(width = 0.9))
plot_shan_genus <- plot_shan_genus + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Treatment_Long, color = Treatment_Long), size = 3,
                                              position = position_dodge(width = 0.9))
plot_shan_genus <- plot_shan_genus + geom_point(aes(color = Treatment_Long), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold"))
plot_shan_genus <- plot_shan_genus + xlab("Time (Days)") + ylab("Shannon") 
plot_shan_genus <- plot_shan_genus + guides(color = guide_legend(title = "Treatment"))

## Plot InvSimpson ====

### Genus
plot_invSimp_genus <- ggplot(alphaDiv.genus, aes(x = as.factor(Time), y = InvSimpson, color = Treatment_Long)) + 
  geom_boxplot(lwd = 1.1, outlier.color = "NA") + stat_summary(fun = mean, geom = "line", mapping = aes(group = Treatment_Long, color = Treatment_Long),
                                                               linewidth = 1.25, position = position_dodge(width = 0.9))
plot_invSimp_genus <- plot_invSimp_genus + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Treatment_Long, color = Treatment_Long), size = 3,
                                                  position = position_dodge(width = 0.9))
plot_invSimp_genus <- plot_invSimp_genus + geom_point(aes(color = Treatment_Long), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold"))
plot_invSimp_genus <- plot_invSimp_genus + xlab("Time (Days)") + ylab("Inverse Simpson") 
plot_invSimp_genus <- plot_invSimp_genus + guides(color = guide_legend(title = "Treatment"))

## Plot - Observed
### Genus
plot_obs_genus <- ggplot(alphaDiv.genus, aes(x = as.factor(Time), y = Observed, color = Treatment_Long)) + 
  geom_boxplot(lwd = 1.1, outlier.color = "NA") + stat_summary(fun = mean, geom = "line", mapping = aes(group = Treatment_Long, color = Treatment_Long),
                                                               linewidth = 1.25, position = position_dodge(width = 0.9))
plot_obs_genus <- plot_obs_genus + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Treatment_Long, color = Treatment_Long), size = 3,
                                                        position = position_dodge(width = 0.9))
plot_obs_genus <- plot_obs_genus + geom_point(aes(color = Treatment_Long), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold"))
plot_obs_genus <- plot_obs_genus + xlab("Time (Days)") + ylab("Observed") 
plot_obs_genus <- plot_obs_genus + guides(color = guide_legend(title = "Treatment"))

## Plot - Chao1
plot_chao1_genus <- ggplot(alphaDiv.genus, aes(x = as.factor(Time), y = Chao1, color = Treatment_Long)) + 
  geom_boxplot(lwd = 1.1, outlier.color = "NA") + stat_summary(fun = mean, geom = "line", mapping = aes(group = Treatment_Long, color = Treatment_Long),
                                                               linewidth = 1.25, position = position_dodge(width = 0.9))
plot_chao1_genus <- plot_chao1_genus + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Treatment_Long, color = Treatment_Long), size = 3,
                                                  position = position_dodge(width = 0.9))
plot_chao1_genus <- plot_chao1_genus + geom_point(aes(color = Treatment_Long), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold"))
plot_chao1_genus <- plot_chao1_genus + xlab("Time (Days)") + ylab("Chao1") 
plot_chao1_genus <- plot_chao1_genus + guides(color = guide_legend(title = "Treatment"))








ggplot(shannon.order.stat, aes(x=Time, y=Shannon, colour=Treatment_Long, group=Treatment_Long)) + 
  scale_color_manual(values = c("black", "#5386B3", "#EF8737", "darkred")) +
  geom_errorbar(aes(ymin=Shannon-ci, ymax=Shannon+ci), width=0.5, size = 1, alpha = 0.5, position=pd) +
  geom_line(position=pd, linewidth = 2) +
  geom_point(position=pd, size=4) +
  xlab("Time (Days)") +
  ggtitle("Shannon Diversity - Order") +
  theme_bw() +  theme(legend.position = c(0.80, 0.1)) + labs(color = "Treatment") +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 12))


shan_genus <- ggplot(shannon.genus.stat, aes(x=Time_TotalDays, y=Shannon, colour=Treatment_Long, group=Treatment_Long)) + 
  scale_color_manual(values = c("black", "#5386B3", "#EF8737", "darkred")) +
  geom_errorbar(aes(ymin=(Shannon-ci), ymax=(Shannon+ci)), width=0.5, size = 1, alpha = 0.5, position=pd) +
  geom_line(position=pd, linewidth = 2) +
  geom_point(position=pd, size=4) +
  xlab("Time (Days)") +
  ggtitle("Shannon Diversity - Genus") +
  theme_bw() +  theme(legend.position = c(0.80, 0.1)) + labs(color = "Treatment") +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 12))

ggplot2::ggsave(here::here("Output Files/06 - Alpha Diversity - Output/01 - plot_shan_genus.png"), shan_genus,
                height = 450, width = 700, units = "mm",
                scale = 0.5, dpi = 1000)

obs_genus <- ggplot(obs.genus.stat, aes(x=Time_TotalDays, y=Observed, colour=Treatment_Long, group=Treatment_Long)) + 
  scale_color_manual(values = c("black", "#5386B3", "#EF8737", "darkred")) +
  geom_errorbar(aes(ymin=Observed-ci, ymax=Observed+ci), width=0.5, size = 1, alpha = 0.5, position=pd) +
  geom_line(position=pd, linewidth = 2) +
  geom_point(position=pd, size=4) +
  xlab("Time (Days)") +
  ggtitle("Observed - Genus") +
  theme_bw() +  theme(legend.position = c(0.80, 0.1)) + labs(color = "Treatment") +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 12))

ggplot2::ggsave(here::here("Output Files/06 - Alpha Diversity - Output/02 - plot_obs_genus.png"), obs_genus,
                height = 450, width = 700, units = "mm",
                scale = 0.5, dpi = 1000)

subset(alphaDiv.genus, Treatment_Long == "Blank") -> alphaDiv.genus.blank

# Blank over time
mod_shan_blank <- lmer(Shannon ~ Time + (1|TankID), data = alphaDiv.genus.blank) # boundary (singular) fit: see help('isSingular')
anova(mod_shan) -> anova_shan # Time and interaction of treatment and time is sig (p = 0.001991, p = 0.001571)
broom::tidy(anova_shan) -> anova_shan
anova_shan$Metric <- "Shannon"
# seems like there is a significant effect of diff in shannon diversity over time in blank
# but that may not matter if I'm using within time point comparisons 

# Between group comparisons, T0
alphaDiv.genus$Time <- as.numeric(alphaDiv.genus$Time)
subset(alphaDiv.genus, Time == "0") -> alphaDiv.genus.T0

mod_shan <- lmer(Shannon ~ Treatment_Long + (1|TankID), data = alphaDiv.genus.T0)
anova(mod_shan) -> anova_shan # Time and interaction of treatment and time is sig (p = 0.001991, p = 0.001571)
broom::tidy(anova_shan) -> anova_shan
anova_shan$Metric <- "Shannon"

## Maybe for making comparisons within a time point, I should do just a linear model?
model <- lm(Shannon ~ Treatment, data = alphaDiv.genus.T0)
summary(model) 
# coefficient for diff treats indicates the effect of that treatment on shannon
# intercept indicates estimated shannon when treatment is 0??
# Adjusted R-squared:  0.03622 (not good at all)

head(mtcars)
mod <- lm(mpg ~ wt, data = mtcars)
summary(mod)
# Adj R squared: 0.7446

# What if i just do anova without including the random effect of tank
data("ToothGrowth")
df <- ToothGrowth

# Perform the ANOVA test
aov.shan.T0 <- aov(Shannon ~ Treatment_Long, data = alphaDiv.genus.T0)
summary(aov.shan.T0) # p = 0.356

# T2
alphaDiv.genus$Time_TotalDays <- as.numeric(alphaDiv.genus$Time_TotalDays)
subset(alphaDiv.genus, Time_TotalDays == "2") -> alphaDiv.genus.T2
aov.shan.T2 <- aov(Shannon ~ Treatment_Long, data = alphaDiv.genus.T2)
summary(aov.shan.T2) # p = 0.00623
TukeyHSD(aov.shan.T2) # sig between blank - abx, blank - abxtemp, trending sig between temp (i.e, should be essentially control) and abxtemp

# T4
subset(alphaDiv.genus, Time_TotalDays == "4") -> alphaDiv.genus.T4
aov.shan.T4 <- aov(Shannon ~ Treatment_Long, data = alphaDiv.genus.T4)
summary(aov.shan.T4) # p = 0.0153
TukeyHSD(aov.shan.T4) # sig between Temp-Antibiotic, temp-abxTemp

# T9
subset(alphaDiv.genus, Time_TotalDays == "9") -> alphaDiv.genus.T9
aov.shan.T9 <- aov(Shannon ~ Treatment_Long, data = alphaDiv.genus.T9)
summary(aov.shan.T9) # p = 0.000326
TukeyHSD(aov.shan.T9) # sig between blank-antibiotic, temp-antibiotic, blank-antibiotictemp, temp-antibiotic-temp

# T14
subset(alphaDiv.genus, Time_TotalDays == "14") -> alphaDiv.genus.T14
aov.shan.T14 <- aov(Shannon ~ Treatment_Long, data = alphaDiv.genus.T14)
summary(aov.shan.T14) # p = 0.00527
TukeyHSD(aov.shan.T14) # sig between blank-antibiotic, temp-antibiotic

# T19
subset(alphaDiv.genus, Time_TotalDays == "19") -> alphaDiv.genus.T19
aov.shan.T19 <- aov(Shannon ~ Treatment_Long, data = alphaDiv.genus.T19)
summary(aov.shan.T19) # p = 0.0581

# T24
subset(alphaDiv.genus, Time_TotalDays == "24") -> alphaDiv.genus.T24
aov.shan.T24 <- aov(Shannon ~ Treatment_Long, data = alphaDiv.genus.T24)
summary(aov.shan.T24) # p = 0.0193
TukeyHSD(aov.shan.T24) # sig between Blank-antibiotic, trending sig between blank-abxtemp

# T29
subset(alphaDiv.genus, Time_TotalDays == "29") -> alphaDiv.genus.T29
aov.shan.T29 <- aov(Shannon ~ Treatment_Long, data = alphaDiv.genus.T29)
summary(aov.shan.T29) # p = 0.00372
TukeyHSD(aov.shan.T29) # sig between blank-antibiotic, blank-abxtemp, blank-temp

# T34
subset(alphaDiv.genus, Time_TotalDays == "34") -> alphaDiv.genus.T34
aov.shan.T34 <- aov(Shannon ~ Treatment_Long, data = alphaDiv.genus.T34)
summary(aov.shan.T34) # p = 0.228
TukeyHSD(aov.shan.T34) 

# Shannon
## T34 NS
## T29 sig between blank-abx, blank-abxtemp, blank-temp
## T24 sig between blank-abx, trending sig between blank-abxtemp
## T19 p value = 0.0581
## T14 sig between blank-abx, temp-abx
## T9 sig between blank-abx, temp-abx, blank-abxtemp, temp-abxtemp
## T4 sig between temp-abx, temp-abxtemp
## T2 sig between blank-abx, blank-abxtemp
## T0 p value = 0.356

# InvSimp
## T34 NS
## T29 sig between blank-abx, blank-abxtemp, blank-temp
## T24 sig between blank-abx, trending sig between blank-abxtemp
## T19 p value = 0.0589
## T14 sig between blank-abx, temp-abx, trending sig between temp-abxtemp
## T9 sig between blank-abx, blank-abxtemp, temp-abxtemp
## T4 sig between temp-abx, temp-abxtemp
## T2 sig between blank-abxtemp, trending sig between blank-abx
## T0 p value = 0.295



## Plotting shannon diversity over time, faceted by treatment
alphaDiv.genus$Treatment_Long <- factor(alphaDiv.genus$Treatment_Long, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
alphaDiv.genus$Time_TotalDays <- factor(alphaDiv.genus$Time_TotalDays, c("0", "2", "4", "9", "14", "19", "24", "29", "34"), ordered = TRUE)
## All plotted together
plot_shan_1 <- alphaDiv.genus %>%
  ggplot(aes(x = Time_TotalDays, y = Shannon, color = Treatment_Long)) +
  geom_boxplot(lwd = 1.1, outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  geom_smooth(aes(x = as.numeric(Time_TotalDays), y = Shannon), method = "lm") +
  theme_bw() +
  scale_x_discrete() +
  scale_color_manual(values = c("black", "#672998","#CB4255","darkred" )) +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 12))+ labs(x = "Time (Days)", y = "Shannon Diversity") + 
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.5),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5)) +
  facet_wrap(~Treatment_Long) + theme(legend.position = "none")

## Plotting overall shannon diversity by group (no time)
alphaDiv.genus$Treatment_Long <- factor(alphaDiv.genus$Treatment_Long, c("Antibiotics + Temperature", "Temperature", "Antibiotics", "No Treatment"), ordered = TRUE)

g_ridges <- 
  ggplot(alphaDiv.genus, aes(Shannon, Treatment_Long, color = Treatment_Long, fill = Treatment_Long)) + 
  theme_bw() + scale_color_manual(values = c("darkred","#CB4255", "#672998","darkgrey")) +
  scale_fill_manual(values = c("darkred","#CB4255", "#672998","darkgrey")) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 12)) + scale_y_discrete(expand = c(.01,.01), position = "right") + labs(x = "Shannon Diversity") + 
  theme(axis.title.y = element_blank()) + theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.5)) 

plot_shan_2 <- g_ridges +
  ggridges::geom_density_ridges(
    alpha = .7, lwd = 1.1
  )

x <- ggarrange(plot_shan_1, plot_shan_2)
ggplot2::ggsave(here::here("Output Files/06 - Alpha Diversity - Output/plot_shannon_combo.svg"), x,
                height = 400, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)
