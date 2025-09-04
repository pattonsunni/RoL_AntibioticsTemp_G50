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
library(EnvStats)

## Load data ====
readRDS(here::here("Output Files/05 - Relative Abundance - Output/ps.rare.rds")) -> ps.rare

## Calculate alpha diversity ====
# Agglomerate to genus level
ps.rare.genus <- tax_glom(ps.rare, taxrank = "Genus", NArm=FALSE)

# Calculate
estimate_richness(ps.rare.genus, measures = c("Observed", "Simpson", "Shannon")) -> alphaDiv.genus

### Add Samples column
alphaDiv.genus$Samples <- rownames(alphaDiv.genus)

### Combine alpha metrics with sam_data
merge(as.matrix(ps.rare@sam_data), alphaDiv.genus, by = "Samples") -> alphaDiv.genus

### Add phase column
alphaDiv.genus$Phase <- paste0(alphaDiv.genus$Time_TotalDays)
alphaDiv.genus$Phase[alphaDiv.genus$Phase == "0"] <- "Antibiotic challenge" 
alphaDiv.genus$Phase[alphaDiv.genus$Phase == "2"] <- "Antibiotic challenge"  
alphaDiv.genus$Phase[alphaDiv.genus$Phase == "4"] <- "Antibiotic challenge"  
alphaDiv.genus$Phase[alphaDiv.genus$Phase == "9"] <- "Temperature ramp"
alphaDiv.genus$Phase[alphaDiv.genus$Phase == "14"] <- "Temperature ramp"  
alphaDiv.genus$Phase[alphaDiv.genus$Phase == "19"] <- "Temperature challenge"  
alphaDiv.genus$Phase[alphaDiv.genus$Phase == "24"] <- "Temperature challenge"  
alphaDiv.genus$Phase[alphaDiv.genus$Phase == "29"] <- "Temperature challenge"  
alphaDiv.genus$Phase[alphaDiv.genus$Phase == "34"] <- "Temperature challenge" 

write.csv(alphaDiv.genus, here::here("alphaDiv.genus.csv"))

#read.csv(here::here("Output files/06 - Alpha Diversity - Output/alphaDiv.genus.csv")) -> alphaDiv.genus

## Look into data distributions ====
shapiro.test(alphaDiv.genus$Shannon) # p = 2.2x10^-9 
shapiro.test(alphaDiv.genus$Simpson) # p = 1.3x10^-8
shapiro.test(alphaDiv.genus$Observed) # p = 1.4x10^-6

# See what distribution is closest to 
descdist(alphaDiv.genus$Shannon, boot=500)   ## Shannon seems to have beta distribution (but isn't because values aren't from 0-1)

# Create plots to visualize different distribution fits
shan_wb <- fitdist(alphaDiv.genus$Shannon, "weibull")
shan_gam <- fitdist(alphaDiv.genus$Shannon, "gamma")
shan_exp <- fitdist(alphaDiv.genus$Shannon, "exp")
shan_lognorm <- fitdist(alphaDiv.genus$Shannon, "lnorm")


par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "gamma", "expo", "lognorm")
denscomp(list(shan_wb, shan_gam, shan_exp, shan_lognorm), legendtext = plot.legend)
qqcomp(list(shan_wb, shan_gam, shan_exp, shan_lognorm), legendtext = plot.legend)
cdfcomp(list(shan_wb, shan_gam, shan_exp, shan_lognorm), legendtext = plot.legend)
ppcomp(list(shan_wb, shan_gam, shan_exp, shan_lognorm), legendtext = plot.legend)
## likely not lognormal or exponential

## Goodness of fit tests to compare models
gofstat(list(shan_wb, shan_gam, shan_lognorm))

gof = gofTest(alphaDiv.genus$Shannon,distribution = "lnorm", test = "sw")
print(gof) 
gof = gofTest(alphaDiv.genus$Shannon,distribution = "weibull", test = "sw")
print(gof) 
gof = gofTest(alphaDiv.genus$Shannon,distribution = "gamma", test = "sw")
print(gof) 
## Looks to be gamma or weibull

 
## Generalized linear model ====
alphaDiv.genus$Time_TotalDays <- factor(alphaDiv.genus$Time_TotalDays)
mod <- glmer(Shannon ~ Pretreatment + Phase + Pretreatment*Phase + (1|TankID) + (1|Time_TotalDays),
             family = Gamma(link = "log"),
             data = alphaDiv.genus) 
shapiro.test(residuals(mod)) # residuals normal
emmeans(mod, pairwise ~ Pretreatment | Phase, adjust = "tukey") 
### Antibiotic challenge
#### Antibiotics vs No treatment sig

### Temperature ramp
#### Antibiotics vs No treatment sig
#### Antibiotics vs Temperature sig
#### Antibiotics + Temperature vs No treatment sig
#### Antibiotics + Temperature vs Temperature sig

### Temperature challenge
#### Antibiotics vs No treatment sig
#### Antibiotics + Temperature vs No treatment sig
#### Antibiotics + Temperature vs Temperature sig
#### No treatment vs Temperature sig

## Figure 3 ====
## alpha diversity boxplot over time  
alphaDiv.genus$Treatment_Long <- factor(alphaDiv.genus$Treatment_Long, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)
alphaDiv.genus$Time_TotalDays <- factor(alphaDiv.genus$Time_TotalDays, c("0", "2", "4", "9", "14", "19", "24", "29", "34"), ordered = TRUE)

boxplot <- alphaDiv.genus %>%
  ggplot(aes(x = Time_TotalDays, y = Shannon, color = Treatment_Long)) +
  geom_boxplot(lwd = 1.1, outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  geom_smooth(aes(x = as.numeric(Time_TotalDays), y = Shannon), method = "lm") +
  theme_bw() +
  scale_x_discrete() +
  scale_color_manual(values = c("black", "#672998","#EF8737","darkred" )) +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 12))+ labs(x = "Time (Days)", y = "Shannon Diversity") + 
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.5),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5)) +
  facet_wrap(~Treatment_Long) + theme(legend.position = "none")

## ridgeline plot by phase
### T0-T4 (antibiotic challenge)
alphaDiv.genus$Pretreatment <- factor(alphaDiv.genus$Pretreatment, c("No Treatment", "Antibiotics"), ordered = TRUE)

plot_shan_abx <- alphaDiv.genus %>%
  filter(Phase == "Antibiotic challenge") %>%
  ggplot(aes(x = Pretreatment, y = Shannon, color = Pretreatment)) +
  geom_boxplot(lwd = 1.1, outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  theme_bw() +
  scale_x_discrete() +
  scale_color_manual(values = c("black", "#672998")) +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 12)) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Shannon Diversity", title = "Antibiotic challenge", color = "Treatment") + 
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.5),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5))  + 
  ylim(0,1.5)

### T9-T14 (temp ramp and antibiotic recovery)
alphaDiv.genus$Pretreatment <- factor(alphaDiv.genus$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)

plot_shan_ramp <- alphaDiv.genus %>%
  filter(Phase == "Temperature ramp") %>%
  ggplot(aes(x = Pretreatment, y = Shannon, color = Pretreatment)) +
  geom_boxplot(lwd = 1.1, outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  theme_bw() +
  scale_x_discrete() +
  scale_color_manual(values = c("black", "#672998", "#EF8737", "darkred")) +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 12)) + 
  labs(title = "Temperature ramp", color = "Treatment", y = "") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.5),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5))+ ylim(0, 1.5)

### T19-T34 Full temperature
alphaDiv.genus$Pretreatment <- factor(alphaDiv.genus$Pretreatment, c("No Treatment", "Antibiotics", "Temperature", "Antibiotics + Temperature"), ordered = TRUE)

plot_shan_temp <- alphaDiv.genus %>%
  filter(Phase == "Temperature challenge") %>%
  ggplot(aes(x = Pretreatment, y = Shannon, color = Pretreatment)) +
  geom_boxplot(lwd = 1.1, outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  theme_bw() +
  scale_x_discrete() +
  scale_color_manual(values = c("black", "#672998", "#EF8737", "darkred")) +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 12)) + 
  labs(title = "Temperature challenge", color = "Treatment", y = "") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.5),
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5)) + ylim(0, 1.5)

y <- ggarrange(plot_shan_abx,plot_shan_ramp, plot_shan_temp, ncol = 3, legend = "none")

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
  scale_color_manual(values = c("black", "#672998","#EF8737","darkred" )) +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.5),
        legend.position = "none",
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5)) + 
  labs(x = "Time (Days)", y = "Shannon Diversity") + 
  facet_wrap(~Treatment_Long) 

## Plotting overall shannon diversity by group (no time)
alphaDiv.genus$Treatment_Long <- factor(alphaDiv.genus$Treatment_Long, c("Antibiotics + Temperature", "Temperature", "Antibiotics", "No Treatment"), ordered = TRUE)

g_ridges <- 
  ggplot(alphaDiv.genus, aes(Simpson, Treatment_Long, color = Treatment_Long, fill = Treatment_Long)) + 
  theme_bw() + scale_color_manual(values = c("darkred","#EF8737", "#672998","darkgrey")) +
  scale_fill_manual(values = c("darkred","#EF8737", "#672998","darkgrey")) +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.5),
        legend.position = "none",
        strip.background = element_rect(
          color = "black", # Border color
          fill = "white", # Background fill color of the strip
          size = 1.5)) + 
  scale_y_discrete(expand = c(.01,.01), position = "right") + labs(x = "Shannon Diversity") 

plot_shan_2 <- g_ridges +
  ggridges::geom_density_ridges(
    alpha = .7, lwd = 1.1
  )

x <- ggarrange(plot_shan_1, plot_shan_2)
xx <- ggarrange(x, y, widths = 1:2)
ggplot2::ggsave(here::here("Output Files/06 - Alpha Diversity - Output/plot_shannon_combo.svg"), xx,
                height = 400, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)

























## Simpson diversity ====
descdist(alphaDiv.genus$Simpson, boot=500)   ## Shannon seems to have beta distribution (but isn't because values aren't from 0-1)

# Create plots to visualize different distribution fits
simp_wb <- fitdist(alphaDiv.genus$Simpson, "weibull")
simp_gam <- fitdist(alphaDiv.genus$Simpson, "gamma")
simp_exp <- fitdist(alphaDiv.genus$Simpson, "exp")
simp_lognorm <- fitdist(alphaDiv.genus$Simpson, "lnorm")


par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "gamma", "expo", "lognorm")
denscomp(list(simp_wb, simp_gam, simp_exp, simp_lognorm), legendtext = plot.legend)
qqcomp(list(simp_wb, simp_gam, simp_exp, simp_lognorm), legendtext = plot.legend)
cdfcomp(list(simp_wb, simp_gam, simp_exp, simp_lognorm), legendtext = plot.legend)
ppcomp(list(simp_wb, simp_gam, simp_exp, simp_lognorm), legendtext = plot.legend)

## GLM
alphaDiv.genus$Time_TotalDays <- factor(alphaDiv.genus$Time_TotalDays)
## running the same model as with shannon is giving the boundary is.singular warning
## since simpson diversity is bound between 0-1, it's likely beta dist

mod <- glmmTMB(Simpson ~ Pretreatment + Phase + Pretreatment*Phase + (1|TankID) + (1|Time_TotalDays),
        family = beta_family(link = "logit"),
        data = alphaDiv.genus)
shapiro.test(residuals(mod)) # residuals not normal
emmeans(mod, pairwise ~ Pretreatment | Phase, adjust = "tukey") 