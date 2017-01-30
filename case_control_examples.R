# ------------------------------------------------------------------------------
# Title: Simulating and Visualizing Case Control Studies
# Author: Ryan Gan
# Date: Jan 28th 2017
# ------------------------------------------------------------------------------

# loading libraries -----
library(tidyverse)
library(broom)
library(gganimate)

# simulation population data ----
pop_data <- data_frame(x = rep(seq(from = 1, to = 100, by =1), 100)) %>% 
  arrange(x) %>%  # sort x by ascending values
  cbind(y = rep(seq(from = 1, to = 100, by =1), 100)) %>% 
  # create randomly assigned proportion of disease
  mutate(expose = rbinom(10000, size = 1, prob = 0.5),
         exp_yn = ifelse(expose == 1, "Yes", "No"))

# finding the formula of the baseline disease probability I want
#1/(1+exp(-(-3.6 + 0.69))) 
# relationship between dis and exp
logit_form = -3.25 + 0.69*pop_data$expose # linear combination with a bias
pr = 1/(1+exp(-logit_form))   

pop_data <- pop_data %>% 
  cbind(disease = rbinom(10000, size = 1, prob = pr)) %>% 
  mutate(dis_yn = ifelse(disease == 1, "Yes", "No"),
         # make 4 category exposure/disease variable
         exp_dis = as.factor(
                   ifelse(expose == 0 & disease == 0, "Exp = N, Dis = N",
                   ifelse(expose == 0 & disease == 1, "Exp = N, Dis = Y",
                   ifelse(expose == 1 & disease == 0, "Exp = Y, Dis = N",
                   ifelse(expose == 1 & disease == 1, "Exp = Y, Dis = Y", 
                          NA))))))

summary(pop_data$exp_dis)

# plot of population at time 0 ----
ggplot(pop_data, aes(x=x, y=y)) +
  geom_point(color = "grey") +
  # title at time 0
  ggtitle("Population - Time 0") +
  theme(panel.background = element_rect("white"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

# plot of exposure at time 0 ----
ggplot(pop_data, aes(x=x, y=y, color = exp_yn)) +
  geom_point() +
  scale_color_manual(guide = guide_legend("Exposure"), 
                     values = c("grey", "blue")) +
  # title at time 0
  ggtitle("Population Exposure - Time 0") +
  theme(panel.background = element_rect("white"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

# proportion of exposed in population
mean(pop_data$expose)
summary(as.factor(pop_data$exp_dis))
# plot of population exposure/disease -----
ggplot(pop_data, aes(x=x, y=y, color = exp_dis, shape = exp_dis)) +
  geom_point() +
  scale_color_manual(guide = guide_legend("Exposure/Disease"), 
                     values = c("grey", "grey", "blue", "blue")) +
  # custom shape
  scale_shape_manual(guide = guide_legend("Exposure/Disease"), 
                     values = c(20, 4, 20, 4)) + 
  # title at time 0
  ggtitle("Population Exposure/Disease - Time 0") +
  theme(panel.background = element_rect("white"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

# cross tabs (2x2 table)
xtabs(~ expose + disease, pop_data)
# hand calc
rel_risk <- (330/(330+4645))/(176/(176+4849))
rel_risk
# calculate relative risk
rel_mod <- tidy(glm(disease ~ expose, data = pop_data, 
                    family = "poisson"(link="log")))

# OR hand calc
odds_ratio <- (375/4664)/(177/4784)
odds_ratio

# test logistic model
or_mod <- tidy(glm(disease ~ expose, 
               data = pop_data, family = "binomial"))

# case control study from source population ----
# sample 100 cases and controls from the population
case_control_data <- pop_data %>% 
  # create selection variable
  group_by(disease) %>% 
  do(sample_n(.,200))

xtabs(~disease + expose, case_control_data)

# plot case control sample ----
ggplot(case_control_data, aes(x=x, y=y, color = exp_dis, shape = exp_dis)) +
  geom_point() +
  scale_color_manual(guide = guide_legend("Exposure/Disease"), 
                     values = c("grey", "grey", "blue", "blue")) +
  # custom shape
  scale_shape_manual(guide = guide_legend("Exposure/Disease"), 
                     values = c(20, 4, 20, 4)) + 
  # title at time 0
  ggtitle("Case Control Study Sample") +
  theme(panel.background = element_rect("white"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

or_mod <- tidy(glm(disease ~ expose, 
               data = case_control_data,
               family = "binomial"(link="logit")))
or_mod
exp(0.98)

rr_mod <- tidy(glm(disease ~ expose,
                   data = case_control_data,
                   family = "poisson"(link="log")))
rr_mod
# plot of population selected for cohort study ----
ggplot(pop_data, aes(x=x, y=y)) +
  geom_point(color = "grey") +
  # sample box
  geom_rect(aes(xmin = 50, xmax = 100, ymin = 50, ymax = 100), 
            fill = NA, color = "red") +
  ggtitle("Cohort Sample of Population") +
  theme(panel.background = element_rect("white"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

# limit sample df to cohort df by bounding box
cohort_data <- pop_data %>% 
  filter(x >= 51 & x <= 100) %>%  # limit by x
  filter(y >= 51 & y <= 100) # limit by y

# notice exposure and disease proportion of the cohort is roughly equal to 
# sample population ----
summary(cohort_data)

# Risk Ratio and OR calculation of cohorot design
rel_mod <- tidy(glm(disease ~ expose, data = cohort_data, 
                    family = "poisson"(link="log")))
rel_mod
exp(0.56)
