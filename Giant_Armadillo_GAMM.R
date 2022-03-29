
## Run and viz GAMMs for giant armadillo behavioral states and land cover classes ##


library(tidyverse)
library(viridis)
library(mgcv)


dat2<- read.csv("Input Data for GAMM.csv")




#######################################################
### Analyze relationships between states and covars ###
#######################################################


#remove observations w/ "Unclassified" state (based on z.post.thresh)
dat3<- dat2 %>% 
  filter(z.post.thresh != "Unclassified") %>%  #remove Unclassified obs
  mutate(across(id, factor))  #convert ID to factor


state.pal<- viridis(n=4, option = 'inferno')



## Compare VE vs all others
dat.ve<- dat3
dat.ve$state<- ifelse(dat.ve$z.post.thresh == "VE", 1, 0)
table(dat.ve$id, dat.ve$state)  #check n per combo
ve.mod<- gam(state ~ s(Forest, k = 3, bs = "cs") + #s(Closed_Savanna, k = 5, bs = "cs") +
               s(Open_Savanna, k = 3, bs = "cs") + s(Floodable, k = 3, bs = "cs") +
               s(id, bs = "re"),
             data = dat.ve, family = binomial, method = "REML")
summary(ve.mod)
plot.gam(ve.mod, pages = 1, shade = TRUE)




## define theme for GAMM ggplots
theme_gam <- function(){ 
  font <- "Arial"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 18),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 14),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10))
      
      
    )
}


## Manually create partial effects plots
x.seq<- seq(0, 1, len = 100)

new.dat<- expand.grid(Forest = x.seq, Open_Savanna = 0, Floodable = 0, id = "Blanca")
new.dat<- rbind(new.dat,
                expand.grid(Forest = 0, Open_Savanna = x.seq, Floodable = 0, id = "Blanca"))
new.dat<- rbind(new.dat,
                expand.grid(Forest = 0, Open_Savanna = 0, Floodable = x.seq, id = "Blanca"))


pred.ve.mod<- predict(ve.mod, newdata = new.dat, type = "link", se.fit = TRUE,
                      exclude = "s(id)")


pred.ve.df<- cbind(pred.ve.mod$fit, pred.ve.mod$se.fit) %>% 
  data.frame() %>% 
  rename(fit = X1, se = X2) %>% 
  mutate(lower = plogis(fit - 1.96*se),
         upper = plogis(fit + 1.96*se),
         x = rep(x.seq, 3),
         LULC = factor(rep(c("Forest", "Open Savanna", "Floodable"), each = 100),
                       levels = c("Forest", "Open Savanna", "Floodable"))
  ) %>% 
  mutate(fit = plogis(fit),
         mod = "Pr(VE)")






## Compare Local Search vs all others
dat.ls<- dat3
dat.ls$state<- ifelse(dat.ls$z.post.thresh == "Local Search", 1, 0)
table(dat.ls$id, dat.ls$state)  #check n per combo
ls.mod<- gam(state ~ s(Forest, k = 3, bs = "cs") + #s(Closed_Savanna, k = 5, bs = "cs") +
               s(Open_Savanna, k = 3, bs = "cs") + s(Floodable, k = 3, bs = "cs") +
               s(id, bs = "re"),
             data = dat.ls, family = binomial, method = "REML")
summary(ls.mod)
plot.gam(ls.mod, pages = 1, shade = TRUE)


## Manually create partial effects plots
pred.ls.mod<- predict(ls.mod, newdata = new.dat, type = "link", se.fit = TRUE,
                      exclude = "s(id)")


pred.ls.df<- cbind(pred.ls.mod$fit, pred.ls.mod$se.fit) %>% 
  data.frame() %>% 
  rename(fit = X1, se = X2) %>% 
  mutate(lower = plogis(fit - 1.96*se),
         upper = plogis(fit + 1.96*se),
         x = rep(x.seq, 3),
         LULC = factor(rep(c("Forest", "Open Savanna", "Floodable"), each = 100),
                       levels = c("Forest", "Open Savanna", "Floodable"))
  ) %>% 
  mutate(fit = plogis(fit),
         mod = "Pr(Local Search)")







## Compare Exploratory vs all others
dat.exp<- dat3
dat.exp$state<- ifelse(dat.exp$z.post.thresh == "Exploratory", 1, 0)
table(dat.exp$id, dat.exp$state)  #check n per combo
exp.mod<- gam(state ~ s(Forest, k = 3, bs = "cs") + #s(Closed_Savanna, k = 5, bs = "cs") +
               s(Open_Savanna, k = 3, bs = "cs") + s(Floodable, k = 3, bs = "cs") +
               s(id, bs = "re"),
             data = dat.exp, family = binomial, method = "REML")
summary(exp.mod)
plot.gam(exp.mod, pages = 1, shade = TRUE)


## Manually create partial effects plots
pred.exp.mod<- predict(exp.mod, newdata = new.dat, type = "link", se.fit = TRUE,
                       exclude = "s(id)")


pred.exp.df<- cbind(pred.exp.mod$fit, pred.exp.mod$se.fit) %>% 
  data.frame() %>% 
  rename(fit = X1, se = X2) %>% 
  mutate(lower = plogis(fit - 1.96*se),
         upper = plogis(fit + 1.96*se),
         x = rep(x.seq, 3),
         LULC = factor(rep(c("Forest", "Open Savanna", "Floodable"), each = 100),
                       levels = c("Forest", "Open Savanna", "Floodable"))
  ) %>% 
  mutate(fit = plogis(fit),
         mod = "Pr(Exploratory)")






## Compare Transit vs all others
dat.t<- dat3
dat.t$state<- ifelse(dat.t$z.post.thresh == "Transit", 1, 0)
table(dat.t$id, dat.t$state)  #check n per combo
t.mod<- gam(state ~ s(Forest, k = 3, bs = "cs") + #s(Closed_Savanna, k = 5, bs = "cs") +
               s(Open_Savanna, k = 3, bs = "cs") + s(Floodable, k = 3, bs = "cs") +
               s(id, bs = "re"),
             data = dat.t, family = binomial, method = "REML")
summary(t.mod)
plot.gam(t.mod, pages = 1, shade = TRUE)


## Manually create partial effects plots
pred.t.mod<- predict(t.mod, newdata = new.dat, type = "link", se.fit = TRUE,
                     exclude = "s(id)")


pred.t.df<- cbind(pred.t.mod$fit, pred.t.mod$se.fit) %>% 
  data.frame() %>% 
  rename(fit = X1, se = X2) %>% 
  mutate(lower = plogis(fit - 1.96*se),
         upper = plogis(fit + 1.96*se),
         x = rep(x.seq, 3),
         LULC = factor(rep(c("Forest", "Open Savanna", "Floodable"), each = 100),
                       levels = c("Forest", "Open Savanna", "Floodable"))
  ) %>% 
  mutate(fit = plogis(fit),
         mod = "Pr(Transit)")




## Combine all results for single facet plot
pred.df<- rbind(pred.ve.df, pred.ls.df, pred.exp.df, pred.t.df)
pred.df$mod<- factor(pred.df$mod, levels = c("Pr(VE)", "Pr(Local Search)", "Pr(Exploratory)",
                                             "Pr(Transit)"))

# Remove predictions for Forest at proportions > 0.75 due to limited data
# apply(dat3[,c('Forest','Open_Savanna','Floodable')],2,quantile,0.95)
pred.df2<- pred.df %>% 
  filter(!(LULC == "Forest" & x > 0.75))

ggplot(pred.df2, aes(x, fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = mod), alpha = 0.5) +
  geom_line() +
  scale_fill_manual(values = state.pal, guide = "none") +
  theme_gam() +
  theme(strip.text = element_text(size = 18),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title = element_blank()) +
  facet_grid(mod~LULC, scales = "free_y", switch = "both")
