
#remove observations w/ "Unclassified" state (based on z.post.thresh)
dat3<- dat2 %>% 
  filter(z.post.thresh != "Unclassified") %>% 
  mutate(across(c(id,season2), factor))


test.mod <-
  brm(data = dat3,
      family = categorical(link = logit, refcat = "Slow-Turn"),
      bf(z.post.thresh ~ 1,
         nlf(mu2 ~ a2 + b2 * lulc),
         nlf(mu3 ~ a3 + b3 * lulc),
         nlf(mu4 ~ a4 + b4 * lulc),
         a2 + a3 + a4 + b2 + b3 + b4 ~ 1),
      prior = c(prior(normal(0, 1.5), class = b, nlpar = a2),
                prior(normal(0, 1.5), class = b, nlpar = a3),
                prior(normal(0, 1.5), class = b, nlpar = a4),
                prior(normal(0, 1), class = b, nlpar = b2),
                prior(normal(0, 1), class = b, nlpar = b3),
                prior(normal(0, 1), class = b, nlpar = b4)),
      iter = 2000, warmup = 1000, cores = 4, chains = 4,
      seed = 11)


fit1 <-
  brm(data = dat3, 
      family = categorical(link = logit, refcat = 'Slow-Turn'),
      z.post.thresh ~ 0 + Intercept + lulc,
      prior(normal(0, 5), class = b),
      iter = 2000, warmup = 1000, cores = 4, chains = 4,
      seed = 22)

print(fit1)


fit1_loo <- add_criterion(fit1, "loo")
loo(fit1_loo)


#make predictions
nd <- tibble(Forest = seq(from = 0, to = 1, length.out = 60),
             Closed_Savanna = seq(from = 0, to = 1, length.out = 60),
             Open_Savanna = seq(from = 0, to = 1, length.out = 60),
             Floodable = seq(from = 0, to = 1, length.out = 60))
f <- fitted(fit1, newdata = nd)


#Plot probabilities per state for each land class

# wrangle
rbind(f[, , 1],
      f[, , 2],
      f[, , 3],
      f[, , 4]) %>% 
  data.frame() %>% 
  bind_cols(nd %>% expand(z.post.thresh = 1:4, Forest)) %>% 
  data.frame() %>% 
  mutate(Closed_Savanna = rep(nd$Closed_Savanna, 4), Open_Savanna = rep(nd$Open_Savanna, 4),
         Floodable = rep(nd$Floodable, 4)) %>% 
  mutate_at("z.post.thresh", factor) %>% 
  # mutate(career = str_c("career: ", career)) %>% 
  
  # plot
  ggplot(aes(x = Floodable, y = Estimate,
             ymin = Q2.5, ymax = Q97.5,
             fill = z.post.thresh, color = z.post.thresh)) +
  geom_ribbon(alpha = 2/3, size = 0) +
  geom_line(size = 3/4) +
  scale_fill_manual(values = wes_palette("Moonrise2")[c(4, 3, 2, 1)]) +
  scale_color_manual(values = wes_palette("Moonrise2")[c(4, 3, 2, 1)]) +
  scale_x_continuous(breaks = 0:2 / 2) +
  scale_y_continuous("probability", limits = c(0, 1),
                     breaks = 0:3 / 3, labels = c("0", ".33", ".67", "1")) +
  theme(axis.text.y = element_text(hjust = 0),
        legend.position = "none") +
  facet_wrap(~ z.post.thresh)





# define the model
code_lulc <- "
data{
  int N; // number of observations
  int K; // number of possible states 
  int state[N]; // outcome
  vector[K] lulc;
}
parameters{
  vector[K - 1] a; // intercepts
  real<lower=0> b; // association of lulc with choice
}
model{
  vector[K] p;
  vector[K] s;
  a ~ normal(0, 1);
  b ~ normal(0, 0.5);
  s[1] = 0; // pivot
  s[2] = a[1] + b * lulc[1]; 
  s[3] = a[2] + b * lulc[2];
  s[4] = a[3] + b * lulc[3];
  p = softmax(s);
  state ~ categorical(p);
} 
"

# wrangle the data
dat_list <- 
  list(N = nrow(dat3), 
       K = 4, 
       state = dat3$z.post.thresh, 
       lulc = dat3$lulc)

# fit the model
new.mod <- 
  stan(data = dat_list,
       model_code = code_lulc,
       chains = 4)
