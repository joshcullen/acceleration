
library(brms)

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



## Based on Solomon Kurz code from: https://solomonkurz.netlify.app/post/2021-11-17-conditional-logistic-models-with-brms-rough-draft/

#create custom distr for conditional logistic regression
cond_log_1 <- custom_family(
  name     = "cond_log_1", 
  dpars    = c("mu", "mub", "muc"), 
  links    = "identity", 
  type     = "int",
  vars     = c("n_cat"),
  specials = "categorical"
)

#define custom PMF for predictions
stan_lpmf_1 <- stanvar(block = "functions", 
                       scode = "
real cond_log_1_lpmf(int y, real mu, real mu_b, real mu_c, int n_cat) {
  real p_mu  = inv_logit(mu);
  real p_mub = inv_logit(mu_b);
  real p_muc = inv_logit(mu_c);
  vector[n_cat] prob;
  prob[1] = p_mu;
  prob[2] = p_mub * (1 - p_mu);
  prob[3] = p_muc * (1 - p_mub) * (1 - p_mu);
  prob[4] = (1 - p_mu) * (1 - p_mub) * (1 - p_muc);
  return(categorical_lpmf(y | prob));
}

vector cond_log_1_pred(int y, real mu, real mu_b, real mu_c, int n_cat) {
  real p_mu  = inv_logit(mu);
  real p_mub = inv_logit(mu_b);
  real p_muc = inv_logit(mu_c);
  vector[n_cat] prob;
  prob[1] = p_mu;
  prob[2] = p_mub * (1 - p_mu);
  prob[3] = p_muc * (1 - p_mub) * (1 - p_mu);
  prob[4] = (1 - p_mu) * (1 - p_mub) * (1 - p_muc);
  return(prob);
}
")


#additional info
stanvars <- stanvar(x = 4, name = "n_cat", scode = "  int n_cat;")


## fit the model
fit2 <-
  brm(data = dat3, 
      family = cond_log_1,
      z.post.thresh ~ Forest + I(Forest^2) + Open_Savanna + I(Open_Savanna^2) + 
        Floodable + I(Floodable^2) + (1|id),
      prior = c(prior(normal(0, 5), class = Intercept, dpar = muLocalSearch),
                prior(normal(0, 5), class = Intercept, dpar = muExploratory),
                prior(normal(0, 5), class = Intercept, dpar = muTransit),
                prior(normal(0, 5), class = b, dpar = muLocalSearch),
                prior(normal(0, 5), class = b, dpar = muExploratory),
                prior(normal(0, 5), class = b, dpar = muTransit),
                prior(exponential(1), class = sd, dpar = muLocalSearch),
                prior(exponential(1), class = sd, dpar = muExploratory),
                prior(exponential(1), class = sd, dpar = muTransit)),
      iter = 2000, warmup = 1000, cores = 4, chains = 4,
      seed = 123,
      stanvars = stan_lpmf_1 + stanvars)

#for post-processing w/ brms
expose_functions(fit2, vectorize = TRUE)

#plot results
posterior_epred_cond_log_1 <- function(prep) {
  mu   <- brms:::get_dpar(prep, "muLocalSearch")
  mu_b <- brms:::get_dpar(prep, "muExploratory")
  mu_c <- brms:::get_dpar(prep, "muTransit")
  n_cat <- prep$data$n_cat
  y <- prep$data$Y
  prob <- cond_log_1_pred(y = y, mu = mu, mu_b = mu_b, mu_c = mu_c, n_cat = n_cat)
  dim(prob) <- c(dim(prob)[1], dim(mu))
  prob <- aperm(prob, c(2,3,1))
  dimnames(prob) <- list(
    as.character(seq_len(dim(prob)[1])), 
    NULL, 
    as.character(seq_len(dim(prob)[3]))
  )
  prob
}

ce <- conditional_effects(
  fit2, 
  categorical = T,
  effects = "Forest")

plot(ce, plot = FALSE)[[1]] + 
  scale_fill_viridis_d(option = "F", begin = .15, end = .85) +
  scale_color_viridis_d(option = "F", begin = .15, end = .85)


## posterior predictive check
posterior_predict_cond_log_1 <- function(i, prep, ...) {
  mu   <- brms::get_dpar(prep, "muLocalSearch", i = i)
  mu_b <- brms::get_dpar(prep, "muExploratory", i = i)
  mu_c <- brms::get_dpar(prep, "muTransit", i = i)
  n_cat <- prep$data$n_cat
  y <- prep$data$Y[i]
  prob <- cond_log_1_pred(y, mu, mu_b, mu_c, n_cat)
  # make sure you have the extraDistr package
  extraDistr::rcat(length(mu), t(prob))
}

pp_check(fit2, 
         type = "bars", 
         ndraws = 100, 
         size = 1/2, 
         fatten = 2)



## prior predictive check
fit0 <-
  brm(data = dat3, 
      family = categorical(link = logit, refcat = 'VE'),
      z.post.thresh ~ Forest + I(Forest^2) + Open_Savanna + I(Open_Savanna^2) + 
        Floodable + I(Floodable^2) + (1|id),
      prior = c(prior(normal(0, 5), class = Intercept, dpar = muLocalSearch),
                prior(normal(0, 5), class = Intercept, dpar = muExploratory),
                prior(normal(0, 5), class = Intercept, dpar = muTransit),
                prior(normal(0, 5), class = b, dpar = muLocalSearch),
                prior(normal(0, 5), class = b, dpar = muExploratory),
                prior(normal(0, 5), class = b, dpar = muTransit),
                prior(exponential(1), class = sd, dpar = muLocalSearch),
                prior(exponential(1), class = sd, dpar = muExploratory),
                prior(exponential(1), class = sd, dpar = muTransit)),
      iter = 2000, warmup = 1000, cores = 4, chains = 4,
      seed = 123, sample_prior = "only")

pp_check(fit0, 
         type = "bars", 
         ndraws = 200, 
         size = 1/2, 
         fatten = 2)

## run model
fit1 <-
  brm(data = dat3, 
      family = categorical(link = logit, refcat = 'VE'),
      z.post.thresh ~ Forest + I(Forest^2) + Open_Savanna + I(Open_Savanna^2) + 
        Floodable + I(Floodable^2) + (1|id),
      prior = c(prior(normal(0, 5), class = Intercept, dpar = muLocalSearch),
                prior(normal(0, 5), class = Intercept, dpar = muExploratory),
                prior(normal(0, 5), class = Intercept, dpar = muTransit),
                prior(normal(0, 5), class = b, dpar = muLocalSearch),
                prior(normal(0, 5), class = b, dpar = muExploratory),
                prior(normal(0, 5), class = b, dpar = muTransit),
                prior(exponential(1), class = sd, dpar = muLocalSearch),
                prior(exponential(1), class = sd, dpar = muExploratory),
                prior(exponential(1), class = sd, dpar = muTransit)),
      iter = 2000, warmup = 1000, cores = 4, chains = 4,
      seed = 123)

print(fit1)
plot(fit1)  #traceplots and density plots for each param

ce <- conditional_effects(
  fit1, 
  categorical = TRUE,
  effects = c("Forest", "Open_Savanna", "Floodable"))

ce.df <- ce %>% 
  map({. %>% dplyr::select(effect1__:upper__)}) %>% 
  bind_rows(.id = "lulc") %>% 
  mutate(lulc = gsub(":.*", "", lulc)) %>% 
  mutate(lulc = gsub("Open_Savanna", "Open Savanna", lulc)) %>% 
  mutate(across(lulc, factor, levels = c('Forest', 'Open Savanna', 'Floodable')))

ggplot(ce.df, aes(effect1__, estimate__, group = effect2__)) +
  geom_line(aes(color = effect2__), size = 1) +
  scale_color_viridis_d("", option = "inferno", begin = 0, end = 0.8) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = effect2__), alpha = 0.6) +
  scale_fill_viridis_d("", option = "inferno") +
  theme_bw() +
  ylab("Probability") +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1)) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 20),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.text = element_text(size = 16),
        legend.position = "top",
        panel.spacing.x = unit(1, "lines")) +
  facet_wrap(~ lulc, nrow = 1, strip.position = "bottom")

# ggsave('Figure 5_new.png', width = 12, height = 9, units = "in", dpi = 330)


## posterior predictive check
pp_check(fit1, 
         type = "bars", 
         ndraws = 200, 
         size = 1/2, 
         fatten = 2)


#make predictions
nd <- tibble(Forest = seq(from = 0, to = 1, length.out = 60),
             # Closed_Savanna = seq(from = 0, to = 1, length.out = 60),
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
  bind_cols(nd %>% 
              expand(z.post.thresh = 1:4, Forest)
            ) %>% 
  data.frame() %>% 
  mutate(Open_Savanna = rep(nd$Open_Savanna, 4),
         Floodable = rep(nd$Floodable, 4)) %>% 
  mutate_at("z.post.thresh", factor) %>% 
  # mutate(career = str_c("career: ", career)) %>% 
  
  # plot
  ggplot(aes(x = Forest, y = Estimate,
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
