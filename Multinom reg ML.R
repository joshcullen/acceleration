 library(nnet)
 
 dry.ind<- which(dat$season2 == "Dry")
 
 # Proportions LULC w/in 30 m buffer (by season)
 
 #Flood
 flood.extr.lulc<- extract(flood.lulc, dat[,c('x','y')])  #extract lulc from single cell
 dry.extr.lulc<- extract(dry.lulc, dat[,c('x','y')])
 
 # extr.lulc2<- purrr::map(flood.extr.lulc, ~{prop.table(table(.))}) %>%  #convert to proportions lulc
 #   purrr::map(., ~{
 #     mat<- matrix(0, 1, 4)
 #     mat[1, as.numeric(names(.x))]<- .x
 #     
 #     df<- data.frame(mat)
 #     df
 #   }) %>% 
 #   bind_rows()
 # dry.lulc2<- purrr::map(dry.extr.lulc, ~{prop.table(table(.))}) %>%  #convert to proportions lulc
 #   purrr::map(., ~{
 #     mat<- matrix(0, 1, 4)
 #     mat[1, as.numeric(names(.x))]<- .x
 #     
 #     df<- data.frame(mat)
 #     df
 #   }) %>% 
 #   bind_rows()
 
 extr.lulc2<- flood.extr.lulc
 extr.lulc2[dry.ind]<- dry.extr.lulc[dry.ind]
 extr.lulc3<- factor(extr.lulc2)
 levels(extr.lulc3)<- c("Forest", "Closed_Savanna", "Open_Savanna", "Floodable")
 

 dat2$lulc<- extr.lulc3
 dat3<- dat2 %>% 
   filter(z.post.thresh != "Unclassified")
 boop<- dat3
 boop$z.post.thresh<- factor(as.character(boop$z.post.thresh), levels = unique(boop$z.post.thresh))
 boop$z.post.thresh2 <- relevel(boop$z.post.thresh, ref = "Slow-Unif")
 test <- multinom(z.post.thresh2 ~ lulc, data = boop) 

 summary(test) 

 ## extract the coefficients from the model and exponentiate
 exp(coef(test)) 

 head(pp <- fitted(test)) 
 
 
 
 d.lulc <- data.frame(lulc = c("Forest", "Closed_Savanna", "Open_Savanna", "Floodable"))
 predict(test, newdata = d.lulc, "probs")
 