
merge_ac_gps = function(data.AC, data.GPS, tol) {
  #data.AC is a data.frame per ID for Activity.Count
  #data.GPS is a data frame per ID for GPS coords
  #tol is the tolerance (in minutes) for which to include Activity Counts

### Only select AC and GPS values when recorded within 1 min of each other

tmp.AC<- data.AC
tmp.GPS<- data.GPS


## First, create `start.time` and `end.time` columns for each datetime of the primary DF (tmp.GPS) on which to later perform a left_join; this should be done using the time threshold/tolerance for retaining Activity.Counts

tmp.GPS$start.time<- tmp.GPS$date - lubridate::minutes(tol)
tmp.GPS$end.time<- tmp.GPS$date + lubridate::minutes(tol)


## Find Activity Counts w/in 1 min of recorded GPS locs and interpolate

for (i in 1:nrow(tmp.GPS)) {
  # print(i)
  ind<- which(tmp.AC$Acquisition.Time >= tmp.GPS$start.time[i] &
          tmp.AC$Acquisition.Time <= tmp.GPS$end.time[i])
  
  if (length(ind) == 0) {
    tmp.GPS$Activity.Count[i]<- NA
  } else {
  tmp.GPS$Activity.Count[i]<- round(tmp.AC$Activity.Count[ind] * (7/5))
  }
}

tmp.GPS
}

#----------------------------------------

fill_NA = function(data, int, units = "min") {
  # dt must be in minutes until I make updates to allow other units
  # units must be "min" for now
  
  dat.list<- bayesmove::df_to_list(data, "id")
  
  dat.list<- purrr::map(dat.list, ~{
    
    dat<- data.frame(.)
    ind<- which(!is.na(dat$dt) & dat$dt > int)
    ind2<- ind - 1
    
    for (i in 1:length(ind)) {
      if (dat$dt[ind[i]] >= 2*int & dat$dt[ind[i]] < 7*24*60) {
        
        # cond<- dat$dt[ind[i]] %% int == 0  #check if dt is multiple of int
        vec.length<- floor(dat$dt[ind[i]]/int)  #find multiple of int in dt to determine seq length
        seq.dates<- seq(dat$date[ind2[i]], dat$date[ind[i]], by = paste(int, units))[1:vec.length]
        tmp1<- as.data.frame(lapply(dat[ind2[i],], rep, length(seq.dates)))
        tmp1$date<- seq.dates
        tmp1$dt<- int
        
        dat[seq(ind2[i]+vec.length, nrow(dat)+vec.length-1),]<- dat[ind[i]:nrow(dat),]
        dat[ind2[i]:(ind2[i]+vec.length-1),]<- tmp1
        
        ind<- ind + (vec.length - 1)  #update index
        ind2<- ind2 + (vec.length - 1)  #update index
        
      } else {
        dat<- dat
      }
    }
    
    
    #update dt column and remove any duplicate rows
    dat$dt<- as.numeric(difftime(dat$date, lag(dat$date), units = "min"))
    dat<- dplyr::distinct(dat, date, .keep_all = TRUE)
    
    dat
    })
  
  dat.out<- dplyr::bind_rows(dat.list)
  
  dat.out
}



