
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
