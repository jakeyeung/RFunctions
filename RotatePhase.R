# Jake Yeung
# RotatePhase.R
#  
# 2017-10-30

RotatePhase <- Vectorize(function(phase, rotate.hr=0){
  # rotate phase, in hours. If negative add 24. If > 24 minus 24.
  if (is.na(phase)) return(phase)
  phase <- phase + rotate.hr
  if (phase < 0) phase <- phase + 24
  if (phase > 24) phase <- phase - 24
  return(phase)
}, "phase")
