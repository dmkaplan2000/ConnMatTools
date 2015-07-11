# Function used to fix problem with quantile estimation when you have many exact 0's and
# 1's at beginning and end, respectively or probabilty distribution
rel.conn.approxfun <- function(x,y,xmin=0,xmax=1,...) {
  if (any(x>xmax) || any(x<xmin))
    stop(paste("Values outside of supposed bounds: [",xmin,",",xmax,"]"))
  
  I = which(x > xmin & x < xmax)
  
  I1 = max(1,min(I)-1)
  I2 = min(length(x),max(I)+1)
  
  return(approxfun(x[I1:I2],y[I1:I2],...))
}
