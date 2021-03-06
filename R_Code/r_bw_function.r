# Code for the r.bw function.
#
# This function originates from the "asbio" pacakge and full credit should be given to the respected authors of said package.
# This function has been included within this code due to some issues encountered when installing and loading the asbio package
#
# https://cran.r-project.org/web/packages/asbio/index.html
# https://cran.r-project.org/web/packages/asbio/asbio.pdf

r.bw<-function(X,Y=NULL) {
  
  U.i<-(X-median(X))/(9*qnorm(.75)*mad(X))
  a.i<-ifelse(U.i<=-1|U.i>=1,0,1)
  n<-nrow(as.matrix(X))
  nx<-sqrt(n)*sqrt(sum((a.i*((X-median(X))^2))*((1-U.i^2)^4)))
  dx<-abs(sum(a.i*(1-U.i^2)*(1-5*U.i^2)))
  S.xx<-(nx/dx)^2
  if(!is.null(Y)){
    V.i<-(Y-median(Y))/(9*qnorm(.75)*mad(Y))
    b.i<-ifelse(V.i<=-1|V.i>=1,0,1)
    ny<-sqrt(n)*sqrt(sum((b.i*((Y-median(Y))^2))*((1-V.i^2)^4)))
    dy<-abs(sum(b.i*(1-V.i^2)*(1-5*V.i^2)))
    S.yy<-(ny/dy)^2
    S.xy<-n*sum((a.i*(X-median(X)))*((1-U.i^2)^2)*(b.i*(Y-median(Y)))*((1-V.i^2)^2))/((sum((a.i*(1-U.i^2))*(1-5*U.i^2)))*(sum((b.i*(1-V.i^2))*(1-5*V.i^2))))
    R.xy<-S.xy/(sqrt(S.xx*S.yy))}
  if(is.null(Y))res<-data.frame(S.xx=S.xx)
  if(!is.null(Y))res<-data.frame(s.xx=S.xx,s.yy=S.yy,s.xy=S.xy,r.xy=R.xy)
  res
}
