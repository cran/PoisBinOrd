gen.PoisBinOrd <-
function(n, n.P, n.B, n.O, lambda.vec=NULL, prop.vec=NULL, prop.list=NULL, final.corr.mat){

   if(missing(n)==TRUE)          stop("n was not specified! \n")
   if(missing(final.corr.mat))   stop("Final correlation matrix was not specified! \n")

   validation.bin(n.B, prop.vec)
   validation.ord(n.O, prop.list)


   if(ncol(final.corr.mat) != (n.P+n.B+n.O)) {
     stop("Dimension of final correlation matrix does not match the number of variables! \n")
   } #if

   myz<-rmvnorm(n, mean=rep(0,(n.P+n.B+n.O)),final.corr.mat)
   
   if (n.P!=0) {
     PP=matrix(0, n, n.P)
     for (i in 1: length(lambda.vec)) {
       PP[,i]=qpois(pnorm(myz[,i]), lambda.vec[i])}
   } else PP=NULL
   
   
   if (n.B!=0){
     BB=matrix(0,n, n.B) 
     for (j in (n.P+1):(n.P+n.B)){
       for (i in 1:n){
         if (1*myz[i,j]>qnorm(1-prop.vec[j-n.P])) BB[i,j-n.P]=1 else BB[i,j-n.P]=0}
     }} else BB=NULL
   
   
   if (n.O!=0){
     OO=matrix(0, n, n.O)
     for (i in 1:length(prop.list)){
       y=myz[,(n.P+n.B+i)]
       pvec=prop.list[[i]]
       yord = numeric(n)
       for (r in 1:length(pvec)){
         if (r !=length(pvec)) {
           t1 = qnorm(pvec[r])
           t2 = qnorm(pvec[r+1] )
           yord[(t1<y)&(y<=t2)]= r
         } else {
           yord[y>qnorm(pvec[r])]= r
         }
       }
       OO[,i]=yord}
   } else OO=NULL
   
   
   mydata=cbind(PP, BB, OO)

return(mydata)
}
