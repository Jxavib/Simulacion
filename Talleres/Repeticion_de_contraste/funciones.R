potencia.n<-function(theta,n,alpha){                              #
  1-((1-alpha))/theta^n   # 1-((1- alpha))/theta^n                #
}                                                                 #

# Gráfico de potencia
#n.v=tamaños muestrales
#alpha=nivel de significancia

potencia.g<-function(n.v,alpha){
  x<-seq(1,2,0.1)                                                   #
  a<-expression(pi(theta)==1-over((1-alpha),theta^n))               #
  plot(c(0,0),type="n",ylim=c(0.05,1),       #
       xlim=c(0.8,2),ylab="Potencia",xlab="Theta",main=a,cex.main=1,#
       xaxs="i",bty="n")                                            #
for (i in seq(n.v))
  {
  curve(potencia.n(x,n.v[i],alpha),from=1,to=2,add=T,col=i,lty=1,     #
        lwd=2)                                                      #
  }

  legend("right",,paste("n",n.v,sep="="),lty=1:3,col=seq(n.v),   #
         title=paste("Alpha=",alpha*100,"%"),bty="n"  )    
  
}

# función para construir los estadísticos

estadistico<-function(n,theta){  
  set.seed(7654321)
  theta.null=1
  uniforme<-function(n,theta){
    nsim<<-1000 
    runif(n*nsim,0,theta)
  }
  simulacion<-mapply(uniforme,n,theta)
  maxi<-matrix(,nrow=1000,ncol=length(n))
  p.valor<-matrix(,nrow=1000,ncol=length(n))
  for (i in seq(n)){    

    simulacion.m<-matrix(simulacion[[i]],nsim,n[i])
    maxi[,i]<-apply(simulacion.m,1,max) # estadístico
    p.valor[,i]<-1-pmin(1,(maxi[,i]/theta.null)^n[i]) #p-valor
  }
  return(list(maxi,p.valor))
}

EDA2<-function (x, trim = 0.05) 
{
  require(e1071)
  Altblue <- "#A9E2FF"
  Adkblue <- "#0080FF"
  Ared <- "#C51111"
  # varname <- deparse(substitute(x))
  varname<-"p-valores"
  N <- length(x)
  UM <- sum(is.na(x))
  n <- N - UM
  x <- x[!(is.na(x) > 0)]
  LQ1 <- (n + 1)/4
  LQ3 <- (3 * (n + 1))/4
  Sort <- sort(x)
  V1 <- floor(LQ1)
  V2 <- floor(LQ3)
  V3 <- V1 + 1
  V4 <- V2 + 1
  Q1 <- round(Sort[V1] + (LQ1 - V1) * (Sort[V3] - Sort[V1]), 
              3)
  Q3 <- round(Sort[V2] + (LQ3 - V2) * (Sort[V4] - Sort[V2]), 
              3)
  IQR <- round(Q3 - Q1, 3)
  Min <- round(min(x), 3)
  Max <- round(max(x), 3)
  Stdev <- round(sd(x, na.rm = TRUE), 3)
  Mean <- round(mean(x, na.rm = TRUE), 3)
  Median <- round(median(x, na.rm = TRUE), 3)
  TrMean <- round(mean(x, trim = trim), 3)
  Var <- round(var(x, na.rm = TRUE), 3)
  SE <- round(Stdev/sqrt(n), 3)
  Range <- round(Max - Min, 3)
  par(mfrow=c(2,1))
  hist(x, probability = TRUE, col = Adkblue, xlab = "", ylab = "", 
       axes = FALSE, main = paste("Histogram of", varname))

  l.out <- x[x < (Q1 - 1.5 * IQR)]
  r.out <- x[x > (Q3 + 1.5 * IQR)]
  outliers <- c(l.out, r.out)
  rest <- x[x > (Q1 - 1.5 * IQR) & x < (Q3 + 1.5 * IQR)]
  Minrest <- min(rest)
  Maxrest <- max(rest)
  plot(x, x, main = paste("Boxplot of", varname), xlab = "", 
       ylab = "", axes = FALSE, type = "n", xlim = c(min(x), 
                                                     max(x)), ylim = c(0, 1))
  box()
  polygon(c(Q1, Q1, Q3, Q3), c(0.3, 0.7, 0.7, 0.3), density = -1, 
          col = Altblue)
  points(outliers, c(rep(0.5, length(outliers))), col = Ared)
  lines(c(min(rest), Q1), c(0.5, 0.5), lty = 1)
  lines(c(Q3, max(rest)), c(0.5, 0.5), lty = 1)
  lines(c(min(rest), min(rest)), c(0.4, 0.6))
  lines(c(max(rest), max(rest)), c(0.4, 0.6))
  lines(c(Q1, Q1), c(0.3, 0.7))
  lines(c(Q3, Q3), c(0.3, 0.7))
  lines(c(Median, Median), c(0.3, 0.7))
  lines(c(Q1, Q3), c(0.3, 0.3))
  lines(c(Q1, Q3), c(0.7, 0.7))
  points(Mean, 0.5, pch = 16, col = "black")
#   qqnorm(x, col = "black", main = paste("Q-Q Plot of", varname), 
#          xlab = "", ylab = "", axes = FALSE)
#   qqline(x, col = Ared)
#   box()
  mtext("EXPLORATORY  DATA  ANALYSIS", side = 3, outer = TRUE, 
        cex = 1.5, col = Adkblue, line = 1)
  par(oma = c(0, 0, 0, 0))
  par(mfrow = c(1, 1))
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(omi = c(0, 0, 0, 0))
  par(pty = "m")
  #SW <- shapiro.test(x)
  K <- round(kurtosis(x), 3)
  S <- round(skewness(x), 3)
  #SWpval <- round(SW$p.value, 3)
  TOT <- c(n, UM, Min, Q1, Mean, Median, TrMean, Q3, Max, Stdev, 
           Var, SE, IQR, Range, K, S)
  names(TOT) <- c("Size (n)", "Missing", "Minimum", " 1st Qu", 
                  "   Mean", " Median", "TrMean", " 3rd Qu", "   Max.", 
                  " Stdev.", "   Var.", "SE Mean", " I.Q.R.", "  Range", 
                  "Kurtosis", "Skewness")
  par(mfrow=c(1,1))
  return(TOT)
}
