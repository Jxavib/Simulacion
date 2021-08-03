# ==================================================================
# Práctica 9: Reducción de la varianza
# ==================================================================

# ------------------------------------------------------------------
# Ejercicio 1
# ------------------------------------------------------------------
# Variables antitéticas en integración Monte Carlo 

# Función objetivo
a <- 0; b <- 2
ftn <- function(x) return(exp(x)/(b-a))

curve(ftn, a, b, ylim=c(0,4))
abline(h=0,lty=2)
abline(v=c(a,b),lty=2)

# Se trata de calcular la media de exp(U(a,b))
res <- (exp(b)-exp(a))/(b-a)
res

# -------------------------------
mc.integral <- function(ftn, a, b, n, plot=TRUE) {
# Integración Monte Carlo
# -------------------------------
  fx <- sapply(runif(n, a, b), ftn)*(b-a)
  if (plot) {
    estint <- cumsum(fx)/(1:n)
    esterr <- sqrt(cumsum((fx-estint)^2))/(1:n)
    plot(estint, ylab="Media y rango de error", type="l", lwd= 2, ylim=mean(fx)+2*c(-esterr[1],esterr[1]), xlab="Iteraciones")
    lines(estint+2*esterr, col="gold", lwd=2)
    lines(estint-2*esterr, col="gold", lwd=2)
    return(list(valor=estint[n], error=2*esterr[n]))  
  } else return(list(valor=mean(fx), error=2*sd(fx)/sqrt(n)))
}  

set.seed(54321)
mc.integral(ftn, a, b, 500)
abline(h=res, lty=2)

# -------------------------------
mc.integrala <- function(ftn, a, b, n, plot=TRUE,...) {
# Integración Monte Carlo con variables antitéticas
# (Solo se genera la mitad)
# Se compara con el mismo nº de evaluaciones
# -------------------------------
  x <- runif(n%/%2, a, b)
  # La siguiente línea solo para representar alternando
  x <- as.numeric(matrix(c(x,a+b-x),nrow=2,byrow=TRUE))
  # bastaría con emplear p.e. c(x,a+b-x)
  fx <- sapply(x, ftn)*(b-a)
  
  if (plot) {
    estint <- cumsum(fx)/(1:n)
    esterr <- sqrt(cumsum((fx-estint)^2))/(1:n)
#    plot(estint, ylab="Media y rango de error",type="l",lwd= 2,ylim=mean(fx)+2*c(-esterr[1],esterr[1]),xlab="Iteraciones",...)
    lines(estint,lwd= 2,col="green")
    lines(estint+2*esterr,col="green",lwd=2)
    lines(estint-2*esterr,col="green",lwd=2)
    return(list(valor=estint[n],error=2*esterr[n]))
  } else return(list(valor=mean(fx),error=2*sd(fx)/sqrt(n)))
}

set.seed(54321)
mc.integrala(ftn, a, b, 500)

# Porcentaje de reducción del error ESTIMADO
100*(0.1619886-0.1610066)/0.1619886
# La reducción teórica de la varianza es del 96.7%
# Problema: en el segundo caso se está estimando mal la varianza...


# -------------------------------
mc.integrala2 <- function(ftn, a, b, n, plot = TRUE,...) {
# Integración Monte Carlo con variables antitéticas
# (Solo se genera la mitad)
# Se compara con el mismo nº de evaluaciones
# -------------------------------
  x <- runif(n%/%2, a, b)
  # La siguiente línea solo para representar alternando
  x <- matrix(c(x,a+b-x),nrow=2,byrow=TRUE)
  # bastaría con emplear p.e. c(x,a+b-x)
  fx <- apply(x, 1,  ftn)*(b-a)
  corr <- cor(fx[,1], fx[,2])
  fx <- as.numeric(fx)
  return(list(valor=mean(fx), error=2*sd(fx)/sqrt(n)*sqrt(1+corr)))
}

set.seed(54321)
mc.integrala2(ftn, a, b, 500)

# Porcentaje de reducción del error ESTIMADO
100*(0.1603748-0.05700069)/0.1603748


# ------------------------------------------------------------------
# Ejercicio 2
# ------------------------------------------------------------------
# Integración Monte Carlo con estratificación

# -------------------------------
mc.integrale <- function(ftn, a, b, n, k) {
# Integración Monte Carlo con estratificación
# -------------------------------
  l <- n%/%k
  int <- seq(a, b, len=k+1)
  x <- runif(l*k, rep(int[-(k+1)], each=l), rep(int[-1], each=l))
  # l uniformes en cada uno de los intervalos [(j-1)/k , j/k]
  fx <- sapply(x, ftn)*(b-a)
  return(list(valor=mean(fx), error=2*sd(fx)/sqrt(n)))   # error mal calculado
}


set.seed(54321)
mc.integral(ftn, a, b, 500)
abline(h=res,lty=2)

set.seed(54321)
mc.integrale(ftn, a, b, 500, 50)

set.seed(54321)
mc.integrale(ftn, a, b, 500, 100)
# No se tiene en cuenta la variabilidad en el estrato
# El tamaño de las submuestras debería incrementarse hacia el extremo superior


# ------------------------------------------------------------------
# Ejercicio 3
# ------------------------------------------------------------------
# Integración Monte Carlo con variables de control

# Se trata de calcular la media de exp(U(a,b))
a <- 0; b <- 2
res <- (exp(b)-exp(a))/(b-a)
res

set.seed(54321)
nsim <- 1000
u <- runif(nsim, a, b)
expu <- exp(u)
mean(expu) # aproximación de la característica de interés

# Con variable control
plot(u,expu)
reg <- lm(expu~u)$coef
abline(reg,col='blue')

reg[1]+reg[2] # Conincidirá con la solución mean(expuc)
# Lo siguiente ya no sería necesario
expuc <- expu-reg[2]*(u-1)
mean(expuc)  

# Estimación del porcentaje de reducción en la varianza
100*(var(expu)-var(expuc))/var(expu)


# ------------------------------------------------------------------
# Ejercicio 4
# ------------------------------------------------------------------

set.seed(54321)
nsim <- 1000
x <- 2.5
thet <- rnorm(nsim,mean=x)
delt <- thet/(1+thet^2)
mean(delt) # aproximación de la característica de interés

# Con variable control
# plot(thet,delt)
moms <- thet - x # variable control de media cero
reg <- lm(delt~moms)$coef
reg[1] # La solución mean(delta) coincidirá con reg[1]
# Lo siguiente ya no sería necesario
delta <- delt-reg[2]*moms
mean(delta)   

# Gráfico evolución
plot(cumsum(delt)/(1:nsim), type="l", lwd=2, lty=2)
lines(cumsum(delta)/(1:nsim), lwd=2)
#Compararlo con integrate
integrando <- function(x) {x/(x^2+1)*dnorm(x,2.5)}
abline(h=integrate(integrando, -Inf, Inf)$value)

# Estimación del porcentaje de reducción en la varianza
100*(var(delt)-var(delta))/var(delt)







