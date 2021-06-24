# ==================================================================
# Práctica 2: Métodos universales para la simulación de variables continuas
# ==================================================================


# ------------------------------------------------------------------
# Ejercicio 1
# ------------------------------------------------------------------
# Distribución doble exponencial por inversión

# -------------------------------
# Apartado a

ddexp <- function(x, lambda = 1){
# Densidad doble exponencial
  lambda*exp(-lambda*abs(x))/2
}

rdexp <- function(lambda = 1){
# Simulación por inversión
# Doble exponencial
  U <- runif(1)
  if (U<0.5) {
    return(log(2*U)/lambda)
  } else {
    return(-log(2*(1-U))/lambda)
  }
}

rdexpn <- function(n = 1000, lambda = 1) {
# Simulación n valores de doble exponencial
    x <- numeric(n)
    for(i in 1:n) x[i]<-rdexp(lambda)
    return(x)
}


# -------------------------------
# Apartado b

set.seed(54321)
system.time(x <- rdexpn(10^4, 2))

# -------------------------------
# Apartado c

hist(x, breaks = "FD", freq = FALSE)
curve(ddexp(x, 2), add = TRUE)



# ------------------------------------------------------------------
# Ejercicio 2
# ------------------------------------------------------------------
# Simulación por aceptación-rechazo de normal estandar a partir de doble exponencial

# EJECUTAR CÓDIGO DEL APARTADO A DEL EJERCICIO 1

rnormAR <- function() {
# Simulación por aceptación-rechazo
# Normal estandar a partir de doble exponencial
  while (TRUE) {
    U <- runif(1)
    X <- rdexp(1)
    ngen <<- ngen + 1 # Comentar esta línea para uso normal
    if (U*exp((X^2+1)*0.5-abs(X)) <= 1) return(X)
  }
}

rnormARn <- function(n=1000) {
# Simulación n valores N(0,1)
    x <- numeric(n)
    ngen <<- 0 # Comentar esta línea para uso normal
    for(i in 1:n) x[i]<-rnormAR()
    return(x)
}


# Grafico
curve(sqrt(2*exp(1)/pi)*ddexp(x), ylab = "c*g(x)", xlim=c(-4,4), lty = 2)
curve(dnorm(x), add = TRUE)


# -------------------------------
# Apartado a

set.seed(54321)
nsim <- 10^4

# ngen <- 0
system.time(x <- rnormARn(nsim))
# Nº generaciones
{
cat("\nNº de generaciones = ", ngen)
cat("\nNº medio de generaciones = ", ngen/nsim)
cat("\nProporción de rechazos = ", 1-nsim/ngen, "\n")
}

# -------------------------------
# Apartado b

hist(x, breaks="FD", freq=FALSE)
curve(dnorm(x), add=TRUE)

# -------------------------------
# Apartado c

# Obtención de un valor c óptimo aproximado
# NOTA: a partir del gráfico del encabezado se deduce que en x=+-1 se obtiene el máximo
optimize(f=function(x){dnorm(x)/ddexp(x)}, maximum=TRUE, interval=c(0.5,1.5))$objective
# Valor óptimo real
sqrt(2*exp(1)/pi)


# -------------------------------
# Apartado d

# Obtención de valores c y lambda óptimos aproximados
fopt <- function(lambda) {
# Obtiene c fijado lambda
  optimize(f=function(x){dnorm(x)/ddexp(x,lambda)}, maximum=TRUE, interval=c(0.5,1.5))$objective
}

# Encontar lambda que minimiza
optimize(f=function(x){fopt(x)}, interval=c(0.5,2))


# ------------------------------------------------------------------
# Ejercicio 3
# ------------------------------------------------------------------
# Aceptación-Rechazo en Inferencia Bayesiana

# -------------------------------
# Apartado a

mu0 <- 1
n <- 10
nsim <- 10^3
set.seed(54321)
x <- rnorm(n, mean = mu0)

# Función de verosimilitud
lik <- function(mu){prod(dnorm(x, mean = mu))}

# Cota óptima
# Estimación por máxima verosimilitud
emv <- optimize(f = lik, int = range(x), maximum = TRUE)
emv
c <- emv$objective

mean(x)
c <- lik(mean(x))
c

# Simulación distribución a posteriori
# Simulación por aceptación-rechazo a partir de Cauchy(0,1)
ngen <- nsim
Y <- rcauchy(nsim)
ind <- (c*runif(nsim) > sapply(Y, lik)) # TRUE si no verifica condición
# Volver a generar si no verifica condición
while (sum(ind)>0){
  le <- sum(ind)
  ngen <- ngen + le
  Y[ind] <- rcauchy(le)
  ind[ind] <- (c*runif(le) > sapply(Y[ind], lik)) # TRUE si no verifica condición
}


# Nº generaciones
{
cat("\nNº de generaciones = ", ngen)
cat("\nNº medio de generaciones = ", ngen/nsim)
cat("\nProporción de rechazos = ", 1-nsim/ngen,"\n")
}

# Intervalo de probabilidad al 95% (IC Bayes)
q <- quantile(Y, c(0.025, 0.975))

# Representar estimador e IC Bayes
hist(Y, freq=FALSE, main="Distribución a posteriori")
abline(v = mean(x), lty = 2, lwd = 2)
abline(v = q, lty = 2)

