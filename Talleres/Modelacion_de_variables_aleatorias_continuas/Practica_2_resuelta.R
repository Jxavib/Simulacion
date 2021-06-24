# ==================================================================
# A entregar por los alumnos
# ==================================================================

# ------------------------------------------------------------------
# Ejercicio 4
# ------------------------------------------------------------------
# Distribución de Cauchy por inversión

rcauchyn <- function(n = 1000) {
# Simulación n valores de distribución de Cauchy
    return(tan(pi*runif(n)))
}


# -------------------------------
# Apartado a

set.seed(1)
nsim <- 10^4
system.time(x <- rcauchyn(nsim))


# -------------------------------
# Apartado b

hist(x, breaks = "FD", freq = FALSE, xlim=c(-10,10))
curve(dcauchy(x), add = TRUE)


# -------------------------------
# Apartado c

plot(1:nsim, cumsum(x)/(1:nsim), type="l", ylab="Media muestral", xlab="Nº de simulaciones")


# ------------------------------------------------------------------
# Ejercicio 5
# ------------------------------------------------------------------

dtiempo <- function(x) ifelse(x > 0, x * exp(-x), 0)
curve(dtiempo, 0, 10)

lambda.opt <- 0.5
c.opt <- 4/exp(1)

rtiempoAR <- function() {
# Simulación por aceptación-rechazo
# dtiempo a partir de exponencial
  while (TRUE) {
    U <- runif(1)
    X <- rexp(1, lambda.opt)
    ngen <<- ngen + 1 # Comentar esta línea para uso normal
    if (c.opt * U * dexp(X, lambda.opt) <= dtiempo(X)) return(X)
  }
}

rtiempoARn <- function(n=1000) {
# Simulación n valores dtiempo
    x <- numeric(n)
    ngen <<- 0 # Comentar esta línea para uso normal
    for(i in 1:n) x[i]<-rtiempoAR()
    return(x)
}


# Grafico
curve(c.opt * dexp(x, lambda.opt), ylab = "c*g(x)", xlim=c(0, 10), lty = 2)
curve(dtiempo(x), add = TRUE)


# -------------------------------
# Apartado a

# Obtención de un valor c óptimo aproximado
# NOTA: a partir del gráfico del encabezado se deduce que en x=+-2 se obtiene el máximo
optimize(f=function(x){dtiempo(x)/dexp(x, lambda.opt)}, maximum=TRUE, interval=c(1,3))$objective
# Valor óptimo real
c.opt


# Obtención de valores c y lambda óptimos aproximados
fopt <- function(lambda) {
# Obtiene c fijado lambda
  optimize(f=function(x){dtiempo(x)/dexp(x,lambda)}, maximum=TRUE, interval=c(1,3))$objective
}

# Encontar lambda que minimiza
optimize(f = function(x){fopt(x)}, interval = c(0.1,1))


# -------------------------------
# Apartado b

set.seed(54321)
nsim <- 10^4

ngen <- 0
system.time(x <- rtiempoARn(nsim))
# Nº generaciones
{
cat("\nNº de generaciones = ", ngen)
cat("\nNº medio de generaciones = ", ngen/nsim)
cat("\nProporción de rechazos = ", 1-nsim/ngen, "\n")
}

# -------------------------------
# Apartado c

hist(x, breaks="FD", freq=FALSE)
curve(dtiempo(x), add=TRUE)
