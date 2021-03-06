# ==================================================================
# A entregar por los alumnos
# ==================================================================

# ------------------------------------------------------------------
# Ejercicio 4
# ------------------------------------------------------------------
# Distribuci�n de Cauchy por inversi�n

rcauchyn <- function(n = 1000) {
# Simulaci�n n valores de distribuci�n de Cauchy
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

plot(1:nsim, cumsum(x)/(1:nsim), type="l", ylab="Media muestral", xlab="N� de simulaciones")


# ------------------------------------------------------------------
# Ejercicio 5
# ------------------------------------------------------------------

dtiempo <- function(x) ifelse(x > 0, x * exp(-x), 0)
curve(dtiempo, 0, 10)

lambda.opt <- 0.5
c.opt <- 4/exp(1)

rtiempoAR <- function() {
# Simulaci�n por aceptaci�n-rechazo
# dtiempo a partir de exponencial
  while (TRUE) {
    U <- runif(1)
    X <- rexp(1, lambda.opt)
    ngen <<- ngen + 1 # Comentar esta l�nea para uso normal
    if (c.opt * U * dexp(X, lambda.opt) <= dtiempo(X)) return(X)
  }
}

rtiempoARn <- function(n=1000) {
# Simulaci�n n valores dtiempo
    x <- numeric(n)
    ngen <<- 0 # Comentar esta l�nea para uso normal
    for(i in 1:n) x[i]<-rtiempoAR()
    return(x)
}


# Grafico
curve(c.opt * dexp(x, lambda.opt), ylab = "c*g(x)", xlim=c(0, 10), lty = 2)
curve(dtiempo(x), add = TRUE)


# -------------------------------
# Apartado a

# Obtenci�n de un valor c �ptimo aproximado
# NOTA: a partir del gr�fico del encabezado se deduce que en x=+-2 se obtiene el m�ximo
optimize(f=function(x){dtiempo(x)/dexp(x, lambda.opt)}, maximum=TRUE, interval=c(1,3))$objective
# Valor �ptimo real
c.opt


# Obtenci�n de valores c y lambda �ptimos aproximados
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
# N� generaciones
{
cat("\nN� de generaciones = ", ngen)
cat("\nN� medio de generaciones = ", ngen/nsim)
cat("\nProporci�n de rechazos = ", 1-nsim/ngen, "\n")
}

# -------------------------------
# Apartado c

hist(x, breaks="FD", freq=FALSE)
curve(dtiempo(x), add=TRUE)
