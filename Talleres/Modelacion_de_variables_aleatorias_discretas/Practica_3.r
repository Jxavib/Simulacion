# ==================================================================
# Métodos universales para variables discretas
# ==================================================================


# ------------------------------------------------------------------
# Ejercicio 1
# ------------------------------------------------------------------
# Distribución binomial

# -------------------------------
# Apartado a

rfmp <- function(x, prob = 1/length(x), nsim = 1000) {
# Simulación nsim v.a. discreta a partir de fmp
# por inversión generalizada (transformación cuantil)
  # Inicializar FD
  Fx <- cumsum(prob)
  # Simular
  X <- numeric(nsim)
  U <- runif(nsim)
  ncomp <<- 0 # Comentar esta línea para uso normal
  for(j in 1:nsim) {
    i <- 1
    while (Fx[i] < U[j]) i <- i + 1
    X[j] <- x[i]
    ncomp <<- ncomp + i # Comentar esta línea para uso normal
  }
  return(X)
}


# Simular binomial
set.seed(54321)
n <- 10
p <- 0.5
nsim <- 10^5

x <- 0:n
fmp <- dbinom(x, n, p)

# ncomp <- 0
system.time(
rx <- rfmp(x, fmp, nsim)
)

# Número de comparaciones
ncomp/nsim
sum((1:length(x))*fmp)

# Aproximación media
mean(rx)
n*p

# Análisis resultados
res <- as.data.frame(table(rx)/nsim)
names(res) <- c("x", "psim")

# Comparación teórica
plot(as.matrix(res), type="h")
points(x, fmp, pch=4)

res$pteor <- fmp
res

# Errores
max(abs(res$psim - res$pteor))
max(abs(res$psim - res$pteor) / res$pteor)

# NOTA: Puede ocurrir que no todos los valores sean generados en la simulación
# Si length(x) > length(psim) el código anterior (res$pteor <- fmp) producirá un error
# Alternativamente:
psim <- rep(0, length(x))
names(psim) <- as.character(x)
psim[as.character(res$x)] <- res$psim
psim
res <- data.frame(x=x, pteor=fmp, psim =psim)
res


# -------------------------------
# Apartado b
# Ordenar de forma decreciente las probabilidades

tini <- proc.time()

ind <- order(fmp, decreasing=TRUE)
rx <- rfmp(x[ind], fmp[ind], nsim)

tiempo <- proc.time() - tini
tiempo

# -------------------------------
# Comandos R: Función sample
# sample(x, size, replace = FALSE, prob = NULL)

system.time(
rx <- sample(x, nsim, replace = TRUE, prob = fmp)
)
# Método recomendado en R


# -------------------------------
# Si n >> 30
# Aproximación normal con corrección continuidad
mean <- n*p
sd <- sqrt(n*p*(1-p))

res$pnorm <- pnorm(x+0.5,mean,sd)-pnorm(x-0.5,mean,sd)
res
max(abs(res$pnorm - res$pteor))
max(abs(res$pnorm - res$pteor) / res$pteor)

system.time(
rx <- pmax(0, pmin(n, round( rnorm(nsim, mean, sd) )))
)

table(rx)/nsim

# Realmente las prob de los extremos serían mayores...
res$pnorm[1] <- pnorm(0.5, mean, sd)
res$pnorm[n+1] <- 1 - pnorm(n-0.5, mean, sd)


# -------------------------------
# Apartado c

rfmp.tabla <- function(x, prob = 1/length(x), m, nsim = 1000) {
# Simulación v.a. discreta a partir de función de masa de probabilidad
# por tabla guia de tamaño m
  # Inicializar tabla y FD
  Fx <- cumsum(prob)
  g <- rep(1,m)
  i <- 1
  for(j in 2:m) {
    while (Fx[i] < (j-1)/m) i <- i + 1
    g[j] <- i
  }
  # Generar valores
  X <- numeric(nsim)
  U <- runif(nsim)
  # ncomp <<- 0 # Comentar esta línea para uso normal
  for(j in 1:nsim) {
    i <- g[floor(U[j]*m) + 1]
    # ncomp <<- ncomp - i + 1 # Comentar esta línea para uso normal
    while (Fx[i] < U[j]) i <- i + 1
    # ncomp <<- ncomp + i     # Comentar esta línea para uso normal
    X[j] <- x[i]
  }
  return(X)
}

system.time(
rx <- rfmp.tabla(x, fmp, n-1, nsim)
)

# Análisis resultados
res <- as.data.frame(table(rx)/nsim)
names(res) <- c("x", "psim")

# Comparación teórica
plot(as.matrix(res), type="h")
points(x, fmp, pch=4)

