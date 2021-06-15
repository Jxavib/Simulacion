# ==================================================================
# Generador congruencial de n?meros pseudoaleatorios
# ==================================================================

# ------------------------------------------------------------------
# initRANDC(semilla,a,c,m)
#   Selecciona el generador congruencial
#   Por defecto RANDU de IBM con semilla del reloj
#   OJO: No se hace ninguna verificaci?n de los par?metros
# ------------------------------------------------------------------
initRANDC <- function(semilla=as.numeric(Sys.time()), a=2^16+3, c=0, m=2^31) {
  .semilla <<- as.double(semilla) %% m  #C?lculos en doble precisi?n
  .a <<- a
  .c <<- c
  .m <<- m
  return(invisible(list(semilla=.semilla,a=.a,c=.c,m=.m))) #print(initRANDC())
}

# ------------------------------------------------------------------
# RANDC()
#   Genera un valor pseudoaleatorio con el generador congruencial
#   Actualiza la semilla (si no existe llama a initRANDC)
# ------------------------------------------------------------------
RANDC <- function() {
    if (!exists(".semilla", envir=globalenv())) initRANDC()
    .semilla <<- (.a * .semilla + .c) %% .m
    return(.semilla/.m)
}

# ------------------------------------------------------------------
# RANDCN(n)
#   Genera un vector de valores pseudoaleatorios con el generador congruencial
#   (por defecto de dimensi?n 1000)
#   Actualiza la semilla (si no existe llama a initRANDC)
# ------------------------------------------------------------------
RANDCN <- function(n=1000) {
    x <- numeric(n)
    for(i in 1:n) x[i]<-RANDC()
    return(x)
    # return(replicate(n,RANDC()))  # Alternativa m?s r?pida    
}

n <- 50
nsim <- 500
estadistico <- numeric(nsim)
pvalor <- numeric(nsim)

for(isim in 1:nsim) {
  u <- RANDCN(n) # Generar
  # u <- runif(n)
  tmp <- Box.test(u,lag=10, type="Ljung")
  estadistico[isim] <- tmp$statistic
  pvalor[isim] <- tmp$p.value
}
plot(x=estadistico,y=pvalor,main = "Grafico Estadistico vs P-valor",xlab ="estadistico",
     ylab = "p-valor")