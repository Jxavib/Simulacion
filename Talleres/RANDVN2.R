# ==================================================================
# Generador Von Neumann de n?meros pseudoaleatorios
# ==================================================================

# Borrar todos los objetos (problemas semilla)
rm(list=ls(all=TRUE))

# ------------------------------------------------------------------
# initRANDVN(semilla,n)
#   Inicia el generador 
#   n n?mero de digitos centrales, 4 por defecto (debe ser un n? par)
#   Por defecto semilla del reloj
#   OJO: No se hace ninguna verificaci?n de los par?metros
# ------------------------------------------------------------------
initRANDVN <- function(semilla=as.numeric(Sys.time()), n=4) {
  .semilla <<- as.double(semilla) %% 10^n  #C?lculos en doble precisi?n
  .n <<- n
  .aux <<- 10^(2*n-n/2)
  .aux2 <<- 10^(n/2)
  return(invisible(list(semilla=.semilla,n=.n)))
}

# ------------------------------------------------------------------
# RANDVN()
#   Genera un valor pseudoaleatorio con el generador de Von Neumann
#   Actualiza la semilla (si no existe llama a initRANDVN)
# ------------------------------------------------------------------
RANDVN <- function() {
  if (!exists(".semilla", envir=globalenv())) initRANDVN(1750,4)
  z <- .semilla^2
  .semilla <<- trunc((z-trunc(z/.aux)*.aux)/.aux2)
  return(.semilla/10^.n)
}


# ------------------------------------------------------------------
# RANDVNN(n)
#   Genera un vector de valores pseudoaleatorios con el generador congruencial
#   (por defecto de dimensi?n 1000)
#   Actualiza la semilla (si no existe llama a initRANDVN)
# ------------------------------------------------------------------
RANDVNN <- function(n=500) {
  x <- numeric(n)
  for(i in 1:n) x[i] <- RANDVN()
  return(x)
  # return(replicate(n,RANDVN()))  # Alternativa m?s r?pida
}

# x <- RANDVNN()
# x
# plot(x,xlab = "Indice",ylab = "Aleatorios",main ="Grafica con semilla=1750")
