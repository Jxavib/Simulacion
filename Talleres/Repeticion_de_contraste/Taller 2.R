# Ejercicio 0.1
###################################################################
##    P-valor es una variable aleatoria Uniforme[0,1]             #              
###################################################################
# Para verificar mediante simulación que el p-valor tiene una     #
# distribución Uniforme [0,1] bajo hipótesis nula H0:theta=theta0#
# El escenario de simulación es el siguiente:                     #
# -  Sistema de hipótesis: Ho: mu=0   VS.    Ha: mu <> 0.          #
# -  Poblaciones con distribución normal (mu,1). (mu=0, mu=2)      #
# -  Muestra aleatoria de tamaño n=20 tomada de una población      #
#   normal.                                                       #
# -	Aplicación de una prueba t, es decir el estadístico de prueba # 
#   sigue una distribución t-student con (n-1) grados de libertad.# 
# -	Se realizan nsim=10000 simulaciones                           #
#-----------------------------------------------------------------#
# Bajo la hipótesis nula Ho: mu=0, se considera los siguientes    #
# casos:                                                          #
#-----------------------------------------------------------------#
set.seed(7654321)                                                 #
nsim<-10000                                                       #
n<-20                                                             #
#
estadistico<-function(x,mu){                                      #
  n<-length(x)                                                    #
  xbar<-mean(x)                                                   #
  s<-sqrt(var(x))                                                 #
  T<-abs((xbar-mu)/(s/sqrt(n)))                                   # 
  p<-2*(1-pt(T,n-1))                                              #
}                                                                 #
#
par(mfrow=c(1,2))                                                 #
###################### H0 mu=0  vs. H0 mu<>0 ######################                     
# Generamos nsim muestras distribuidas normalmente de tamaño n    #
# con media mu diferente a cero                                   #
#-----------------------------------------------------------------#
# CASO 1: N(mu=0,1)                                               #
# Descripción resultados:                                         # 
# El histograma de los p-valores toma la forma de una             #
# distribución plana y uniforme sobre el intervalo [0, 1].        #
# Por tanto podemos verificar que un p-valor no es la probabilidad#
# de que Ho sea cierto, ya que al considerar una muestra de una   #
# población normal (0,1), claramente Ho es cierta, sin embargo el #
# p-valor está uniformemente distribuido entre cero y uno.        #
#-----------------------------------------------------------------#         
mu<-0                                                             #
X<-matrix(rnorm(n*nsim,mu,1),ncol=nsim,nrow=n)                    #
p<-apply(X,2,estadistico,mu=0)                                    #
hist(p,main=c("Caso1:mu = 0"))                                    #
#-----------------------------------------------------------------#   
# CASO 2: N(mu=0.5,1)                                             #   
#-----------------------------------------------------------------#
# Descripción resultados:                                         #
# El histograma de los p-valores ya no es uniforme.               #
# Puede observarse que existe más chance de obtener p-valores     #
# menores al nivel de significación, será más alto estos valores  #
# bajo la hipótesis alternativa que bajo la hipótesis nula y ese  #
# efecto es más claro a medida que mu incrementa su valor.        #
#-----------------------------------------------------------------#
mu<-0.5                                                           #
X<-matrix(rnorm(n*nsim,mu,1),ncol=nsim,nrow=n)                    #
p<-apply(X,2,estadistico,mu=0)                                    #
hist(p,main=c("Caso2:mu = 0.5"))                                  #
#-----------------------------------------------------------------#

# Ejercicio 0.2                                                   
###################################################################
##         Gráficos de Potencia (Ho=1 vs Ho>1)                    #
###################################################################
potencia.n<-function(theta,n,alpha){                              #
  1-((1-alpha))/theta^n   # 1-((1- alpha))/theta^n                        #
}                                                                 #
#-----------------------------------------------------------------#
# alpha=0.05                                                      #
#-----------------------------------------------------------------#
alpha<-0.05                                                       #
# Seteamos el área del gráfico                                    #  
x<-seq(1,2,0.1)                                                   #
a<-expression(pi(theta)==1-over((1-alpha),theta^n))               #
plot(potencia.n(x,alpha=0.05,n=10),type="n",ylim=c(0.05,1),       #
     xlim=c(0.8,2),ylab="Potencia",xlab="Theta",main=a,cex.main=1,#
     xaxs="i",bty="n")                                            #
# gráficos las diferentes curvas de potencia                      #
curve(potencia.n(x,alpha,n=10),from=1,to=2,add=T,col=1,lty=1,     #
      lwd=2)                                                      #
curve(potencia.n(x,alpha,n=20),from=1,to=2,add=T,col=2,lty=2,     #
      lwd=2)                                                      #
curve(potencia.n(x,alpha,n=50),from=1,to=2,add=T,col=3,lty=3,     #
      lwd=2)                                                      #
# 
legend("right",,paste("n",c(10,20,50),sep="="),lty=1:3,col=1:3,   #
       title="Alpha=0.05",,bty="n")                               # 
#-----------------------------------------------------------------#
# alpha=0.1                                                       #
#-----------------------------------------------------------------#
alpha<-0.1                                                        #
# Seteamos el área del gráfico                                    #  
x<-seq(1,2,0.1)                                                   #
a<-expression(pi(theta)==1-over((1-alpha),theta^n))               #
plot(potencia.n(x,alpha=0.1,n=10),type="n",ylim=c(0.05,1),        #
     xlim=c(0.8,2),ylab="Potencia",xlab="Theta",main=a,cex.main=1,#
     xaxs="i",bty="n")                                            #
# gráficos las diferentes curvas de potencia                      #
curve(potencia.n(x,alpha,n=10),from=1,to=2,add=T,col=1,lty=1,     #
      lwd=2)                                                      #
curve(potencia.n(x,alpha,n=20),from=1,to=2,add=T,col=2,lty=2,     #
      lwd=2)                                                      #
curve(potencia.n(x,alpha,n=50),from=1,to=2,add=T,col=3,lty=3,     #
      lwd=2)                                                      #
# 
legend("right",,paste("n",c(10,20,50),sep="="),lty=1:3,col=1:3,   #
       title="Alpha=0.05",,bty="n")                               # 
#-----------------------------------------------------------------#
# n=10,20 y 50                                                            #
#-----------------------------------------------------------------#
# Seteamos el área del gráfico                                    #  
x<-seq(1,2,0.1)                                                   #
a<-expression(pi(theta)==1-over((1-alpha),theta^n))               #
plot(potencia.n(x,alpha=0.1,n),type="n",ylim=c(0.05,1),           #
     xlim=c(0.8,2),ylab="Potencia",xlab="Theta",main=a,cex.main=1,#
     xaxs="i",bty="n")                                            #
# gráficos las diferentes curvas de potencia                      #
n=10
curve(potencia.n(x,alpha=0.05,n),from=1,to=2,add=T,col=1,lty=1,   #
      lwd=2)                                                      #
curve(potencia.n(x,alpha=0.1,n),from=1,to=2,add=T,col=2,lty=2,    #
      lwd=2)                                                      #
n<-20                                                             # 
curve(potencia.n(x,alpha=0.05,n),from=1,to=2,add=T,col=1,lty=1,   #
      lwd=2)                                                      #
curve(potencia.n(x,alpha=0.1,n),from=1,to=2,add=T,col=2,lty=2,    #
      lwd=2)                                                      #
#
n<-50                                                             # 
curve(potencia.n(x,alpha=0.05,n),from=1,to=2,add=T,col=1,lty=1,   #
      lwd=2)                                                      #
curve(potencia.n(x,alpha=0.1,n),from=1,to=2,add=T,col=2,lty=2,    #
      lwd=2)                                                      #
#
legend("right",,paste("alpha",c(0.05,0.1),sep="="),lty=1:2,col=1:2#
       ,title="Alpha=0.05",,bty="n")                              # 
#-----------------------------------------------------------------#
# Ejercicio Práctico                                                   
###################################################################
##                  Comportamiento p-valor                        #
###################################################################

#d)  Prepare algunas simulaciones para estudiar el comportamiento del p-valor:
#  - Simular 1000 muestras de tamaño n (por ejemplo, n = 10, 20, 50) de 
#una Uniforme [0, theta] distribución (por ejemplo, theta = 1, 1,05, 1,10).
#- Para cada muestra, realizar la prueba y calcular el correspondiente valor de p.
#- Utilizar las estadísticas de resumen, histogramas y / o box-plot para describir 
#el comportamiento de los p-valores obtenidos. Comentar los resultados.


n<-c(10,20,50)
theta<-c(1,1.05,1.10)

n<-10
theta<-4
nsim<-10000
set.seed(7654321)
muestra<-matrix(runif(n*nsim,0,theta),ncol=nsim,nrow=n)
estadistico<-apply(muestra,2,max)
hist(estadistico)
valores<-quantile(estadistico,c(0.90,0.95,0.99));valores
evaluar<-(valores/theta)^n;evaluar
p.valor<-(estadistico/theta)^n
hist(p.valor)
boxplot(p.valor)
valores2<-quantile(p.valor,c(0.90,0.95,0.99));valores2

alpha<-c(0.90,0.95,0.99)
(1-alpha)^(1/n)

