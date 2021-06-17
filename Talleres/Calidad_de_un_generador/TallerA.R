#' ------
#' 
#' Ejercicios de fin de práctica
#' -----------------------------
#' 
#' ------
#' 
#' 6.  Uno de los primeros generadores fue el denominado método de los
#'     cuadrados medios propuesto por Von Neumann (1946). Con este
#'     procedimiento se generan números pseudoaleatorios de 4 dígitos de la
#'     siguiente forma:
#' 
#'     i.  Se escoge un número de cuatro dígitos $x_{0}$ (semilla).
#' 
#'     ii.  Se eleva al cuadrado ($x_{0}^{2}$) y se toman los cuatro dígitos
#'         centrales ($x_{1}$).
#' 
#'     iii.  Se genera el número pseudo-aleatorio
#'         como$$u_{1}=\frac{x_{1}}{10^{4}}.$$
#' 
#'     iv.  Volver al paso ii y repetir el proceso.
#' 
#'     Para obtener los $k$ (número par) dígitos centrales de $x_{i}^{2}$
#'     se puede utilizar
#'     que:$$x_{i+1}=\left\lfloor \left(  x_{i}^{2}-\left\lfloor \dfrac{x_{i}^{2}}{10^{(2k-\frac{k}{2})}}\right\rfloor 10^{(2k-\frac{k}{2})}\right)
#'     /10^{\frac{k}{2}}\right\rfloor$$ El algoritmo está **implementado
#'     en el fichero RANDVN.R**. Estudiar las características del
#'     generador de cuadrados medios a partir de una secuencia de 500
#'     valores. Emplear únicamente métodos gráficos.
#' 
#' 7.  Considerando el generador congruencial multiplicativo de parámetros
#'     $a=7^{5}=16807$, $c=0$ y $m=2^{31}-1$. ¿Se observan los mismos problemas 
#'     que con el algoritmo RANDU al considerar las tripletas $(x_{k},x_{k+1},x_{k+2})$?
#' 
#' 
#' ------
#' 
#' Ejercicio para entregar
#' -----------------------
#' 
#' 8.  Considera el generador congruencial definido por: $$\begin{aligned}
#'     x_{n+1}  & =(65x_{n}+1)\ \operatorname{mod}\ 2048,\nonumber\\
#'     u_{n+1}  & =\frac{x_{n+1}}{2048},\ n=0,1,\dots\nonumber\end{aligned}$$
#'     
#'     a)  Indicar razonadamente si es de ciclo máximo.
#'     
#'     b)  Generar 1000 valores tomando como semilla inicial el n de grupo 
#'         multiplicado por 100 y obtener el tiempo de CPU. Representar
#'         gráficamente el ajuste a la densidad teórica y realizar el
#'         correspondiente contraste de Kolmogorov-Smirnov.
#' 
#'     c)  Representar los pares de datos $\left( u_{i}, u_{i+1} \right)$, 
#'         ¿se observa algún problema?.
#' 
#'     d)  Estudiar la aleatoriedad de este generador empleando
#'         repetidamente el test de Ljung-Box, considerando 500 pruebas con
#'         muestras de tamaño 50 y hasta el salto 10
#'         (`Box.test(u,lag=10, type=Ljung)`). Comparar el
#'         ajuste de las distribuciones del estadístico y $p$-valor a las
#'         de referencia.
#' 
