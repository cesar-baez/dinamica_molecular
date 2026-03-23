# dinamica_molecular Ejemplo de un programa de dinámica molecular basada en el algoritmo *Verlet*


Este programa se creó para la materia de simulaciones avanzadas del semestre septiembre 2025 -enero 2026 del Programa de la Maestría en Ciencias Físicas de la Unidad Académica de Física de la Universidad Autónoma de Zacatecas.

La intención de este programa es mostrar el uso de un programa de dinámica molecular, así como algunos de los algoritmos básicos para su funcionamiento.

Este código no contiene optimizaciones avanzadas. Además sólo se implementa el algoritmo directo para el cálculo de fuerzas.  Por estas dos razones este código  es principalmente de uso educativo.

En esta simulación se modela un fluido con interacciones Lennard-Jones bidimensional.

En este archivo se encuentra la siguiente información:

- Instrucciones de compilación.
- Uso de la simulación.
- Descripción de las variables del archivo variables_entrada.dat.
- Descripción de los archivos de salida.
- Descripción de los archivos que se compilan.

En este repostitorio se incluye una versión modificada del archivo tauss88.f90 de Allan Miller. Tiene 2 modificaciones. La primera es que solo se incluye el generador de números aleatorios sin el programa que realiza una prueba estadística del generador. La segunda modificación es la declaración de la variable dp al inicio del código.

El código original del archivo tauss88.f90 pueden encontrarse en el repositorio de generadores de números aleatorios para Fortran de Allan Miller: [https://wp.csiro.au/alanmiller/random.html](https://wp.csiro.au/alanmiller/random.html)



## Instrucciones de compilación

Los siguientes comandos son para el compilador gfortran, algunas flags cambian para cada compilador en específico.

El orden para compilar los archivos es el siguiente

```
gfortran mod_parametros.f90 taus88.f90 mod_aleatorio.f90 mod_inicializacion.f90 mod_fuerza_pasos.f90 main_simulacion.f90  -o dinamica_molecular
```

La flag -o dinamica_molecular (con o minúscula) es una opción que indica el nombre del programa ya compilado, en este caso el nombre del programa es **dinamica_molecular**.

Es buena idea usar la optimización -O2 y símbolos de backtrace (la optimización es con «O» mayúcsculo, los símbolos de backtrace perimiten indicar la línea del código que se ejecutaba cuando el programa falla). Esto se hace con el siguiente comando

```
gfortran -fbacktrace -O2 mod_parametros.f90 taus88.f90 mod_aleatorio.f90 mod_inicializacion.f90 mod_fuerza_pasos.f90 main_simulacion.f90  -o dinamica_molecular.
```

Es recomendable revisar las [flags de gfortran](https://gcc.gnu.org/onlinedocs/gfortran/Option-Summary.html) para usar otras opciones de compilación.


Link de las flags de Fortran


[https://gcc.gnu.org/onlinedocs/gfortran/Option-Summary.html](https://gcc.gnu.org/onlinedocs/gfortran/Option-Summary.html)



## Uso de la simulación

Una vez que se compila el código se recomienda crear una carpeta para correr el programa.

En la carpeta en la que se corre el programa, debe existir un archivo de nombre «**variables_entrada.dat**». Este archivo contiene algunas de las variables de entrada del programa. En la siguiente sección se encuentra la escripción de las variables que se encuentran en este archivo.


Al modificar este archivo, se modifican las variables de la simulación. Se puede usar cualquier editor de texto para modificar este archivo.


Una vez que se modifica el archivo, simplemente hay que correr el programa **dinamica_molecular** en la misma carpeta que el archivo **variables_entrada.dat**.


Para correr el programa en linux, se usa el siguiente código en la terminal en la misma carpeta donde se encuentra el programa,


```
./dinamica_molecular
```



## Descripción de las variables del archivo variables_entrada.dat

El archivo **variables_entrada.dat** es un archivo de texto simple que permite modificar las variables del programa de simulación sin necesitar de compilar cada que se requiera cambiar alguna de las variables.


El contenido de un archivo **variables_entrada.dat** es el siguiente,


```
&LISTA_FISICA
 NUM_PARTICULAS=900        ,
 LADO_CAJA=  30.000000000000000     ,
 PHI= 0.10000000000000001     ,
 EPSILON_LJ=  1.0000000000000000     ,
 DELTA_T=  1.0000000000000000E-003,
 MASA=  1.0000000000000000     ,
 TEMPERATURA_INI=  1.0000000000000000     ,
 /
&LISTA_SIMULACION
 INICIO_ORDENADO=T,
 SALIDA_ANIMACION=F,
 NUM_INTENTOS_EMPALMES=2000       ,
 NUM_PASOS_TERMALIZACION=100000     ,
 NUM_PASOS_PROMEDIOS=100000     ,
 DELTA_N_TERM=10         ,
 DELTA_N_PROM=10         ,
 INT_RADIO_CORTE=4          ,
 SEMILLA_1=-1         ,
 SEMILLA_2=-123       ,
 SEMILLA_3=-12345     ,
 ARCHIVO_TERMAL="observables_termalizacion.dat                                                                                                                         ",
 ARCHIVO_PROMEDIOS="observables_promedios.dat                                                                                                                             ",
 ARCHIVO_LOG="log_simulacion.dat                                                                                                                                    ",
 ARCHIVO_GR="g_r.dat                                                                                                                                               ",
 PREFIJO_CFG="./animacion/cfg_                                                                                                                                      ",
 /
&LISTA_ESTADISTICA
 NUM_R=1000       ,
 NUM_PASOS_SALIDA=300        ,
 DELTA_R=  1.0000000000000000E-002,
 /
```

A continuación se muestra la descripción y el tipo de las variables de este archivo

- **NUM_PARTICULAS**. Tipo: **integer 32**. 
	
	El número de partículas en la simulación. 
- **LADO_CAJA**. Tipo: **real 64**.	
	
	Las partículas se colocan en una caja cuadrada, este parámetro modifica la longitud del lado de la caja. 
- **PHI**. Tipo: **real 64**.
	
	Fracción de llenado.

	Para un sistema monodisperso en 2D de $N$ partículas de diámetro $\sigma$ en una caja cuadrada de lado $l$, la fracción de llenado es $\phi =\frac{ N \pi \sigma^{2}}{4 L}$.
- **EPSILON_LJ**. Tipo: **real 64**. 
	
	Profundidad del pozo de la interacción Lennard-Jones.
- **DELTA_T**. Tipo: **real 64**. 
	
	Tamaño del paso.
- **MASA**. Tipo: **real 64**. 
	
	El valor de la masa de todas las partículas.
- **TEMPERATURA_INI**. Tipo: **real 64**. 
	
	Valor al cual se escala la temperatura inicial.
- **INICIO_ORDENADO**. Tipo: **logical**. 
	
	Variable lógica, si es verdadera inicia con una estructura cuadrada cristalina y si es falsa el sistema inicia con partículas colocadas al azar.
- **SALIDA_ANIMACION**. Tipo: **logical**.
	
	Variable lógica, si es verdadera el programa escribe una secuencia de archivos con las posiciones de todas las partículas para crear una animación.
- **NUM_INTENTOS_EMPALMES**. Tipo: **integer 32**. 
	
	El número de intentos que se realizan para insertar a una partícula con una posición al azar al seleccionarse un inicio aleatorio.
- **NUM_PASOS_TERMALIZACION**. Tipo: **integer 32**. 
	
	El número de pasos durante la etapa de termalización.
- **NUM_PASOS_PROMEDIOS**. Tipo: **integer 32**. 
	
	El número de pasos para el registro de datos del cálculo de promedios.
- **DELTA_N_TERM**. Tipo: **integer 32**. 
	
	Número de pasos entre cada registro de datos en la etapa de termalización.
- **DELTA_N_PROM**. Tipo: **integer 32**. 
	
	Número de pasos entre cada registro de datos en la etapa del cálculo de promedios.
- **INT_RADIO_CORTE**. Tipo: **integer 32**. 
	
	Distancia del radio de corte de la interacción en unidades enteras del diámetro de las partículas.
- **SEMILLA_1**. Tipo: **integer 32**. 
	
	Semilla 1 para el generador de números aleatorios.
- **SEMILLA_2**.  Tipo: **integer 32**. 
	
	Semilla 2 para el generador de números aleatorios.
- **SEMILLA_3**.  Tipo: **integer 32**. 
	
	Semilla 3 para el generador de números aleatorios.
- **ARCHIVO_TERMAL** . Tipo: **character**. 

	Nombre del archivo con una tabla de observables escalares durante la etapa de termalización.
- **ARCHIVO_PROMEDIOS**. Tipo: **character**. 
	
	 Nombre del archivo con una tabla de observables escalares durante la etapa de termalización.
- **ARCHIVO_LOG**. Tipo: **character**. 
	
	Nombre del archivo con la salida del programa.
- **ARCHIVO_GR**. Tipo: **character**. 
	
	Nombre del archivo de la función de distribución por pares, $g(r)$.
- **PREFIJO_CFG**. Tipo: **character**.  
	
	Prefijo para la secuencia de los archivos de animación.
- **NUM_R**.   Tipo: **integer 32**. 
	
	Número de bines de la $g(r)$.
- **NUM_PASOS_SALIDA**.  Tipo: **integer 32**.  
	
	Número de archivos para la animación.
- **DELTA_R**. Tipo: **real 64**.
	
	Tamaño del $\Delta r$ para los bines de la $g(r)$.


## Descripción de los archivos de salida


Al terminar la simulación se crean los siguientes archivos. Los nombres de los archivos que se indican en mayúscula se modifican con el archivo de entrada.

- **01_variables_info.dat** Un archivo de ejemplo para el archivo **archivo_entrada.dat** con los valores de esas variables al momento de la compilación. Se crea al inicio de la simulación.
- **configuracion_inicial.dat** Un archivo XYZ con las posiciones iniciales de las partículas.
- **ARCHIVO_GR** El archivo de salida de la $g(r)$.
- **ARCHIVO_LOG** Un archivo con el registro de la simulación el cual es idéntico a la salida en pantalla.
- **ARCHIVO_PROMEDIOS**
- **ARCHIVO_TERMAL**





Ambos archivos, **ARCHIVO_TERMAL** y **ARCHIVO_PROMEDIOS**  contienen el mismo tipo de observables. Solo que uno es contiene observables de la etapa de termalización y la otra para la etapa de medición.

Ambos arhivos son de 4 columnas.

La primer columna es el tiempo y las otras 3 columnas son cantidades instantáneas (que sean instantáneas singnifica que solo se calculan para el paso en el que se reportan).

1. El tiempo en unidades de simulación.
2. Energía cinética media instantánea (En la simulación se indica como temperatura).
3. Energía de interacción media instantánea (Sin la corrección debido al radio de corte).
4. Virial medio instantáneo.


Para encontrar las energías y viriales de la simulación, se toma el promedio temporal de estas columnas. Es decir, se suman y se divide entre el número de pasos.



## Descripción de los archivos que se compilan

A continuación se describe el contenido de cada archivo que se compila, así como una descripción de las subrutinas que se encuentran en cada módulo.


- **mod_parametros.f90** Módulo en el que se definen las variables globales de la simulación
- **taus88.f90** Generador de números aleatorio Tausworth de ciclo $2^{88}$ del repositorio de [Allan Miller](https://wp.csiro.au/alanmiller/random.html). 
- **mod_aleatorio.f90** Contiene el algoritmo de Numerical Recipes, no se usa en esta simulación.
- **mod_inicializacion.f90** Módulo con las rutinas de inicialización de la simulación
	- Lectura de archivo externo
	- Escritura del registro de las variables de la simulación
	- Inicialización de variables
	- Inicialización de velocidades y posiciones de las partículas
	- Salida de la configuración inicial
- **mod_fuerza_pasos.f90**
	- Cálculo de fuerzas.
	- Algoritmo de Verlet.
	- Ciclo de los pasos de termalización.
	- Ciclo de los pasos para el registro de promedios.
	- Cálculo de la función de correlación por pares.
	- Sálida de configuraciones para la animación.
- **main_simulacion.f90**
	- Programa princial
Simulación de dinámica molecular con el algoritmo verlet para Fortran (En español)
