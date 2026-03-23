module modulo_parametros
implicit none
!Variables globales

integer, parameter                     ::  dp = kind(1d0)
integer, parameter                     ::  num_dimensiones = 2
integer, parameter                     ::  num_char = 150

real(dp),parameter                     ::  pi = 3.141592653589793238462643383279     !pi

logical, parameter                     ::  escritura_log = .True.

integer, parameter                     ::  archivo_ejemplo = 1

!Parámetros Físicos

integer                                ::  num_particulas = 900
real(dp)                               ::  lado_caja = 30.0_dp
real(dp)                               ::  phi = 0.1_dp            !Fracción de llenado
real(dp)                               ::  epsilon_lj = 1.0_dp     !Pozo lennard Jones


real(dp)                               ::  delta_t = 0.001_dp

real(dp)                               ::  masa = 1.0_dp
real(dp)                               ::  temperatura_ini = 1.0_dp


!Simulación

logical                                ::  inicio_ordenado = .True.
logical                                ::  salida_animacion = .False.

integer                                ::  num_intentos_empalmes = 2000
integer                                ::  num_pasos_termalizacion = 100000
integer                                ::  num_pasos_promedios = 100000

integer                                ::  delta_n_term = 10
integer                                ::  delta_n_prom = 10

integer                                ::  int_radio_corte = 4


integer                                ::  semilla_1 = -1
integer                                ::  semilla_2 = -123
integer                                ::  semilla_3 = -12345

character(num_char)                    ::  archivo_cfg_inicial = 'configuracion_inicial.dat'
character(num_char)                    ::  archivo_termal = 'observables_termalizacion.dat'
character(num_char)                    ::  archivo_promedios = 'observables_promedios.dat'
character(num_char)                    ::  archivo_log = 'log_simulacion.dat'
character(num_char)                    ::  archivo_gr= 'g_r.dat'
character(num_char)                    ::  prefijo_cfg = './animacion/cfg_'


!Estadística

integer                                ::  num_r = 1000
integer                                ::  num_pasos_salida = 300 ! 10 s, a 30 fps
real(dp)                               ::  delta_r = 0.01_dp



character(num_char)                    ::  variables_momento_compilacion='01_variables_info.dat'
character(num_char)                    ::  archivo_entrada_variables='variables_entrada.dat'


namelist /lista_fisica/                num_particulas, lado_caja, phi,epsilon_lj,delta_t, masa, &
                                    &  temperatura_ini

namelist /lista_simulacion/            inicio_ordenado, salida_animacion, num_intentos_empalmes,  &
                                    &  num_pasos_termalizacion,num_pasos_promedios, delta_n_term, &
                                    &  delta_n_prom, int_radio_corte, &
                                    &  semilla_1, semilla_2, semilla_3, archivo_termal, &
                                    &  archivo_promedios, archivo_log, archivo_gr,  &
                                    &  prefijo_cfg 

namelist /lista_estadistica/           num_r, num_pasos_salida, delta_r


end module modulo_parametros