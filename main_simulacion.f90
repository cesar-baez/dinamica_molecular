program ejemplo_simulacion_fortran
use modulo_parametros
use modulo_inicializacion_programa
use pasos
implicit none

integer                                ::  n_prom

real(dp)                               ::  radio, radio_corte, densidad

real(dp), allocatable, dimension(:)    ::  x, y, vel_x, vel_y, x_m1, y_m1, g_r
real(dp), allocatable, dimension(:,:)  ::  cfg_x,cfg_y


!Inicialización del programa y de la simulación
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call inicializacion(n_prom, radio, radio_corte, densidad, &
               &  x, y, vel_x, vel_y, x_m1, y_m1, g_r, cfg_x, cfg_y, &
               &  num_particulas, phi, epsilon_lj,delta_t,masa,temperatura_ini, inicio_ordenado, &
               &  salida_animacion, num_intentos_empalmes, num_pasos_termalizacion, num_pasos_promedios, &
               &  delta_n_term, delta_n_prom, int_radio_corte, lado_caja, &
               &  semilla_1, semilla_2, semilla_3, archivo_cfg_inicial, archivo_termal, &
               &  archivo_promedios, archivo_log, archivo_gr, prefijo_cfg, num_r, &
               &  num_pasos_salida, delta_r)




!Termalización
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call termalizacion(x_m1, y_m1, num_particulas, num_pasos_termalizacion, delta_n_term, epsilon_lj, &
                &  delta_t, masa, lado_caja, radio, radio_corte, x, y, vel_x, vel_y, &
                &  archivo_termal, archivo_log)


!Recolección de mediciones
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call configuraciones_para_promedios(cfg_x, cfg_y, &!temperatura_t, u_t, virial_t &
                        &   num_particulas, num_pasos_promedios, delta_n_prom, n_prom, &
                        &   epsilon_lj, delta_t, masa, lado_caja, radio, radio_corte, x, y, &
                        &   vel_x, vel_y, x_m1 ,y_m1, archivo_promedios, archivo_log)



!Análisis estadístico
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call correlacion_por_pares(num_particulas, n_prom, num_r, radio_corte, lado_caja, cfg_x, cfg_y, &
                        &  archivo_gr)


!Salida para realizar animaciones
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (salida_animacion) then
    call salida_cfgx(num_particulas, n_prom, num_pasos_salida, lado_caja, radio, &
                            & cfg_x, cfg_y, prefijo_cfg)
end if

end program ejemplo_simulacion_fortran