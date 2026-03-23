module modulo_inicializacion_programa
use modulo_parametros, only: dp, pi, num_dimensiones, num_char,escritura_log,archivo_ejemplo
use Ecuyer_random, only: init_seeds, taus88
implicit none

private

public inicializacion, temperatura





integer, parameter                     ::  v = 777                  !Unidad para archivo de salida

!El cambio en el número de espacios de las strings de formato
!sirve para compenzar que un acento equivale a 2 caracteres
!en el conteo de espacios para el formato
                                                            !10+39+33+8=90
character(num_char), parameter         ::  formato_entero = '(10x,a39,33x,i8)'  !39+24=63
character(num_char), parameter         ::  formato_entero1 = '(10x,a40,33x,i8)' !40+24=64
character(num_char), parameter         ::  formato_entero2 = '(10x,a41,33x,i8)' !41+24=65

character(num_char), parameter         ::  formato_real = '(10x,a39,30x,f11.5)' !39+21=60v
character(num_char), parameter         ::  formato_real1 = '(10x,a40,30x,f11.5)'!40+21=61
character(num_char), parameter         ::  formato_real2 = '(10x,a40,31x,f11.5)'!40+22=62

character(num_char), parameter         ::  formato_string = '(10x,a39,2x,a39)'
character(num_char), parameter         ::  formato_string1 = '(10x,a40,2x,a39)'

character(num_char), parameter         ::  formato_logical = '(10x,a39,33x,l8)' 
character(num_char), parameter         ::  formato_logical1 = '(10x,a40,33x,l8)'  

!write(*,'(10x,a,39x,F11.5)')  'Densidad de phi (=1):',sqrt(variables%phi_t/pi)

contains

   subroutine inicializacion(n_prom, radio, radio_corte, densidad, &
               &  x, y, vel_x, vel_y, x_m1, y_m1, g_r, cfg_x, cfg_y, &
               &  num_particulas, phi, epsilon_lj,delta_t,masa,temperatura, inicio_ordenado, &
               &  salida_animacion, num_intentos_empalmes, num_pasos_termalizacion, num_pasos_promedios, &
               &  delta_n_term, delta_n_prom, int_radio_corte, lado_caja, &
               &  semilla_1, semilla_2, semilla_3, archivo_cfg_inicial, archivo_termal,&
               &  archivo_promedios, archivo_log, archivo_gr, prefijo_cfg, num_r, &
               &  num_pasos_salida, delta_r)
   implicit none
   
   logical, intent(inout)              ::  inicio_ordenado
   logical, intent(inout)              ::  salida_animacion

   integer, intent(inout)              ::  num_particulas
   
   integer, intent(inout)              ::  num_intentos_empalmes
   integer, intent(inout)              ::  num_pasos_termalizacion
   integer, intent(inout)              ::  num_pasos_promedios

   integer, intent(inout)              ::  delta_n_term
   integer, intent(inout)              ::  delta_n_prom

   integer, intent(inout)              ::  int_radio_corte
   integer, intent(inout)              ::  semilla_1, semilla_2, semilla_3
   
   integer, intent(inout)              ::  num_r
   integer, intent(inout)              ::  num_pasos_salida

   real(dp), intent(inout)             ::  phi
   real(dp), intent(inout)             ::  epsilon_lj
   real(dp), intent(inout)             ::  delta_t
   real(dp), intent(inout)             ::  masa
   real(dp), intent(inout)             ::  temperatura

   real(dp), intent(inout)             ::  lado_caja

   real(dp), intent(inout)             ::  delta_r

   character(num_char), intent(inout)  ::  archivo_cfg_inicial, archivo_termal, archivo_promedios
   character(num_char), intent(inout)  ::  archivo_log, archivo_gr, prefijo_cfg

   integer, intent(out)                ::  n_prom
   real(dp), intent(out)               ::  radio, radio_corte, densidad
   
   real(dp), intent(out), &
   & allocatable, dimension(:)         ::  x, y, vel_x, vel_y, x_m1, y_m1, g_r

   real(dp), intent(out), &
   & allocatable, dimension(:,:)       ::  cfg_x, cfg_y

   integer                             ::  total_intentos




   call lectura_inicializacion()


   call info_lectura(num_particulas, phi, epsilon_lj,delta_t,masa,temperatura, inicio_ordenado, &
                  &  salida_animacion, num_intentos_empalmes, num_pasos_termalizacion, &
                  &  num_pasos_promedios, delta_n_term, delta_n_prom, int_radio_corte, lado_caja, &
                  &  semilla_1, semilla_2, semilla_3, archivo_termal, archivo_promedios, &
                  &  archivo_log, archivo_gr, prefijo_cfg, num_r, num_pasos_salida, delta_r)




   call  inicializacion_variables_simulacion(n_prom, radio, radio_corte, densidad, x, y, &
                                           & vel_x, vel_y, x_m1, y_m1, g_r, cfg_x, cfg_y, &
                                           & num_particulas, num_pasos_promedios, delta_n_prom, &
                                           & int_radio_corte, semilla_1, semilla_2, semilla_3, &
                                           & num_r, phi, lado_caja)

   call inicio_simulacion(total_intentos, x, y, vel_x, vel_y, inicio_ordenado, num_particulas,  &
                        & num_intentos_empalmes, masa, lado_caja, radio, archivo_cfg_inicial)



   call info_inicializacion(num_particulas, n_prom, total_intentos, masa, radio, radio_corte, densidad, &
                           &vel_x, vel_y, archivo_log)

   end subroutine inicializacion

   subroutine lectura_inicializacion()
   use modulo_parametros, only: lista_fisica,lista_simulacion,lista_estadistica, &
                                & variables_momento_compilacion, archivo_entrada_variables
   implicit none



   if (archivo_ejemplo==1) then
      open(7001,file=variables_momento_compilacion,status='unknown',delim='quote')
      write(7001,nml=lista_fisica)
      write(7001,nml=lista_simulacion)
      write(7001,nml=lista_estadistica)
      close(7001)
      write(*,*) 'Se escribió un archivo de entrada de ejemplo'
      write(*,*) 'Este archivo puede usarse en vez del archivo de entrada'


      write(*,*) 'Nombre del archivo de ejemplo:'
      write(*,*) trim(variables_momento_compilacion)            
      write(*,*)

   end if




   open(7002,file=archivo_entrada_variables,status='old',delim='quote')
   read(7002,nml=lista_fisica)
   read(7002,nml=lista_simulacion)
   read(7002,nml=lista_estadistica)
   close(7002)

   end subroutine lectura_inicializacion

! Esta rutina rompe la convención enteros, flotantes, vectores, etc.
! para modificarse simultáneamente con el módulo parámetros

   subroutine info_lectura(num_particulas, phi, epsilon_lj,delta_t,masa,temperatura, inicio_ordenado, &
                        &  salida_animacion, num_intentos_empalmes, num_pasos_termalizacion, &
                        &  num_pasos_promedios, delta_n_term, delta_n_prom, int_radio_corte, &
                        &  lado_caja, semilla_1, semilla_2, semilla_3, archivo_termal, &
                        &  archivo_promedios, archivo_log, archivo_gr, prefijo_cfg, &
                        &  num_r, num_pasos_salida, delta_r)
   implicit none

   logical, intent(in)                 ::  inicio_ordenado
   logical, intent(in)                 ::  salida_animacion
   

   integer, intent(in)                 ::  num_particulas
   
   integer, intent(in)                 ::  num_intentos_empalmes
   integer, intent(in)                 ::  num_pasos_termalizacion
   integer, intent(in)                 ::  num_pasos_promedios

   integer, intent(in)                 ::  delta_n_term
   integer, intent(in)                 ::  delta_n_prom

   integer, intent(in)                 ::  int_radio_corte
   integer, intent(in)                 ::  semilla_1, semilla_2, semilla_3
   
   integer, intent(in)                 ::  num_r
   integer, intent(in)                 ::  num_pasos_salida

   real(dp), intent(in)                ::  phi
   real(dp), intent(in)                ::  epsilon_lj
   real(dp), intent(in)                ::  delta_t
   real(dp), intent(in)                ::  masa
   real(dp), intent(in)                ::  temperatura

   real(dp), intent(in)                ::  lado_caja

   real(dp), intent(in)                ::  delta_r

   character(num_char), intent(in)     ::  archivo_termal, archivo_promedios, archivo_log
   character(num_char), intent(in)     ::  archivo_gr, prefijo_cfg

   integer                             ::  horas(8)
   

   character(8)                        ::  fechas
   character(5)                        ::  zona
   character(10)                       ::  hora

   character(num_char)                 ::  string_temp, string_hora




   

   call  date_and_time(fechas,hora,zona,horas)
   write(string_hora,'(10x,5(I2,1x,a,1x),I4)') horas(5), 'h', horas(6), 'm', horas(7), 's del', &
                                             & horas(3), 'de', horas(2), 'de', horas(1) 


   write(*,*) trim(string_hora)


   string_temp = 'Inicio de la simulación'
   write(*,*) trim(string_temp)


   



   if (escritura_log) then
      open(v,file=archivo_log,status='unknown')
      write(v,*) trim(string_temp)  
      write(v,*) trim(string_hora)
      write(v,*)  
      close(v)
   end if



   write(*,formato_real  ) 
   
   write(*,*)
   write(*,*)
   write(*,*) '** Parámetros físicos del sistema **'
   write(*,*)
   write(*,formato_entero2) 'Número de partículas:', num_particulas
   write(*,formato_real1  ) 'Lado de la caja de simulación:', lado_caja
   write(*,formato_real1  ) 'Fracción de llenado:', phi
   write(*,formato_real   ) 'Epsilon, pozo Lennard-Jones:', epsilon_lj
   write(*,formato_real1  ) 'Tamaño del paso:', delta_t
   write(*,formato_real1  ) 'Masa de las partículas:', masa
   write(*,formato_real   ) 'Temperatura:', temperatura

   !------------------------------------------------------

   write(*,*)
   write(*,*)
   write(*,*) '** Parámetros de simulación **'
   write(*,*)

   write(*,formato_logical )'Inicio ordenado:', inicio_ordenado
   write(*,formato_logical1)'Archivos para animación:', salida_animacion
   write(*,formato_entero2) 'Máximo de intentos al colocar partículas:', num_intentos_empalmes
   write(*,formato_entero2) 'Número de pasos de termalización:', num_pasos_termalizacion
   write(*,formato_entero1) 'Número de pasos para toma de promedios:', num_pasos_promedios
   write(*,*)
   write(*,formato_entero1) 'Pasos entre toma de datos termalización:', delta_n_term
   write(*,formato_entero ) 'Pasos entre toma de datos para promedios:', delta_n_prom
   write(*,formato_entero1) 'Número de vecinos para radio de corte:', int_radio_corte
   write(*,*)
   write(*,formato_entero ) 'Semilla 1:', semilla_1
   write(*,formato_entero ) 'Semilla 2:', semilla_2
   write(*,formato_entero ) 'Semilla 3:', semilla_3
   write(*,*)
   write(*,formato_string1) 'Archivo observables de termalización:', trim(archivo_termal)
   write(*,formato_string ) 'Archivo de observables para promedios:', trim(archivo_promedios)
   write(*,formato_string ) 'Archivo para g(r):', trim(archivo_gr)
   write(*,formato_string ) 'Prefijo archivos de configuraciones:', trim(prefijo_cfg)
   write(*,formato_string ) 'Archivo del log:', trim(archivo_log)

   !------------------------------------------------------

   write(*,*)
   write(*,*)
   write(*,*) '**Parámetros para gráficas y estadística **'
   write(*,*)

   write(*,formato_entero1) 'Número de bines histogramas r:', num_r
   write(*,formato_entero1) 'Número de archivos de configuraciones:', num_pasos_salida
   write(*,formato_real1  ) 'Delta r para el tamaño del bin:',delta_r
   write(*,*)
   write(*,*)


   if (escritura_log) then
      open(v,file=archivo_log,status='unknown',position='append')

   
      write(v,*)
      write(v,*)
      write(v,*) '** Parámetros físicos del sistema **'
      write(v,*)
      write(v,formato_entero2) 'Número de partículas:', num_particulas
      write(v,formato_real1   ) 'Lado de la caja de simulación:', lado_caja
      write(v,formato_real1  ) 'Fracción de llenado:', phi
      write(v,formato_real  ) 'Epsilon, pozo Lennard-Jones:', epsilon_lj
      write(v,formato_real1  ) 'Tamaño del paso:', delta_t
      write(v,formato_real1  ) 'Masa de las partículas:', masa
      write(v,formato_real  ) 'Temperatura:', temperatura

      !------------------------------------------------------

      write(v,*)
      write(v,*)
      write(v,*) '** Parámetros de simulación **'
      write(v,*)

      write(v,formato_logical) 'Inicio ordenado:', inicio_ordenado
      write(v,formato_logical1)'Archivos para animación:', salida_animacion
      write(v,formato_entero2) 'Intentos máximos al colocar partículas:', num_intentos_empalmes
      write(v,formato_entero2) 'Número de pasos de termalización:', num_pasos_termalizacion
      write(v,formato_entero1) 'Número de pasos para toma de promedios:', num_pasos_promedios
      write(v,*)
      write(v,formato_entero1) 'Pasos entre toma de datos termalización:', delta_n_term
      write(v,formato_entero ) 'Pasos entre toma de datos para promedios:', delta_n_prom
      write(v,formato_entero1) 'Número de vecinos para radio de corte:', int_radio_corte
      
      write(v,*)
      write(v,formato_entero ) 'Semilla 1:', semilla_1
      write(v,formato_entero ) 'Semilla 2:', semilla_2
      write(v,formato_entero ) 'Semilla 3:', semilla_3
      write(v,*)
      write(v,formato_string1) 'Archivo observables de termalización:', trim(archivo_termal)
      write(v,formato_string ) 'Archivo de observables para promedios:', trim(archivo_promedios)
      write(v,formato_string ) 'Archivo para g(r):', trim(archivo_gr)
      write(v,formato_string ) 'Prefijo archivos de configuraciones:', trim(prefijo_cfg)
      write(v,formato_string ) 'Archivo del log:', trim(archivo_log)

      !------------------------------------------------------

      write(v,*)
      write(v,*)
      write(v,*) '** Parámetros para gráficas y estadística **'
      write(v,*)

      write(v,formato_entero1) 'Número de bines histogramas r:', num_r
      write(v,formato_entero1) 'Número de archivos de configuraciones:', num_pasos_salida
      write(v,formato_real1  ) 'Delta r para el tamaño del bin:',delta_r
      write(v,*)
      write(v,*)

      close(v)
   end if 


   end subroutine info_lectura

   subroutine inicializacion_variables_simulacion(n_prom, radio, radio_corte, densidad, &
                                                & x, y, vel_x, vel_y, x_m1, y_m1, &
                                                & g_r, cfg_x, cfg_y, &
                                                & num_particulas, num_pasos_promedios, &
                                                & delta_n_prom, int_radio_corte, &
                                                & semilla_1, semilla_2, semilla_3, num_r ,&
                                                & phi, lado_caja)
   implicit none

   integer, intent(in)                 ::  num_particulas, num_pasos_promedios, delta_n_prom
   integer, intent(in)                 ::  int_radio_corte
   integer, intent(in)                 ::  semilla_1, semilla_2, semilla_3
   integer, intent(in)                 ::  num_r

   real(dp), intent(in)                ::  phi, lado_caja

   !Salida

   integer, intent(out)                ::  n_prom
   
   real(dp), intent(out)               ::  radio, radio_corte, densidad
   
   real(dp), intent(out) &
   &, allocatable, dimension(:)        ::  x, y , vel_x, vel_y, x_m1, y_m1, g_r

   real(dp), intent(out) &
   &, allocatable, dimension(:,:)      ::  cfg_x, cfg_y

   !Variables locales

   integer                             ::  estado_allocate

   real(dp)                            ::  area




   call init_seeds(semilla_1, semilla_2, semilla_3)

   n_prom = int( real(num_pasos_promedios, dp)/real(delta_n_prom, dp) )

   area = lado_caja * lado_caja

   densidad = real(num_particulas,dp)/area

   radio = sqrt(phi/(densidad*pi) )

   radio_corte = int_radio_corte*radio*2.0_dp


   allocate(x(num_particulas),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('x')
      stop
   end if 

   allocate(y(num_particulas),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('y')
      stop
   end if 


   allocate(vel_x(num_particulas),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('vel_x')
      stop
   end if 


   allocate(vel_y(num_particulas),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('vel_y')
      stop
   end if 

   allocate(x_m1(num_particulas),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('x_m1')
      stop
   end if 

   allocate(y_m1(num_particulas),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('y_m1')
      stop
   end if 

   allocate(g_r(num_r),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('g_r')
      stop
   end if 

   allocate(cfg_x(num_particulas, n_prom),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('cfg_x')
      stop
   end if 

   allocate(cfg_y(num_particulas, n_prom),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('cfg_y')
      stop
   end if 

   x = 0.0_dp
   y = 0.0_dp
   vel_x = 0.0_dp
   vel_y = 0.0_dp

   end subroutine inicializacion_variables_simulacion

   
   function hay_empalmes(num_particulas, lado_caja, radio, x, y)
   implicit none

   integer, intent(in)                 ::  num_particulas

   real(dp), intent(in)                ::  lado_caja, radio
   real(dp), intent(in), &
   & dimension(num_particulas)         ::  x, y

   logical                             ::  hay_empalmes


   integer                             ::  i,j

   real(dp)                            ::  sigma, sigma_2, delta_x, delta_y, r_ij_2
   real(dp)                            ::  x_temp, y_temp, inv_lado
   logical                             ::  empalme

   inv_lado = 1.0_dp / lado_caja
   sigma = 2.0_dp * radio
   sigma_2 = sigma * sigma
   empalme = .false.

   ciclo_i: do i = 1, num_particulas-1
      x_temp = x(i)
      y_temp = y(i)
      ciclo_j: do j = i+1, num_particulas
         delta_x = x_temp - x(j)
         delta_y = y_temp - y(j)

         delta_x = delta_x - lado_caja * idnint(delta_x*inv_lado)
         delta_y = delta_y - lado_caja * idnint(delta_y*inv_lado) 

         r_ij_2 = delta_x*delta_x+ delta_y*delta_y

         if (r_ij_2 .le. sigma_2) then
            empalme = .true.
            exit ciclo_j
         end if

      end do ciclo_j

      if (empalme) exit ciclo_i 
   end do ciclo_i


   hay_empalmes = empalme

   end function hay_empalmes

   function temperatura(num_particulas,masa, vel_x,vel_y)
   implicit none

   integer, intent(in)                 ::  num_particulas

   real(dp), intent(in)                ::  masa

   real(dp), intent(in), &
   & dimension(num_particulas)         ::  vel_x, vel_y


   real(dp)                            ::  temperatura

   real(dp)                            ::  inv_num, cinetica_media

   real(dp), parameter                 ::  inv_dim = 1.0_dp/real(num_dimensiones,dp)

   inv_num = 1.0_dp/real(num_particulas*num_dimensiones,dp)

   cinetica_media = 0.5_dp *inv_num * masa* sum(vel_x*vel_x + vel_y*vel_y)

   temperatura = 2.0_dp* inv_dim* cinetica_media


   end function temperatura


   subroutine remover_empalmes(total_intentos, x, y, num_particulas, num_intentos_empalmes, &
                              &lado_caja, radio)
   implicit none

   integer, intent(in)                 ::  num_particulas, num_intentos_empalmes

   real(dp), intent(in)                ::  lado_caja, radio

   integer, intent(out)                ::  total_intentos

   real(dp), intent(out), &
   & dimension(num_particulas)         ::  x, y

   integer                             ::  i, j, k, intento,resto_intentos

   real(dp)                            ::  sigma, sigma_2, delta_x, delta_y, r_ij_2
   real(dp)                            ::  x_temp, y_temp
   real(dp)                            ::  medio_lado, ran_x_temp, ran_y_temp

   real(dp)                            ::  inv_lado

   logical                             ::  empalme

   inv_lado = 1.0_dp / lado_caja
   
   medio_lado = lado_caja*0.5_dp 

   sigma = 2.0_dp * radio
   sigma_2 = sigma * sigma
   empalme = .false.

   total_intentos = 0

   ciclo_i: do i = 1, num_particulas-1
      x_temp = x(i)
      y_temp = y(i)
      ciclo_j: do j = i+1, num_particulas
         delta_x = x_temp - x(j)
         delta_y = y_temp - y(j)

         delta_x = delta_x - lado_caja * idnint(delta_x*inv_lado)
         delta_y = delta_y - lado_caja * idnint(delta_y*inv_lado) 

         r_ij_2 = delta_x*delta_x+ delta_y*delta_y

         if (r_ij_2 .le. sigma_2) then

            

            resto_intentos = num_intentos_empalmes- total_intentos

            if ( resto_intentos .le. 0 ) then
               write(*,*) 'Se acabaron los intentos al remover empalmes.'
               stop
            end if

            ciclo_intentos: do intento = 1, resto_intentos
               ran_x_temp = taus88()*lado_caja - medio_lado
               x(j) = ran_x_temp

               ran_y_temp = taus88()*lado_caja - medio_lado
               y(j) = ran_y_temp

               total_intentos = total_intentos + 1

               ciclo_k: do k = 1, j-1

                  delta_x = ran_x_temp - x(k)
                  delta_y = ran_y_temp - y(k)

                  delta_x = delta_x - lado_caja * idnint(delta_x*inv_lado)
                  delta_y = delta_y - lado_caja * idnint(delta_y*inv_lado)

                  r_ij_2 = delta_x*delta_x+ delta_y*delta_y

                  if (r_ij_2 .le. sigma_2) then
                     cycle ciclo_intentos
                  end if

               end do ciclo_k

               exit ciclo_intentos

            end do ciclo_intentos

            resto_intentos = num_intentos_empalmes- total_intentos

            if ( resto_intentos .le. 0 ) then
               write(*,*) 'Se acabaron los intentos al remover empalmes.'
               stop
            end if
            !empalme = .true.
            
         end if

      end do ciclo_j

      !if (empalme) exit ciclo_i 
   end do ciclo_i

   write(*,*) 'Se removieron los empalmes con',total_intentos , ' intentos'

   end subroutine remover_empalmes

   subroutine inicio_simulacion(total_intentos, x, y, vel_x, vel_y, inicio_ordenado, num_particulas, &
                              & num_intentos_empalmes, masa, lado_caja, radio,archivo_cfg_inicial)   
   implicit none

   logical, intent(in)                 ::  inicio_ordenado

   integer, intent(in)                 ::  num_particulas, num_intentos_empalmes

   real(dp), intent(in)                ::  masa, lado_caja

   real(dp), intent(in)                ::  radio

   character(num_char), intent(in)     ::  archivo_cfg_inicial

   integer, intent(out)                ::  total_intentos

   real(dp), intent(out), &
   &dimension(num_particulas)          ::  x, y, vel_x, vel_y

   integer                             ::  i, j, int_lado_caja, i_particula

   real(dp)                            ::  medio_lado, r_temp, inv_num,temperatura_ini
   real(dp)                            ::  factor_escala_vel, v_cm_x, v_cm_y

   logical                             ::  empalmes

   inv_num = 1.0_dp/real(num_particulas,dp)!*num_dimensiones,dp)

   medio_lado = lado_caja*0.5_dp 

   int_lado_caja = int(lado_caja)


   if (inicio_ordenado) then

      do i_particula = 1, num_particulas

         i = mod(i_particula, int_lado_caja)
         j = floor(real(i_particula)/real(int_lado_caja))

         x(i_particula) = real(i,dp)
         y(i_particula) = real(j,dp)



      end do

      call salida_cfg_inicial(num_particulas, lado_caja, radio, x, y, archivo_cfg_inicial)
      empalmes = hay_empalmes(num_particulas, lado_caja, radio, x, y)

      if (empalmes) then
         write(*,*) 'Hay empalmes en el inicio ordenado, rutina: inicio_simulacion'
         stop
      else
         write(*,*) 'Se iniciaron las partículas en un arreglo ordenado'
      end if

   else
      do i = 1, num_particulas
         r_temp = taus88()*lado_caja - medio_lado
         x(i) = r_temp
      end do

      do i = 1, num_particulas
         r_temp = taus88()*lado_caja - medio_lado
         y(i) = r_temp
      end do

      !call histograma()
      call salida_cfg_inicial(num_particulas, lado_caja, radio, x, y, archivo_cfg_inicial)
      empalmes = hay_empalmes(num_particulas, lado_caja, radio, x, y)

      if (empalmes) then
         call remover_empalmes(total_intentos, x, y, num_particulas, num_intentos_empalmes, &
                              &lado_caja,radio)
      end if

      empalmes = hay_empalmes(num_particulas, lado_caja, radio, x, y)

      if (empalmes) then
         write(*,*) 'No se removieron los empalmes, rutina: inicio_simulacion'
         stop
      else
         write(*,*) 'Se removieron los empalmes de forma exitosa'
      end if
   
   end if

   

   !call histograma()

   do i = 1, num_particulas
      r_temp = taus88()*2.0_dp - 1.0_dp
      vel_x(i) = r_temp
   end do

   do i = 1, num_particulas
      r_temp = taus88()*2.0_dp - 1.0_dp
      vel_y(i) = r_temp
   end do

   v_cm_x = sum(vel_x)*inv_num
   v_cm_y = sum(vel_y)*inv_num

   vel_x = (vel_x - v_cm_x)
   vel_y = (vel_y - v_cm_y)!*factor_escala_vel

   temperatura_ini = temperatura(num_particulas,masa,vel_x,vel_y)
   factor_escala_vel = sqrt(1.0_dp/temperatura_ini)

   vel_x = vel_x*factor_escala_vel
   vel_y = vel_y*factor_escala_vel

   end subroutine inicio_simulacion


   subroutine info_inicializacion(num_particulas, n_prom, total_intentos, masa,radio,radio_corte,densidad, &
                                 &vel_x, vel_y, archivo_log)
   implicit none

   integer, intent(in)                 ::  num_particulas, n_prom, total_intentos

   real(dp), intent(in)                ::  masa, radio, radio_corte, densidad
   real(dp), intent(in), &
   & dimension(num_particulas)         ::  vel_x,vel_y

   character(num_char), intent(in)     ::  archivo_log

   integer                             ::  bytes_memoria_cfg_dp

   real(dp)                            ::  temperatura_ini, inv_num, vel_cm_x, vel_cm_y
   real(dp)                            ::  mega_bytes_memoria_cfg_dp
   character(num_char)                 ::  formato_local = '(10x,a39,31x,i10)'

   inv_num = 1.0_dp/real(num_particulas,dp)

   temperatura_ini = temperatura(num_particulas,masa,vel_x,vel_y)

   vel_cm_x = sum(vel_x)*inv_num
   vel_cm_y = sum(vel_x)*inv_num

   bytes_memoria_cfg_dp = 16*num_particulas*n_prom

   mega_bytes_memoria_cfg_dp = real(bytes_memoria_cfg_dp,dp)/real(1024*1024,dp)

   write(*,*) '** Variables de la simulación **'
   write(*,*)
   write(*,*)

   write(*,formato_entero1) 'Número de configuraciones en promedios:', n_prom
   write(*,formato_local) 'Bytes de configuraciones promedios:', bytes_memoria_cfg_dp
   write(*,formato_real   ) 'MB las configuraciones para promedios:', mega_bytes_memoria_cfg_dp
   write(*,*)   
   write(*,formato_real1  ) 'Radio de las partículas:', radio
   write(*,formato_real   ) 'Radio de corte:', radio_corte
   write(*,formato_real1  ) 'Densidad numérica:', densidad
   

   write(*,*) '** Velocidades al iniciar **'
   write(*,*)
   write(*,*)

   write(*,formato_real   ) 'Temperatura:', temperatura_ini

   write(*,formato_real   ) 'Velocidad del centro de masa, x:', vel_cm_x
   write(*,formato_real   ) 'Velocidad del centro de masa, y:', vel_cm_y
   write(*,*)
   write(*,*)

   if (escritura_log) then
      open(v,file=archivo_log,status='unknown', position='append')
      write(v,*)
      write(v,*) 

      write(v,*) '              Se removieron los empalmes con',total_intentos , ' intentos'
      write(v,*)
      write(v,*)

      write(v,*) '** Variables de la simulación **'
      write(v,*)
      write(v,*)

      write(v,formato_entero1) 'Número de configuraciones en promedios:', n_prom
      write(v,formato_local) 'Bytes de configuraciones promedios:', bytes_memoria_cfg_dp
      write(v,formato_real   ) 'MB las configuraciones para promedios:', mega_bytes_memoria_cfg_dp
      write(v,*)
      write(v,formato_real1  ) 'Radio de las partículas:', radio
      write(v,formato_real   ) 'Radio de corte:', radio_corte
      write(v,formato_real1  ) 'Densidad numérica:', densidad

      write(v,*) '** Velocidades al iniciar **'
      write(v,*)
      write(v,*)

      write(v,formato_real   ) 'Temperatura:', temperatura_ini

      write(v,formato_real   ) 'Velocidad del centro de masa, x:', vel_cm_x
      write(v,formato_real   ) 'Velocidad del centro de masa, y:', vel_cm_y
      write(v,*)
      write(v,*)
      close(v)
   end if

   end subroutine info_inicializacion


   subroutine salida_cfg_inicial(num_particulas, lado_caja, radio, x, y, archivo_cfg_inicial)
   implicit none

   integer, intent(in)                 ::  num_particulas

   real(dp), intent(in)                ::  lado_caja, radio

   real(dp), intent(in), &
   & dimension(num_particulas)         ::  x, y

   character(num_char), intent(in)     ::  archivo_cfg_inicial

   integer                             ::  i, u

   real(dp)                            :: inv_lado, x_i, y_i


   u = 819

   inv_lado = 1.0_dp/lado_caja

   


   open(u,file=archivo_cfg_inicial,status='replace')
   do i = 1, num_particulas
      
      x_i = x(i)- lado_caja * idnint(x(i)*inv_lado)
      y_i = y(i)- lado_caja * idnint(y(i)*inv_lado)

      write(u,'(2(f11.5,4x), f11.5)') x_i, y_i, radio
   end do
   close(u)

   

   end subroutine salida_cfg_inicial


end module modulo_inicializacion_programa