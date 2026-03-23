module pasos
use modulo_parametros, only: dp, pi, num_dimensiones, num_char,escritura_log
use modulo_inicializacion_programa, only: temperatura
implicit none

private
public  termalizacion, configuraciones_para_promedios, correlacion_por_pares, salida_cfgx




integer, parameter                     ::  v = 888                  !Unidad para archivo de salida
!character(num_char)                    ::  formato_banner = '(i8,2x,5(f11.5,5x),a25)'
character(num_char)                    ::  formato_banner  = '(6(f11.5,5x),a20)'
character(num_char)                    ::  formato_banner2 = '(4(f11.5,5x),a20)'


contains

   subroutine fuerza(u_total, virial, f_x, f_y, num_particulas, epsilon_lj, lado_caja, radio, &
                  & radio_corte, x, y)
   implicit none

   integer, intent(in)                 ::  num_particulas

   real(dp), intent(in)                ::  epsilon_lj, lado_caja, radio, radio_corte

   real(dp), intent(in), &
   & dimension(num_particulas)         ::  x , y

   real(dp), intent(out)               ::  u_total, virial

   real(dp), intent(out), &
   & dimension(num_particulas)         ::  f_x, f_y


   integer                             ::  i, j

   real(dp)                            ::  inv_lado, r_corte_2, x_i, y_i, delta_x, delta_y, r_ij_2
   real(dp)                            ::  sigma_2, r_2_inv, sigma_radio_2, sigma_radio_6
   real(dp)                            ::  epsilon_4, epsilon_24, f_ij_x, f_ij_y,u_lj, f_ij
   real(dp)                            ::  virial_temp, inv_num


   inv_lado = 1.0_dp/lado_caja
   r_corte_2 = radio_corte * radio_corte

   sigma_2= 4.0_dp * radio * radio

   epsilon_4 = 4.0_dp * epsilon_lj
   epsilon_24 = 24.0_dp * epsilon_lj

   inv_num = 1.0_dp /real(num_particulas,dp)

   u_total = 0.0_dp
   virial = 0.0_dp


   f_x = 0.0_dp
   f_y = 0.0_dp

   do i = 1, num_particulas-1
      x_i = x(i)
      y_i = y(i)

      do j = i+1, num_particulas
         delta_x = x_i - x(j)
         delta_y = y_i - y(j)

         delta_x = delta_x - lado_caja * idnint(delta_x *inv_lado )
         delta_y = delta_y - lado_caja * idnint(delta_y *inv_lado )

         r_ij_2 = delta_x*delta_x +delta_y*delta_y

         if (r_ij_2.lt.r_corte_2) then
            r_2_inv = 1.0_dp /r_ij_2

            sigma_radio_2 = sigma_2*r_2_inv
            sigma_radio_6 = sigma_radio_2* sigma_radio_2 *sigma_radio_2
            
            u_lj = sigma_radio_6 * (sigma_radio_6 - 1.0_dp)! * epsilon_4
            u_total = u_total + u_lj

            virial_temp = sigma_radio_6 * (2.0_dp * sigma_radio_6 - 1.0_dp)! * epsilon_24 

            virial = virial + virial_temp

            f_ij = r_2_inv * virial_temp

            f_ij_x = f_ij * delta_x
            f_ij_y = f_ij * delta_y

            f_x(i) =f_x(i) + f_ij_x
            f_y(i) =f_y(i) + f_ij_y

            f_x(j) =f_x(j) - f_ij_x
            f_y(j) =f_y(j) - f_ij_y

         end if

      end do

   end do

   virial = virial * epsilon_24* inv_num   
   u_total = u_total* epsilon_4*inv_num

   f_x = f_x * epsilon_24
   f_y = f_y * epsilon_24


   end subroutine fuerza

   subroutine paso_verlet(x_nueva, y_nueva ,vel_x_nueva, vel_y_nueva,&
                        & num_particulas, delta_t, masa, x, y, f_x, f_y, x_m1, y_m1)
   implicit none

   integer, intent(in)                 ::  num_particulas

   real(dp), intent(in)                ::  delta_t, masa

   real(dp), intent(in), &
   & dimension(num_particulas)         ::  x, y, f_x, f_y, x_m1, y_m1

   real(dp), intent(out), &
   & dimension(num_particulas)         ::  x_nueva,y_nueva, vel_x_nueva,vel_y_nueva
   
   real(dp)                            ::  m_inv, factor

   m_inv = 1.0_dp /masa

   factor = 1.0_dp /(2.0_dp * delta_t)

   x_nueva = 2.0_dp* x - x_m1 + f_x* m_inv* delta_t * delta_t

   y_nueva = 2.0_dp* y - y_m1 + f_y* m_inv* delta_t * delta_t



   vel_x_nueva = (x_nueva - x_m1)* factor 
   vel_y_nueva = (y_nueva - y_m1)* factor 

   end subroutine

   subroutine termalizacion(x_m1, y_m1, num_particulas, num_pasos_termalizacion, delta_n_term, &
                        &   epsilon_lj, delta_t, masa, lado_caja, radio, radio_corte, x, y, &
                        &   vel_x, vel_y, archivo_termal, archivo_log)

   implicit none

   integer, intent(in)                 ::  num_particulas, num_pasos_termalizacion,delta_n_term

   real(dp), intent(in)                ::  epsilon_lj, delta_t, masa, lado_caja, radio, radio_corte

   real(dp), intent(inout), &
   & dimension(num_particulas)         ::  x, y, vel_x, vel_y

   character(num_char), intent(in)     ::  archivo_termal, archivo_log

   real(dp), intent(out), &
   & dimension(num_particulas)         ::  x_m1 ,y_m1


   integer                             ::  diez_porciento_nt, i, nt_total, nt_t, estado_allocate

   real(dp)                            ::  num_p_inv, temperatura_i, v_cm_x, v_cm_y, u_total, virial
   real(dp)                            ::  porcentaje,inv_nt_total
   real(dp)                            ::  temperatura_promedio, virial_promedio, u_promedio

   real(dp), dimension(num_particulas) ::  f_x, f_y, x_nueva, y_nueva
   real(dp), dimension(num_particulas) ::  vel_x_nueva,vel_y_nueva

   real(dp), allocatable, &
   & dimension(:)                      ::  tiempo, temperatura_t, virial_t, u_t   

   integer                             ::  horas(8)
   

   character(8)                        ::  fechas
   character(5)                        ::  zona
   character(10)                       ::  hora

   character(num_char)                 ::  string_hora,string_banner


   nt_t = ceiling(real(num_pasos_termalizacion,dp)/real(delta_n_term,dp))

   allocate(tiempo(nt_t),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('tiempo')
      stop
   end if 

   allocate(temperatura_t(nt_t),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('temperatura_t')
      stop
   end if 

   allocate(virial_t(nt_t),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('virial_t')
      stop
   end if 

   allocate(u_t(nt_t),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('u_t')
      stop
   end if 


   nt_total = 0

   diez_porciento_nt = int ( real(num_pasos_termalizacion,dp)*0.1)

   num_p_inv = 1.0_dp/real(num_particulas,dp)

   x_m1 = x - vel_x*delta_t
   y_m1 = y - vel_y*delta_t

   temperatura_promedio = 0.0_dp
   virial_promedio = 0.0_dp
   u_promedio = 0.0_dp





   string_banner = '    Avance     /          velocidad CM       /    Temperatura  /  U promedio  /     Virial    /            Ho'
   string_banner = trim(string_banner) // trim('ra        ')
                  !   10.00000         0.00000        -0.00000       390.71663        -0.65106        -0.09477     12 h 30 m 27 s del 16 de  7 de
   write(*,*) string_banner

   if (escritura_log) then
      open(v,file=archivo_log,status='unknown', position='append')
      write(v,*)
      write(v,*)
      write(v,*) string_banner
      write(v,*)
      write(v,*)
      close(v)
   end if

   !if (escritura_log) then
      !open(v,file=archivo_termal,status='replace')
      !close(v)
   !end if


   do i = 1, num_pasos_termalizacion
      call fuerza(u_total, virial, f_x, f_y, num_particulas, epsilon_lj, lado_caja, radio, &
               &  radio_corte, x, y)

      call paso_verlet(x_nueva, y_nueva ,vel_x_nueva ,vel_y_nueva, &
                     & num_particulas, delta_t, masa,&
                     & x, y, f_x, f_y, x_m1, y_m1)

      x_m1 = x
      y_m1 = y

      x = x_nueva
      y = y_nueva

      if (mod(i, delta_n_term).eq.0) then
         nt_total = nt_total + 1

         temperatura_i = temperatura(num_particulas,masa, vel_x_nueva,vel_y_nueva)
         !temperatura_promedio = temperatura_promedio + temperatura_i
         !u_promedio = u_promedio + u_total
         !virial_promedio = virial_promedio + virial

         tiempo(nt_total) = i * delta_t
         temperatura_t(nt_total) = temperatura_i
         u_t(nt_total) = u_total
         virial_t(nt_total) = virial

         !v_cm_x = sum(vel_x_nueva)*num_p_inv
         !v_cm_y = sum(vel_y_nueva)*num_p_inv
         !if (escritura_log) then
            !open(v,file=archivo_termal,status='unknown', position='append')                        
            !write(v,*) i, temperatura_i, u_total, virial                        
            !close(v)
         !end if
      end if 


      if (mod(i, diez_porciento_nt).eq. 0 ) then
         
         call  date_and_time(fechas,hora,zona,horas)
         
         !write(string_hora,'(5(I2,1x,a,1x),I4)') horas(5), 'h', horas(6), 'm', horas(7), &
         !                                          & 's del',  horas(3), 'de', horas(2), 'de', &
         !                                          & horas(1) 

         write(string_hora,'(5(I0.2,a),I4)') horas(5),':', horas(6), ':', horas(7), ' ',&
                                       &     horas(3), '/',horas(2), '/', horas(1) 

         temperatura_i = temperatura(num_particulas,masa, vel_x_nueva,vel_y_nueva)
         v_cm_x = sum(vel_x_nueva)*num_p_inv
         v_cm_y = sum(vel_y_nueva)*num_p_inv
         porcentaje = real(i,dp)/real(num_pasos_termalizacion,dp) * 100.0_dp

         write(*,formato_banner) porcentaje, v_cm_x,v_cm_y,temperatura_i, u_total, virial, trim(string_hora)

         if (escritura_log) then
            open(v,file=archivo_log,status='unknown', position='append')
            write(v,formato_banner) porcentaje, v_cm_x,v_cm_y,temperatura_i, u_total, virial, trim(string_hora)
            close(v)
         end if

      end if

   end do

   inv_nt_total = 1.0_dp /real(nt_total, dp)

   temperatura_promedio = sum(temperatura_t)*inv_nt_total
   u_promedio = sum(u_t)*inv_nt_total
   virial_promedio = sum(virial_t)*inv_nt_total

   if (escritura_log) then
      open(v,file=archivo_termal,status='replace')
      do i = 1, nt_total
         write(v,*) i*delta_n_term, temperatura_t(i), u_t(i), virial_t(i)
      end do
      close(v)
   end if


   end subroutine termalizacion

   subroutine configuraciones_para_promedios(cfg_x, cfg_y, & !temperatura_t, u_t, virial_t, &
                        &   num_particulas, num_pasos_promedios, delta_n_prom, n_prom, &
                        &   epsilon_lj, delta_t, masa, lado_caja, radio, radio_corte, x, y, &
                        &   vel_x, vel_y, x_m1 ,y_m1, archivo_promedios, archivo_log)

   implicit none

   integer, intent(in)                 ::  num_particulas, num_pasos_promedios,delta_n_prom, n_prom

   real(dp), intent(in)                ::  epsilon_lj, delta_t, masa, lado_caja, radio, radio_corte

   real(dp), intent(inout), &
   & dimension(num_particulas)         ::  x, y, vel_x, vel_y, x_m1 ,y_m1

   character(num_char), intent(in)     ::  archivo_promedios, archivo_log

   real(dp), intent(out), &
   & dimension(num_particulas,n_prom)  ::  cfg_x, cfg_y

   integer                             ::  diez_porciento_nt, i, nt_total, nt_t, estado_allocate

   real(dp)                            ::  num_p_inv, temperatura_i, v_cm_x, v_cm_y, u_total, virial
   real(dp)                            ::  porcentaje,inv_nt_total
   real(dp)                            ::  temperatura_promedio, virial_promedio, u_promedio

   real(dp), dimension(num_particulas) ::  f_x, f_y, x_nueva, y_nueva
   real(dp), dimension(num_particulas) ::  vel_x_nueva, vel_y_nueva

   real(dp), allocatable, &
   & dimension(:)                      ::  tiempo, temperatura_t, virial_t, u_t   

   integer                             ::  horas(8)
   

   character(8)                        ::  fechas
   character(5)                        ::  zona
   character(10)                       ::  hora

   character(num_char)                 ::  string_hora,string_banner


   nt_t = ceiling(real(num_pasos_promedios,dp)/real(delta_n_prom,dp))

   allocate(tiempo(nt_t),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('tiempo')
      stop
   end if 

   allocate(temperatura_t(nt_t),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('temperatura_t')
      stop
   end if 

   allocate(virial_t(nt_t),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('virial_t')
      stop
   end if 

   allocate(u_t(nt_t),stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo montar' // trim('u_t')
      stop
   end if 


   nt_total = 0

   diez_porciento_nt = int ( real(num_pasos_promedios,dp)*0.1)

   num_p_inv = 1.0_dp/real(num_particulas,dp)

   !x_m1 = x - vel_x*delta_t
   !y_m1 = y - vel_y*delta_t

   temperatura_promedio = 0.0_dp
   virial_promedio = 0.0_dp
   u_promedio = 0.0_dp





   string_banner = '    Avance     /    Temperatura  /  U promedio  /     Virial    /            Hora        '
                  !   10.00000         0.00000        -0.00000       390.71663        -0.65106        -0.09477     12 h 30 m 27 s del 16 de  7 de
   write(*,*) string_banner

   if (escritura_log) then
      open(v,file=archivo_log,status='unknown', position='append')
      write(v,*)
      write(v,*)
      write(v,*) string_banner
      write(v,*)
      write(v,*)
      close(v)
   end if

   !if (escritura_log) then
      !open(v,file=archivo_promedios,status='replace')
      !close(v)
   !end if


   do i = 1, num_pasos_promedios
      call fuerza(u_total, virial, f_x, f_y, num_particulas, epsilon_lj, lado_caja, radio, &
               &  radio_corte, x, y)

      call paso_verlet(x_nueva, y_nueva ,vel_x_nueva ,vel_y_nueva, &
                     & num_particulas, delta_t, masa,&
                     & x, y, f_x, f_y, x_m1, y_m1)

      x_m1 = x
      y_m1 = y

      x = x_nueva
      y = y_nueva

      if (mod(i, delta_n_prom).eq.0) then
         nt_total = nt_total + 1

         temperatura_i = temperatura(num_particulas,masa, vel_x_nueva,vel_y_nueva)
         !temperatura_promedio = temperatura_promedio + temperatura_i
         !u_promedio = u_promedio + u_total
         !virial_promedio = virial_promedio + virial

         tiempo(nt_total) = i * delta_t
         temperatura_t(nt_total) = temperatura_i
         u_t(nt_total) = u_total
         virial_t(nt_total) = virial

         cfg_x(:, nt_total) = x
         cfg_y(:, nt_total) = y

         !v_cm_x = sum(vel_x_nueva)*num_p_inv
         !v_cm_y = sum(vel_y_nueva)*num_p_inv
         !if (escritura_log) then
            !open(v,file=archivo_promedios,status='unknown', position='append')                        
            !write(v,*) i, temperatura_i, u_total, virial                        
            !close(v)
         !end if
      end if 


      if (mod(i, diez_porciento_nt).eq. 0 ) then
         
         call  date_and_time(fechas,hora,zona,horas)
         
         !write(string_hora,'(5(I2,1x,a,1x),I4)') horas(5), 'h', horas(6), 'm', horas(7), &
         !                                          & 's del',  horas(3), 'de', horas(2), 'de', &
         !                                          & horas(1) 

         write(string_hora,'(5(I0.2,a),I4)') horas(5),':', horas(6), ':', horas(7), ' ',&
                                       &     horas(3), '/',horas(2), '/', horas(1) 

         temperatura_i = temperatura(num_particulas,masa, vel_x_nueva,vel_y_nueva)
         v_cm_x = sum(vel_x_nueva)*num_p_inv
         v_cm_y = sum(vel_y_nueva)*num_p_inv
         porcentaje = real(i,dp)/real(num_pasos_promedios,dp) * 100.0_dp

         write(*,formato_banner2) porcentaje, temperatura_i, u_total, virial, trim(string_hora)

         if (escritura_log) then
            open(v,file=archivo_log,status='unknown', position='append')
            write(v,formato_banner2) porcentaje, temperatura_i, u_total, virial, trim(string_hora)
            close(v)
         end if

      end if

   end do

   inv_nt_total = 1.0_dp /real(nt_total, dp)

   temperatura_promedio = sum(temperatura_t)*inv_nt_total
   u_promedio = sum(u_t)*inv_nt_total
   virial_promedio = sum(virial_t)*inv_nt_total

   if (escritura_log) then
      open(v,file=archivo_promedios,status='replace')
      do i = 1, nt_total
         write(v,*) i*delta_n_prom, temperatura_t(i), u_t(i), virial_t(i)
      end do
      close(v)
   end if

   
   deallocate(tiempo, stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo desmontar' // trim('tiempo')
      stop
   end if 
   
   deallocate(temperatura_t, stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo desmontar' // trim('temperatura_t')
      stop
   end if 

   deallocate(u_t, stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo desmontar' // trim('u_t')
      stop
   end if 

   deallocate(virial_t, stat=estado_allocate)
   if (estado_allocate.ne.0) then
      write(*,*) 'No se pudo desmontar' // trim('virial_t')
      stop
   end if 

   end subroutine configuraciones_para_promedios

   subroutine correlacion_por_pares(num_particulas, n_prom, num_r, radio_corte, lado_caja, cfg_x, cfg_y, &
                                 &  archivo_gr)
   implicit none

   integer, intent(in)                 ::  num_particulas, n_prom, num_r

   real(dp), intent(in)                ::  radio_corte, lado_caja

   real(dp), intent(in), &
   & dimension(num_particulas, n_prom) ::  cfg_x, cfg_y

   character(num_char), intent(in)     ::  archivo_gr

   !real(dp), intent(out) &
   !& dimension(num_r)                  ::  r, g_r


   integer                             ::  i, j, i_t, bin

   real(dp)                            ::  factor, delta_x, delta_y, x_i, y_i, inv_lado
   real(dp)                            ::  delta_r, inv_delta, r_ij_2, r_ij,r_corte_2, da

   real(dp), dimension(num_particulas) ::  x_temp, y_temp
   real(dp), dimension(num_r)          ::  r, g_r

   inv_lado = 1.0_dp/lado_caja


   r = 0.0_dp

   r_corte_2 = radio_corte * radio_corte

   delta_r = radio_corte / real(num_r ,dp)
   inv_delta = 1.0_dp/delta_r

   g_r = 0.0_dp


   

   do i_t = 1, n_prom

      x_temp = cfg_x(:, i_t)
      y_temp = cfg_y(:, i_t)

      do i = 1, num_particulas - 1
         x_i = x_temp(i)
         y_i = y_temp(i)

         do j= i+1, num_particulas
            delta_x = x_i - x_temp(j)
            delta_y = y_i - y_temp(j)

            delta_x = delta_x - lado_caja * idnint(delta_x*inv_lado)
            delta_y = delta_y - lado_caja * idnint(delta_y*inv_lado)

            r_ij_2 = delta_x * delta_x + delta_y * delta_y

            if (r_ij_2.lt.r_corte_2) then
               r_ij = sqrt(r_ij_2)
               bin = int(r_ij *inv_delta) + 1
               g_r(bin) =g_r(bin) + 2.0_dp
            end if


         end do
      end do

   end do


   

   do i = 1, num_r
      r(i) = real(i,dp)  * delta_r
   end do

   do i = 1, num_r
      da = 2.0_dp*pi * r(i) * delta_r
      g_r(i) = g_r(i) /da
   end do

   factor = 1.0_dp/(real(num_particulas*n_prom ,dp))


   g_r = g_r * factor

   do i = 1, num_r
      r(i) = ( real(i,dp) -0.5_dp) * delta_r
   end do


   if (escritura_log) then
      open(v,file=archivo_gr,status='replace')
      do i = 1, num_r
         write(v,*) r(i), g_r(i)
      end do
      close(v)
   end if


   end subroutine correlacion_por_pares

   subroutine salida_cfgx(num_particulas, n_prom, num_pasos_salida, lado_caja, radio, &
                        & cfg_x, cfg_y, prefijo_cfg)
   implicit none

   integer, intent(in)                 ::  num_particulas, n_prom, num_pasos_salida

   real(dp), intent(in)                ::  lado_caja, radio

   real(dp), intent(in), &
   & dimension(num_particulas, n_prom) ::  cfg_x, cfg_y

   character(num_char), intent(in)     ::  prefijo_cfg

   integer                             ::  i,j

   real(dp)                            :: inv_lado, x_ji, y_ji

   character(num_char)                 ::  archivo_temp


   inv_lado = 1.0_dp/lado_caja

   

   do i = 1, num_pasos_salida
      write(archivo_temp, '(a, i0.4,a)') trim(prefijo_cfg), i, '.dat'
      open(v,file=archivo_temp,status='replace')
      do j = 1, num_particulas
         
         x_ji = cfg_x(j,i)- lado_caja * idnint(cfg_x(j,i)*inv_lado)
         y_ji = cfg_y(j,i)- lado_caja * idnint(cfg_y(j,i)*inv_lado)

         write(v,'(2(f11.5,4x), f11.5)') x_ji, y_ji, radio
      end do
      close(v)
   end do

   

   end subroutine salida_cfgx


end module pasos