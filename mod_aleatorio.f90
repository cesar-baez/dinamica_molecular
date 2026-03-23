module numeros_aleatorios_gaussianos
use modulo_parametros, only: dp
use Ecuyer_random, only: taus88

private
public gasdev_s!, gasdev

contains

   FUNCTION gasdev()
   implicit none

   REAL(dp)                               ::  gasdev

   !INTEGER                                ::  idum

   !Rutina del numerical recipes 
   !C USES ran1
   !Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
   !as the source of uniform deviates.
   INTEGER,save                           ::  iset=0
   integer                                ::  i
   integer, parameter                     ::  num_intentos = 100
   REAL(dp)                               ::  fac, rsq, v1, v2
   real(dp), save                         ::  gset
   !SAVE iset,gset
    



   !if (idum.lt.0) iset=0                    !Reinitialize.
   if (iset.eq.0) then                       !We don’t have an extra deviate handy, so
      do i= 1, num_intentos
         v1=2.0_dp*taus88()-1.0_dp           ! pick two uniform numbers in the square extend-
         v2=2.0_dp*taus88()-1.0_dp           !ing from -1 to +1 in each direction,
      
         rsq=v1*v1+v2*v2                     !see if they are in the unit circle,
         !if(rsq.ge.1..or.rsq.eq.0.)goto 1   !and if they are not, try again.
         
         if( (rsq .gt. 0.0_dp) .and. (rsq .lt. 1.0_dp)  ) exit
         
         if(i==num_intentos) then
            write(*,*) 'Error números aleatorios gaussianos, se alcanzó el número máximo de intentos'
            write(*,*) 'num_intentos = ', num_intentos
            stop
         end if

      end do 
      
      fac=sqrt(-2.0*log(rsq)/rsq)       !Now make the Box-Muller transformation to get
      gset=v1*fac
                                       !two normal deviates. Return one and save
      gasdev=v2*fac
                                       !the other for next time.
      iset=1                           !Set flag.
   else                                ! We have an extra deviate handy,
      gasdev=gset                      !so return it,
      iset=0                           !and unset the flag.
   endif

   return
   END FUNCTION gasdev


   SUBROUTINE gasdev_s(harvest)
   !USE nrtype
   !USE nr, ONLY : ran1
   IMPLICIT NONE
   REAL(dp), INTENT(OUT)               :: harvest
   !Returns in harvest a normally distributed deviate with zero mean and unit variance, using
   !ran1 as the source of uniform deviates.
   REAL(dp)                            :: rsq,v1,v2
   REAL(dp), SAVE                      :: g
   LOGICAL, SAVE                       :: gaus_stored=.false.
   
   if (gaus_stored) then                  !We have an extra deviate handy,
      harvest=g                           !so return it,
      gaus_stored=.false.                 !and unset the flag.
   else                                   !We don’t have an extra deviate handy, so
   
   do
                                          !pick two uniform numbers in the square ex-
                                          !tending from -1 to +1 in each direction,
                                          
      v1=2.0_dp*taus88()-1.0_dp
      v2=2.0_dp*taus88()-1.0_dp

      rsq=v1*v1+v2*v2                     !see if they are in the unit circle,
      
      if (rsq > 0.0_dp .and. rsq < 1.0_dp) exit
   end do                                 !otherwise try again.
   
   rsq=sqrt(-2.0_dp*log(rsq)/rsq)         !Now make the Box-Muller transformation to
   harvest=v1*rsq
                                          !get two normal deviates. Return one and
   g=v2*rsq                               !save the other for next time.
   
   gaus_stored=.true. !Set flag.
   end if
   
   END SUBROUTINE gasdev_s


end module numeros_aleatorios_gaussianos

