      module timit
!$      use omp_lib
        implicit none
        private
        real(8)::wt0_=0d0
!       ----------------------------------------------------------------
        public::timit_init,timit_timit
!       ----------------------------------------------------------------

      contains

        subroutine timit_init()
          implicit none
!$        wt0_=omp_get_wtime()
        end subroutine

        subroutine timit_timit(mes)
          implicit none
          character(*),intent(in)::mes
          real(8)::wt
!$        wt=omp_get_wtime()
!$        write(*,"('timit:',a10,f15.6,' s')") mes,wt-wt0_
!$        wt0_=wt
        end subroutine

      end module
