      module test_uga
        use uga
        implicit none
        private
        public::test_main

      contains

        subroutine run_(mes,ne,nb,mult,ltab0)
          implicit none
          character(*),intent(in)::mes
          integer,intent(in)::ne,nb,mult,ltab0(:,:,:)
          logical::ok
          integer,allocatable::ltab(:,:,:)
          call uga_pltab(ne,nb,mult,ltab)
!         call pout_(ltab,ltab0)
          ok=sum((ltab(:,:,:)-ltab0(:,:,:))**2)==0
          call asserttrue_(mes,ok)
        end subroutine

!       subroutine pout_(ltab,ltab0)
!         implicit none
!         integer,intent(in)::ltab(:,:,:),ltab0(:,:,:)
!         integer::ic,nc,im,nm
!         nc=SIZE(ltab(1,1,:))
!         nm=SIZE(ltab(1,:,1))
!         do ic=1,nc
!           write(*,"('-',i2)") ic
!           do im=1,nm
!             write(*,"(3i1,x,3i1,x,2i1)")
!    &            ltab(:,im,ic),ltab0(:,im,ic),
!    &            sum((ltab(:,im,ic)-ltab0(:,im,ic))**2)
!           end do
!         end do
!       end subroutine

        subroutine test_uga_2e2o1_()
          implicit none
          integer,parameter::NC=3,NE=2,NB=2,MULT=1
          integer::ltab0(3,NB,NC)
          ltab0(:,:,:)=RESHAPE((/
     &      1,0,1, 1,0,0,
     &      1,0,1, 0,1,0,
     &      1,0,1, 0,0,1
     &      /),(/3,NB,NC/))
          call run_("2e2o1",NE,NB,MULT,ltab0)
        end subroutine

        subroutine test_uga_3e3o2_()
          implicit none
          integer,parameter::NC=8,NE=3,NB=3,MULT=2
          integer::ltab0(3,NB,NC)
          ltab0(:,:,:)=RESHAPE((/
     &      1,1,1, 1,1,0, 1,0,0,
     &      1,1,1, 1,1,0, 0,1,0,
     &      1,1,1, 1,0,1, 1,0,0,
     &      1,1,1, 1,0,1, 0,1,0,
     &      1,1,1, 1,0,1, 0,0,1,
     &      1,1,1, 0,2,0, 0,1,0,
     &      1,1,1, 0,1,1, 0,1,0,
     &      1,1,1, 0,1,1, 0,0,1
     &      /),(/3,NB,NC/))
          call run_("3e3o2",NE,NB,MULT,ltab0)
        end subroutine

        subroutine test_uga_4e4o1_()
          implicit none
          integer,parameter::NC=20,NE=4,NB=4,MULT=1
          integer::ltab0(3,NB,NC)
          ltab0(:,:,:)=0
          ltab0(:,:,:NC)=RESHAPE((/
     &      2,0,2, 2,0,1, 2,0,0, 1,0,0,
     &      2,0,2, 2,0,1, 1,1,0, 1,0,0,
     &      2,0,2, 2,0,1, 1,1,0, 0,1,0,
     &      2,0,2, 2,0,1, 1,0,1, 1,0,0,
     &      2,0,2, 2,0,1, 1,0,1, 0,1,0,
     &      2,0,2, 2,0,1, 1,0,1, 0,0,1,
     &      2,0,2, 1,1,1, 1,1,0, 1,0,0,
     &      2,0,2, 1,1,1, 1,1,0, 0,1,0,
     &      2,0,2, 1,1,1, 1,0,1, 1,0,0,
     &      2,0,2, 1,1,1, 1,0,1, 0,1,0,
     &      2,0,2, 1,1,1, 1,0,1, 0,0,1,
     &      2,0,2, 1,1,1, 0,2,0, 0,1,0,
     &      2,0,2, 1,1,1, 0,1,1, 0,1,0,
     &      2,0,2, 1,1,1, 0,1,1, 0,0,1,
     &      2,0,2, 1,0,2, 1,0,1, 1,0,0,
     &      2,0,2, 1,0,2, 1,0,1, 0,1,0,
     &      2,0,2, 1,0,2, 1,0,1, 0,0,1,
     &      2,0,2, 1,0,2, 0,1,1, 0,1,0,
     &      2,0,2, 1,0,2, 0,1,1, 0,0,1,
     &      2,0,2, 1,0,2, 0,0,2, 0,0,1
     &      /),(/3,NB,NC/))
          call run_("4e4o2",NE,NB,MULT,ltab0)
        end subroutine

        subroutine asserttrue_(mes,ok)
          implicit none
          character(*),intent(in)::mes
          logical,intent(in)::ok
          if (ok) then
            write(*,'(a,"...ok")') TRIM(mes)
          else
            write(*,'(a,"...failure")') TRIM(mes)
            stop 1
          end if
        end subroutine

        subroutine test_main()
          implicit none
          write(*,'("--- UGA ---")')
          call test_uga_2e2o1_()
          call test_uga_3e3o2_()
          call test_uga_4e4o1_()
        end subroutine

      end module

      program main
        use test_uga,only: test_main
        implicit none
        call test_main()
      end program
