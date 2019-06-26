      module test_basis
        use basis
        implicit none
        private
        public::test_main

      contains

        subroutine run_(ian,r,hb0,gbasis)
          implicit none
          integer,intent(in)::ian(:)
          real(8),intent(in)::r(:,:),hb0(:)
          character(*),intent(in)::gbasis
          real(8),parameter::TOL_=1d-7
          logical::ok
          integer::i
          real(8),allocatable::hb(:)
!         allocate(hb(SIZE(hb0)))
          call basis_init()
          do i=1,SIZE(ian)
            call basis_set(gbasis,ian(i),i)
          end do
          call basis_hij(r,ian,hb)
          ok=SUM((hb(:)-hb0(:))**2)<TOL_
          call asserttrue_(gbasis,ok)
        end subroutine

        subroutine test_basis_sto3g_()
!       --------------------------------------------------------
!        $CONTRL RUNTYP=ENERGY NPRINT=3 $END
!        $BASIS GBASIS=STO NGAUSS=3 $END
!        $DATA
!         TEST FOR STO3G
!         C1
!       H           1.0      0.0000      0.0000      0.000
!       HE          2.0      2.0000      0.0000      0.000
!       LI          3.0     -1.0000      1.0000      0.000
!        $END
!       --------------------------------------------------------
          implicit none
          integer,parameter::NA=3,NB2=28
          integer::ian(NA)
          real(8)::r(3,NA),hb0(NB2)
          ian=(/1,2,3/)
          r(:,:)=RESHAPE((/
     &      0.0000000000d0,0.0000000000d0,0.0000000000d0,
     &      3.7794519754d0,0.0000000000d0,0.0000000000d0,
     &     -1.8897259877d0,1.8897259877d0,0.0000000000d0
     &      /),(/3,NA/))
          hb0(:)=(/
     &      -2.111575d0,
     &      -0.188360d0,-2.698357d0,
     &      -0.488586d0,-0.000369d0,-5.114961d0,
     &      -1.025026d0,-0.167019d0,-1.155485d0,-1.749834d0,
     &      -0.868942d0,-0.258517d0,-0.022365d0,-0.186216d0,-1.573836d0,
     &       0.837231d0, 0.087396d0, 0.017149d0, 0.124114d0, 0.070907d0,
     &       -1.525424d0,
     &       0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &       0.000000d0, -1.466620d0
     &      /)
          call run_(ian,r,hb0,'STO3G')
        end subroutine

        subroutine test_basis_321g_()
!       --------------------------------------------------------
!        $CONTRL RUNTYP=ENERGY NPRINT=3 $END
!        $BASIS GBASIS=N21 NGAUSS=3 $END
!        $DATA
!         TEST FOR 321G
!         C1
!       H           1.0      0.0000      0.0000      0.000
!       HE          2.0      2.0000      0.0000      0.000
!       LI          3.0     -1.0000      1.0000      0.000
!        $END
!       --------------------------------------------------------
          implicit none
          integer,parameter::NA=3,NB2=91
          integer::ian(NA)
          real(8)::r(3,NA),hb0(NB2)
          ian=(/1,2,3/)
          r(:,:)=RESHAPE((/
     &      0.0000000000d0,0.0000000000d0,0.0000000000d0,
     &      3.7794519754d0,0.0000000000d0,0.0000000000d0,
     &     -1.8897259877d0,1.8897259877d0,0.0000000000d0
     &      /),(/3,NA/))
          hb0(:)=(/
     &     -1.823244d0,
     &     -1.556259d0,-2.034434d0,
     &     -0.001820d0,-0.155512d0,-2.340175d0,
     &     -0.071079d0,-0.369353d0,-2.112641d0,-2.167249d0,
     &     -0.127917d0,-0.769775d0,-0.000000d0,-0.000175d0,-5.13553d0,
     &     -0.678473d0,-1.259441d0,-0.041911d0,-0.116768d0,-0.89084d0,
     &     -1.73009d0,
     &     -0.682664d0,-0.924375d0,-0.112259d0,-0.255367d0,-0.01801d0,
     &     -0.18755d0, -1.55407d0,
     &      0.673217d0, 0.872627d0, 0.037452d0, 0.086745d0, 0.01382d0,
     &      0.12869d0,  0.07282d0, -1.50588d0,
     &      0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.00000d0,
     &      0.00000d0,  0.00000d0,  0.00000d0, -1.44510d0,
     &     -0.454073d0,-0.977676d0,-0.163588d0,-0.311215d0,-0.70212d0,
     &     -1.30537d0, -0.16254d0,  0.10235d0,  0.00000d0, -1.32352d0,
     &     -0.270624d0,-0.504947d0,-0.305868d0,-0.540128d0,-0.00277d0,
     &     -0.09805d0, -0.82017d0,  0.04060d0,  0.00000d0, -0.16264d0,
     &     -1.01203d0,
     &      0.266365d0, 0.464935d0, 0.102034d0, 0.181863d0, 0.00212d0,
     &      0.06054d0,  0.04060d0, -0.77562d0,  0.00000d0,  0.08435d0,
     &      0.05062d0, -0.92107d0,
     &      0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.00000d0,
     &      0.00000d0,  0.00000d0,  0.00000d0, -0.74615d0,  0.00000d0,
     &      0.00000d0,  0.00000d0, -0.89319d0
     &      /)
          call run_(ian,r,hb0,'321G')
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
          call test_basis_sto3g_()
          call test_basis_321g_()
        end subroutine

      end module

      program main
        use test_basis,only: test_main
        implicit none
        call test_main()
      end program
