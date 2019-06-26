      module test_uhf
        use basis
        use uhf
        implicit none
        private
        real(8),parameter::TOL=1d-7
        public::test_main

      contains

        subroutine run_(mes,nea,neb,ian,r,v0)
          implicit none
          character(*),intent(in)::mes
          integer,intent(in)::nea,neb,ian(:)
          real(8),intent(in)::r(:,:),v0
          integer::i,nb
          real(8)::en,ee
          real(8),allocatable::
     &      ema(:),cma(:,:),emb(:),cmb(:,:),hb(:),sb(:),tb(:)
          logical::ok
          ! Prepare basis set.
          call basis_init()
          do i=1,SIZE(ian)
            call basis_set('STO3G',ian(i),i)
          end do
          ! Get memory.
          call basis_dim(nb)
          allocate(ema(nb))
          allocate(emb(nb))
          allocate(cma(nb,nb))
          allocate(cmb(nb,nb))
          ! Prepare AO integral matrices.
          call basis_enuc(r,ian,en)
          call basis_sij(r,sb)
          call basis_hij(r,ian,hb)
          call basis_ijkl(r,tb)
          cma(:,:)=0d0
          cmb(:,:)=0d0
          call uhf_scf(nea,neb,sb,hb,tb,cma,cmb,ema,emb)
          call uhf_eelc(nea,neb,hb,tb,cma,cmb,ee)
          ok=(ee+en-v0)**2<TOL
          call asserttrue_(mes,ok)
        end subroutine

      subroutine test_uhf_h2_()
!     --------------------------------------------------
!      $CONTRL SCFTYP=UHF RUNTYP=ENERGY MAXIT=30 MULT=1 $END
!      $SYSTEM MWORDS=1 $END
!      $BASIS GBASIS=STO NGAUSS=3 $END
!      $SCF DIRSCF=.TRUE. $END
!      $DATA
!     TITLE
!     C1
!      H           1.0  -0.3561016195   0.0000000000   0.0000000000
!      H           1.0   0.3561016195   0.0000000000   0.0000000000
!      $END
!      --------------------------------------------------
        implicit none
        character(LEN=15)::mes
        integer,parameter::NA=2,MULT=1
        integer::ian(NA),ne,nea,neb
        real(8):: r(3,NA),v0
        mes='H2'
        v0=-1.1175058844d0
        ian(:)=(/1,1/)
        ne=SUM(ian)
        nea=(ne+MULT-1)/2
        neb=(ne-MULT+1)/2
        r(:,:)=RESHAPE((/
     &      -0.6729344846d0,0.0000000000d0,0.0000000000d0,
     &       0.6729344846d0,0.0000000000d0,0.0000000000d0
     &      /),(/3,NA/))
        call run_(mes,nea,neb,ian,r,v0)
      end subroutine

        subroutine test_uhf_ch3_()
!       --------------------------------------------------
!        $CONTRL SCFTYP=UHF RUNTYP=ENERGY MAXIT=30 MULT=2 $END
!        $SYSTEM MWORDS=1 $END
!        $BASIS GBASIS=STO NGAUSS=3 $END
!        $SCF DIRSCF=.TRUE. $END
!        $DATA
!       Title
!        c1
!        C           6.0  -0.0000000013  -0.0000000000  -0.0000000000
!        H           1.0   0.5390728116  -0.9337014944   0.0000000000
!        H           1.0   0.5390728116   0.9337014944   0.0000000000
!        H           1.0  -1.0781456218   0.0000000000   0.0000000000
!        $END
!       --------------------------------------------------
          implicit none
          character(LEN=15)::mes
          integer,parameter::NA=4,MULT=2
          integer::ian(NA),ne,nea,neb
          real(8):: r(3,NA),v0
          mes='CH3'
          v0=-39.0767107618d0
          ian(:)=(/6,1,1,1/)
          ne=SUM(ian)
          nea=(ne+MULT-1)/2
          neb=(ne-MULT+1)/2
          r(:,:)=RESHAPE((/
     &        -0.0000000025d0,-0.0000000000d0,-0.0000000000d0,
     &         1.0186999014d0,-1.7644399787d0, 0.0000000000d0,
     &         1.0186999014d0, 1.7644399787d0, 0.0000000000d0,
     &        -2.0373998001d0, 0.0000000000d0, 0.0000000000d0
     &        /),(/3,NA/))
          call run_(mes,nea,neb,ian,r,v0)
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
          write(*,'("--- UHF ---")')
          call test_uhf_h2_()
          call test_uhf_ch3_()
        end subroutine

      end module

      program main
        use test_uhf,only: test_main
        implicit none
        call test_main()
      end program
