      module test_rhf
        use basis
        use rhf
        implicit none
        private
        public::test_main
        real(8),parameter::TOL_=1d-7

      contains

        subroutine run_(mes,ne,ian,r,v0)
          implicit none
          character(*),intent(in)::mes
          integer,intent(in)::ne,ian(:)
          real(8),intent(in)::r(:,:),v0
          logical::ok
          integer::i,nb
          real(8)::en,ee
          real(8),allocatable::em(:),cm(:,:),hb(:),sb(:),tb(:)
          ! Prepare basis set.
          call basis_init()
          do i=1,SIZE(ian)
            call basis_set('STO3G',ian(i),i)
          end do
          ! Get memory.
          call basis_dim(nb)
          allocate(em(nb))
          allocate(cm(nb,nb))
          ! Prepare AO integral matrices.
          call basis_enuc(r,ian,en)
          call basis_sij(r,sb)
          call basis_hij(r,ian,hb)
          call basis_ijkl(r,tb)
          ! SCF calculation.
          cm(:,:)=0d0
          call rhf_scf(ne,sb,hb,tb,cm,em)
          call rhf_eelc(ne,hb,tb,cm,ee)
          ! Check the result.
          ok=(ee+en-v0)**2<TOL_
          call asserttrue_(mes,ok)
        end subroutine

        subroutine test_rhf_lih()
!       --------------------------------------------------
!        $CONTRL SCFTYP=RHF RUNTYP=ENERGY $END
!        $SYSTEM MWORDS=1 $END
!        $BASIS GBASIS=STO NGAUSS=3 $END
!        $SCF DIRSCF=.T. $END
!        $DATA
!        HCP
!       C1
!        LI          3.0   1.5107517  0.0000000   0.0000000
!        H           1.0   0.0000000  0.0000000   0.0000000
!        $END
!       --------------------------------------------------
          implicit none
          character(LEN=15)::mes
          integer,parameter::NA=2
          integer::ne,ian(NA)
          real(8):: r(3,NA),v0
          mes='LiH'
          v0=-7.8633821393d0
          ian(:)=(/3,1/)
          r(:,:)=0d0
          r(1,1)=2.8549067485d0
          ne=SUM(ian(:na))
          call run_(mes,ne,ian,r,v0)
        end subroutine

        subroutine test_rhf_h2o()
!       --------------------------------------------------
!        $CONTRL SCFTYP=RHF RUNTYP=ENERGY $END
!        $SYSTEM MWORDS=1 $END
!        $BASIS GBASIS=STO NGAUSS=3 $END
!        $SCF DIRSCF=.T. $END
!        $DATA
!        HCP
!       C1
!        H           1.0   0.7581836745  -0.0380971262  -0.0000000000
!        H           1.0  -0.7581836745  -0.0380971263  -0.0000000000
!        O           8.0  -0.0000000000  -0.6738057475   0.0000000000
!        $END
!        --------------------------------------------------
          implicit none
          character(LEN=15)::mes
          integer,parameter::NA=3
          integer::ne,ian(NA)
          real(8):: r(3,NA),v0
          mes='H2O'
          v0=-74.9659002670d0
          ian(:)=(/1,1,8/)
          r(:,:)=RESHAPE((/
     &         1.4327593932d0,-0.0719931294d0,-0.0000000000d0,
     &        -1.4327593932d0,-0.0719931296d0,-0.0000000000d0,
     &        -0.0000000000d0,-1.2733082317d0, 0.0000000000d0
     &        /),(/3,NA/))
          ne=SUM(ian(:))
          call run_(mes,ne,ian,r,v0)
        end subroutine

        subroutine test_rhf_nh3()
!       --------------------------------------------------
!        $CONTRL SCFTYP=RHF RUNTYP=ENERGY $END
!        $SYSTEM MWORDS=1 $END
!        $BASIS GBASIS=STO NGAUSS=3 $END
!        $SCF DIRSCF=.T. $END
!        $DATA
!        HCP
!       C1
!        H           1.0   0.8145530804  -0.0904105230  -0.4702453377
!        H           1.0  -0.8145530804  -0.0904105231  -0.4702453377
!        H           1.0   0.0000000000  -0.0907625405   0.9405791527
!        N           7.0   0.0000000000  -0.5165135398  -0.0000884773
!        $END
!        --------------------------------------------------
          implicit none
          character(LEN=15)::mes
          integer,parameter::NA=4
          integer::ne,ian(NA)
          real(8):: r(3,NA),v0
          mes='NH3'
          v0=-55.4554197963d0
          ian(:)=(/1,1,1,7/)
          r(:,:)=RESHAPE((/
     &         1.5392821244d0,-0.1708511149d0,-0.8886348353d0,
     &        -1.5392821244d0,-0.1708511151d0,-0.8886348353d0,
     &         0.0000000000d0,-0.1715163315d0, 1.7774368684d0,
     &         0.0000000000d0,-0.9760690592d0,-0.0001671979d0
     &        /),(/3,NA/))
          ne=SUM(ian(:))
          call run_(mes,ne,ian,r,v0)
        end subroutine

        subroutine test_rhf_bef2()
!       --------------------------------------------------
!        $CONTRL SCFTYP=RHF RUNTYP=ENERGY $END
!        $SYSTEM MWORDS=1 $END
!        $BASIS GBASIS=STO NGAUSS=3 $END
!        $SCF DIRSCF=.T. $END
!        $DATA
!        HCP
!       C1
!        BE          4.0   0.0000000000   0.0000000000   0.0000000000
!        F           9.0   1.3287209000   0.0000000000   0.0000000000
!        F           9.0  -1.3287209000   0.0000000000   0.0000000000
!        $END
!        --------------------------------------------------
          implicit none
          character(LEN=15)::mes
          integer,parameter::NA=3
          integer::ne,ian(NA)
          real(8):: r(3,NA),v0
          mes='BeF2'
          v0=-210.6444730935d0
          ian(:)=(/4,9,9/)
          r(:,:)=RESHAPE((/
     &         0.0000000000d0, 0.0000000000d0, 0.0000000000d0,
     &         2.5109184152d0, 0.0000000000d0, 0.0000000000d0,
     &        -2.5109184152d0, 0.0000000000d0, 0.0000000000d0
     &        /),(/3,NA/))
          ne=SUM(ian(:))
          call run_(mes,ne,ian,r,v0)
        end subroutine

        subroutine test_rhf_c3h8()
!       --------------------------------------------------
!        $CONTRL SCFTYP=RHF RUNTYP=ENERGY $END
!        $SYSTEM MWORDS=1 $END
!        $BASIS GBASIS=STO NGAUSS=3 $END
!        $SCF DIRSCF=.T. $END
!        $DATA
!        HCP
!       C1
!        C           6.0  -0.0969810098   0.5473580347  -0.0801205231
!        C           6.0  -0.0969842682   1.8287574600  -0.9364263397
!        C           6.0  -0.0969871575  -0.7339672629  -0.9364876570
!        H           1.0   0.7788500676   0.5470433671   0.5663371674
!        H           1.0  -0.9728256336   0.5474134912   0.5663281929
!        H           1.0   0.7826809940   1.8659390647  -1.5721146546
!        H           1.0  -0.9765928064   1.8656907154  -1.5724247880
!        H           1.0  -0.9766089490  -0.7713886447  -1.5722372967
!        H           1.0   0.7825578104  -0.7711457550  -1.5725634068
!        H           1.0  -0.0970130335   2.7115364301  -0.3038663574
!        H           1.0  -0.0969860240  -1.6166069007  -0.3037816658
!        $END
!        --------------------------------------------------
          implicit none
          character(LEN=15)::mes
          integer,parameter::NA=11
          integer::ne,ian(NA)
          real(8):: r(3,NA),v0
          mes='C3H8'
          v0=-116.8864221176d0
          ian(:)=(/6,6,6,1,1,1,1,1,1,1,1/)
          r(:,:)=RESHAPE((/
     &        -0.1832675345d0, 1.0343567028d0,-0.1514058347d0,
     &        -0.1832736920d0, 3.4558504974d0,-1.7695891897d0,
     &        -0.1832791520d0,-1.3869970108d0,-1.7697050626d0,
     &         1.4718132133d0, 1.0337620672d0, 1.0702220630d0,
     &        -1.8383738813d0, 1.0344615004d0, 1.0702051037d0,
     &         1.4790526145d0, 3.5261135421d0,-2.9708659185d0,
     &        -1.8454928057d0, 3.5256442299d0,-2.9714519856d0,
     &        -1.8455233108d0,-1.4577131685d0,-2.9710976784d0,
     &         1.4788198312d0,-1.4572541735d0,-2.9717139372d0,
     &        -0.1833280506d0, 5.1240608586d0,-0.5742241524d0,
     &        -0.1832770100d0,-3.0549440722d0,-0.5740641085d0
     &        /),(/3,NA/))
          ne=SUM(ian(:))
          call run_(mes,ne,ian,r,v0)
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
          write(*,'("--- RHF ---")')
          call test_rhf_lih()
          call test_rhf_h2o()
          call test_rhf_nh3()
          call test_rhf_bef2()
          call test_rhf_c3h8()
        end subroutine

      end module

      program main
        use test_rhf,only: test_main
        implicit none
        call test_main()
      end program
