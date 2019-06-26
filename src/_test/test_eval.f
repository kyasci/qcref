      module test_eval
!     --- Molecular integeral evaluation test ---
!     The refarence integral values are taken from the
!     output of GAMESS with the following input data.
!     NPRINT=3 is for 1e and NPRINT=4 is for 2e integrals.
!     --------------------------------------------------
!      $CONTRL SCFTYP=RHF RUNTYP=ENERGY NPRINT=? $END
!      $SYSTEM MWORDS=1 $END
!      $BASIS GBASIS=STO NGAUSS=3 NDFUNC=1 $END
!      $DATA
!      HCP
!     C1
!      H           1.0  -1.078134   0.000000  -0.000000
!      C           6.0   0.000000   0.000000   0.000000
!      P          15.0   1.478509   0.000000   0.000000
!      $END
!     --------------------------------------------------
        implicit none
        public::test_main

      contains

        subroutine setup_sto3g_(key,lp,ex,cf)
          implicit none
          integer,parameter::NC=3
          character(*),intent(in)::key
          integer,intent(out)::lp(NC)
          real(8),intent(out)::ex(NC),cf(NC)
          lp(:)=(/0,0,0/)
          if (key=='Hs') then
            ex(:)=(/3.4252509d0,0.6239137d0,0.1688554d0/)
            cf(:)=(/0.154328967295d0,0.535328142282d0,0.444634542185d0/)
          else if (key=='Cs') then
            ex(:)=(/2.9412494d0,0.6834831d0,0.2222899d0/)
            cf(:)=(/-0.099967229187d0,0.399512826089d0,0.70011546888d0/)
          else if (key(1:2)=='Cp') then
            ex(:)=(/2.9412494d0,0.6834831d0,0.2222899d0/)
            cf(:)=(/0.155916274999d0,0.607683718598d0,0.391957393099d0/)
          else if (key(1:2)=='Pp') then
            ex(:)=(/1.4618890d0,0.4078633d0,0.1596349d0/)
            cf(:)=(/0.010587604289d0,0.595167005266d0,0.462001011973d0/)
          else if (key(1:2)=='Pd') then
            ex(:)=(/0.55d0,1d0,1d0/)
            cf(:)=(/1d0,0d0,0d0/)
          else
            write(*,'("error: no such key: ",a)') key
            stop 1
          end if
          if (key(2:2)=='p') then
            if (key(3:3)=='x') lp(1)=1
            if (key(3:3)=='y') lp(2)=1
            if (key(3:3)=='z') lp(3)=1
          end if
          if (key(2:2)=='d') then
            if (key(3:3)=='x') lp(1)=1
            if (key(3:3)=='y') lp(2)=1
            if (key(3:3)=='z') lp(3)=1
            if (key(4:4)=='x') lp(1)=lp(1)+1
            if (key(4:4)=='y') lp(2)=lp(2)+1
            if (key(4:4)=='z') lp(3)=lp(3)+1
          end if
        end subroutine

        subroutine testdrv_overlap(mes,ic,basis,v0)
          ! Test driver for overlap integrals.
          use eval
          implicit none
          integer,parameter::NC=3
          real(8),parameter::TOL=1.d-8
          character(LEN=31),intent(in)::mes,basis(2)
          real(8),intent(in)::v0
          integer,intent(in)::ic(2)
          integer::la(3),lb(3)
          real(8)::val,r(3,3),c1(NC),c2(NC),e1(NC),e2(NC)
          logical::ok
          ! Default settings.
          r(:,:)=0d0
          r(1,1)=-2.0373778380d0
          r(1,3)=2.7939768804d0
          ! Set Atomic orbital basis set.
          call setup_sto3g_(basis(1),la,e1,c1)
          call setup_sto3g_(basis(2),lb,e2,c2)
          ! Test run.
          call eval_overlap(
     &      NC,NC,c1,c2,e1,e2,r(:,ic(1)),r(:,ic(2)),la,lb,val)
          ok=(val-v0)**2<TOL
          if (.not.ok) then
            call asserttrue_(mes,ok)
          end if
        end subroutine

        subroutine test_overlap_hshs()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_HsHs'
          ic(:)=(/1,1/)
          basis(:)=(/'Hs','Hs'/)
          v0=1d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_cpxcpx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_CpxCpx'
          ic(:)=(/2,2/)
          basis(:)=(/'Cpx','Cpx'/)
          v0=1d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_pdxxpdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_PdxxPdxx'
          ic(:)=(/3,3/)
          basis(:)=(/'Pdxx','Pdxx'/)
          v0=1d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_hscs()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_HsCs'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs','Cs'/)
          v0=0.497602d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_hscpx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_HsCpx'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs ','Cpx'/)
          v0=-0.471316d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_hscpy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_HsCpy'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs ','Cpy'/)
          v0=0d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_hspdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_HsPdxx'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdxx'/)
          v0=0.048076d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_hspdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_HsPdyy'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdyy'/)
          v0=0.015347d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_hspdxy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_HsPdxy'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdxy'/)
          v0=0d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_cpxppx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_CpxPpx'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx','Ppx'/)
          v0=-0.329314d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_cpxpdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_CpxPdxx'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx ','Pdxx'/)
          v0=0.463055d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_cpxpdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_CpxPdyy'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx ','Pdyy'/)
          v0=0.211329d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine test_overlap_pdxxpdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='overlap_PdxxPdyy'
          ic(:)=(/3,3/)
          basis(:)=(/'Pdxx','Pdyy'/)
          v0=0.3333333d0
          call testdrv_overlap(mes,ic,basis,v0)
        end subroutine

        subroutine testdrv_kinetic(mes,ic,basis,v0)
          ! Test driver for kinetic energy integrals.
          use eval
          implicit none
          integer,parameter::NC=3
          real(8),parameter::TOL=1.d-8
          character(LEN=31),intent(in)::mes,basis(2)
          real(8),intent(in)::v0
          integer,intent(in)::ic(2)
          integer::la(3),lb(3)
          real(8)::val,r(3,3),c1(NC),c2(NC),e1(NC),e2(NC)
          logical::ok
          ! Default settings.
          r(:,:)=0d0
          r(1,1)=-2.0373778380d0
          r(1,3)=2.7939768804d0
          ! Set Atomic orbital basis set.
          call setup_sto3g_(basis(1),la,e1,c1)
          call setup_sto3g_(basis(2),lb,e2,c2)
          ! Test run.
          call eval_kinetic(
     &      NC,NC,c1,c2,e1,e2,r(:,ic(1)),r(:,ic(2)),la,lb,val)
          ok=(val-v0)**2<TOL
          if (.not.ok) then
            call asserttrue_(mes,ok)
          end if
        end subroutine

        subroutine test_kinetic_hshs()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_HsHs'
          ic(:)=(/1,1/)
          basis(:)=(/'Hs','Hs'/)
          v0=0.760032d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_cpxcpx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_CpxCpx'
          ic(:)=(/2,2/)
          basis(:)=(/'Cpx','Cpx'/)
          v0=1.477728d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_pdxxpdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_PdxxPdxx'
          ic(:)=(/3,3/)
          basis(:)=(/'Pdxx','Pdxx'/)
          v0=1.191667d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_hscs()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_HsCs'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs','Cs'/)
          v0=0.110776d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_hscpx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_HsCpx'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs ','Cpx'/)
          v0=-0.268529d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_hscpy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_HsCpy'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs ','Cpy'/)
          v0=0d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_hspdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_HsPdxx'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdxx'/)
          v0=-0.014455d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_hspdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_HsPdyy'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdyy'/)
          v0=-0.008269d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_hspdxy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_HsPdxy'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdxy'/)
          v0=0d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_cpxppx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_CpxPpx'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx','Ppx'/)
          v0=-0.301168d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_cpxpdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_CpxPdxx'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx ','Pdxx'/)
          v0=0.440422d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_cpxpdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_CpxPdyy'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx ','Pdyy'/)
          v0=-0.003513d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine test_kinetic_pdxxpdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='kinetic_PdxxPdxy'
          ic(:)=(/3,3/)
          basis(:)=(/'Pdxx','Pdyy'/)
          v0=-0.091667d0
          call testdrv_kinetic(mes,ic,basis,v0)
        end subroutine

        subroutine testdrv_hcore(mes,ic,basis,v0)
          ! Test driver for kinetic energy integrals.
          use eval
          implicit none
          integer,parameter::NC=3
          real(8),parameter::TOL=1.d-8
          character(LEN=31),intent(in)::mes,basis(2)
          real(8),intent(in)::v0
          integer,intent(in)::ic(2)
          integer::la(3),lb(3),i,ian(3)
          real(8)::val,r(3,3),c1(NC),c2(NC),e1(NC),e2(NC),vt
          logical::ok
          ! Default settings.
          ian(:)=(/1,6,15/)
          r(:,:)=0d0
          r(1,1)=-2.0373778380d0
          r(1,3)=2.7939768804d0
          ! Set Atomic orbital basis set.
          call setup_sto3g_(basis(1),la,e1,c1)
          call setup_sto3g_(basis(2),lb,e2,c2)
          ! Test run.
          call eval_kinetic(
     &      NC,NC,c1,c2,e1,e2,r(:,ic(1)),r(:,ic(2)),la,lb,val)
          do i=1,3
            call eval_attract(
     &        NC,NC,c1,c2,e1,e2,r(:,ic(1)),r(:,ic(2)),la,lb,r(:,i),vt)
            val=val-ian(i)*vt
          end do
          ok=(val-v0)**2<TOL
          if (.not.ok) then
            call asserttrue_(mes,ok)
          end if
        end subroutine

        subroutine test_hcore_hshs()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_HsHs'
          ic(:)=(/1,1/)
          basis(:)=(/'Hs','Hs'/)
          v0=-6.449668d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_cpxcpx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_CpxCpx'
          ic(:)=(/2,2/)
          basis(:)=(/'Cpx','Cpx'/)
          v0=-10.217213d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_pdxxpdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_PdxxPdxx'
          ic(:)=(/3,3/)
          basis(:)=(/'Pdxx','Pdxx'/)
          v0=-11.194280
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_hscs()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_HsCs'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs','Cs'/)
          v0=-4.285071d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_hscpx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_HsCpx'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs ','Cpx'/)
          v0=3.368105d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_hscpy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_HsCpy'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs ','Cpy'/)
          v0=0d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_hspdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_HsPdxx'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdxx'/)
          v0=-0.618205d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_hspdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_HsPdyy'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdyy'/)
          v0=-0.178918d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_hspdxy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_HsPdxy'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdxy'/)
          v0=0d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_cpxppx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_CpxPpx'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx','Ppx'/)
          v0=4.238650d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_cpxpdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_CpxPdxx'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx ','Pdxx'/)
          v0=-5.814851d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_cpxpdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_CpxPdyy'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx ','Pdyy'/)
          v0=-2.602259d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine test_hcore_pdxxpdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          real(8)::v0
          mes='hcore_PdxxPdyy'
          ic(:)=(/3,3/)
          basis(:)=(/'Pdxx','Pdyy'/)
          v0=-4.046020d0
          call testdrv_hcore(mes,ic,basis,v0)
        end subroutine

        subroutine testdrv_repulse(mes,ic,basis,v0)
          ! Test driver for overlap integrals.
          use eval
          implicit none
          integer,parameter::NC=3
          real(8),parameter::TOL=1.d-8
          character(LEN=31),intent(in)::mes,basis(4)
          real(8),intent(in)::v0
          integer,intent(in)::ic(4)
          integer::la(3),lb(3),l3(3),l4(3)
          real(8)::val,r(3,3),
     &        c1(NC),c2(NC),c3(NC),c4(NC),e1(NC),e2(NC),e3(NC),e4(NC)
          logical::ok
          ! Default settings.
          r(:,:)=0d0
          r(1,1)=-2.0373778380d0
          r(1,3)=2.7939768804d0
          ! Set Atomic orbital basis set.
          call setup_sto3g_(basis(1),la,e1,c1)
          call setup_sto3g_(basis(2),lb,e2,c2)
          call setup_sto3g_(basis(3),l3,e3,c3)
          call setup_sto3g_(basis(4),l4,e4,c4)
          ! Test run.
          call eval_repulse(
     &        NC,NC,NC,NC,c1,c2,c3,c4,e1,e2,e3,e4,
     &        r(:,ic(1)),r(:,ic(2)),r(:,ic(3)),r(:,ic(4)),
     &        la,lb,l3,l4,val)
          ok=(val-v0)**2<TOL
          if (.not.ok) then
            call asserttrue_(mes,ok)
          end if
        end subroutine

        subroutine test_repulse_hshshshs()
          implicit none
          character(LEN=31)::mes,basis(4)
          integer::ic(4)
          real(8)::v0
          mes='repulse_HsHsHsHs'
          ic(:)=(/1,1,1,1/)
          basis(:)=(/'Hs','Hs','Hs','Hs'/)
          v0=0.774605944d0
          call testdrv_repulse(mes,ic,basis,v0)
        end subroutine

        subroutine test_repulse_cpxcpxcpxcpx()
          implicit none
          character(LEN=31)::mes,basis(4)
          integer::ic(4)
          real(8)::v0
          mes='repulse_CpxCpxCpxCpx'
          ic(:)=(/2,2,2,2/)
          basis(:)=(/'Cpx','Cpx','Cpx','Cpx'/)
          v0=0.672832726d0
          call testdrv_repulse(mes,ic,basis,v0)
        end subroutine

        subroutine test_repulse_pdxxpdxxpdxxpdxx()
          implicit none
          character(LEN=31)::mes,basis(4)
          integer::ic(4)
          real(8)::v0
          mes='repulse_PdxxPdxxPdxxPdxx'
          ic(:)=(/3,3,3,3/)
          basis(:)=(/'Pdxx','Pdxx','Pdxx','Pdxx'/)
          v0=0.633653982d0
          call testdrv_repulse(mes,ic,basis,v0)
        end subroutine

        subroutine test_repulse_cpxhshshs()
          implicit none
          character(LEN=31)::mes,basis(4)
          integer::ic(4)
          real(8)::v0
          mes='repulse_CpxHsHsHs'
          ic(:)=(/2,1,1,1/)
          basis(:)=(/'Cpx','Hs ','Hs ','Hs '/)
          v0=-0.361043291d0
          call testdrv_repulse(mes,ic,basis,v0)
        end subroutine

        subroutine test_repulse_pdxxhshshs()
          implicit none
          character(LEN=31)::mes,basis(4)
          integer::ic(4)
          real(8)::v0
          mes='repulse_PdxxHsHsHs'
          ic(:)=(/3,1,1,1/)
          basis(:)=(/'Pdxx','Hs  ','Hs  ','Hs  '/)
          v0=0.016862462d0
          call testdrv_repulse(mes,ic,basis,v0)
        end subroutine

        subroutine test_repulse_pdxxcpxhshs()
          implicit none
          character(LEN=31)::mes,basis(4)
          integer::ic(4)
          real(8)::v0
          mes='repulse_PdxxCpxHsHs'
          ic(:)=(/3,2,1,1/)
          basis(:)=(/'Pdxx','Cpx ','Hs  ','Hs  '/)
          v0=0.133233339d0
          call testdrv_repulse(mes,ic,basis,v0)
        end subroutine

        subroutine testdrv_velocity(mes,ic,basis)
          ! Test driver for velocity integrals (numerical testing).
          use eval
          implicit none
          integer,parameter::NC=3,IX=1
          real(8),parameter::TOL=1.d-8,DX=1d-3
          character(LEN=31),intent(in)::mes,basis(2)
          integer,intent(in)::ic(2)
          integer::la(3),lb(3)
          real(8)::val,r(3,3),c1(NC),c2(NC),e1(NC),e2(NC),vf,vb,rt(3,3)
          logical::ok
          ! Default settings.
          r(:,:)=0d0
          r(1,1)=-2.0373778380d0
          r(1,3)=2.7939768804d0
          ! Set Atomic orbital basis set.
          call setup_sto3g_(basis(1),la,e1,c1)
          call setup_sto3g_(basis(2),lb,e2,c2)
          ! Test run.
          call eval_velocity(
     &      NC,NC,c1,c2,e1,e2,r(:,ic(1)),r(:,ic(2)),la,lb,IX,val)
          ! Forward shift.
          rt(:,:)=r(:,:)
          rt(IX,ic(2))=rt(IX,ic(2))-DX
          call eval_overlap(
     &      NC,NC,c1,c2,e1,e2,r(:,ic(1)),rt(:,ic(2)),la,lb,vf)
          ! Backward shift.
          rt(:,:)=r(:,:)
          rt(IX,ic(2))=rt(IX,ic(2))+DX
          call eval_overlap(
     &      NC,NC,c1,c2,e1,e2,r(:,ic(1)),rt(:,ic(2)),la,lb,vb)
          ok=(val-(vf-vb)/(2d0*DX))**2<TOL
          if (.not.ok) then
            call asserttrue_(mes,ok)
          end if
        end subroutine

        subroutine test_velocity_hshs()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_HsHs'
          ic(:)=(/1,1/)
          basis(:)=(/'Hs','Hs'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_cpxcpx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_CpxCpx'
          ic(:)=(/2,2/)
          basis(:)=(/'Cpx','Cpx'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_pdxxpdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_PdxxPdxx'
          ic(:)=(/3,3/)
          basis(:)=(/'Pdxx','Pdxx'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_hscs()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_HsCs'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs','Cs'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_hscpx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_HsCpx'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs ','Cpx'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_hscpy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_HsCpy'
          ic(:)=(/1,2/)
          basis(:)=(/'Hs ','Cpy'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_hspdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_HsPdxx'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdxx'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_hspdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_HsPdyy'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdyy'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_hspdxy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_HsPdxy'
          ic(:)=(/1,3/)
          basis(:)=(/'Hs  ','Pdxy'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_cpxppx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_CpxPpx'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx','Ppx'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_cpxpdxx()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_CpxPdxx'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx ','Pdxx'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_cpxpdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_CpxPdyy'
          ic(:)=(/2,3/)
          basis(:)=(/'Cpx ','Pdyy'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine

        subroutine test_velocity_pdxxpdyy()
          implicit none
          character(LEN=31)::mes,basis(2)
          integer::ic(2)
          mes='velocity_PdxxPdyy'
          ic(:)=(/3,3/)
          basis(:)=(/'Pdxx','Pdyy'/)
          call testdrv_velocity(mes,ic,basis)
        end subroutine
        subroutine test_eval_overlap()
          implicit none
          call test_overlap_hshs()
          call test_overlap_cpxcpx()
          call test_overlap_pdxxpdxx()
          call test_overlap_hscs()
          call test_overlap_hscpx()
          call test_overlap_hscpy()
          call test_overlap_hspdxx()
          call test_overlap_hspdyy()
          call test_overlap_hspdxy()
          call test_overlap_cpxppx()
          call test_overlap_cpxpdxx()
          call test_overlap_cpxpdyy()
          call test_overlap_pdxxpdyy()
          call asserttrue_('overlap',.true.)
        end subroutine

        subroutine test_eval_kinetic()
          implicit none
          call test_kinetic_hshs()
          call test_kinetic_cpxcpx()
          call test_kinetic_pdxxpdxx()
          call test_kinetic_hscs()
          call test_kinetic_hscpx()
          call test_kinetic_hscpy()
          call test_kinetic_hspdxx()
          call test_kinetic_hspdyy()
          call test_kinetic_hspdxy()
          call test_kinetic_cpxppx()
          call test_kinetic_cpxpdxx()
          call test_kinetic_cpxpdyy()
          call test_kinetic_pdxxpdyy()
          call asserttrue_('kinetic',.true.)
        end subroutine

        subroutine test_eval_hcore()
          implicit none
          call test_hcore_hshs()
          call test_hcore_cpxcpx()
          call test_hcore_pdxxpdxx()
          call test_hcore_hscs()
          call test_hcore_hscpx()
          call test_hcore_hscpy()
          call test_hcore_hspdxx()
          call test_hcore_hspdyy()
          call test_hcore_hspdxy()
          call test_hcore_cpxppx()
          call test_hcore_cpxpdxx()
          call test_hcore_cpxpdyy()
          call test_hcore_pdxxpdyy()
          call asserttrue_('hcore',.true.)
        end subroutine

        subroutine test_eval_repulse()
          implicit none
          call test_repulse_hshshshs()
          call test_repulse_cpxcpxcpxcpx()
          call test_repulse_pdxxpdxxpdxxpdxx()
          call test_repulse_cpxhshshs()
          call test_repulse_pdxxhshshs()
          call test_repulse_pdxxcpxhshs()
          call asserttrue_('repulse',.true.)
        end subroutine

        subroutine test_eval_velocity()
          implicit none
          call test_velocity_hshs()
          call test_velocity_cpxcpx()
          call test_velocity_pdxxpdxx()
          call test_velocity_hscs()
          call test_velocity_hscpx()
          call test_velocity_hscpy()
          call test_velocity_hspdxx()
          call test_velocity_hspdyy()
          call test_velocity_hspdxy()
          call test_velocity_cpxppx()
          call test_velocity_cpxpdxx()
          call test_velocity_cpxpdyy()
          call test_velocity_pdxxpdyy()
          call asserttrue_('velocity',.true.)
        end subroutine

        subroutine test_main()
          implicit none
          write(*,'("--- EVAL ---")')
          call test_eval_overlap()
          call test_eval_kinetic()
          call test_eval_hcore()
          call test_eval_repulse()
          call test_eval_velocity()
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

      end module

      program main
        use test_eval,only: test_main
        implicit none
        call test_main()
      end program
