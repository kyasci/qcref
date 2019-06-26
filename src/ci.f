      module ci
        ! Configuration Intaraction (CI) calculation.
        !
        !   Requirements:
        !
        !     * 1e and 2e integrals in AO basis
        !         (see eval.f and basis.f)
        !
        !     * MO or NO coefficients
        !         (see rhf.f and rohf.f)
        !
        !     * 1e and 2e coupling constants
        !         (see uga.f)
        !
        implicit none
        private
        !---------------------------------------------------------------
        public ::
     &    ci_hcsf,
     &    ci_int1,
     &    ci_int2,
     &    ci_dm1,
     &    ci_dm2,
     &    ci_lgr
        !---------------------------------------------------------------

      contains

        subroutine ci_hcsf(lc1,vc1,lc2,vc2,hm,tm,val)
          ! Evaluate CI Hamiltonian matrix.
          implicit none
          integer,intent(in)::lc1(:,:),lc2(:,:)
          real(8),intent(in)::vc1(:),vc2(:),hm(:),tm(:)
          real(8),intent(out),allocatable::val(:)
          real(8),allocatable::wk(:)
          integer::nc,nc2
          ! Get memory.
          nc=MAXVAL(lc1(1,:))
          nc2=ip_(nc,nc)
          allocate(val(nc2),wk(nc2))
          val(:)=0d0
          ! 1e integral contribution.
          call ci_int1(lc1,vc1,hm,wk)
          val(:)=val(:)+wk(:)
          ! 2e integral contribution.
          call ci_int2(lc1,vc1,lc2,vc2,tm,wk)
          val(:)=val(:)+wk(:)
        end subroutine

        subroutine ci_int2(lc1,vc1,lc2,vc2,tm,val)
          ! Calculate CI 2e integrals in the CSF basis.
          implicit none
          integer,intent(in)::lc1(:,:),lc2(:,:)
          real(8),intent(in)::vc1(:),vc2(:),tm(:)
          real(8),intent(out)::val(:)
          real(8),allocatable::vt(:)
          integer::nfc,nc
          allocate(vt(SIZE(val)))
          val(:)=0d0
          ! Core-Core contribution.
          nfc=MINVAL(lc1(3,:))-1
          nc=MAXVAL(lc1(1,:))
          call int2cc_(nfc,nc,tm,vt)
          val(:)=val(:)+vt(:)
          ! Core-Active contribution.
          call int2ca_(nfc,lc1,vc1,tm,vt)
          val(:)=val(:)+vt(:)
          ! Active-Active contribution.
          call int2aa_(lc2,vc2,tm,vt)
          val(:)=val(:)+vt(:)
        end subroutine

        subroutine ci_int1(lc1,vc1,hm,val)
          ! Calculate CI 1e integrals in the CSF basis.
          implicit none
          integer,intent(in)::lc1(:,:)
          real(8),intent(in)::vc1(:),hm(:)
          real(8),intent(out)::val(:)
          integer::nfc,nc
          real(8),allocatable::vt(:)
          allocate(vt(SIZE(val)))
          val(:)=0d0
          ! Core contribution.
          nfc=MINVAL(lc1(3,:))-1
          nc=MAXVAL(lc1(1,:))
          call int1c_(nfc,nc,hm,vt)
          val(:)=val(:)+vt(:)
          ! Active contribution.
          call int1a_(lc1,vc1,hm,vt)
          val(:)=val(:)+vt(:)
        end subroutine

        subroutine ci_dm1(lc1,vc1,c,ca,val)
          ! Reduced 1e density matrix in the MO basis.
          implicit none
          integer,intent(in)::lc1(:,:)
          real(8),intent(in)::vc1(:),c(:),ca(:,:)
          real(8),intent(out)::val(:)
          integer::nc,nc2,ic,jc,ijc,im,jm,ijm,inz,nfc
          real(8)::a1
          real(8),allocatable::pc(:)
          nc=SIZE(ca(:,1))
          nfc=MINVAL(lc1(3,:))-1
          ! Calculate the CI density matrix.
          nc2=ip_(nc,nc)
          allocate(pc(nc2))
          call evalcsfdm_(c,ca,pc)
          ! Calculate the reduced 1e density matrix in the MO basis.
          val(:)=0d0
          ! -- core orbitals contribution --
          do im=1,nfc
            ijm=ip_(im,im)
            val(ijm)=val(ijm)+2d0
          end do
          ! -- active orbitals contribution --
          do inz=1,SIZE(vc1)
            call loadc1nz_(inz,lc1,vc1,ic,jc,im,jm,a1)
            ijc=ip_(ic,jc)
            ijm=ip_(im,jm)
            val(ijm)=val(ijm)+pc(ijc)*a1
          end do
        end subroutine

        subroutine ci_lgr(nfzc,nacc,dm1,dm2,hm,tm,val)
          ! CI Lagrangian matrix.
          implicit none
          integer,intent(in)::nfzc,nacc
          real(8),intent(in)::dm1(:),dm2(:),hm(:),tm(:)
          real(8),intent(out)::val(:,:)
          real(8),allocatable::fa(:),fc(:),tmp(:,:)
          integer::nb,nb2
          integer::i,mi,mq,mu,mv,miq,mvu,muq,my,mx,myuvx,mquvx
          ! Get memory.
          nb2=SIZE(dm1)
          call iup_(nb2,nb,i)
          allocate(fa(nb2),fc(nb2),tmp(nb,nb))
          ! Core (inactive) and valence (active) Fock matrices.
          call evalfi_(nfzc,hm,tm,fc)
          call evalfa_(nfzc,nacc,dm1,tm,fa)
          ! Pure Fock matrix contribution.
          val(:,:)=0d0
          do mq=1,nb
            do mi=1,nfzc
              miq=ip_(mi,mq)
              val(mi,mq)=2d0*(fc(miq)+fa(miq))
            end do
          end do
          ! One-particle density contribution.
          do mq=1,nb
            do mv=nfzc+1,nfzc+nacc
              do mu=nfzc+1,nfzc+nacc
                mvu=ip_(mv,mu)
                muq=ip_(mq,mu)
                val(mv,mq)=val(mv,mq)+dm1(mvu)*fc(muq)
              end do
            end do
          end do
          ! Two-particl density contribution.
          tmp(:,:)=0d0
          do mq=1,nb
            do my=nfzc+1,nfzc+nacc
              do mu=nfzc+1,nfzc+nacc
                do mv=nfzc+1,nfzc+nacc
                  do mx=nfzc+1,nfzc+nacc
                    myuvx=ip2_(my,mu,mv,mx)
                    mquvx=ip2_(mq,mu,mv,mx)
                    tmp(my,mq)=tmp(my,mq)+2d0*dm2(myuvx)*tm(mquvx)
                  end do
                end do
              end do
            end do
          end do
          do mq=1,nb
            val(:,mq)=val(:,mq)+tmp(:,mq)
          end do
        end subroutine

        subroutine evalfa_(nfzc,nacc,dm1,tm,fa)
          ! Generalized Fock matrix (active orbital contribution).
          implicit none
          integer,intent(in)::nfzc,nacc
          real(8),intent(in)::dm1(:),tm(:)
          real(8),intent(out)::fa(:)
          integer::i,j,k,l,ij,kl,ijkl,nb
          fa(:)=0d0
          call iup_(SIZE(fa),nb,i)
          do ij=1,SIZE(fa)
            call iup_(ij,i,j)
            do k=nfzc+1,nfzc+nacc
              do l=nfzc+1,nfzc+nacc
                kl=ip_(k,l)
                ijkl=ip2_(i,j,k,l)
                fa(ij)=fa(ij)+1.0d0*dm1(kl)*tm(ijkl)
                ijkl=ip2_(i,l,j,k)
                fa(ij)=fa(ij)-0.5d0*dm1(kl)*tm(ijkl)
              end do
            end do
          end do
        end subroutine

        subroutine evalfi_(nfzc,hm,tm,fi)
          ! Generalized Fock matrix (inactive orbital contribution).
          implicit none
          integer,intent(in)::nfzc
          real(8),intent(in)::hm(:),tm(:)
          real(8),intent(out)::fi(:)
          integer::ij,i,j,k,ijkl
          ! 1e contribution.
          fi(:)=hm(:)
          if (nfzc<=0) then
            return
          end if
          ! 2e contribution.
          do ij=1,SIZE(fi)
            call iup_(ij,i,j)
            do k=1,nfzc
              ijkl=ip2_(i,j,k,k)
              fi(ij)=fi(ij)+2d0*tm(ijkl)
              ijkl=ip2_(i,k,k,j)
              fi(ij)=fi(ij)-1d0*tm(ijkl)
            end do
          end do
        end subroutine

        subroutine loadc1nz_(inz,lnz,vnz,ic,jc,im,jm,a1)
          ! Return non-zero values of 1e coupling const.
          implicit none
          integer,intent(in)::inz,lnz(:,:)
          real(8),intent(in)::vnz(:)
          integer,intent(out)::ic,jc,im,jm
          real(8),intent(out)::a1
          ic=lnz(1,inz)
          jc=lnz(2,inz)
          im=lnz(3,inz)
          jm=lnz(4,inz)
          a1=vnz(inz)
        end subroutine

        subroutine evalcsfdm_(c,ca,val)
          ! Evaluate the CI density matrix.
          implicit none
          real(8),intent(in)::c(:),ca(:,:)
          real(8),intent(out)::val(:)
          integer::na,ncr2,ijcr,icr,jcr,ia,ja
          real(8)::cc
          na=SIZE(c)
          ncr2=SIZE(val)
          val(:)=0d0
          do ijcr=1,ncr2
            call iup_(ijcr,icr,jcr)
            do ia=1,na
              do ja=1,na
                cc=c(ia)*c(ja)*ca(icr,ia)*ca(jcr,ja)
                val(ijcr)=val(ijcr)+cc
              end do
            end do
          end do
        end subroutine

        subroutine loadc2nz_(inz,lc2,vc2,icr,jcr,im,jm,km,lm,val)
          ! Load the non-zero value of 2e coupling const.
          implicit none
          integer,intent(in)::inz,lc2(:,:)
          real(8),intent(in)::vc2(:)
          integer,intent(out)::icr,jcr,im,jm,km,lm
          real(8),intent(out)::val
          icr=lc2(1,inz)
          jcr=lc2(2,inz)
          im=lc2(3,inz)
          jm=lc2(4,inz)
          km=lc2(5,inz)
          lm=lc2(6,inz)
          val=vc2(inz)
        end subroutine

        subroutine ci_dm2(lc1,vc1,lc2,vc2,c,ca,val)
          ! Reduced 2e density matrix in the MO basis.
          implicit none
          integer,intent(in)::lc1(:,:),lc2(:,:)
          real(8),intent(in)::vc1(:),vc2(:),c(:),ca(:,:)
          real(8),intent(out)::val(:)
          integer::
     &      nc,nc2,nb3,inz,ic,jc,ijc,im,jm,km,lm,ijklm,ijm,klm,nfc
          real(8)::a1,coef
          real(8),allocatable::pc(:)
          ! Calculate the CI density matrix.
          nc=SIZE(ca(:,1))
          nfc=MINVAL(lc1(3,:))-1
          nc2=ip_(nc,nc)
          allocate(pc(nc2))
          call evalcsfdm_(c,ca,pc)
          ! Calculate the reduced 2e density matrix in the MO basis.
          nb3=SIZE(val)
          val(:)=0d0
          ! -- core-core contribution --
          do im=1,nfc
            do jm=1,nfc
              ijklm=ip2_(im,im,jm,jm)
              val(ijklm)=val(ijklm)+2d0
              ijklm=ip2_(im,jm,jm,im)
              val(ijklm)=val(ijklm)-1d0
            end do
          end do
          ! -- core-active contribution --
          do inz=1,SIZE(vc1)
            call loadc1nz_(inz,lc1,vc1,ic,jc,im,jm,a1)
            ijc=ip_(ic,jc)
            coef=pc(ijc)*a1
            if (ic/=jc) then
              coef=coef*2d0
            end if
            do km=1,nfc
              ijklm=ip2_(im,jm,km,km)
              val(ijklm)=val(ijklm)+coef*2d0
              ijklm=ip2_(im,km,km,jm)
              val(ijklm)=val(ijklm)-coef
            end do
          end do
          ! -- active-active contribution --
          do inz=1,SIZE(vc2)
            call loadc2nz_(inz,lc2,vc2,ic,jc,im,jm,km,lm,a1)
            ijc=ip_(ic,jc)
            ijklm=ip2_(im,jm,km,lm)
            coef=pc(ijc)*a1
            if (ic/=jc) then
              coef=coef*2d0
            end if
            val(ijklm)=val(ijklm)+coef
          end do
          ! Restore the weight factors.
          do ijklm=1,nb3
            call iup_(ijklm,ijm,klm)
            call iup_(ijm,im,jm)
            call iup_(klm,km,lm)
            coef=8d0
            if (im==jm) coef=coef/2d0
            if (km==lm) coef=coef/2d0
            if (ijm==klm) coef=coef/2d0
            val(ijklm)=val(ijklm)/coef
          end do
        end subroutine

        subroutine int1c_(nfc,nc,hm,val)
          ! Core contribution to the CI 1e integrals.
          implicit none
          integer,intent(in)::nfc,nc
          real(8),intent(in)::hm(:)
          real(8),intent(out)::val(:)
          integer::im,iim,icr,ijcr
          real(8)::vnn
          val(:)=0d0
          if (nfc<=0) then
            return
!           ******
          end if
          vnn=0d0
          do im=1,nfc
            iim=ip_(im,im)
            vnn=vnn+2d0*hm(iim)
          end do
          do icr=1,nc
            ijcr=ip_(icr,icr)
            val(ijcr)=vnn
          end do
        end subroutine

        subroutine int1a_(lc1,vc1,hm,val)
          ! Active-space contribution to the CI 1e integrals.
          implicit none
          integer,intent(in)::lc1(:,:)
          real(8),intent(in)::vc1(:),hm(:)
          real(8),intent(inout)::val(:)
          integer::ic,jc,im,jm,ijm,inz,ijc
          real(8)::a1
          val(:)=0d0
          do inz=1,SIZE(vc1)
            call loadc1nz_(inz,lc1,vc1,ic,jc,im,jm,a1)
            ijm=ip_(im,jm)
            ijc=ip_(ic,jc)
            val(ijc)=val(ijc)+a1*hm(ijm)
          end do
        end subroutine

        subroutine int2cc_(nfc,nc,tm,val)
          ! Core contribution to the CI 2e integrals.
          implicit none
          integer,intent(in)::nfc,nc
          real(8),intent(in)::tm(:)
          real(8),intent(out)::val(:)
          integer::im,jm,ijklm,icr,ijcr
          real(8)::vnn
          val(:)=0d0
          if (nfc<=0) then
            return
!           ******
          end if
          vnn=0d0
          do im=1,nfc
            do jm=1,nfc
              ijklm=ip2_(im,im,jm,jm)
              vnn=vnn+2d0*tm(ijklm)
              ijklm=ip2_(im,jm,jm,im)
              vnn=vnn-tm(ijklm)
            end do
          end do
          do icr=1,nc
            ijcr=ip_(icr,icr)
            val(ijcr)=vnn
          end do
        end subroutine

        subroutine int2ca_(nfc,lc1,vc1,tm,val)
          ! Core-active contribution to the 2e integrals.
          implicit none
          integer,intent(in)::nfc,lc1(:,:)
          real(8),intent(in)::vc1(:),tm(:)
          real(8),intent(out)::val(:)
          integer::ic,jc,ijc,im,jm,km,ijklm,inz
          real(8)::a1
          val(:)=0d0
          do inz=1,SIZE(vc1)
            call loadc1nz_(inz,lc1,vc1,ic,jc,im,jm,a1)
            ijc=ip_(ic,jc)
            do km=1,nfc
              ijklm=ip2_(im,jm,km,km)
              val(ijc)=val(ijc)+2d0*a1*tm(ijklm)
              ijklm=ip2_(im,km,km,jm)
              val(ijc)=val(ijc)-a1*tm(ijklm)
            end do
          end do
        end subroutine

        subroutine int2aa_(lc2,vc2,tm,val)
          ! Active-space contribution to the CI 2e integrals.
          implicit none
          integer,intent(in)::lc2(:,:)
          real(8),intent(in)::vc2(:),tm(:)
          real(8),intent(inout)::val(:)
          integer::ic,jc,im,jm,km,lm,ijklm,inz,ijc
          real(8)::a1
          val(:)=0d0
          do inz=1,SIZE(vc2)
            call loadc2nz_(inz,lc2,vc2,ic,jc,im,jm,km,lm,a1)
            ijklm=ip2_(im,jm,km,lm)
            ijc=ip_(ic,jc)
            val(ijc)=val(ijc)+a1*tm(ijklm)
          end do
        end subroutine

        integer function ip_(i,j) result(ij)
          ! Canonical index packing for a symmetric matrix.
          implicit none
          integer,intent(in)::i,j
          integer::ii,jj
          ii=max(i,j)
          jj=min(i,j)
          ij=(ii-1)*ii/2+jj
        end function

        integer function ip2_(i,j,k,l) result(ijkl)
          ! Canonical index packing for a symmetric matrix (2e).
          implicit none
          integer,intent(in)::i,j,k,l
          integer::ij,kl
          ij=ip_(i,j)
          kl=ip_(k,l)
          ijkl=ip_(ij,kl)
        end function

        subroutine iup_(ij,i,j)
          ! Unpack the canonical index for a symmetric matrix (1e).
          implicit none
          integer,intent(in)::ij
          integer,intent(out)::i,j
          real(8)::vt
          vt=0.5d0*(-1d0+SQRT(1d0+8*ij))
          i=CEILING(vt)
          j=ij-(i-1)*i/2
        end subroutine

      end module
