      module uhf
        ! Unrestricted Hatree-Fock.
        use mtrx
        implicit none
        private
        integer::maxit_=20
        real(8)::conv_=1d-4
        !---------------------------------------------------------------
        public::
     &      uhf_scf,
     &      uhf_eelc,
     &      uhf_dens,
     &      uhf_fock,
     &      uhf_seti,
     &      uhf_setf
        !---------------------------------------------------------------

      contains

        subroutine uhf_scf(nea,neb,sb,hb,tb,cma,cmb,ema,emb)
          ! UHF-SCF calculation to obtain molecular orbitals.
          implicit none
          integer,intent(in)::nea,neb
          real(8),intent(in)::sb(:),hb(:),tb(:)
          real(8),intent(inout)::cma(:,:),cmb(:,:)
          real(8),intent(out)::ema(:),emb(:)
          integer::i,nb,nb2
          real(8)::delta,ee
          real(8),allocatable::
     &      xb(:,:),pba(:),pbb(:),fba(:),fbap(:),fbb(:),fbbp(:),
     &      pba0(:),pbb0(:)
          nb=SIZE(ema)
          nb2=size(sb)
          allocate(
     &      xb(nb,nb),pba(nb2),pbb(nb2),fba(nb2),fbb(nb2),
     &      fbap(nb2),fbbp(nb2),pba0(nb2),pbb0(nb2)
     &      )
          call mtrx_xcanon(sb,xb)
          call uhf_dens(nea,cma,pba)
          call uhf_dens(neb,cmb,pbb)
          do i=1,maxit_
            pba0(:)=pba(:)
            pbb0(:)=pbb(:)
            call uhf_fock(pba,pbb,hb,tb,fba)
            call uhf_fock(pbb,pba,hb,tb,fbb)
            fbap(:)=fba(:)
            fbbp(:)=fbb(:)
            call mtrx_utau(xb,fbap)
            call mtrx_utau(xb,fbbp)
            call mtrx_eigsp(fbap,ema,cma)
            call mtrx_eigsp(fbbp,emb,cmb)
            call mtrx_uc(xb,cma)
            call mtrx_uc(xb,cmb)
            call uhf_dens(nea,cma,pba)
            call uhf_dens(neb,cmb,pbb)
            delta=SQRT((SUM(pba-pba0)**2+SUM(pbb-pbb0)**2)/4)
            if (delta<conv_) then
              exit
            end if
            if (i==maxit_) then
              call uhf_eelc(nea,neb,hb,tb,cma,cmb,ee)
              write(*,'("SCF did not converge")')
              write(*,'("i,ee,delta=",i4,2f13.6)') i,ee,delta
              stop 1
            end if
          end do
        end subroutine

        subroutine uhf_eelc(nea,neb,hb,tb,cma,cmb,val)
          ! Electronic energy (naive implementation).
          implicit none
          integer::nea,neb
          real(8),intent(in)::hb(:),tb(:),cma(:,:),cmb(:,:)
          real(8),intent(out)::val
          real(8)::vt
          val=0d0
          ! Alpha-MO contributions.
          call eelcaa_(nea,hb,tb,cma,vt)
          val=val+vt
          ! Beta-MO contributions.
          call eelcaa_(neb,hb,tb,cmb,vt)
          val=val+vt
          ! Alpha-Beta MOs contribution.
          call eelcab_(nea,neb,hb,tb,cma,cmb,vt)
          val=val+vt
        end subroutine

        subroutine eelcaa_(ne,hb,tb,cm,val)
          ! Alpha- or beta-MO contribution to the electronic energy.
          implicit none
          integer,intent(in)::ne
          real(8),intent(in)::hb(:),tb(:),cm(:,:)
          real(8),intent(out)::val
          integer::i,j
          real(8),allocatable::hm(:),tm(:)
          ! Get memory.
          allocate(hm(SIZE(hb)),tm(SIZE(tb)))
          ! Transform integrals into the MO basis.
          hm(:)=hb(:)
          tm(:)=tb(:)
          call mtrx_utau(cm,hm)
          call mtrx_utau2(cm,tm)
          val=0d0
          ! 1e integral contributions.
          do i=1,ne
            val=val+hm(ip_(i,i))
          end do
          ! 2e integral contributions.
          do i=1,ne
            do j=1,i
              val=val+tm(ip2_(i,i,j,j))
              val=val-tm(ip2_(i,j,j,i))
            end do
          end do
        end subroutine

        subroutine eelcab_(nea,neb,hb,tb,cma,cmb,val)
          ! Alpha- and beta-MO coupling term for the electronic energy.
          implicit none
          integer,intent(in)::nea,neb
          real(8),intent(in)::hb(:),tb(:),cma(:,:),cmb(:,:)
          real(8),intent(out)::val
          integer::i,j,ij,kl,ijkl,n2
          real(8),allocatable::tmu(:,:)
          ! Get memory.
          n2=SIZE(hb)
          allocate(tmu(n2,n2))
          ! Transform integrals into the MO basis.
          do ij=1,n2
            do kl=1,n2
              ijkl=ip_(ij,kl)
              tmu(ij,kl)=tb(ijkl)
            end do
          end do
          call mtrx_utau2u(cma,cmb,tmu)
          ! 2e integral contributions.
          val=0d0
          do i=1,nea
            do j=1,neb
              ij=ip_(i,i)
              kl=ip_(j,j)
              val=val+tmu(ij,kl)
            end do
          end do
        end subroutine

        subroutine uhf_dens(ne,cm,pb)
          ! AO density matrix from MO coefficients.
          implicit none
          integer,intent(in)::ne
          real(8),intent(in)::cm(:,:)
          real(8),intent(out)::pb(:)
          integer::ib,jb,ij,im,nb
          nb=SIZE(cm(:,1))
          pb(:)=0d0
          do ib=1,nb
            do jb=1,ib
              ij=ip_(ib,jb)
              do im=1,ne
                pb(ij)=pb(ij)+cm(ib,im)*cm(jb,im)
              end do
            end do
          end do
        end subroutine

        subroutine uhf_fock(pba,pbb,hb,tb,val)
          ! Alpha (or beta) Fock matrix F in AO basis.
          implicit none
          real(8),intent(in)::pba(:),pbb(:),hb(:),tb(:)
          real(8),intent(out)::val(:)
          integer :: ijkl,ilkj,ij,kl,i,j,k,l,nb,nb2
          val(:)=hb(:)
          nb2=SIZE(pba)
          call iup_(nb2,i,nb)
          do ij=1,nb2
            call iup_(ij,i,j)
            do k=1,nb
              do l=1,nb
                kl=ip_(k,l)
                ijkl=ip_(ij,kl)
                val(ij)=val(ij)+(pba(kl)+pbb(kl))*tb(ijkl)
                ilkj=ip2_(i,l,k,j)
                val(ij)=val(ij)-pba(kl)*tb(ilkj)
              end do
            end do
          end do
        end subroutine

        subroutine uhf_seti(key,n)
          ! Set an integer parameter.
          implicit none
          character(*),intent(in)::key
          integer,intent(in)::n
          if (key=='maxit') then
            maxit_=n
          else
            write(*,'("error: no such key: ",a)') key
            stop 1
          end if
        end subroutine

        subroutine uhf_setf(key,v)
          ! Set a float parameter.
          implicit none
          character(*),intent(in)::key
          real(8),intent(in)::v
          if (key=='conv') then
            conv_=v
          else
            write(*,'("error: no such key: ",a)') key
            stop 1
          end if
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

        integer function ip2_(i,j,k,l) result(ijkl)
          ! Canonical index packing for a symmetric matrix.
          implicit none
          integer,intent(in)::i,j,k,l
          integer::ij,kl
          ij=ip_(i,j)
          kl=ip_(k,l)
          ijkl=ip_(ij,kl)
        end function

      end module
