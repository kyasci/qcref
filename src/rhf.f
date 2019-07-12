      module rhf
        ! Restricted closed-shell Hartree-Fock.
        use mtrx
        implicit none
        private
        integer::maxit_=20
        real(8)::conv_=1d-4
        !---------------------------------------------------------------
        public::
     &    rhf_scf,
     &    rhf_eelc,
     &    rhf_dens,
     &    rhf_seti,
     &    rhf_setf,
     &    rhf_fock
        !---------------------------------------------------------------

      contains

        subroutine rhf_scf(ne,sb,hb,tb,cm,em)
          ! RHF-SCF calculation to obtain molecular orbitals.
          implicit none
          integer,intent(in)::ne
          real(8),intent(in)::sb(:),hb(:),tb(:)
          real(8),intent(inout)::cm(:,:)
          real(8),intent(out)::em(:)
          integer::i,nb,nb2
          real(8)::ee,delta
          real(8),allocatable::xb(:,:),pb(:),fb(:),pb0(:)
          nb=SIZE(em)
          nb2=SIZE(sb)
          allocate(xb(nb,nb),pb(nb2),fb(nb2),pb0(nb2))
          call mtrx_xcanon(sb,xb)
          call rhf_dens(ne,cm,pb)
          do i=1,maxit_
            pb0(:)=pb(:)
            call rhf_fock(pb,hb,tb,fb)
            call mtrx_utau(xb,fb)
            call mtrx_eigsp(fb,em,cm)
            call mtrx_uc(xb,cm)
            call rhf_dens(ne,cm,pb)
            delta=SQRT(SUM((pb-pb0)**2)/4)
            if (delta<conv_) then
              exit
            end if
            if (i==maxit_) then
              call rhf_eelc(ne,hb,tb,cm,ee)
              write(*,'("SCF did not converge")')
              write(*,'("i,ee,delta=",i4,2f13.6)') i,ee,delta
              stop 1
            end if
          end do
        end subroutine

        subroutine rhf_eelc(ne,hb,tb,cm,val)
          ! Electronic energy (naive implementation).
          implicit none
          integer,intent(in)::ne
          real(8),intent(in)::hb(:),tb(:),cm(:,:)
          real(8),intent(out)::val
          integer::i,j,ndoc
          real(8),allocatable::hm(:),tm(:)
          ! Transform integrals into the MO basis.
          ndoc=ne/2
          call mtrx_utaup(cm(:,:ndoc),hb,hm)
          call mtrx_utau2p(cm(:,:ndoc),tb,tm)
          val=0d0
          ! 1e integral contributions.
          do i=1,ndoc
            val=val+2d0*hm(ip_(i,i))
          end do
          ! 2e integral contributions.
          do i=1,ndoc
            do j=1,ndoc
              val=val+2d0*tm(ip2_(i,i,j,j))
              val=val-tm(ip2_(i,j,j,i))
            end do
          end do
        end subroutine

        subroutine rhf_dens(ne,cm,pb)
          ! AO density matrix for RHF.
          implicit none
          integer,intent(in)::ne
          real(8),intent(in)::cm(:,:)
          real(8),intent(out)::pb(:)
          integer::ib,jb,ij,im,nb,ndoc
          nb=SIZE(cm(:,1))
          ndoc=ne/2
          pb(:)=0d0
          do ib=1,nb
            do jb=1,ib
              ij=ip_(ib,jb)
              do im=1,ndoc
                pb(ij)=pb(ij)+2d0*cm(ib,im)*cm(jb,im)
              end do
            end do
          end do
        end subroutine

        subroutine rhf_seti(key,n)
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

        subroutine rhf_setf(key,v)
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

        subroutine rhf_fock(pb,hb,tb,val)
          ! Fock matrix F in AO basis.
          implicit none
          real(8),intent(in)::pb(:),hb(:),tb(:)
          real(8),intent(out)::val(:)
          integer :: ijkl,ilkj,ij,kl,i,j,k,l,nb,nb2
          val(:)=hb(:)
          nb2=SIZE(pb)
          call iup_(nb2,nb,i)
          do ij=1,nb2
            call iup_(ij,i,j)
            do k=1,nb
              do l=1,nb
                kl=ip_(k,l)
                ijkl=ip_(ij,kl)
                ilkj=ip2_(i,l,k,j)
                val(ij)=val(ij)+pb(kl)*(tb(ijkl)-0.5d0*tb(ilkj))
              end do
            end do
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
