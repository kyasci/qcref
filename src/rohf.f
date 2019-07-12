      module rohf
        ! Restricted open-shell Hartree-Fock.
        use mtrx
        use uhf
        implicit none
        private
        integer::maxit_=20,nfock_=0
        real(8)::conv_=1d-4
        !---------------------------------------------------------------
        public::
     &      rohf_scf,
     &      rohf_eelc,
     &      rohf_seti,
     &      rohf_setf
        !---------------------------------------------------------------

      contains

        subroutine rohf_scf(nea,neb,sb,hb,tb,cm,em)
          implicit none
          integer,intent(in)::nea,neb
          real(8),intent(in)::sb(:),hb(:),tb(:)
          real(8),intent(inout)::cm(:,:)
          real(8),intent(out)::em(:)
          integer::i,nb,nb2
          real(8)::delta,ee
          real(8),allocatable::
     &      xb(:,:),pba(:),pbb(:),fb(:),pba0(:),fba(:),fbb(:)
          nb=SIZE(em)
          nb2=nb*(nb+1)/2
          allocate(
     &      xb(nb,nb),pba(nb2),pbb(nb2),fb(nb2),
     &      pba0(nb2),fba(nb2),fbb(nb2)
     &      )
          call mtrx_xcanon(sb,xb)
          call uhf_dens(nea,cm,pba)
          call uhf_dens(neb,cm,pbb)
          do i=1,maxit_
            pba0(:)=pba(:)
            call uhf_fock(pba,pbb,hb,tb,fba)
            call uhf_fock(pbb,pba,hb,tb,fbb)
            call fock_(nea,neb,cm,sb,fba,fbb,fb)
            call mtrx_utau(xb,fb)
            call mtrx_eigsp(fb,em,cm)
            call mtrx_uc(xb,cm)
            call uhf_dens(nea,cm,pba)
            call uhf_dens(neb,cm,pbb)
            delta=SQRT((SUM(pba-pba0)**2)/4)
            if (delta<conv_) then
              exit
            end if
            if (i==maxit_) then
              call rohf_eelc(nea,neb,hb,tb,cm,ee)
              write(*,'("SCF did not converge")')
              write(*,'("i,ee,delta=",i4,2f13.6)') i,ee,delta
              stop 1
            end if
          end do
        end subroutine

        subroutine rohf_eelc(nea,neb,hb,tb,cm,val)
          ! Electronic energy.
          implicit none
          integer,intent(in)::nea,neb
          real(8),intent(in)::hb(:),tb(:),cm(:,:)
          real(8),intent(out)::val
          call uhf_eelc(nea,neb,hb,tb,cm,cm,val)
        end subroutine

        subroutine rohf_seti(key,n)
          ! Set an integer parameter.
          implicit none
          character(*),intent(in)::key
          integer,intent(in)::n
          if (key=='maxit') then
            maxit_=n
          else if (key=='nfock') then
            nfock_=n
          else
            write(*,'("error: no such key: ",a)') key
            stop 1
          end if
        end subroutine

        subroutine rohf_setf(key,v)
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

        subroutine fockm_(nea,neb,cm,fba,fbb,val)
          ! Combined Fock matrix in the MO basis.
          implicit none
          integer,intent(in)::nea,neb
          real(8),intent(in)::cm(:,:),fba(:),fbb(:)
          real(8),intent(out)::val(:)
          integer::i,j,ij,nb
          real(8)::acc,bcc,aoo,boo,avv,bvv
          real(8),allocatable::fa(:),fb(:)
          call fockparam_(acc,bcc,aoo,boo,avv,bvv)
          call iup_(SIZE(fba),i,nb)
          ! Fock matrix in the MO basis.
          call mtrx_utaup(cm,fba,fa)
          call mtrx_utaup(cm,fbb,fb)
          ! Form combined Fock matrix in the MO basis.
          val(:)=0d0
          do i=1,nb
            do j=1,i
              ij=ip_(i,j)
              if (i<=neb.and.j<=neb) then
                val(ij)=acc*fa(ij)+bcc*fb(ij)
              else if (neb<i.and.i<=nea.and.neb<j.and.j<=nea) then
                val(ij)=aoo*fa(ij)+boo*fb(ij)
              else if (i>nea.and.j>nea) then
                val(ij)=avv*fa(ij)+bvv*fb(ij)
              else if (neb<i.and.i<=nea.and.j<=neb) then
                val(ij)=fb(ij)
              else if (nea<i.and.j<=neb) then
                val(ij)=(fa(ij)+fb(ij))/2
              else if (nea<i.and.neb<j.and.j<=nea) then
                val(ij)=fa(ij)
              else
                write(*,'("error: out of range: fock()")')
              end if
            end do
          end do
        end subroutine

        subroutine fock_(nea,neb,cm,sb,fba,fbb,val)
          ! Combined Fock matrix in the AO basis.
          implicit none
          integer,intent(in)::nea,neb
          real(8),intent(in)::cm(:,:),sb(:),fba(:),fbb(:)
          real(8),intent(out)::val(:)
          integer::i,j,ij,nb
          real(8),allocatable::fa(:),fb(:),ct(:,:)
          ! Get memory.
          call iup_(SIZE(fba),i,nb)
          allocate(ct(nb,nb))
          ! The inverse of MO coefficients.
          do ij=1,SIZE(sb)
            call iup_(ij,i,j)
            ct(i,j)=sb(ij)
            ct(j,i)=sb(ij)
          end do
          call mtrx_uc(TRANSPOSE(cm),ct)
          ! Fock matrix in the MO basis.
          call mtrx_utaup(cm,fba,fa)
          call mtrx_utaup(cm,fbb,fb)
          ! Form combined Fock matrix in the MO basis.
          call fockm_(nea,neb,cm,fba,fbb,val)
          ! Transform the Fock matrix to the AO basis.
          call mtrx_utau(ct,val)
        end subroutine

        subroutine fockparam_(acc,bcc,aoo,boo,avv,bvv)
          ! Canonicalization coefficients.
          implicit none
          real(8),intent(out)::acc,bcc,aoo,boo,avv,bvv
          ! Roothaan single matrix
          if (nfock_==0) then
            acc=-0.5d0
            bcc=1.5d0
            aoo=0.5d0
            boo=0.5d0
            avv=1.5d0
            bvv=-0.5d0
          else
            write(*,'("error: no such param: fockparam_()")')
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

      end module
