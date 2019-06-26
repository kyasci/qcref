      module mtrx
        implicit none
        private
        real(8),parameter::TOL_=1d-10
        !---------------------------------------------------------------
        public::
     &      mtrx_xcanon,
     &      mtrx_uc,
     &      mtrx_utau,
     &      mtrx_utau2,
     &      mtrx_utau2u,
     &      mtrx_eigsp,
     &      mtrx_expv,
     &      mtrx_expv2
        !---------------------------------------------------------------

      contains

        subroutine mtrx_xcanon(sb,val)
          ! Canonical transformation matrix X in AO basis.
          implicit none
          real(8),intent(in)::sb(:)
          real(8),intent(out)::val(:,:)
          integer::i,nb
          real(8),allocatable::e(:),sbt(:)
          nb=SIZE(val(:,1))
          allocate(e(nb))
          allocate(sbt(SIZE(sb)))
          sbt(:)=sb(:)
          e(:)=0d0
          call mtrx_eigsp(sbt,e,val)
          do i=1,nb
            val(:,i)=val(:,i)*e(i)**(-0.5d0)
          end do
        end subroutine

        subroutine mtrx_utau(u,a)
          ! Unitary transformation of a 1e matrix (2-index).
          implicit none
          real(8),intent(in)::u(:,:)
          real(8),intent(out)::a(:)
          integer::i,j,k,ij,n
          real(8),allocatable::vt(:,:)
          n=SIZE(u(:,1))
          allocate(vt(n,n))
          vt(:,:)=0d0
          do i=1,n
            do j=1,n
              do k=1,n
                vt(i,j)=vt(i,j)+a(ip_(i,k))*u(k,j)
              end do
            end do
          end do
          a(:)=0d0
          do i=1,n
            do j=1,i
              ij=ip_(i,j)
              do k=1,n
                a(ij)=a(ij)+u(k,i)*vt(k,j)
              end do
            end do
          end do
        end subroutine

        subroutine mtrx_utau2(cm,ab)
          ! Unitary transformation of a 2e matrix (4-index).
          implicit none
          real(8),intent(in)::cm(:,:)
          real(8),intent(inout)::ab(:)
          integer::n,n3,i,j,k,l,ijkl,im,jm,km,lm,ijklm
          real(8),allocatable::wk1(:,:,:,:),wk2(:,:,:,:)
          ! Get memory.
          n=SIZE(cm(:,1))
          n3=SIZE(ab)
          allocate(wk1(n,n,n,n))
          allocate(wk2(n,n,n,n))
          wk1(:,:,:,:)=0d0
          ! Transf.: a(i,j,k,l) -> a(i,j,k,lm)
!$omp     parallel private(i,j,k,l,lm,ijkl)
!$omp     do reduction(+:wk1)
          do i=1,n
            do j=1,n
              do k=1,n
                do lm=1,n
                  do l=1,n
                    ijkl=ip2_(i,j,k,l)
                    wk1(i,j,k,lm)=wk1(i,j,k,lm)+cm(l,lm)*ab(ijkl)
                  end do
                end do
              end do
            end do
          end do
!$omp     end do
!$omp     end parallel
          wk2(:,:,:,:)=0d0
          ! Transf.: a(i,j,k,lm) -> a(i,j,km,lm)
!$omp     parallel private(i,j,km,lm,k)
!$omp     do reduction(+:wk2)
          do i=1,n
            do j=1,n
              do km=1,n
                do lm=1,n
                  do k=1,n
                    wk2(i,j,km,lm)=wk2(i,j,km,lm)+cm(k,km)*wk1(i,j,k,lm)
                  end do
                end do
              end do
            end do
          end do
!$omp     end do
!$omp     end parallel
          wk1(:,:,:,:)=wk2(:,:,:,:)
          wk2(:,:,:,:)=0d0
          ! Transf.:(i,j,km,lm) -> (i,jm,km,lm)
!$omp     parallel private(i,jm,km,lm,j)
!$omp     do reduction(+:wk2)
          do i=1,n
            do jm=1,n
              do km=1,n
                do lm=1,n
                  do j=1,n
                    wk2(i,jm,km,lm)=wk2(i,jm,km,lm)
     &                +cm(j,jm)*wk1(i,j,km,lm)
                  end do
                end do
              end do
            end do
          end do
!$omp     end do
!$omp     end parallel
          ! Transf.:(i,jm,km,lm) -> (im,jm,km,lm)
          ab(:)=0d0
          do ijklm=1,n3
            call iup2_(ijklm,im,jm,km,lm)
            do i=1,n
              ab(ijklm)=ab(ijklm)+cm(i,im)*wk2(i,jm,km,lm)
            end do
          end do
        end subroutine

        subroutine mtrx_utau2u(cma,cmb,ab)
          ! Unitary transformation of a 2e matrix (4-index).
          implicit none
          real(8),intent(in)::cma(:,:),cmb(:,:)
          real(8),intent(inout)::ab(:,:)
          integer::n,n2,i,j,k,l,ij,kl,ijkl,im,jm,km,lm,ijm,klm
          real(8),allocatable::wk1(:,:,:,:),wk2(:,:,:,:)
          ! Get memory.
          n=SIZE(cma(:,1))
          n2=SIZE(ab(:,1))
          allocate(wk1(n,n,n,n))
          allocate(wk2(n,n,n,n))
          wk1(:,:,:,:)=0d0
          ! Transf.: a(i,j,k,l) -> a(i,j,k,lm)
!$omp     parallel private(i,j,k,l,lm,ijkl)
!$omp     do reduction(+:wk1)
          do i=1,n
            do j=1,n
              ij=ip_(i,j)
              do k=1,n
                do lm=1,n
                  do l=1,n
                    kl=ip_(k,l)
                    wk1(i,j,k,lm)=wk1(i,j,k,lm)
     &                +cmb(l,lm)*ab(ij,kl)
                  end do
                end do
              end do
            end do
          end do
!$omp     end do
!$omp     end parallel
          wk2(:,:,:,:)=0d0
          ! Transf.: a(i,j,k,lm) -> a(i,j,km,lm)
!$omp     parallel private(i,j,km,lm,k)
!$omp     do reduction(+:wk2)
          do i=1,n
            do j=1,n
              do km=1,n
                do lm=1,n
                  do k=1,n
                    wk2(i,j,km,lm)=wk2(i,j,km,lm)
     &                +cmb(k,km)*wk1(i,j,k,lm)
                  end do
                end do
              end do
            end do
          end do
!$omp     end do
!$omp     end parallel
          wk1(:,:,:,:)=wk2(:,:,:,:)
          wk2(:,:,:,:)=0d0
          ! Transf.:(i,j,km,lm) -> (i,jm,km,lm)
!$omp     parallel private(i,jm,km,lm,j)
!$omp     do reduction(+:wk2)
          do i=1,n
            do jm=1,n
              do km=1,n
                do lm=1,n
                  do j=1,n
                    wk2(i,jm,km,lm)=wk2(i,jm,km,lm)
     &                +cma(j,jm)*wk1(i,j,km,lm)
                  end do
                end do
              end do
            end do
          end do
!$omp     end do
!$omp     end parallel
          ! Transf.:(i,jm,km,lm) -> (im,jm,km,lm)
          ab(:,:)=0d0
!         do ijklm=1,n3
!           call iup2_(ijklm,im,jm,km,lm)
!           do i=1,n
!             ab(ijklm)=ab(ijklm)+cm(i,im)*wk2(i,jm,km,lm)
!           end do
!         end do
          do ijm=1,n2
            call iup_(ijm,im,jm)
            do klm=1,n2
              call iup_(klm,km,lm)
              do i=1,n
                ab(ijm,klm)=ab(ijm,klm)+cma(i,im)*wk2(i,jm,km,lm)
              end do
            end do
          end do
        end subroutine

        subroutine mtrx_eigsp(ap,e,c)
          ! Real symmetric matrix diagonalization (packed strage).
          implicit none
          real(8),intent(inout)::ap(:)
          real(8),intent(out)::e(:),c(:,:)
          integer::n,lwork,liwork,info
          integer,allocatable::iwork(:)
          real(8),allocatable::work(:)
          n=size(e)
          lwork=1+6*n+n**2
          liwork=3+5*n
          allocate(work(lwork),iwork(liwork))
          call DSPEVD('V','U',n,ap,e,c,n,work,lwork,iwork,liwork,info)
        end subroutine

        subroutine mtrx_uc(u,c)
          ! Unitary transformation of orbitals.
          implicit none
          real(8),intent(in)::u(:,:)
          real(8),intent(inout)::c(:,:)
          real(8),allocatable::vt(:,:)
          integer::n
          n=SIZE(c(:,1))
          allocate(vt(n,n))
          vt=MATMUL(u,c)
          c(:,:)=vt(:,:)
        end subroutine

        subroutine mtrx_expv(pb,hb,val)
          ! Expectation value of the 1e operator.
          implicit none
          real(8),intent(in)::pb(:),hb(:)
          real(8),intent(out)::val
          integer::i,j,ij,nb2
          real(8),allocatable::vt(:)
          nb2=SIZE(pb)
          allocate(vt(nb2))
          do ij=1,nb2
            vt(ij)=pb(ij)*2d0
            call iup_(ij,i,j)
            if (i==j) vt(ij)=vt(ij)/2d0
          end do
          val=SUM(vt(:)*hb(:))
        end subroutine

        subroutine mtrx_expv2(pb,tb,val)
          ! Expectation value of the 2e operator.
          implicit none
          real(8),intent(in)::pb(:),tb(:)
          real(8),intent(out)::val
          integer::i,j,k,l,ij,kl,ijkl,nb3
          real(8),allocatable::vt(:)
          nb3=SIZE(tb)
          allocate(vt(nb3))
          do ijkl=1,nb3
            vt(ijkl)=pb(ijkl)*8d0
            call iup_(ijkl,ij,kl)
            call iup_(ij,i,j)
            call iup_(kl,k,l)
            if (i==j) vt(ijkl)=vt(ijkl)/2d0
            if (k==l) vt(ijkl)=vt(ijkl)/2d0
            if (ij==kl) vt(ijkl)=vt(ijkl)/2d0
          end do
          val=SUM(vt(:)*tb(:))
        end subroutine

        integer function ip_(i,j) result(ij)
          ! Canonical index packing for a symmetric matrix (1e).
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

        subroutine iup2_(ijkl,i,j,k,l)
          ! Unpack the canonical index for a symmetric matrix (2e).
          implicit none
          integer,intent(in)::ijkl
          integer,intent(out)::i,j,k,l
          integer::ij,kl
          call iup_(ijkl,ij,kl)
          call iup_(ij,i,j)
          call iup_(kl,k,l)
        end subroutine

      end module
