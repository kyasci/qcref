      module uga
        ! Unitary group approach configuration interaction (UGA-CI).
        ! Reference:
        !   J. Paldus, J. Chem. Phys. 61, 5321 (1974).
        implicit none
        private
        real(8),parameter::TOL_=1d-10,TOL_CC_=1d-5
        logical::init_=.false.
        integer::lbfnz_=1000000,nnz_
        integer,allocatable::lnz1_(:,:),idxnz_(:),mcr_(:),mrc_(:)
        real(8),allocatable::vnz1_(:)
        !---------------------------------------------------------------
        public::
     &    uga_init,
     &    uga_del,
     &    uga_pltab,
     &    uga_dim,
     &    uga_conf,
     &    uga_cc1,
     &    uga_cc2,
     &    uga_seti
        !---------------------------------------------------------------

      contains

        subroutine uga_pltab(nea,nac,mult,lpaldus)
          ! Make a list of Paldus tables.
          implicit none
          integer,intent(in)::nea,nac,mult
          integer,intent(out),allocatable::lpaldus(:,:,:)
          integer::ma,mb,mc,np,lt(3,4),ip,im,ih,it,nc,nt
          integer,allocatable::mwork(:,:,:)
          ! Generate Paldus table for a complete active space.
          mb=mult-1
          ma=(nea-mb)/2
          mc=nac-ma-mb
          call dimcas_(nac,ma,mb,mc,nc)
          allocate(lpaldus(3,nac,nc))
          allocate(mwork(3,nac,nc))
          lpaldus(:,:,:)=0
          mwork(:,:,:)=0
          lpaldus(:,1,1)=(/ma,mb,mc/)
          np=1
          do im=1,nac-1
            ih=0
            do ip=1,np
              ma=lpaldus(1,im,ip)
              mb=lpaldus(2,im,ip)
              mc=lpaldus(3,im,ip)
              call table1_(ma,mb,mc,nt,lt)
              do it=1,nt
                ih=ih+1
                mwork(:,:,ih)=lpaldus(:,:,ip)
                mwork(:,im+1,ih)=(/ma,mb,mc/)-lt(:,it)
              end do
              np=ih
            end do
            lpaldus(:,:,:)=mwork(:,:,:)
          end do
        end subroutine

        subroutine uga_init(ne,nfc,nac,mult,nx)
          ! Prepare the 1e coupling constants for CAS.
          implicit none
          integer,intent(in)::ne,nfc,nac,mult,nx
          integer::nea
          integer,allocatable::lpaldus(:,:,:)
          ! Clean up the data for the previous calculation.
          if (init_) then
            call uga_del()
          end if
          init_=.true.
          ! Make a list of Paldus tables.
          nea=ne-2*nfc
          call uga_pltab(nea,nac,mult,lpaldus)
          ! Make a list of index for RAS and CAS.
          call mkras_(nea,mult,nx,lpaldus)
          ! Make nonzero values of the 1e coupling constants (CAS).
          call mknz1_(nfc,nac,lpaldus)
          ! Index table for the quick access to 1e coupling constant.
          call mkidx_()
        end subroutine

        subroutine uga_del()
          ! Destructor.
          implicit none
          deallocate(lnz1_)
          deallocate(vnz1_)
          deallocate(mcr_)
          deallocate(mrc_)
          deallocate(idxnz_)
          init_=.false.
        end subroutine

        subroutine uga_dim(nr)
          ! Return the number of RAS CSFs.
          implicit none
          integer,intent(out)::nr
          nr=SIZE(mrc_)
        end subroutine

        subroutine uga_conf(ir,ls)
          ! Get electronic configuration in the active space.
          integer,intent(in)::ir
          integer,intent(out)::ls(:)
          integer::im,ic
          real(8)::a1
          ls(:)=0
          ic=mrc_(ir)
          do im=1,SIZE(ls)
            call loadcc1_(ic,ic,im,im,a1)
            ls(im)=NINT(a1)
          end do
        end subroutine

        subroutine uga_seti(key,num)
          ! Interface for change the integer parameters.
          implicit none
          character(*),intent(in)::key
          integer,intent(in)::num
          if (key=='lbfnz') then
            lbfnz_=num
          else
            write(*,'("error: unrecognized key: ",a)') key
            stop 1
          end if
        end subroutine

        subroutine mkras_(nea,mult,nx,lpaldus)
          ! Make a list of index for the relationship of RAS and CAS.
          implicit none
          integer,intent(in)::nea,mult,nx,lpaldus(:,:,:)
          integer::ic,nt,im,nc,nr
          integer,allocatable::mwork(:)
          real(8)::vt
          nc=SIZE(lpaldus(1,1,:))
          allocate(mwork(nc))
          allocate(mcr_(nc))
          mcr_(:)=-1
          nr=0
          do ic=1,nc
            nt=0
            do im=1,(nea+(mult-1))/2
              call genweight_(ic,im,lpaldus,vt)
              nt=nt+NINT(vt)
            end do
            if (nt>=(nea-nx)) then
              nr=nr+1
              mwork(nr)=ic
              mcr_(ic)=nr
            end if
          end do
          allocate(mrc_(nr))
          mrc_(:nr)=mwork(:nr)
        end subroutine

        subroutine mknz1_(nfc,nac,lpaldus)
!$        use omp_lib
          ! Nonzero values for 1e coupling constants (CAS).
          implicit none
          integer,intent(in)::nfc,nac,lpaldus(:,:,:)
          integer::ic,jc,im,jm,nc,inz
          real(8)::a1
          integer,allocatable::libuf(:,:,:),lnnz(:)
          real(8),allocatable::vibuf(:,:)
          ! Get memory.
          nc=SIZE(lpaldus(1,1,:))
          allocate(lnnz(nc),libuf(4,lbfnz_,nc),vibuf(lbfnz_,nc))
          lnnz(:)=0
          ! Evaluate non-zero coupling constants.
!$omp     parallel private(ic,jc,im,jm,a1)
!$omp     do
          do ic=1,nc
            do jc=1,ic
              do im=nfc+1,nfc+nac
                do jm=nfc+1,im
                  call evalcc1_(ic,jc,im,jm,nfc,lpaldus,a1)
                  if (ABS(a1)>TOL_CC_) then
                    lnnz(ic)=lnnz(ic)+1
                    if (lnnz(ic)>lbfnz_) then
                      write(*,'("error: nnz_>lbfnz_: ci.mknzcas_()")')
                      stop 1
                    end if
                    libuf(:,lnnz(ic),ic)=(/ic,jc,im,jm/)
                    vibuf(lnnz(ic),ic)=a1
                  end if
                end do
              end do
            end do
          end do
!$omp     end do
!$omp     end parallel
          ! Get memroy for the module variables.
          nnz_=SUM(lnnz)
          allocate(lnz1_(4,nnz_),vnz1_(nnz_))
          ! Set return values.
          inz=0
          do ic=1,nc
            do jc=1,lnnz(ic)
              inz=inz+1
              lnz1_(:,inz)=libuf(:,jc,ic)
              vnz1_(inz)=vibuf(jc,ic)
            end do
          end do
        end subroutine

        subroutine mkidx_()
          ! Index table for the quick access to 1e coupling constant.
          implicit none
          integer::i,ic,jc,ijc,ijc0,nc,ncc2
          nc=SIZE(mcr_)
          ncc2=ip_(nc,nc)
          allocate(idxnz_(ncc2))
          idxnz_(:)=0
          ijc0=0
          do i=1,nnz_
            ic=lnz1_(1,i)
            jc=lnz1_(2,i)
            ijc=ip_(ic,jc)
            if (ijc0/=ijc) then
              idxnz_(ijc)=i
              ijc0=ijc
            end if
          end do
        end subroutine

        subroutine uga_cc1(lc1,vc1)
          ! 1e coupling constatns in the compressed matrix form.
          implicit none
          integer,intent(out),allocatable::lc1(:,:)
          real(8),intent(out),allocatable::vc1(:)
          integer::nnz1,ir,inz
          integer,allocatable::lbuf(:,:)
          real(8),allocatable::vbuf(:)
          ! Get memory for the working space.
          allocate(lbuf(4,lbfnz_),vbuf(lbfnz_))
          ! Pick up nonzero values for RAS from CAS.
          nnz1=0
          do inz=1,nnz_
            ir=mcr_(lnz1_(1,inz))
            if (ir>0) then
              nnz1=nnz1+1
              lbuf(1,nnz1)=ir
              lbuf(2,nnz1)=mcr_(lnz1_(2,inz))
              lbuf(3,nnz1)=lnz1_(3,inz)
              lbuf(4,nnz1)=lnz1_(4,inz)
              vbuf(nnz1)=vnz1_(inz)
            end if
          end do
          ! Get memory and set the picked values.
          allocate(lc1(4,nnz1),vc1(nnz1))
          lc1(:,:nnz1)=lbuf(:,:nnz1)
          vc1(:nnz1)=vbuf(:nnz1)
        end subroutine

        subroutine uga_cc2(lc2,vc2)
          ! 2e coupling constants in the compressed matrix form.
          implicit none
          integer,intent(out),allocatable::lc2(:,:)
          real(8),intent(out),allocatable::vc2(:)
          integer::nnz2
          integer,allocatable::lbuf(:,:)
          real(8),allocatable::vbuf(:)
          ! Get memory for the working space.
          allocate(lbuf(6,lbfnz_),vbuf(lbfnz_))
          ! Evaluate (E_ij E_kl) part with resolution of identity.
          nnz2=0
          call cc2eijekl_(nnz2,lbuf,vbuf)
          ! Evaluate (delta_kj E_il) part.
          call cc2dkjeil_(nnz2,lbuf,vbuf)
          ! Get memory and set the values evaluated above.
          allocate(lc2(6,nnz2),vc2(nnz2))
          lc2(:,:nnz2)=lbuf(:,:nnz2)
          vc2(:nnz2)=vbuf(:nnz2)
        end subroutine

        subroutine cc2eijekl_(nnz2,lbuf,vbuf)
          ! Evaluate (E_ij E_kl) part of 2e coupling constant.
          implicit none
          integer,intent(inout)::nnz2,lbuf(:,:)
          real(8),intent(inout)::vbuf(:)
          integer::istart,istop,ic,jc,kc,nc,ir,jr,nr,im,jm,km,lm,inz
          real(8)::a1,a2
          integer,allocatable::libuf(:,:,:),lnnz(:)
          real(8),allocatable::vibuf(:,:)
          istart=MINVAL(lnz1_(3,:))
          istop=MAXVAL(lnz1_(3,:))
          nc=SIZE(mcr_)
          allocate(lnnz(nc),libuf(6,lbfnz_,nc),vibuf(lbfnz_,nc))
          lnnz(:)=0
          nr=SIZE(mrc_)
!$omp     parallel private(ir,ic,im,jm,kc,jr,km,lm,jc,a1,a2)
!$omp     do
          do ir=1,nr
            ic=mrc_(ir)
            do im=istart,istop
              do jm=istart,istop
                do kc=1,nc
                  call loadcc1_(ic,kc,im,jm,a1)
                  if (ABS(a1)<TOL_CC_) then
                    cycle
                  end if
                  do jr=1,ir
                    do km=istart,istop
                      do lm=istart,istop
                        jc=mrc_(jr)
                        call loadcc1_(kc,jc,km,lm,a2)
                        if (ABS(a2)<TOL_CC_) then
                          cycle
                        end if
                        lnnz(ic)=lnnz(ic)+1
                        if (nnz2>lbfnz_) then
                          write(*,'("error: nnz2>lbfnz_: ci.mknz2_()")')
                          stop 1
                        end if
                        libuf(:,lnnz(ic),ic)=(/ir,jr,im,jm,km,lm/)
                        vibuf(lnnz(ic),ic)=a1*a2*0.5d0
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
!$omp     end do
!$omp     end parallel
          ! Set return values.
          nnz2=SUM(lnnz)
          inz=0
          do ic=1,nc
            do jc=1,lnnz(ic)
              inz=inz+1
              lbuf(:,inz)=libuf(:,jc,ic)
              vbuf(inz)=vibuf(jc,ic)
            end do
          end do
        end subroutine

        subroutine cc2dkjeil_(nnz2,lbuf,vbuf)
          ! Evaluate (delta_kj E_il) part of the 2e coupling constant.
          implicit none
          integer,intent(inout)::nnz2,lbuf(:,:)
          real(8),intent(inout)::vbuf(:)
          integer::istart,istop,ic,jc,nc,ir,jr,nr,im,jm,lm
          real(8)::a1
          istart=MINVAL(lnz1_(3,:))
          istop=MAXVAL(lnz1_(3,:))
          nc=SIZE(mcr_)
          nr=SIZE(mrc_)
          do ir=1,nr
            ic=mrc_(ir)
            do jr=1,ir
              jc=mrc_(jr)
              do im=istart,istop
                do lm=istart,istop
                  call loadcc1_(ic,jc,im,lm,a1)
                  if (ABS(a1)>TOL_CC_) then
                    do jm=istart,istop
                      nnz2=nnz2+1
                      if (nnz2>lbfnz_) then
                        write(*,'("error: nnz2>lbfnz_: ci.mknz2_()")')
                        stop 1
                      end if
                      lbuf(:,nnz2)=(/ir,jr,im,jm,jm,lm/)
                      vbuf(nnz2)=-a1*0.5d0
                    end do
                  end if
                end do
              end do
            end do
          end do
        end subroutine

        recursive subroutine loadcc1_(ic,jc,im,jm,val)
          ! Load the 1e coupling constant using quick access index.
          implicit none
          integer,intent(in)::ic,jc,im,jm
          real(8),intent(out)::val
          logical::match
          integer::i,istart,istop
          val=0d0
          if (ic>=jc) then
            call rangennz1_(ic,jc,istart,istop)
            if (istart>0) then
              do i=istart,istop
                match=lnz1_(3,i)==im.and.lnz1_(4,i)==jm
                if (match) then
                  val=vnz1_(i)
                  exit
                end if
                if (lnz1_(3,i)>im) then
                  exit
                end if
              end do
            end if
          else
            call loadcc1_(jc,ic,jm,im,val)
          end if
        end subroutine

        subroutine rangennz1_(ic,jc,istart,istop)
          ! Returns the range of 1e coupling const. to be scanned.
          implicit none
          integer,intent(in)::ic,jc
          integer,intent(out)::istart,istop
          integer::MAXRUN=10000
          integer::ijc,i
          ijc=ip_(ic,jc)
          if (idxnz_(ijc)==0) then
            istart=0
            istop=0
          else
            istart=idxnz_(ijc)
            do i=1,MAXRUN
              if (ijc+i>SIZE(idxnz_)) then
                istop=SIZE(vnz1_)
                exit
              else if (idxnz_(ijc+i)>0) then
                istop=idxnz_(ijc+i)-1
                exit
              end if
            end do
          end if
        end subroutine

        logical function is_possible_(ic,jc,im,jm,lpaldus) result(res)
          ! Apply the selection rules of 1e constants.
          implicit none
          integer,intent(in)::ic,jc,im,jm,lpaldus(:,:,:)
          real(8)::v1,v2
          res=.false.
          call genweight_(jc,jm,lpaldus,v1)
          call genweight_(ic,jm,lpaldus,v2)
          if ((v1-v2-1d0)**2>TOL_) then
            return
          end if
          call genweight_(jc,im,lpaldus,v1)
          call genweight_(ic,im,lpaldus,v2)
          if ((v1-v2+1d0)**2>TOL_) then
            return
          end if
          res=.true.
        end function

        recursive subroutine evalcc1_(ic,jc,im,jm,nfc,lpaldus,val)
          ! Evaluate CAS 1e constant.
          implicit none
          integer,intent(in)::ic,jc,im,jm,nfc,lpaldus(:,:,:)
          real(8),intent(out)::val
          val=0d0
          ! Weight generator.
          if (ic==jc) then
            if (im==jm) then
              call genweight_(ic,im-nfc,lpaldus,val)
            end if
          ! Raising generator.
          else if (ic<jc) then
            if (im+1==jm) then
              call genbasicr_(ic,jc,im-nfc,lpaldus,val)
            else if (im+1<jm) then
              call genraise_(ic,jc,im-nfc,jm-nfc,lpaldus,val)
            end if
          ! Lowering generator.
          else
            call evalcc1_(jc,ic,jm,im,nfc,lpaldus,val)
          end if
        end subroutine

        subroutine genbasicr_(ic,jc,im,lpaldus,val)
          ! Basic raising generator for the case (im+1=jm).
          integer,intent(in)::ic,jc,im,lpaldus(:,:,:)
          real(8),intent(out)::val
          integer::nac
          integer,allocatable::diff(:,:)
          real(8)::ai,bi,ci,aip,bip,cip,aim,bim,cim
          if (ic>=jc) then
            write(*,'("error: ic>=jc: genbasicr_()")')
            write(*,'(2i5)')ic,jc
            stop 1
          end if
          val=0d0
          if (.not.is_possible_(ic,jc,im,im+1,lpaldus)) then
            return
          end if
          nac=SIZE(lpaldus(1,:,1))
          allocate(diff(3,nac))
          ! Start evaluation.
          ai=lpaldus(1,nac+1-im,jc)
          bi=lpaldus(2,nac+1-im,jc)
          ci=lpaldus(3,nac+1-im,jc)
          aip=lpaldus(1,nac+1-(im+1),jc)
          bip=lpaldus(2,nac+1-(im+1),jc)
          cip=lpaldus(3,nac+1-(im+1),jc)
          aim=0d0
          bim=0d0
          cim=0d0
          if (im>1) then
            aim=lpaldus(1,nac+1-(im-1),jc)
            bim=lpaldus(2,nac+1-(im-1),jc)
            cim=lpaldus(3,nac+1-(im-1),jc)
          end if
          if (ai==aim.and.aim==aip-1.and.ABS(bi)>TOL_) then
            diff(:,:)=lpaldus(:,:,jc)
            diff(:,nac+1-im)=diff(:,nac+1-im)+(/1,-1,0/)
            diff(:,:)=diff(:,:)-lpaldus(:,:,ic)
            if (SUM(diff**2)<TOL_) then
              val=((bi*(bi+1d0))/((bip+1d0)*(bim+1d0)))**0.5d0
            end if
          end if
          if (ci==cip.and.cip==cim+1.and.ABS(ci)>TOL_) then
            diff(:,:)=lpaldus(:,:,jc)
            diff(:,nac+1-im)=diff(:,nac+1-im)+(/0,1,-1/)
            diff(:,:)=diff(:,:)-lpaldus(:,:,ic)
            if (SUM(diff**2)<TOL_) then
              val=(((bi+1d0)*(bi+2d0))/((bip+1d0)*(bim+1d0)))**0.5d0
            end if
          end if
        end subroutine

        recursive subroutine genraise_(ic,jc,im,jm,lpaldus,val)
          ! Raising generator for the case (im+1<jm).
          implicit none
          integer,intent(in)::ic,jc,im,jm,lpaldus(:,:,:)
          real(8),intent(out)::val
          real(8)::v11,v12,v21,v22
          integer::kc
          val=0d0
          if (.not.is_possible_(ic,jc,im,jm,lpaldus)) then
            return
          end if
          do kc=ic+1,jc-1
            call genbasicr_(kc,jc,jm-1,lpaldus,v12)
            call genbasicr_(ic,kc,jm-1,lpaldus,v21)
            if (abs(v12)>TOL_.or.abs(v21)>TOL_) then
              if (im+2==jm) then
                call genbasicr_(ic,kc,im,lpaldus,v11)
                call genbasicr_(kc,jc,im,lpaldus,v22)
              else
                call genraise_(ic,kc,im,jm-1,lpaldus,v11)
                call genraise_(kc,jc,im,jm-1,lpaldus,v22)
              end if
              val=val+v11*v12-v21*v22
            end if
          end do
        end subroutine

        subroutine genweight_(ic,im,lpaldus,val)
          ! Weight generator for the case (im=jm).
          integer,intent(in)::ic,im,lpaldus(:,:,:)
          real(8),intent(out)::val
          real(8)::ai,bi,aim,bim
          integer::nac
          nac=SIZE(lpaldus(1,:,1))
          ai=lpaldus(1,nac+1-im,ic)
          bi=lpaldus(2,nac+1-im,ic)
          aim=0d0
          bim=0d0
          if (im>1) then
            aim=lpaldus(1,nac+1-(im-1),ic)
            bim=lpaldus(2,nac+1-(im-1),ic)
          end if
          val=2d0*(ai-aim)+(bi-bim)
        end subroutine

        subroutine table1_(ma,mb,mc,n,num)
          ! Get the subtration table for the (i-1)'s row.
          implicit none
          integer,intent(in)::ma,mb,mc
          integer,intent(out)::n,num(:,:)
          num(:,:)=0
          if (ma/=0.and.mb/=0.and.mc/=0) then
            n=4
            num(:,1)=(/0,0,1/)
            num(:,2)=(/0,1,0/)
            num(:,3)=(/1,-1,1/)
            num(:,4)=(/1,0,0/)
          else if (ma/=0.and.mb/=0.and.mc==0) then
            n=2
            num(:,1)=(/0,1,0/)
            num(:,2)=(/1,0,0/)
          else if (ma/=0.and.mb==0.and.mc/=0) then
            n=3
            num(:,1)=(/0,0,1/)
            num(:,2)=(/1,-1,1/)
            num(:,3)=(/1,0,0/)
          else if (ma==0.and.mb/=0.and.mc/=0) then
            n=2
            num(:,1)=(/0,0,1/)
            num(:,2)=(/0,1,0/)
          else if (ma/=0.and.mb==0.and.mc==0) then
            n=1
            num(:,1)=(/1,0,0/)
          else if (ma==0.and.mb/=0.and.mc==0) then
            n=1
            num(:,1)=(/0,1,0/)
          else if (ma==0.and.mb==0.and.mc/=0) then
            n=1
            num(:,1)=(/0,0,1/)
          else
            write(*,'("error: unclassified: ci.ndim_()")')
            stop 1
          end if
        end subroutine

        subroutine dimcas_(n,ma,mb,mc,num)
          ! Get the CAS-CI dimension.
          implicit none
          integer,intent(in)::n,ma,mb,mc
          integer,intent(out)::num
          num=(mb+1)
          num=num*binom_(n+1,ma)
          num=num*binom_(n+1,mc)
          num=num/(n+1)
        end subroutine

        integer function binom_(n,k) result(m)
          ! Binomial constant.
          implicit none
          integer,intent(in)::n,k
          integer::ki,i
          ki=min(k,(n-k))
          m=1
          do i=1,ki
            m=m*(n+1-i)
          end do
          m=m/factorial_(ki)
        end function

        recursive integer function factorial_(n) result(m)
          ! Recursive implementation of the factorial calculations.
          implicit none
          integer,intent(in)::n
          if (n<=1) then
            m=1
          else
            m=n*factorial_(n-1)
          end if
        end function

        integer function ip_(i,j) result(ij)
          ! Canonical index packing for a symmetric matrix.
          implicit none
          integer,intent(in)::i,j
          integer::ii,jj
          ii=max(i,j)
          jj=min(i,j)
          ij=(ii-1)*ii/2+jj
        end function

      end module
