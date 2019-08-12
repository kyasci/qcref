      module uga
        ! Unitary group approach configuration interaction (UGA-CI).
        !
        ! Reference:
        !   J. Paldus, J. Chem. Phys. 61, 5321 (1974).
        !
        implicit none
        private
        real(8),parameter::TOL_=1d-10,TOL_CC_=1d-5
        logical::init_=.false.
        integer::lbfnz_=10000000,n1_
        integer,allocatable::l1_(:,:),mcr_(:),mrc_(:),mcf_(:,:)
        real(8),allocatable::v1_(:)
        ! Used for profiling.
        !---------------------------------------------------------------
        public::
     &    uga_init,
     &    uga_del,
     &    uga_pltab,
     &    uga_dim,
     &    uga_cc1,
     &    uga_cc2,
     &    uga_ecf,
     &    uga_seti
        !---------------------------------------------------------------

      contains

        subroutine uga_init(ne,nfc,nac,mult,nx)
          ! Prepare the 1e coupling constants for CAS.
          implicit none
          integer,intent(in)::ne,nfc,nac,mult,nx
          integer::nea
          integer,allocatable::lplds(:,:,:)
          ! Clean up the data for the previous calculation.
          if (init_) then
            call uga_del()
          end if
          init_=.true.
          ! Make a list of Paldus tables.
          nea=ne-2*nfc
          call uga_pltab(nea,nac,mult,lplds)
          ! Make a list of configuration.
          call mkconf_(lplds)
          ! Make a list of index for RAS and CAS.
          call mkras_(nx,lplds)
          ! Make nonzero values of the 1e coupling constants (CAS).
          call mknz1_(nfc,nac,lplds)
          ! Index table for the quick access to 1e coupling constant.
        end subroutine

        subroutine uga_ecf(lecf)
          ! Return the list of electronic configurations.
          implicit none
          integer,intent(out),allocatable::lecf(:,:)
          integer::nac,nc
          nc=SIZE(mcf_(1,:))
          nac=SIZE(mcf_(:,1))
          allocate(lecf(nac,nc))
          lecf(:,:)=mcf_(:,:)
        end subroutine

        subroutine mkconf_(lplds)
          ! Make a list of electronic configurations.
          implicit none
          integer,intent(in)::lplds(:,:,:)
          integer::nc,nac,ic,im,md(3),ldwn(3)
          nc=SIZE(lplds(1,1,:))
          nac=SIZE(lplds(1,:,1))
          ! Get memory.
          allocate(mcf_(nac,nc))
          mcf_(:,:)=0
          ! Electronic configuration:
          !     2: doubly occupied
          !     1: alpha spin occupied.
          !    -1: beta spin occupied.
          !     0: unoccupied.
          do ic=1,nc
            do im=1,nac
              if (im==nac) then
                ldwn(:)=(/0,0,0/)
              else
                ldwn(:)=lplds(:,im+1,ic)
              end if
              md(:)=lplds(:,im,ic)-ldwn(:)
              if (md(1)==0.and.md(2)==0.and.md(3)==1) then
                mcf_(nac-im+1,ic)=0
              else if (md(1)==0.and.md(2)==1.and.md(3)==0) then
                mcf_(nac-im+1,ic)=1
              else if (md(1)==1.and.md(2)==-1.and.md(3)==1) then
                mcf_(nac-im+1,ic)=-1
              else if (md(1)==1.and.md(2)==0.and.md(3)==0) then
                mcf_(nac-im+1,ic)=2
              end if
            end do
          end do
        end subroutine

        subroutine uga_pltab(nea,nac,mult,lplds)
          ! Make a list of Paldus tables.
          implicit none
          integer,intent(in)::nea,nac,mult
          integer,intent(out),allocatable::lplds(:,:,:)
          integer::ma,mb,mc,np,lt(3,4),ipt,im,ih,it,nc,nt
          integer,allocatable::mwork(:,:,:)
          ! Generate Paldus table for a complete active space.
          mb=mult-1
          ma=(nea-mb)/2
          mc=nac-ma-mb
          call dimcas_(nac,ma,mb,mc,nc)
          allocate(lplds(3,nac,nc))
          allocate(mwork(3,nac,nc))
          lplds(:,:,:)=0
          mwork(:,:,:)=0
          lplds(:,1,1)=(/ma,mb,mc/)
          np=1
          do im=1,nac-1
            ih=0
            do ipt=1,np
              ma=lplds(1,im,ipt)
              mb=lplds(2,im,ipt)
              mc=lplds(3,im,ipt)
              call table1_(ma,mb,mc,nt,lt)
              do it=1,nt
                ih=ih+1
                mwork(:,:,ih)=lplds(:,:,ipt)
                mwork(:,im+1,ih)=(/ma,mb,mc/)-lt(:,it)
              end do
              np=ih
            end do
            lplds(:,:,:)=mwork(:,:,:)
          end do
        end subroutine

        subroutine uga_del()
          ! Destructor.
          implicit none
          deallocate(l1_)
          deallocate(v1_)
          deallocate(mcr_)
          deallocate(mrc_)
          deallocate(mcf_)
          init_=.false.
        end subroutine

        subroutine uga_dim(nr)
          ! Return the number of RAS CSFs.
          implicit none
          integer,intent(out)::nr
          nr=SIZE(mrc_)
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

        subroutine mkras_(nx,lplds)
          ! Make a list of index for the relationship of RAS and CAS.
          implicit none
          integer,intent(in)::nx,lplds(:,:,:)
          integer::ic,nt,im,nm,nc,nr,nd
          integer,allocatable::mwork(:)
          ! Get memoery.
          nc=SIZE(lplds(1,1,:))
          nm=SIZE(lplds(1,:,1))
          allocate(mwork(nc))
          allocate(mcr_(nc))
          ! Set -1 to show there is no corresponding index.
          mcr_(:)=-1
          nr=0
          ! Evaluation.
          do ic=1,nc
            nd=0
            do im=1,nm
              nt=ABS(mcf_(im,ic)-mcf_(im,1))
              if (mcf_(im,ic)==-1.and.mcf_(im,1)==2) then
                nt=1
              else if (mcf_(im,ic)==2.and.mcf_(im,1)==-1) then
                nt=1
              end if
              nd=nd+nt
            end do
            nd=nd/2
            if (nd<=nx) then
              nr=nr+1
              mwork(nr)=ic
              mcr_(ic)=nr
            end if
          end do
          allocate(mrc_(nr))
          mrc_(:nr)=mwork(:nr)
        end subroutine

        subroutine uga_cc1(lc1,vc1)
          ! 1e coupling constatns in the compressed matrix form.
          implicit none
          integer,intent(out),allocatable::lc1(:,:)
          real(8),intent(out),allocatable::vc1(:)
          integer::nnz1,ir,inz
          integer,allocatable::lbf(:,:)
          real(8),allocatable::vbf(:)
          ! Get memory for the working space.
          allocate(lbf(4,lbfnz_))
          allocate(vbf(lbfnz_))
          ! Pick up nonzero values for RAS from CAS.
          nnz1=0
          do inz=1,n1_
            ir=mcr_(l1_(1,inz))
            if (ir>0) then
              nnz1=nnz1+1
              lbf(1,nnz1)=ir
              lbf(2,nnz1)=mcr_(l1_(2,inz))
              lbf(3,nnz1)=l1_(3,inz)
              lbf(4,nnz1)=l1_(4,inz)
              vbf(nnz1)=v1_(inz)
            end if
          end do
          ! Get memory and set the picked values.
          allocate(lc1(4,nnz1),vc1(nnz1))
          lc1(:,:nnz1)=lbf(:,:nnz1)
          vc1(:nnz1)=vbf(:nnz1)
        end subroutine

        subroutine uga_cc2(lc2,vc2)
          ! 2e coupling constants in the compressed matrix form.
          implicit none
          integer,intent(out),allocatable::lc2(:,:)
          real(8),intent(out),allocatable::vc2(:)
          integer::n2
          integer,allocatable::lbf(:,:)
          real(8),allocatable::vbf(:)
          ! Get memory for the working space.
          allocate(lbf(6,lbfnz_))
          allocate(vbf(lbfnz_))
          ! Evaluate (E_ij E_kl) part with resolution of identity.
          n2=0
          call cc2eijekl_(n2,lbf,vbf)
          ! Evaluate (delta_kj E_il) part.
          call cc2dkjeil_(n2,lbf,vbf)
          ! Get memory and set the values evaluated above.
          allocate(lc2(6,n2))
          allocate(vc2(n2))
          lc2(:,:n2)=lbf(:,:n2)
          vc2(:n2)=vbf(:n2)
        end subroutine

        subroutine addnz1_(ic,jc,im,jm,val,n,l,v)
          ! Add indies and value of nonzero 1e constants to the list.
          implicit none
          integer,intent(in)::ic,jc,im,jm
          real(8),intent(in)::val
          integer,intent(inout)::n,l(:,:)
          real(8),intent(inout)::v(:)
          n=n+1
          l(:,n)=(/ic,jc,im,jm/)
          v(n)=val 
        end subroutine

        subroutine add1weight_(nfc,nac,lplds,nbf,lbf,vbf)
          ! Weight generator contributions to 1e constants (CAS).
          implicit none
          integer,intent(in)::nfc,nac,lplds(:,:,:)
          integer,intent(inout)::nbf,lbf(:,:)
          real(8),intent(inout)::vbf(:)
          integer::im,ic,nc
          real(8)::a1
          nc=SIZE(lplds(1,1,:))
          do im=1,nac
            do ic=1,nc
              a1=ABS(mcf_(im,ic))
              if (ABS(a1)>TOL_CC_) then
                call addnz1_(ic,ic,nfc+im,nfc+im,a1,nbf,lbf,vbf)
              end if
            end do
          end do
        end subroutine

        subroutine add1basic_(nfc,nac,lplds,nbf,lbf,vbf)
          ! Basic generator contributions to 1e constants (CAS).
          integer,intent(in)::nfc,nac,lplds(:,:,:)
          integer,intent(inout)::nbf,lbf(:,:)
          real(8),intent(inout)::vbf(:)
          integer::im,ic,jc,nc
          real(8)::a1
          nc=SIZE(lplds(1,1,:))
          do im=1,nac-1
            do jc=1,nc
              ! Type A
              call genbasica_(jc,im,lplds,ic,a1)
              if (ABS(a1)>TOL_CC_) then
                call addnz1_(jc,ic,nfc+im+1,nfc+im,a1,nbf,lbf,vbf)
              end if
              call genbasicb_(jc,im,lplds,ic,a1)
              if (ABS(a1)>TOL_CC_) then
                call addnz1_(jc,ic,nfc+im+1,nfc+im,a1,nbf,lbf,vbf)
              end if
            end do
          end do
        end subroutine

        subroutine mountbt_(nac,nfc,nbf,lbf,vbf,nbr,lbr,ibr,vbr)
          ! Mount the nonzero values of the basic generator.
          implicit none
          integer,intent(in)::nac,nfc,nbf,lbf(:,:)
          real(8),intent(in)::vbf(:)
          integer,intent(out)::nbr
          integer,intent(out),allocatable::lbr(:,:),ibr(:)
          real(8),intent(out),allocatable::vbr(:)
          integer::it,ic,jc,im,jm
          real(8)::a1
          allocate(lbr(3,lbfnz_))
          allocate(ibr(nac))
          allocate(vbr(lbfnz_))
          nbr=0
          ibr(:)=0
          do it=1,nbf
            jc=lbf(1,it)
            ic=lbf(2,it)
            jm=lbf(3,it)-nfc
            im=lbf(4,it)-nfc
            a1=vbf(it)
            if (im+1==jm) then
              nbr=nbr+1
              lbr(:,nbr)=(/im,ic,jc/)
              vbr(nbr)=a1
              if (ibr(im)==0) then
                ibr(im)=nbr
              end if
            end if
          end do
        end subroutine

        subroutine cpbrnz_(im,nbr,lbr,vbr,gr)
          ! copy the nonzero values of basic generator matrix elements.
          implicit none
          integer,intent(in)::im,nbr,lbr(:,:)
          real(8),intent(in)::vbr(:)
          real(8),intent(out)::gr(:)
          integer::it,ijc
          gr(:)=0d0
          do it=1,nbr
            if (lbr(1,it)==im) then
              ijc=ip_(lbr(2,it),lbr(3,it))
              gr(ijc)=vbr(it)
            end if
            if (lbr(1,it)>im) then
              exit
            end if
          end do
        end subroutine

        subroutine evalgr_(jm,gr,nbr,ibr,lbr,vbr)
          ! Evaluate raising generator matrix elements for the given jm.
          implicit none
          integer,intent(in)::jm,nbr,ibr(:),lbr(:,:)
          real(8),intent(in)::vbr(:)
          real(8),intent(inout)::gr(:)
          real(8),allocatable::gt(:)
          integer::ic,jc,kc,ijc,nc,it,ngr,nc2
          nc2=SIZE(gr)
          call iup_(nc2,nc,ic)
          allocate(gt(SIZE(gr)))
          gt(:)=0d0
          ngr=0
          do it=ibr(jm),nbr
            if (lbr(1,it)>jm) then
              exit
            end if
            ! Frist term contribution.
            kc=lbr(2,it)
            jc=lbr(3,it)
            do ic=1,jc-1
              ijc=ip_(ic,jc)
              gt(ijc)=gt(ijc)+gr(ip_(ic,kc))*vbr(it)
            end do
            ! Second term contribution.
            ic=lbr(2,it)
            kc=lbr(3,it)
            do jc=ic+1,nc
              ijc=ip_(ic,jc)
              gt(ijc)=gt(ijc)-vbr(it)*gr(ip_(kc,jc))
            end do
          end do
          do ijc=1,nc2
            if (ABS(gt(ijc))>TOL_CC_) then
              call iup_(ijc,jc,ic)
            end if
          end do
          gr(:)=gt(:)
        end subroutine

        subroutine add1raising_(nfc,nac,lplds,nbf,lbf,vbf)
          ! Reising generator contributions to 1e constants (CAS).
          integer,intent(in)::nfc,nac,lplds(:,:,:)
          integer,intent(inout)::nbf,lbf(:,:)
          real(8),intent(inout)::vbf(:)
          integer::im,jm,ic,jc,nc,nc2,ijc
          integer::nbr
          integer,allocatable::lbr(:,:),ibr(:)
          real(8),allocatable::vbr(:),gr(:)
          ! Get memory.
          nc=SIZE(lplds(1,1,:))
          nc2=ip_(nc,nc)
          ! Load and mount the basic generator.
          call mountbt_(nac,nfc,nbf,lbf,vbf,nbr,lbr,ibr,vbr)
          ! Evaluate raising generators.
          allocate(gr(nc2))
          do im=1,nac-1
            ! Set the initial value.
            call cpbrnz_(im,nbr,lbr,vbr,gr)
            ! Evaluate the 1e constants and dump them on the memory.
            do jm=im+1,nac-1
              call evalgr_(jm,gr,nbr,ibr,lbr,vbr)
              ! Damp
              do ijc=1,nc2
                if (ABS(gr(ijc))>TOL_CC_) then
                call iup_(ijc,jc,ic)
                  call addnz1_(
     &              jc,ic,nfc+jm+1,nfc+im,gr(ijc),nbf,lbf,vbf)
                end if
              end do
            end do
          end do
        end subroutine

        subroutine mknz1_(nfc,nac,lplds)
          ! Evaluate 1e coupling constants (CAS) and store nonzeros.
          implicit none
          integer,intent(in)::nfc,nac,lplds(:,:,:)
          integer::nbf
          integer,allocatable::lbf(:,:)
          real(8),allocatable::vbf(:)
          ! Get memory.
          nbf=0
          allocate(lbf(4,lbfnz_))
          allocate(vbf(lbfnz_))
          ! Weight generators.
          call add1weight_(nfc,nac,lplds,nbf,lbf,vbf)
          ! Basic raising generators.
          call add1basic_(nfc,nac,lplds,nbf,lbf,vbf)
          ! Raising genrators.
          call add1raising_(nfc,nac,lplds,nbf,lbf,vbf)
          ! Set the module variables.
          n1_=nbf
          allocate(l1_(4,n1_))
          allocate(v1_(n1_))
          l1_(:,:n1_)=lbf(:,:n1_)
          v1_(:n1_)=vbf(:n1_)
        end subroutine

        subroutine loadcc1_(it,ic,jc,im,jm,a1)
          ! Load the set of indices and value of 1e constants.
          implicit none
          integer,intent(in)::it
          integer,intent(out)::ic,jc,im,jm
          real(8),intent(out)::a1
          ic=l1_(1,it)
          jc=l1_(2,it)
          im=l1_(3,it)
          jm=l1_(4,it)
          a1=v1_(it)
        end subroutine

        subroutine getlcmp(ncmp,lcmp)
          ! Get the list of index to be scanned.
          integer,parameter::MXLCMP=1000
          integer,intent(out),allocatable::ncmp(:),lcmp(:,:)
          integer::nc,it,ic,jc
          nc=MAXVAL(l1_(1,:))
          allocate(ncmp(nc))
          allocate(lcmp(MXLCMP,nc))
          ncmp(:)=0
          do it=1,n1_
            ic=l1_(1,it)
            ncmp(ic)=ncmp(ic)+1
            lcmp(ncmp(ic),ic)=it
            jc=l1_(2,it)
            if (ic/=jc) then
              ncmp(jc)=ncmp(jc)+1
              lcmp(ncmp(jc),jc)=it
            end if
          end do
        end subroutine

        subroutine cc2eijekl_(nbf,lbf,vbf)
          ! Evaluate (delta_kj E_il) part of the 2e coupling constant.
          implicit none
          integer,intent(inout)::nbf,lbf(:,:)
          real(8),intent(inout)::vbf(:)
          integer::it,jt,iw,jw,ic,jc,im,jm
          real(8)::a1
          integer::ipt
          integer,allocatable::lcmp(:,:),ncmp(:)
          logical::go
          ! Get the list of index to be scanned.
          call getlcmp(ncmp,lcmp)
          ! Evaluate 2e constants using nonzero values of 1e constants.
          do it=1,n1_
            call loadcc1_(it,ic,jc,im,jm,a1)
            do iw=1,ncmp(ic)
              jt=lcmp(iw,ic)
              call addifnz2_(it,jt,nbf,lbf,vbf)
            end do
            ipt=1
            do jw=1,ncmp(jc)
              jt=lcmp(jw,jc)
              go=.true.
              do iw=ipt,ncmp(ic)
                if (jt==lcmp(iw,ic)) then
                  go=.false.
                  ipt=iw
                  exit
                end if
              end do
              if (go) then
                call addifnz2_(it,jt,nbf,lbf,vbf)
              end if
            end do
          end do
        end subroutine

        subroutine addifnz2_(it,jt,nbf,lbf,vbf)
          ! Evaluate 2e constants and add the nonzero values.
          implicit none
          integer,intent(in)::it,jt
          integer,intent(inout)::nbf,lbf(:,:)
          real(8),intent(inout)::vbf(:)
          integer::ic,jc,kc,lc,im,jm,km,lm
          real(8)::a1,a2
          ! Load the 1e constants.
          call loadcc1_(it,ic,jc,im,jm,a1)
          call loadcc1_(jt,kc,lc,km,lm,a2)
          ! Evaluate 2e constants and add the nonzero values.
          if (ic/=kc.and.ic/=lc.and.jc/=kc.and.jc/=lc) then
            return
          end if
          if (ic==jc.and.jc==kc.and.kc==lc) then
            call addnz2r_(ic,lc,im,jm,km,lm,a1,a2,nbf,lbf,vbf)
          else
            if (ic==kc.and.jc>=lc) then
              call addnz2r_(jc,lc,jm,im,km,lm,a1,a2,nbf,lbf,vbf)
            end if
            if (ic==lc.and.jc>=kc) then
              call addnz2r_(jc,kc,jm,im,lm,km,a1,a2,nbf,lbf,vbf)
            end if
            if (jc==kc.and.ic>=lc.and.ic/=jc) then
              call addnz2r_(ic,lc,im,jm,km,lm,a1,a2,nbf,lbf,vbf)
            end if
            if (jc==lc.and.ic>=kc.and.lc/=kc) then
              call addnz2r_(ic,kc,im,jm,lm,km,a1,a2,nbf,lbf,vbf)
            end if
          end if
        end subroutine

        subroutine addnz2r_(ic,jc,im,jm,km,lm,a1,a2,n,l,v)
          ! Add the non-zero values of the 2e constants (RAS).
          implicit none
          integer,intent(in)::ic,jc,im,jm,km,lm
          real(8),intent(in)::a1,a2
          integer,intent(inout)::n,l(:,:)
          real(8),intent(inout)::v(:)
          integer::ir,jr
          ir=mcr_(ic)
          jr=mcr_(jc)
          if (ir>0.and.jr>0) then
            n=n+1
            l(:,n)=(/ir,jr,im,jm,km,lm/)
            v(n)=a1*a2/2d0
          end if
        end subroutine

        subroutine cc2dkjeil_(n2,lbf,vbf)
          ! Evaluate (delta_kj E_il) part of the 2e coupling constant.
          implicit none
          integer,intent(inout)::n2,lbf(:,:)
          real(8),intent(inout)::vbf(:)
          integer::istart,istop,ic,jc,im,jm,lm,it
          real(8)::a1
          istart=MINVAL(l1_(3,:))
          istop=MAXVAL(l1_(3,:))
          do it=1,n1_
            ic=l1_(1,it)
            jc=l1_(2,it)
            im=l1_(3,it)
            lm=l1_(4,it)
            a1=v1_(it)
            do jm=istart,istop
              call addnz2r_(ic,jc,im,jm,jm,lm,-1d0,a1,n2,lbf,vbf)
            end do
          end do
        end subroutine

        subroutine genbasica_(jc,im,lplds,ic,val)
          ! Basic raising generator for the case (im+1=jm) (Type A).
          integer,intent(in)::jc,im,lplds(:,:,:)
          integer,intent(out)::ic
          real(8),intent(out)::val
          integer::nac,it
          integer,allocatable::diff(:,:)
          real(8)::ai,bi,aip,bip,aim,bim
          nac=SIZE(lplds(1,:,1))
          allocate(diff(3,nac))
          ! Start evaluation.
          ai=lplds(1,nac+1-im,jc)
          bi=lplds(2,nac+1-im,jc)
          aip=lplds(1,nac+1-(im+1),jc)
          bip=lplds(2,nac+1-(im+1),jc)
          aim=0d0
          bim=0d0
          if (im>1) then
            aim=lplds(1,nac+1-(im-1),jc)
            bim=lplds(2,nac+1-(im-1),jc)
          end if
          ic=0
          val=0d0
          if (ai==aim.and.aim==aip-1.and.ABS(bi)>TOL_) then
            diff(:,:)=lplds(:,:,jc)
            diff(:,nac+1-im)=diff(:,nac+1-im)+(/1,-1,0/)
            val=((bi*(bi+1d0))/((bip+1d0)*(bim+1d0)))**0.5d0
          end if
          if (ABS(val)>TOL_) then
            do it=jc,1,-1
              if (SUM((diff(:,:)-lplds(:,:,it))**2)<TOL_) then
                ic=it
                exit
              end if
            end do
          end if
        end subroutine

        subroutine genbasicb_(jc,im,lplds,ic,val)
          ! Basic raising generator for the case (im+1=jm) (Type B).
          integer,intent(in)::jc,im,lplds(:,:,:)
          integer,intent(out)::ic
          real(8),intent(out)::val
          integer::nac,it
          integer,allocatable::diff(:,:)
          real(8)::bi,ci,bip,cip,bim,cim
          nac=SIZE(lplds(1,:,1))
          allocate(diff(3,nac))
          ! Start evaluation.
          bi=lplds(2,nac+1-im,jc)
          ci=lplds(3,nac+1-im,jc)
          bip=lplds(2,nac+1-(im+1),jc)
          cip=lplds(3,nac+1-(im+1),jc)
          bim=0d0
          cim=0d0
          if (im>1) then
            bim=lplds(2,nac+1-(im-1),jc)
            cim=lplds(3,nac+1-(im-1),jc)
          end if
          ic=0
          val=0d0
          if (ci==cip.and.cip==cim+1.and.ABS(ci)>TOL_) then
            diff(:,:)=lplds(:,:,jc)
            diff(:,nac+1-im)=diff(:,nac+1-im)+(/0,1,-1/)
            val=(((bi+1d0)*(bi+2d0))/((bip+1d0)*(bim+1d0)))**0.5d0
          end if
          if (ABS(val)>TOL_) then
            do it=jc,1,-1
              if (SUM((diff(:,:)-lplds(:,:,it))**2)<TOL_) then
                ic=it
                exit
              end if
            end do
          end if
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
