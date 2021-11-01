      module drt
        ! Distinct Row Table (DRT) and its utilities.
        !
        ! Reference:
        !
        !   [1] https://www.univie.ac.at/columbus/workshops/argonne2005/pdf/Shavitt_Argonne_GUGA.pdf
        !   [2] J. Paldus, J. Chem. Phys. 61, 5321 (1974).
        !
        implicit none
        private
        ! Private variables.
        integer::NOVAL_=-1
        integer,allocatable::ly_(:,:),lk_(:,:),lx_(:),labc_(:,:)
        ! Public methods.
        public::
     &    drt_init,
     &    drt_del,
     &    drt_getld,
     &    drt_getic,
     &    drt_getwgt,
     &    drt_getbra,
     &    drt_getbrb,
     &    drt_showtab

      contains

        subroutine mkdrt_(mabc)
          ! Make the distinct row table (DRT).
          implicit none
          integer,intent(in)::mabc(3)
          integer,parameter::MXNJ=3000
          integer::labct(3,MXNJ),nj
          ! Generate the list of (a,b,c).
          call mkabc_(mabc,nj,labct)
          ! Get memorey.
          allocate(labc_(3,nj))
          allocate(ly_(4,nj))
          allocate(lk_(4,nj))
          allocate(lx_(nj))
          ly_(:,:)=0
          lx_(:)=0
          labc_(:,:nj)=labct(:,:nj)
          call mkkd_(labc_,lk_)
          call mkydx_(labc_,lk_,ly_,lx_)
        end subroutine

        subroutine mkkd_(labc,lkd)
          ! Make link k_dj.
          implicit none
          integer,intent(in)::labc(:,:)
          integer,intent(out)::lkd(:,:)
          integer::j,k,nj,mabcj(3),mabck(3),lsft(3,4),nsft,isft
          nj=size(lkd(1,:))
          lkd(:,:)=NOVAL_
          do j=1,nj-1
            mabcj(:)=labc(:,j)
            call table1_(mabcj,nsft,lsft)
            do isft=1,nsft
              mabck(:)=mabcj(:)-lsft(:,isft)
              do k=j,nj
                if (sum(abs(mabck(:)-labc(:,k)))==0) then
                  lkd(sid_(lsft(:,isft)),j)=k
                end if
              end do
            end do
          end do
        end subroutine

        subroutine mkydx_(labc,lkd,lyd,lx)
          ! Make the list of y_dj and x_j.
          implicit none
          integer,intent(in)::labc(:,:),lkd(:,:)
          integer,intent(out)::lyd(:,:),lx(:)
          integer::j,nj,md,kdjp,mydjp
          logical::is1st
          nj=size(labc(1,:))
          lyd(:,:)=NOVAL_
          lx(:)=0
          lx(nj)=1
          do j=nj-1,1,-1
            is1st=.True.
            do md=1,4
              if (lkd(md,j)==NOVAL_) then
                cycle
              else if (is1st) then
                lyd(md,j)=0
                mydjp=0
                kdjp=lkd(md,j)
                is1st=.False.
              else
                lyd(md,j)=mydjp+lx(kdjp)
                mydjp=lyd(md,j)
                kdjp=lkd(md,j)
              end if
            end do
            lx(j)=mydjp + lx(kdjp)
          end do
        end subroutine

        integer function sid_(msft) result(num)
          integer,intent(IN)::msft(3)
          if (msft(1)==0.and.msft(2)==0.and.msft(3)==1) then
            num=1
          else if (msft(1)==0.and.msft(2)==1.and.msft(3)==0) then
            num=2
          else if (msft(1)==1.and.msft(2)==-1.and.msft(3)==1) then
            num=3
          else if (msft(1)==1.and.msft(2)==0.and.msft(3)==0) then
            num=4
          else
            print *, 'error in drt.sftid()'
            stop 1
          end if
        end function

        subroutine mkabc_(mabc,nj,labct)
          ! Make (a,b,c) rows of DRT.
          implicit none
          integer,intent(in)::mabc(3)
          integer,intent(out)::nj,labct(:,:)
          integer::
     &      mabcj(3),mabck(3),i,j,k,nlevel,jmax,lsft(3,4),nsft,isft
          nlevel=sum(mabc)
          nj=1
          labct(:,:)=0
          labct(:,1)=mabc(:)
          jmax=size(labct(1,:))
          do i=nlevel,2,-1
            do j=1,jmax
              if (sum(labct(:,j))==i) then
                mabcj(:)=labct(:,j)
                call table1_(mabcj,nsft,lsft)
                do isft=1,nsft
                  mabck(:)=mabcj(:)-lsft(:,isft)
                  do k=1,jmax
                    if (sum(abs(mabck(:)-labct(:,k)))==0) then
                      exit
                    else if (sum(labct(:,k))==0) then
                      nj=nj+1
                      labct(:,nj)=mabck(:)
                      exit
                    end if
                  end do
                end do
              end if
            end do
          end do
          nj=nj+1
          labct(:,nj)=(/0,0,0/)
        end subroutine

        subroutine table1_(mabc,n,num)
          ! Get the subtration table for the (i-1)'s row.
          implicit none
          integer,intent(in)::mabc(:)
          integer,intent(out)::n,num(:,:)
          integer::ma,mb,mc
          ma=mabc(1)
          mb=mabc(2)
          mc=mabc(3)
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
            write(*,'("error: unclassified: drt.table1_()")')
            stop 1
          end if
        end subroutine

        subroutine drt_init(ma,mb,mc)
          ! DRT constractor.
          implicit none
          integer,intent(in)::ma,mb,mc
          integer::mabc(3)
          mabc(:)=(/ma,mb,mc/)
          call mkdrt_(mabc)
        end subroutine

        subroutine drt_del()
          ! DRT destoractor.
          implicit none
          deallocate(ly_)
          deallocate(lk_)
          deallocate(lx_)
          deallocate(labc_)
        end subroutine

        subroutine drt_getld(ic,ld)
          ! Get d's to define the walk corresponding to the given index.
          implicit none
          integer,intent(in)::ic
          integer,intent(out)::ld(:)
          integer::i,j,nlv,md,mdt,mw,mwt
          if (ic>lx_(1)) then
            print '("error: ic > N(CSF)")'
            stop 1
          end if
          nlv=sum(labc_(:,1))
          mw = 0
          j=1
          ld(:)=NOVAL_
          do i=1,nlv
            do md=1,4
              if (ly_(md,j)/=NOVAL_) then
                if (mw+ly_(md,j)>=ic) then
                  exit
                end if
                mwt=ly_(md,j)
                mdt=md
              end if
            end do
            mw=mw+mwt
            j=lk_(mdt,j)
            ld(i)=mdt
          end do
        end subroutine

        subroutine drt_getic(ld,ic)
          ! Get the index corresponding to the walk defined by d's.
          implicit none
          integer,intent(in)::ld(:)
          integer,intent(out)::ic
          integer::j,md
          ic=1
          j=1
          do md=1,size(ld)-1
            ic=ic+ly_(ld(md),j)
            j=lk_(ld(md),j)
          end do
        end subroutine

        subroutine drt_getwgt(im,ic,val)
          ! Matrix element of the weight operator.
          implicit none
          integer,intent(in)::im,ic
          real(8),intent(out)::val
          integer::nlv,labc(3,3)
          integer,allocatable::ld(:)
          nlv=sum(labc_(:,1))
          allocate(ld(nlv))
          call get_abcsub_(im,ic,labc,ld)
          val=0d0
          val=val+2*(labc(1,2)-labc(1,3))
          val=val+(labc(2,2)-labc(2,3))
        end subroutine

        subroutine drt_getbra(im,ic,jc,val)
          ! Matrix element and index of the basic raising (Type A).
          implicit none
          integer,intent(in)::im,ic
          integer,intent(out)::jc
          real(8),intent(out)::val
          integer::nlv,labc(3,3)
          integer,allocatable::ld(:)
          nlv=sum(labc_(:,1))
          allocate(ld(nlv))
          call get_abcsub_(im,ic,labc,ld)
          jc=0
          val=0.d0
          if (
     &      labc(1,2)==labc(1,3).and.labc(1,2)==labc(1,1)-1
     &      .and.labc(2,2)/=0) then
            ! Calculte value.
            val=labc(2,2)*(labc(2,2)+1)
            val=val/((labc(2,1)+1)*(labc(2,3)+1))
            val=val**0.5
            ! Calculate index.
            if (ld(nlv-im)==4.and.ld(nlv-im+1)==1) then
              ld(nlv-im)=2
              ld(nlv-im+1)=3
            else if (ld(nlv-im)==4.and.ld(nlv-im+1)==2) then
              ld(nlv-im)=2
              ld(nlv-im+1)=4
            else if (ld(nlv-im)==3.and.ld(nlv-im+1)==2) then
              ld(nlv-im)=1
              ld(nlv-im+1)=4
            else
              print '("error: drt_getbra(): unexpected ld type")'
              stop 1
            end if
            call drt_getic(ld,jc)
          end if
        end subroutine

        subroutine drt_getbrb(im,ic,jc,val)
          ! Matrix element and index of the basic raising (Type B).
          implicit none
          integer,intent(in)::im,ic
          integer,intent(out)::jc
          real(8),intent(out)::val
          integer::nlv,labc(3,3)
          integer,allocatable::ld(:)
          nlv=sum(labc_(:,1))
          allocate(ld(nlv))
          call get_abcsub_(im,ic,labc,ld)
          jc=0
          val=0.d0
          if (
     &      labc(3,2)==labc(3,1).and.labc(3,2)==labc(3,3)+1
     &      .and.labc(3,2)/=0) then
            ! Calculte value.
            val=(labc(2,2)+1)*(labc(2,2)+2)
            val=val/((labc(2,1)+1)*(labc(2,3)+1))
            val=val**0.5
            ! Calculate index.
            if (ld(nlv-im)==4.and.ld(nlv-im+1)==1) then
              ld(nlv-im)=3
              ld(nlv-im+1)=2
            else if (ld(nlv-im)==4.and.ld(nlv-im+1)==3) then
              ld(nlv-im)=3
              ld(nlv-im+1)=4
            else if (ld(nlv-im)==2.and.ld(nlv-im+1)==1) then
              ld(nlv-im)=1
              ld(nlv-im+1)=2
            else if (ld(nlv-im)==2.and.ld(nlv-im+1)==3) then
              ld(nlv-im)=3
              ld(nlv-im+1)=2
            else
              print '("error: drt_getbrb(): unexpected ld type")'
              stop 1
            end if
            call drt_getic(ld,jc)
          end if
        end subroutine

        subroutine get_abcsub_(im,ic,labc,ld)
          ! Get Paldus subtable of (im-1,im,im+1) levels.
          implicit none
          integer,intent(in)::im,ic
          integer,intent(out)::labc(:,:),ld(:)
          integer::jm,it,ilv,nlv
          nlv=sum(labc_(:,1))
          if (im+1>nlv) then
            print '("error: get_abcsub_(): im+1>nlv")'
            stop 1
          end if
          call drt_getld(ic,ld)
          jm=1
          do it=1,nlv+1
            ilv=sum(labc_(:,jm))
            if (ilv-1==im) then
              labc(:,1)=labc_(:,jm)
            else if (ilv==im) then
              labc(:,2)=labc_(:,jm)
            else if (ilv+1==im) then
              labc(:,3)=labc_(:,jm)
            end if
            if (it>=nlv+1) then
              exit
            end if
            jm=lk_(ld(it),jm)
          end do
        end subroutine

        subroutine drt_showtab()
          ! Show the DRT.
          implicit none
          integer::j,nj
          character(len=256)::s
          nj=size(lx_)
          s="  i|   j|  a  b  c | k_0 k_1 k_2 k_3|"
          s=trim(s)//"       y_0       y_1       y_2       y_3|"
          s=trim(s)//"         x"
          print '(a)', trim(s)
          do j=1,nj
            print '(i3,"|",i4,"|",3i3," |",4i4,"|",4i10,"|",i10)',
     &        sum(labc_(:,j)),j,labc_(:,j),lk_(:,j),ly_(:,j),lx_(j)
          end do
        end subroutine

      end module
!
!      program main
!        use drt
!        implicit none
!        integer::ma,mb,mc,ic,im,jc
!        integer,allocatable::ld(:)
!        real(8)::val
!        ma=2
!        mb=1
!        mc=3
!        ic=1
!        im=1
!        allocate(ld(ma+mb+mc))
!        call drt_init(ma,mb,mc)
!        call drt_showtab()
!        call drt_getld(ic,ld)
!        call drt_getic(ld,ic)
!        call drt_getwgt(im,ic,val)
!        print *, "w",ic,im,val
!        call drt_getbra(im,ic,jc,val)
!        if (jc/=0) print *, "a",jc,ic,im,im+1,val
!        call drt_getbrb(im,ic,jc,val)
!        if (jc/=0) print *, "b",jc,ic,im,im+1,val
!        call drt_del()
!      end program
