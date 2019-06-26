      module basis
        ! Atomic basis set and integral matrices as a wrapper of "eval".
        use eval
        implicit none
        private
        integer,parameter::MXB=1000,MXC=20
        integer::nb_,ls_(3,MXB),is_(MXB),ns_(MXB)
        real(8)::cs_(MXC,MXB),es_(MXC,MXB)
        ! --------------------------------------------------------------
        public ::
     &    basis_init,
     &    basis_add,
     &    basis_set,
     &    basis_dim,
     &    basis_sij,
     &    basis_hij,
     &    basis_dijk,
     &    basis_ijkl,
     &    basis_enuc
        ! --------------------------------------------------------------

      contains

        subroutine basis_init()
          implicit none
          nb_=0
          ls_(:,:)=0
          is_(:)=0
          ns_(:)=0
          cs_(:,:)=0d0
          es_(:,:)=0d0
        end subroutine

        subroutine basis_dim(num)
          ! The number of basis set.
          implicit none
          integer,intent(out)::num
          num=nb_
        end subroutine

        subroutine basis_set(key,ia,i1)
          ! Add an atomic basis set.
          implicit none
          character(*),intent(in)::key
          integer,intent(in)::ia,i1
          if (key=='STO3G') then
            call bs_sto3g_(ia,i1)
          else if (key=='321G') then
            call bs_321g_(ia,i1)
          else
            write(*,'("error: no such basis set: ",a)') TRIM(key)
            stop 1
          end if
        end subroutine

        subroutine basis_add(i1,n1,l1,c1,e1)
          ! Add Gaussian basis set to the list.
          implicit none
          integer,intent(in)::i1,n1,l1(3)
          real(8),intent(in)::c1(n1),e1(n1)
          nb_=nb_+1
          if (nb_>MXB) then
            write(*,'("error: out of range: basis_add()")')
            stop 1
          end if
          ns_(nb_)=n1
          is_(nb_)=i1
          ls_(:3,nb_)=l1(:3)
          cs_(:n1,nb_)=c1(:n1)
          es_(:n1,nb_)=e1(:n1)
        end subroutine

        subroutine basis_enuc(r,ian,val)
          ! Nuclear repulsion energy.
          integer,intent(in)::ian(:)
          real(8),intent(in)::r(:,:)
          real(8),intent(out)::val
          real(8)::rij
          integer::i,j
          val=0d0
          do i=1,SIZE(ian)
            do j=1,i-1
              rij=SQRT(SUM((r(:,i)-r(:,j))**2))
              val=val+ian(i)*ian(j)/rij
            end do
          end do
        end subroutine

        subroutine basis_sij(r,val)
          ! AO overlap matrix (i|j).
          implicit none
          real(8),intent(in)::r(:,:)
          real(8),intent(out),allocatable::val(:)
          integer::ij,i,j,nb2
          real(8)::vt
          ! Get memory.
          nb2=ip_(nb_,nb_)
          allocate(val(nb2))
          ! Evaluate all the symmetry-unique matrix elements.
          do ij=1,SIZE(val)
            call iup_(ij,i,j)
            if (i==j) then
              val(ij)=1d0
            else
              call eval_overlap(
     &          ns_(i),ns_(j),cs_(:,i),cs_(:,j),es_(:,i),es_(:,j),
     &          r(:,is_(i)),r(:,is_(j)),ls_(:,i),ls_(:,j),vt)
              val(ij)=vt
            end if
          end do
        end subroutine

        subroutine basis_hij(r,ian,val)
          ! AO core Hamiltonian matrix (i|h|j).
          implicit none
          integer,intent(in)::ian(:)
          real(8),intent(in)::r(:,:)
          real(8),intent(out),allocatable::val(:)
          integer::ij,i,j,ic,nb2
          real(8)::vt
          ! Get memory.
          nb2=ip_(nb_,nb_)
          allocate(val(nb2))
          ! Evaluate all the symmetry-unique matrix elements.
          do ij=1,SIZE(val)
            call iup_(ij,i,j)
            call eval_kinetic(
     &        ns_(i),ns_(j),cs_(:,i),cs_(:,j),es_(:,i),es_(:,j),
     &        r(:,is_(i)),r(:,is_(j)),ls_(:,i),ls_(:,j),vt)
            val(ij)=vt
            do ic=1,SIZE(ian)
              call eval_attract(
     &          ns_(i),ns_(j),cs_(:,i),cs_(:,j),es_(:,i),es_(:,j),
     &          r(:,is_(i)),r(:,is_(j)),ls_(:,i),ls_(:,j),r(:,ic),vt)
              val(ij)=val(ij)-ian(ic)*vt
            end do
          end do
        end subroutine

        subroutine basis_dijk(r,val)
          ! AO derivative coupling matrix (i|d/dRk|j)
          implicit none
          real(8),intent(in)::r(:,:)
          real(8),intent(out)::val(:,:,:,:)
          real(8),allocatable::vb(:,:,:)
          integer::nb,i,j,ix,ia
          real(8)::vt
          nb=SIZE(val(:,1,1,1))
          allocate(vb(nb,nb,3))
          ! AO velocity integral (i|d/dr|j)
          do ix=1,3
            do i=1,nb
              do j=1,nb
                if (i==j) then
                  vt=0d0
                else
                  call eval_velocity(
     &              ns_(i),ns_(j),cs_(:,i),cs_(:,j),es_(:,i),es_(:,j),
     &              r(:,is_(i)),r(:,is_(j)),ls_(:,i),ls_(:,j),ix,vt)
                end if
                vb(i,j,ix)=vt
              end do
            end do
          end do
          val(:,:,:,:)=0d0
          do ia=1,SIZE(r(1,:))
            do ix=1,3
              do i=1,nb
                do j=1,nb
                  if (is_(j)==ia) then
                    val(i,j,ix,ia)=(-1d0)*vb(i,j,ix)
                  end if
                end do
              end do
            end do
          end do
        end subroutine

        subroutine basis_ijkl(r,val)
          ! AO 2e integral matrix (ij|kl).
          implicit none
          real(8),intent(in)::r(:,:)
          real(8),intent(out),allocatable::val(:)
          integer::ijkl,i,j,k,l,nb3
          real(8)::vt
          ! Get memory.
          nb3=ip2_(nb_,nb_,nb_,nb_)
          allocate(val(nb3))
          ! Evaluate all the symmetry-unique matrix elements.
!$omp     parallel private(ijkl,i,j,k,l,vt)
!$omp     do
          do ijkl=1,SIZE(val)
            call iup2_(ijkl,i,j,k,l)
            call eval_repulse(
     &        ns_(i),ns_(j),ns_(k),ns_(l),
     &        cs_(:,i),cs_(:,j),cs_(:,k),cs_(:,l),
     &        es_(:,i),es_(:,j),es_(:,k),es_(:,l),
     &        r(:,is_(i)),r(:,is_(j)),r(:,is_(k)),r(:,is_(l)),
     &        ls_(:,i),ls_(:,j),ls_(:,k),ls_(:,l),vt)
            val(ijkl)=vt
          end do
!$omp     end do
!$omp     end parallel
        end subroutine

        subroutine add_s_(i1,cf,ex)
          implicit none
          integer,intent(in)::i1
          real(8),intent(in)::cf(:),ex(:)
          integer::l1(3)
          l1(:)=0
          call basis_add(i1,SIZE(ex),l1,cf,ex)
        end subroutine

        subroutine add_l_(i1,cfs,cfp,ex)
          implicit none
          integer,intent(in)::i1
          real(8),intent(in)::cfs(:),cfp(:),ex(:)
          integer::l1(3),i
          l1(:)=0
          call basis_add(i1,SIZE(ex),l1,cfs,ex)
          do i=1,3
            l1(:)=0
            l1(i)=1
            call basis_add(i1,SIZE(ex),l1,cfp,ex)
          end do
        end subroutine

        subroutine bs_sto3g_(ia,i1)
          ! STO-3G basis set (H-Ne).
          implicit none
          integer,intent(in)::ia,i1
          integer,parameter::MXNA=10
          integer,parameter::NG=3
          real(8),parameter::
     &      zt1(MXNA)=(/
     &        1.24d0,1.69d0,2.69d0,3.68d0,4.68d0,
     &        5.67d0,6.67d0,7.66d0,8.65d0,9.64d0
     &        /),
     &      zt2(MXNA)=(/
     &        0.00d0,0.00d0,0.80d0,1.15d0,1.45d0,
     &        1.72d0,1.95d0,2.25d0,2.55d0,2.88d0
     &        /),
     &      ex1s(NG)=(/
     &        2.22766d0,0.405771d0,0.109818d0
     &        /),
     &      cf1s(NG)=(/
     &        0.154328967295d0,0.535328142282d0,0.444634542185d0
     &        /),
     &      ex2sp(NG)=(/
     &        0.994203d0,0.231031d0,0.0751386d0
     &        /),
     &      cf2s(NG)=(/
     &        -0.099967229187d0,0.399512826089d0,0.700115468880d0
     &        /),
     &      cf2p(NG)=(/
     &        0.155916274999d0,0.607683718598d0,0.391957393099d0
     &        /)
          real(8)::ex(NG)
          ! H - He.
          if (1<=ia.and.ia<=2) then
            ex(:)=ex1s(:)*zt1(ia)**2
            call add_s_(i1,cf1s,ex)
          ! Li - Ne.
          else if (3<=ia.and.ia<=10) then
            ex(:)=ex1s(:)*zt1(ia)**2
            call add_s_(i1,cf1s,ex)
            ex(:)=ex2sp(:)*zt2(ia)**2
            call add_l_(i1,cf2s,cf2p,ex)
          else
            write(*,'("error: out of range: basis.bs_sto3g_()")')
            stop 1
          end if
        end subroutine

        subroutine bs_321g_(ia,i1)
          ! 3-21G basis set (H-Ne).
          implicit none
          integer,intent(in)::ia,i1
          ! H:
          if (ia==1) then
            call add_s_(i1,
     &        (/0.156284978695d0,0.904690876670d0/),
     &        (/5.4471780d0,0.8245472d0/)
     &        )
            call add_s_(i1,
     &        (/1.0d0/),(/0.1831916d0/) 
     &        )
          ! He:
          else if (ia==2) then
            call add_s_(i1,
     &        (/0.175229871839d0,0.893482346517d0/),
     &        (/13.6267000d0,1.9993500d0/)
     &        )
            call add_s_(i1,
     &        (/1.0d0/),(/0.3829930d0/)
     &        )
          ! Li:
          else if (ia==3) then
            call add_s_(i1,
     &        (/0.069668663812d0,0.381346349289d0,0.681702624397d0/),
     &        (/36.8382000d0,5.4817200d0,1.1132700d0/) 
     &        )
            call add_l_(i1,
     &        (/-0.263126405761d0,1.143387417797d0/),
     &        (/0.161545970837d0,0.915662834700d0/),
     &        (/0.5402050d0,0.1022550d0/)
     &        )
            call add_l_(i1,
     &        (/1.0d0/),(/1.0d0/),(/0.0285645d0/)
     &        )
          ! Be:
          else if (ia==4) then
            call add_s_(i1,
     &        (/0.064426309746d0,0.366096055378d0,0.695934105272d0/),
     &        (/71.8876000d0,10.7289000d0,2.2220500d0/)
     &        )
            call add_l_(i1,
     &        (/-0.421064065879d0,1.224070191516d0/),
     &        (/0.205131923736d0,0.882527671894d0/),
     &        (/1.2954800d0,0.2688810d0/) 
     &        )
            call add_l_(i1,
     &        (/1.0d0/),(/1.0d0/),(/0.0773501d0/)
     &        )
          ! B:
          else if (ia==5) then
            call add_s_(i1,
     &        (/0.062960465890d0,0.363303803175d0,0.697254622252d0/),
     &        (/116.4340000d0,17.4314000d0,3.6801600d0/)
     &        )
            call add_l_(i1,
     &        (/-0.368663477319d0,1.199444806451d0/),
     &        (/0.231151902315d0,0.866763633704d0/),
     &        (/2.2818700d0,0.4652480d0/)
     &        )
            call add_l_(i1,
     &         (/1.0d0/),(/1.0d0/),(/0.1243280d0/)
     &        )
          ! C:
          else if (ia==6) then
            call add_s_(i1,
     &        (/0.061766907377d0,0.358794042852d0,0.700713083689d0/),
     &        (/172.2560000d0,25.9109000d0,5.5333500d0/)
     &        )
            call add_l_(i1,
     &        (/-0.395895162119d0,1.215834355681d0/),
     &        (/0.236459946619d0,0.860618805716d0/),
     &        (/3.6649800d0,0.7705450d0/)
     &        )
            call add_l_(i1,
     &         (/1.0d0/),(/1.0d0/),(/0.1958570d0/)
     &        )
          ! N:
          else if (ia==7) then
            call add_s_(i1,
     &        (/0.059865700508d0,0.352955002994d0,0.706513005993d0/),
     &        (/242.7660000d0,36.4851000d0,7.8144900d0/)
     &        )
            call add_l_(i1,
     &        (/-0.413300077430d0,1.224417266851d0/),
     &        (/0.237972016222d0,0.858953058551d0/),
     &        (/5.4252200d0,1.1491500d0/)
     &        )
            call add_l_(i1,
     &         (/1.0d0/),(/1.0d0/),(/0.2832050d0/)
     &        )
          ! O:
          else if (ia==8) then
            call add_s_(i1,
     &        (/0.059239393389d0,0.351499960776d0,0.707657921031d0/),
     &        (/322.0370000d0,48.4308000d0,10.4206000d0/)
     &        )
            call add_l_(i1,
     &        (/-0.404453583190d0,1.221561761397d0/),
     &        (/0.244586106967d0,0.853955373466d0/),
     &        (/7.4029400d0,1.5762000d0/)
     &        )
            call add_l_(i1,
     &         (/1.0d0/),(/1.0d0/),(/0.3736840d0/)
     &        )
          ! F:
          else if (ia==9) then
            call add_s_(i1,
     &        (/0.058548302929d0,0.349308017477d0,0.709632035505d0/),
     &        (/413.8010000d0,62.2446000d0,13.4340000d0/)
     &        )
            call add_l_(i1,
     &        (/-0.407326277682d0,1.223137830990d0/),
     &        (/0.246680003198d0,0.852321011049d0/),
     &        (/9.7775900d0,2.0861700d0/)
     &        )
            call add_l_(i1,
     &         (/1.0d0/),(/1.0d0/),(/0.4823830d0/)
     &        )
          ! Ne:
          else if (ia==10) then
            call add_s_(i1,
     &        (/0.058143030441d0,0.347951182169d0,0.710714372092d0/),
     &        (/515.7240000d0,77.6538000d0,16.8136000d0/)
     &        )
            call add_l_(i1,
     &        (/-0.409922320837d0,1.224310958240d0/),
     &        (/0.247459983573d0,0.851742943461d0/),
     &        (/12.4830000d0,2.6645100d0/)
     &        )
            call add_l_(i1,
     &         (/1.0d0/),(/1.0d0/),(/0.6062500d0/)
     &        )
          else
            write(*,'("error: out of range: basis.bs_sto3g_()")')
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

        integer function ip2_(i,j,k,l) result(ijkl)
          ! Canonical index packing for a symmetric matrix.
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
