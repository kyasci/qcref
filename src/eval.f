      module eval
        ! Molecular integral evaluation (S, P, and D).
        !
        ! Reference:
        !
        !     H. Taketa et. al., J. Phys. Soc. Jpn. 21, 2313 (1966).
        !
        implicit none
        private
        real(8),parameter::PI=3.1415926535897931d0,TWOETOL=1d-8
        ! --------------------------------------------------------------
        public ::
     &      eval_overlap,
     &      eval_kinetic,
     &      eval_attract,
     &      eval_repulse,
     &      eval_velocity
        ! --------------------------------------------------------------

      contains

        subroutine eval_overlap(n1,n2,c1,c2,e1,e2,r1,r2,l1,l2,val)
          implicit none
          integer,intent(in)::n1,n2,l1(3),l2(3)
          real(8),intent(in)::c1(n1),c2(n2),e1(n1),e2(n2),r1(3),r2(3)
          real(8),intent(out)::val
          integer::i,j
          real(8)::vt
          val=0d0
          do i=1,n1
            do j=1,n2
              call overlap_pg_(e1(i),e2(j),r1,r2,l1,l2,vt)
              vt=vt*pgn_(l1,e1(i))*pgn_(l2,e2(j))
              vt=vt*c1(i)*c2(j)
              val=val+vt
            end do
          end do
        end subroutine

        subroutine overlap_pg_(e1i,e2i,r1,r2,l1,l2,val)
          implicit none
          real(8),intent(in)::e1i,e2i,r1(3),r2(3)
          real(8),intent(out)::val
          integer,intent(in)::l1(3),l2(3)
          integer::ix,nn,i
          real(8)::gm,rp(3),pa,pb,vt,w,vf
          gm=e1i+e2i
          rp(:)=(e1i*r1(:)+e2i*r2(:))/gm
          val=(PI/gm)**1.5
          val=val*EXP(-e1i*e2i*SUM((r1-r2)**2)/gm)
          do ix=1,3
            pa=rp(ix)-r1(ix)
            pb=rp(ix)-r2(ix)
            nn=FLOOR((l1(ix)+l2(ix))/2d0)
            vt=0d0
            do i=0,nn
              w=nfact2_(i)/((2d0*gm)**i)
              call fcoef_(2*i,l1(ix),l2(ix),pa,pb,vf)
              w=w*vf
              vt=vt+w
            end do
            val=val*vt
          end do
        end subroutine

        integer function nfact2_(i)
          implicit none
          integer,intent(in)::i
          integer,parameter::MX=5,NLS(MX)=(/1,1,3,15,105/)
          if (i+1>MX) then
            write(*,'("error: out of range: nfact2_()")')
            stop 1
          end if
          nfact2_=NLS(i+1)
        end function

        subroutine fcoef_(i,m,l,a,b,val)
          implicit none
          integer,intent(in)::i,m,l
          real(8),intent(in)::a,b
          real(8),intent(out)::val
          integer::j,k
          if (i>6) then
            write(*,'("error: out of range: fcoef_()")')
            stop 1
          end if
          val=0d0
          do j=0,i
            k=i-j
            if (m>=j.and.l>=k) then
              val=val+nbin_(m,j)*a**(m-j)*nbin_(l,k)*b**(l-k)
            end if
          end do
        end subroutine

        integer function nbin_(m,l)
          ! Binomial constant.
          implicit none
          integer,intent(in)::m,l
          integer,parameter::MX=7,NLS(MX,MX)=RESHAPE((/
     &      1,0,0,0,0,0,0,
     &      1,1,0,0,0,0,0,
     &      1,2,1,0,0,0,0,
     &      1,3,3,1,0,0,0,
     &      1,4,6,4,1,0,0,
     &      1,5,10,10,5,1,0,
     &      1,6,15,20,15,6,1
     &      /),(/MX,MX/))
          if (l+1>MX.or.m+1>MX) then
            write(*,'("error: out of range: nbin_()")')
            stop 1
          end if
          nbin_=NLS(l+1,m+1)
        end function

        real(8) function pgn_(l1,e1i) result(vn)
          ! Norm of a primitive Gaussian function.
          implicit none
          integer,intent(in)::l1(3)
          real(8),intent(in)::e1i
          integer::i
          vn=2d0**(2d0*SUM(l1)+1.5d0)
          vn=vn*e1i**(SUM(l1)+1.5d0)
          do i=1,3
            vn=vn/nfact2_(l1(i))
          end do
          vn=vn/PI**1.5
          vn=vn**0.5d0
        end function

        subroutine eval_kinetic(n1,n2,c1,c2,e1,e2,r1,r2,l1,l2,val)
          implicit none
          integer,intent(in)::n1,n2,l1(3),l2(3)
          real(8),intent(in)::c1(n1),c2(n2),e1(n1),e2(n2),r1(3),r2(3)
          real(8),intent(out)::val
          integer::i,j
          real(8)::vt
          val=0d0
          do i=1,n1
            do j=1,n2
              call kinetic_pg_(e1(i),e2(j),r1,r2,l1,l2,vt)
              vt=vt*pgn_(l1,e1(i))*pgn_(l2,e2(j))
              vt=vt*c1(i)*c2(j)
              val=val+vt
            end do
          end do
        end subroutine

        subroutine kinetic_pg_(e1i,e2i,r1,r2,l1,l2,val)
          implicit none
          real(8),intent(in)::e1i,e2i,r1(3),r2(3)
          real(8),intent(out)::val
          integer,intent(in)::l1(3),l2(3)
          integer::lps(3),ix
          real(8)::vt,st
          call overlap_pg_(e1i,e2i,r1,r2,l1,l2,st)
          val=e2i*(2*SUM(l2)+3)*st
          vt=0d0
          do ix=1,3
            lps(:)=0
            lps(ix)=lps(ix)+2
            call overlap_pg_(e1i,e2i,r1,r2,l1,l2+lps,st)
            vt=vt+st
          end do
          val=val-2d0*e2i**2d0*vt
          vt=0d0
          do ix=1,3
            if (l2(ix)>=2) then
              lps(:)=0
              lps(ix)=lps(ix)-2
              call overlap_pg_(e1i,e2i,r1,r2,l1,l2+lps,st)
              vt=vt+l2(ix)*(l2(ix)-1)*st
            end if
          end do
          val=val-0.5d0*vt
        end subroutine

        subroutine eval_attract(n1,n2,c1,c2,e1,e2,r1,r2,l1,l2,r3,val)
          implicit none
          integer,intent(in)::n1,n2,l1(3),l2(3)
          real(8),intent(in)::
     &      c1(n1),c2(n2),e1(n1),e2(n2),r1(3),r2(3),r3(3)
          real(8),intent(out)::val
          integer::i,j
          real(8)::vt
          val=0d0
          do i=1,n1
            do j=1,n2
              call attract_pg_(e1(i),e2(j),r1,r2,l1,l2,r3,vt)
              vt=vt*pgn_(l1,e1(i))*pgn_(l2,e2(j))
              vt=vt*c1(i)*c2(j)
              val=val+vt
            end do
          end do
        end subroutine

        subroutine attract_pg_(e1i,e2i,r1,r2,l1,l2,r3,val)
          implicit none
          real(8),intent(in)::e1i,e2i,r1(3),r2(3),r3(3)
          real(8),intent(out)::val
          integer,intent(in)::l1(3),l2(3)
          integer::ix,ii,jj,kk
          real(8) ::
     &      gm,rp(3),gls(5,3),pa(3),pb(3),pc(3),vt,vb,gpc2
          gm=e1i+e2i
          val=2d0*PI/gm*EXP(-e1i*e2i*SUM((r1-r2)**2)/gm)
          rp(:)=(e1i*r1(:)+e2i*r2(:))/gm
          pa(:)=rp(:)-r1(:)
          pb(:)=rp(:)-r2(:)
          pc(:)=rp(:)-r3(:)
          gpc2=gm*SUM(pc**2)
          do ix=1,3
            call glist_(l1(ix),l2(ix),pa(ix),pb(ix),pc(ix),gm,gls(:,ix))
          end do
          vt=0d0
          do ii=0,l1(1)+l2(1)
            do jj=0,l1(2)+l2(2)
              do kk=0,l1(3)+l2(3)
                call boys_(ii+jj+kk,gpc2,vb)
                vt=vt+gls(ii+1,1)*gls(jj+1,2)*gls(kk+1,3)*vb
              end do
            end do
          end do
          val=val*vt
        end subroutine

        subroutine boys_(n,x,val)
          ! Boys function interface.
          implicit none
          integer,parameter::NMAX=10
          integer,intent(in)::n
          real(8),intent(in)::x
          real(8),intent(out)::val
          if (x<1d-6) then
            val=1d0/(2d0*n+1)-x/(2d0*n+3d0)
          else
            val=fboys_(n,x)
          end if
        end subroutine

        recursive real(8) function fboys_(n,x) result(v)
          ! Reculsive implementation of Boys function.
          implicit none
          integer,intent(in)::n
          real(8),intent(in)::x
          real(8)::expmx
          if (n==0) then
            v=SQRT(PI/x)*ERF(SQRT(x))/2d0
          else
            ! Avoiding underflow:
            expmx=0d0
            if (x<1d3) then
              expmx=EXP(-x)
            end if
            v=((2*n-1)*fboys_(n-1,x)-expmx)/(2d0*x)
          end if
        end function

        subroutine glist_(l1w,l2w,paw,pbw,pcw,gm,glsw)
          implicit none
          integer,intent(in)::l1w,l2w
          real(8),intent(in)::paw,pbw,pcw,gm
          real(8),intent(out)::glsw(:)
          real(8)::f0,f1,f2,f3
          glsw(:)=0d0
          if (l1w+l2w==0) then
            glsw(1)=1d0
          else if (l1w+l2w==1) then
            call fcoef_(0,l1w,l2w,paw,pbw,f0)
            glsw(1)=f0
            glsw(2)=-pcw
          else if (l1w+l2w==2) then
            call fcoef_(0,l1w,l2w,paw,pbw,f0)
            call fcoef_(1,l1w,l2w,paw,pbw,f1)
            glsw(1)=f0+1/(2*gm)
            glsw(2)=-f1*pcw-1/(2*gm)
            glsw(3)=pcw**2
          else if (l1w+l2w==3) then
            call fcoef_(0,l1w,l2w,paw,pbw,f0)
            call fcoef_(1,l1w,l2w,paw,pbw,f1)
            call fcoef_(2,l1w,l2w,paw,pbw,f2)
            glsw(1)=f0+f2/(2*gm)
            glsw(2)=-f1*pcw-f2/(2*gm)-3*pcw/(2*gm)
            glsw(3)=f2*pcw**2+3*pcw/(2*gm)
            glsw(4)=-pcw**3
          else if (l1w+l2w==4) then
            call fcoef_(0,l1w,l2w,paw,pbw,f0)
            call fcoef_(1,l1w,l2w,paw,pbw,f1)
            call fcoef_(2,l1w,l2w,paw,pbw,f2)
            call fcoef_(3,l1w,l2w,paw,pbw,f3)
            glsw(1)=f0+f2/(2*gm)+3/(4*gm**2)
            glsw(2)=-f1*pcw-f2/(2*gm)-3*pcw*f3/gm-3/(2*gm**2)
            glsw(3)=f2*pcw**2+3*pcw*f3/(2*gm)+3*pcw**2/gm+3/(4*gm**2)
            glsw(4)=-f3*pcw**3-3*pcw**2/gm
            glsw(5)=pcw**4
          else
            write(*,'("error: out of range: glist_()")')
            stop 1
          end if
        end subroutine

        subroutine eval_repulse(
     &      n1,n2,n3,n4,c1,c2,c3,c4,e1,e2,e3,e4,
     &      r1,r2,r3,r4,l1,l2,l3,l4,val)
          implicit none
          integer,intent(in)::n1,n2,n3,n4,l1(3),l2(3),l3(3),l4(3)
          real(8),intent(in)::
     &      c1(n1),c2(n2),c3(n3),c4(n4),e1(n1),e2(n2),e3(n3),e4(n4),
     &      r1(3),r2(3),r3(3),r4(3)
          real(8),intent(out)::val
          real(8)::vt
          integer::i,j,k,l
          val=0d0
          do i=1,n1
            do j=1,n2
              do k=1,n3
                do l=1,n4
                  call repulse_pg_(
     &              e1(i),e2(j),e3(k),e4(l),r1,r2,r3,r4,l1,l2,l3,l4,vt)
                  vt=vt*pgn_(l1,e1(i))*pgn_(l2,e2(j))
                  vt=vt*pgn_(l3,e3(k))*pgn_(l4,e4(l))
                  vt=vt*c1(i)*c2(j)*c3(k)*c4(l)
                  val=val+vt
                end do
              end do
            end do
          end do
        end subroutine

        subroutine repulse_pg_(
     &      e1i,e2i,e3i,e4i,r1,r2,r3,r4,l1,l2,l3,l4,val)
          implicit none
          integer,intent(in)::l1(3),l2(3),l3(3),l4(3)
          real(8),intent(in) ::
     &      e1i,e2i,e3i,e4i,r1(3),r2(3),r3(3),r4(3)
          real(8),intent(out)::val
          integer::ix,ii,jj,kk
          real(8)::gm1,gm2,rp(3),rq(3),rab2,rcd2,rpq2,delta,
     &      cls(9,3),vt,vb
          gm1=e1i+e2i
          gm2=e3i+e4i
          rp(:)=(e1i*r1(:)+e2i*r2(:))/gm1
          rq(:)=(e3i*r3(:)+e4i*r4(:))/gm2
          rab2=SUM((r1-r2)**2)
          rcd2=SUM((r3-r4)**2)
          rpq2=SUM((rp-rq)**2)
          delta=(1d0/gm1+1d0/gm2)/4d0
          val=1d0
          val=val*2*PI**2.5/(gm1*gm2)/((gm1+gm2)**0.5)
          val=val*EXP(-(e1i*e2i*rab2/gm1)-(e3i*e4i*rcd2/gm2))
          if (ABS(val)>TWOETOL) then
            do ix=1,3
              call clist_(
     &          r1(ix),l1(ix),r2(ix),l2(ix),rp(ix),e1i+e2i,
     &          r3(ix),l3(ix),r4(ix),l4(ix),rq(ix),e3i+e4i,cls(:,ix))
            end do
            vt=0d0
            do ii=0,l1(1)+l2(1)+l3(1)+l4(1)
              do jj=0,l1(2)+l2(2)+l3(2)+l4(2)
                do kk=0,l1(3)+l2(3)+l3(3)+l4(3)
                  call boys_(ii+jj+kk,rpq2/(4*delta),vb)
                  vt=vt+cls(ii+1,1)*cls(jj+1,2)*cls(kk+1,3)*vb
                end do
              end do
            end do
            val=val*vt
          end if
        end subroutine

        subroutine clist_(
     &      raw,l1w,rbw,l2w,rpw,gm1,rcw,l3w,rdw,l4w,rqw,gm2,clsw)
          implicit none
          integer,intent(in)::l1w,l2w,l3w,l4w
          real(8),intent(in)::raw,rbw,rpw,gm1,rcw,rdw,rqw,gm2
          real(8),intent(out)::clsw(:)
          real(8)::pa,pb,qc,qd,delta,hls1(5),hls2(5),f3n,f3d
          integer::i,j,ii,u
          pa=rpw-raw
          pb=rpw-rbw
          qc=rqw-rcw
          qd=rqw-rdw
          delta=0.25d0*(1d0/gm1+1d0/gm2)
          clsw(:)=0d0
          call hlist_(l1w,l2w,pa,pb,gm1,hls1)
          call hlist_(l3w,l4w,qc,qd,gm2,hls2)
          do i=0,l1w+l2w
            do j=0,l3w+l4w
              do u=0,INT(FLOOR((i+j)/2d0))
                f3n=1d0
                f3n=f3n*nfact_(i+j)*(-1)**u
                f3n=f3n*(rqw-rpw)**(i+j-2*u)
                f3d=1d0
                f3d=f3d*nfact_(u)
                f3d=f3d*nfact_(i+j-2*u)
                f3d=f3d*delta**(i+j-u)
                ii=i+j-u+1
                if (ii>SIZE(clsw)) then
                  write(*,'("error: out of range: clist_()")')
                  stop 1
                end if
                clsw(ii)=clsw(ii)+hls1(i+1)*(-1)**j*hls2(j+1)*f3n/f3d
              end do
            end do
          end do
        end subroutine

        integer function nfact_(i)
          implicit none
          integer,intent(in)::i
          integer,parameter::
     &      MX=10,NLS(MX)=(/1,1,2,6,24,120,720,5040,40320,362880/)
          if (i+1>MX) then
            write(*,'("error: out of range: ncact()")')
            stop 1
          end if
          nfact_=NLS(i+1)
        end function

        subroutine hlist_(l1w0,l2w0,pa0,pb0,gm,hls)
          implicit none
          integer,intent(in)::l1w0,l2w0
          real(8),intent(in)::pa0,pb0,gm
          integer::l1w,l2w
          real(8)::pa,pb
          real(8),intent(out)::hls(:)
          if (l1w0>=l2w0) then
            l1w=l1w0
            l2w=l2w0
            pa=pa0
            pb=pb0
          else
            l1w=l2w0
            l2w=l1w0
            pa=pb0
            pb=pa0
          end if
          hls(:)=0d0
          if (l1w==0.and.l2w==0) then
            hls(1)=1d0
          else if (l1w==1.and.l2w==0) then
            hls(1)=pa
            hls(2)=1/(4*gm)
          else if (l1w==1.and.l2w==1) then
            hls(1)=pa*pb+1d0/(2*gm)
            hls(2)=(pa+pb)/(4*gm)
            hls(3)=1/(4*gm)**2
          else if (l1w==2.and.l2w==0) then
            hls(1)=pa**2+1d0/(2*gm)
            hls(2)=pa/(2*gm)
            hls(3)=1/(4*gm)**2
          else if (l1w==2.and.l2w==1) then
            hls(1)=pa**2*pb+(2*pa+pb)/(2*gm)
            hls(2)=(2*pa*pb+pa**2)/(4*gm)+6*(1/(4*gm)**2)
            hls(3)=(2*pa+pb)*(1/(4*gm))**2
            hls(4)=(1/(4*gm))**3
          else if (l1w==2.and.l2w==2) then
            hls(1)=(pa*pb)**2+(pa**2+4*pa*pb+pb**2)/(2*gm)
     &        +12*(1/(4*gm)**2)
            hls(2)=pa*pb*(pa+pb)/(2*gm)+12*(pa+pb)*(1/(4*gm))**2
            hls(3)=(pa**2+4*pa*pb+pb**2)*(1/(4*gm))**2+12*(1/(4*gm))**3
            hls(4)=2*(pa+pb)*(1/(4*gm))**3
            hls(5)=(1/(4*gm))**4
          else
            write(*,'("error: out of range: hlist_()")')
            stop 1
          end if
        end subroutine

        subroutine eval_velocity(n1,n2,c1,c2,e1,e2,r1,r2,l1,l2,ix,val)
          implicit none
          integer,intent(in)::n1,n2,l1(3),l2(3),ix
          real(8),intent(in)::c1(n1),c2(n2),e1(n1),e2(n2),r1(3),r2(3)
          real(8),intent(out)::val
          integer::i,j
          real(8)::vt
          val=0d0
          do i=1,n1
            do j=1,n2
              call velocity_pg_(e1(i),e2(j),r1,r2,l1,l2,ix,vt)
              vt=vt*pgn_(l1,e1(i))*pgn_(l2,e2(j))
              vt=vt*c1(i)*c2(j)
              val=val+vt
            end do
          end do
        end subroutine

        subroutine velocity_pg_(e1i,e2i,r1,r2,l1,l2,ix,val)
          implicit none
          integer,intent(in)::l1(3),l2(3),ix
          real(8),intent(in)::e1i,e2i,r1(3),r2(3)
          real(8),intent(out)::val
          integer::lps(3)
          real(8)::vt
          lps(:)=0
          lps(ix)=lps(ix)+1
          call overlap_pg_(e1i,e2i,r1,r2,l1,l2+lps,vt)
          val=-2d0*e2i*vt
          if (l2(ix)>0) then
            call overlap_pg_(e1i,e2i,r1,r2,l1,l2-lps,vt)
            val=val+l2(ix)*vt
          end if
        end subroutine

      end module
