! This is a new file to contain all subroutines found in Dr. Handler's
! polymer code without cluttering the pre-existing DNS code from
! Goldstein's group. Make sure this code is added to the compile list,
! and everything should work fine. - Ryan 3-1-22

subroutine coriolis1(u,v,w,gn,fn,omz)
!************************************************************************
      use grid_size
      complex u(nyp,nz,nxh),v(nyp,nz,nxh),w(nyp,nz,nxh)
      complex gn(nyp,nz,nxh),fn(nyp,nz,nxh),omz(nyp,nz,nxh)
      common/thermcor/ omearth,rotang,grav,virtual

          twopi = 8.*atan(1.)

          ang =  rotang*(twopi/360.)
          csang = cos(ang)
          siang = sin(ang) 

          do k=1,nxh
              do j=1,nz
                  do i=1,nyp
                       gn(i,j,k) = gn(i,j,k) - 2.0*omearth*(siang*w(i,j,k) + csang*v(i,j,k)) 
                       fn(i,j,k) = fn(i,j,k) + 2.0*omearth*csang*u(i,j,k)    
                       omz(i,j,k)= omz(i,j,k) + 2.0*omearth*siang*u(i,j,k)
                  enddo
              enddo
          enddo

end



subroutine thermal1(fn,gn,scalar,scn,csource)
!************************************************************************
      use grid_size

      integer src_start,src_stop

      complex gn(nyp,nz,nxh),fn(nyp,nz,nxh)
      complex scn(nyp,nz,nxh), csource(nyp,nz,nxh)
      complex scalar(nyp,nz,nxh)
      complex forcebdy1(nyp,nz,nxh),forcebdy2(nyp,nz,nxh)

      common/fbdy/ forcebdy1,forcebdy2
      common/thermcor/ omearth,rotang,grav,virtual
      common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
      common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
      common/iocntrl/irstrt,nsteps,iprnfrq
      common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z
      common/polymer2/ ampbdy1,ampbdy2,c11amp,tfrac,polyrate,src_start,src_stop
      common/pre5f/ it

!         introduce buoyancy
    
          buoy = grav/virtual

!         time parameters

          tau =   tfrac*nsteps*dt       !radius from t0 to step function transition zones

          if (time .gt. tau) then
              tamp = 0.0
              tamp1 = 1.0
          else
              tamp = 1.0
            tamp1 = 0.0
          end if
        ! Adding scalar source delay
        if (it .ge. src_start .and. it .le. src_stop) then
            tamp2 = 1.0
        else
            tamp2 = 0.0
        end if

         do k=1,nxh
             do j=1,nz
                 do i=1,nyp 
                    scn(i,j,k)  = scn(i,j,k) + tamp2*csource(i,j,k)
                     fn(i,j,k)  = fn(i,j,k) +  tamp1*buoy*scalar(i,j,k) !vertical direction force
                     gn(i,j,k)  = gn(i,j,k) +  tamp*forcebdy1(i,j,k) + forcebdy2(i,j,k) !  x-direction force
                 end do
             end do
         end do
end



subroutine derivscji(u,v,w,c11,c12,c13,c21,c22,c23,c31,c32,c33, &
        u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw,           &
        dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313,  &
        dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323,  &
        dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333)
!***************************************************************
!  calculate velocity gradient tensor and derivs of conformation tensor  
!  tj and rah 6/16/2014
!
!************************************************************************
      use grid_size
      complex u(nyp,nz,nxh),v(nyp,nz,nxh),w(nyp,nz,nxh)
!        conformation tensor
      complex c11(nyp,nz,nxh),c12(nyp,nz,nxh),c13(nyp,nz,nxh)
      complex c21(nyp,nz,nxh),c22(nyp,nz,nxh),c23(nyp,nz,nxh)
      complex c31(nyp,nz,nxh),c32(nyp,nz,nxh),c33(nyp,nz,nxh)      
!        velocity gradient tensor 
      complex u11(nyp,nz,nxh),u12(nyp,nz,nxh),u13(nyp,nz,nxh)
      complex u21(nyp,nz,nxh),u22(nyp,nz,nxh),u23(nyp,nz,nxh) 
      complex u31(nyp,nz,nxh),u32(nyp,nz,nxh),u33(nyp,nz,nxh)
!        Laplacian
      complex,dimension(nyp,nz,nxh) :: Lu,Lv,Lw
!        derivatives of conformation tensor
      complex dc111(nyp,nz,nxh),dc112(nyp,nz,nxh),dc113(nyp,nz,nxh) 
      complex dc211(nyp,nz,nxh),dc212(nyp,nz,nxh),dc213(nyp,nz,nxh)
      complex dc311(nyp,nz,nxh),dc312(nyp,nz,nxh),dc313(nyp,nz,nxh)  

      complex dc121(nyp,nz,nxh),dc122(nyp,nz,nxh),dc123(nyp,nz,nxh)
      complex dc221(nyp,nz,nxh),dc222(nyp,nz,nxh),dc223(nyp,nz,nxh)
      complex dc321(nyp,nz,nxh),dc322(nyp,nz,nxh),dc323(nyp,nz,nxh)

      complex dc131(nyp,nz,nxh),dc132(nyp,nz,nxh),dc133(nyp,nz,nxh)  
      complex dc231(nyp,nz,nxh),dc232(nyp,nz,nxh),dc233(nyp,nz,nxh)  
      complex dc331(nyp,nz,nxh),dc332(nyp,nz,nxh),dc333(nyp,nz,nxh)  

      complex wrk1(nyp,nz,nxh),wrk2(nyp,nz,nxh),wrk3(nyp,nz,nxh)   
      complex wrk4(nyp,nz,nxh),wrk5(nyp,nz,nxh),wrk6(nyp,nz,nxh)
      complex wrk7(nyp,nz,nxh),wrk8(nyp,nz,nxh),wrk9(nyp,nz,nxh)

      complex im
      common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
   
         im = (0.0,1.0)

!        calculate the velocity gradient tensor

        call cderiv(u,u12) !du/dy
        call cderiv(v,u22) !dv/dy
        call cderiv(w,u32) !dw/dy

        ! Laplacian terms
        call cderiv(u12,wrk1) ! d^2 u / dy^2
        call cderiv(u22,wrk2) ! d^2 v / dy^2
        call cderiv(u32,wrk3) ! d^2 w / dy^2

      do k=1,nxh
          do j=1,nz
              do i=1,nyp
                  u11(i,j,k) = im*wavx(k)*u(i,j,k)
                  u13(i,j,k) = im*wavz(j)*u(i,j,k)
                  u21(i,j,k) = im*wavx(k)*v(i,j,k)
                  u23(i,j,k) = im*wavz(j)*v(i,j,k)
                  u31(i,j,k) = im*wavx(k)*w(i,j,k)
                  u33(i,j,k) = im*wavz(j)*w(i,j,k)

                ! Laplacian terms
                Lu(i,j,k) = im*wavx(k)*u11(i,j,k) + wrk1(i,j,k) + im*wavz(j)*u13(i,j,k)
                Lv(i,j,k) = im*wavx(k)*u21(i,j,k) + wrk2(i,j,k) + im*wavz(j)*u23(i,j,k)
                Lw(i,j,k) = im*wavx(k)*u31(i,j,k) + wrk3(i,j,k) + im*wavz(j)*u33(i,j,k)
              end do
          end do
      end do


!    Calculate derivs of conformation tensor 

        call cderiv(c11,wrk1) !dc11/dy
        call cderiv(c12,wrk2) !dc12/dy
        call cderiv(c13,wrk3) !dc13/dy
        call cderiv(c21,wrk4) !dc21/dy
        call cderiv(c22,wrk5) !dc22/dy
        call cderiv(c23,wrk6) !dc23/dy
        call cderiv(c31,wrk7) !dc31/dy
        call cderiv(c32,wrk8) !dc32/dy
        call cderiv(c33,wrk9) !dc33/dy

      do k=1,nxh
          do j=1,nz
              do i=1,nyp
                   dc111(i,j,k) = im*wavx(k)*c11(i,j,k)
                   dc112(i,j,k) = wrk1(i,j,k)
                   dc113(i,j,k) = im*wavz(j)*c11(i,j,k)
            
                   dc121(i,j,k) = im*wavx(k)*c12(i,j,k)
                   dc122(i,j,k) = wrk2(i,j,k)
                   dc123(i,j,k) = im*wavz(j)*c12(i,j,k)
            
                   dc131(i,j,k) = im*wavx(k)*c13(i,j,k)
                   dc132(i,j,k) = wrk3(i,j,k)
                   dc133(i,j,k) = im*wavz(j)*c13(i,j,k)
            
                   dc211(i,j,k) = im*wavx(k)*c21(i,j,k)
                   dc212(i,j,k) = wrk4(i,j,k)
                   dc213(i,j,k) = im*wavz(j)*c21(i,j,k)
            
                   dc221(i,j,k) = im*wavx(k)*c22(i,j,k)
                   dc222(i,j,k) = wrk5(i,j,k)
                   dc223(i,j,k) = im*wavz(j)*c22(i,j,k)
            
                   dc231(i,j,k) = im*wavx(k)*c23(i,j,k)
                   dc232(i,j,k) = wrk6(i,j,k)
                   dc233(i,j,k) = im*wavz(j)*c23(i,j,k)
            
                   dc311(i,j,k) = im*wavx(k)*c31(i,j,k)
                   dc312(i,j,k) = wrk7(i,j,k)
                   dc313(i,j,k) = im*wavz(j)*c31(i,j,k)
            
                   dc321(i,j,k) = im*wavx(k)*c32(i,j,k)
                   dc322(i,j,k) = wrk8(i,j,k)
                   dc323(i,j,k) = im*wavz(j)*c32(i,j,k)
            
                   dc331(i,j,k) = im*wavx(k)*c33(i,j,k)
                   dc332(i,j,k) = wrk9(i,j,k)
                   dc333(i,j,k) = im*wavz(j)*c33(i,j,k)
              end do
          end do
      end do    

end



subroutine polyforce(str11n,str12n,str13n,str22n,str23n,str33n,t1,t2,t3)
      ! Calculates divergence of stress tensor
      use grid_size
      complex str11n(nyp,nz,nxh),str12n(nyp,nz,nxh),str13n(nyp,nz,nxh)
      complex str22n(nyp,nz,nxh),str23n(nyp,nz,nxh),str33n(nyp,nz,nxh) 
      complex t1(nyp,nz,nxh),t2(nyp,nz,nxh),t3(nyp,nz,nxh)
      complex wrk1(nyp,nz,nxh),wrk2(nyp,nz,nxh),wrk3(nyp,nz,nxh)
      complex im
      common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)

   
      im = (0.0,1.0)

      call cderiv(str12n,wrk1) !ds12/dy
      call cderiv(str22n,wrk2) !ds22/dy
      call cderiv(str23n,wrk3) !ds23/dy


      do k=1,nxh
        do j=1,nz
          do i=1,nyp
            t1(i,j,k)=im*(wavx(k)*str11n(i,j,k)+wavz(j)*str13n(i,j,k))+wrk1(i,j,k)
            t2(i,j,k)=im*(wavx(k)*str12n(i,j,k)+wavz(j)*str23n(i,j,k))+wrk2(i,j,k)
            t3(i,j,k)=im*(wavx(k)*str13n(i,j,k)+wavz(j)*str33n(i,j,k))+wrk3(i,j,k)
          end do
        end do
      end do

end


subroutine subforce(gn,fn,omz,t1,t2,t3)

    use grid_size

    complex gn(nyp,nz,nxh),fn(nyp,nz,nxh),omz(nyp,nz,nxh)
    complex t1(nyp,nz,nxh),t2(nyp,nz,nxh),t3(nyp,nz,nxh)

    do k=1,nxh
        do j=1,nz
            do i=1,nyp
                gn(i,j,k)  = gn(i,j,k)  + t1(i,j,k)
                fn(i,j,k)  = fn(i,j,k)  + t2(i,j,k)
                omz(i,j,k) = omz(i,j,k) + t3(i,j,k)
            enddo
        enddo
    enddo
end

subroutine temperbc(bcscltop,bcsclbot,w1,w2,wfft1,wfft2,wfft3,wfft4)
      use grid_size
      complex bcscltop(nz,nxh),bcsclbot(nz,nxh) 
      real tempreal(nyp,nz,nx)
      complex tempcomp(nyp,nz,nxh),wrk1(nyp,nz,nxh)
      real wfft1(1),wfft2(1),wfft3(1),wfft4(1),w1(1),w2(1)
      common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
      common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
      common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z

      common/sclgeom/ xshift,yshift,zshift,sigmax,sigmay,sigmaz

      
     do  k = 1,nyp
     !   ycor = yl/2.0*cos(pi*float(k-1)/float(ny))
     !   ysq = (ycor - (0.0 + yshift))**2
     !   betay = ysq/(2.*sigmay**2) 
        do  j = 1,nz       
           zcor = zl*(float(j-1)/(float(nz)))
           zc1 = zl/2.0 + zshift
           zsq = (zcor - zc1)**2
           betaz = zsq/(2.*sigmaz**2)
           do  i = 1,nx
              xcor = xl*(float(i-1)/(float(nx)))
              xc1 = xl/2.0 + xshift
              xsq = (xcor - xc1)**2
              betax = xsq/(2.*sigmax**2)
              
              tempreal(k,j,i) = deltat*exp(-(betax + betay))
           end do
        end do
     end do

!        xc = xl/2.0! -1.0  !xl/10.0
!        zc = zl/2.0 !zl/2.0-zl/10.0
!        xc1 = xc+xl/10.0
!        zc1 = zl/2.0+zl/10.0
!        rnx = float(nx)
!        rnz = float(nz)
!        sigmax = xl/40.0
!        sigmaz = sigmax !zl/10.0
!        gamma1 = 2.0
!
!      do k=1,nx
!          xcor1 = xl*(float(k-1)/rnx)
!          xsq = (xcor1 - xc)**2
!          beta1 = xsq/(2.*sigmax**2)
!          xsq1 = (xcor1 - xc1)**2      !new
!          beta11 = xsq1/(2.*sigmax**2) !new
!
!          do j=1,nz
!              zcor1 = zl*(float(j-1)/rnz)
!              zsq = (zcor1 - zc)**2
!              beta2 = zsq/(2.*sigmaz**2)
!              zsq1 = (zcor1 - zc1)**2      !new
!              beta21 = zsq1/(2.*sigmaz**2) !new
!        
!              do i=1,nyp
!                  tempreal(i,j,k) = deltat*exp( -(beta1+beta2) ) 
!!                  tempreal(i,j,k) = deltat*exp( -(beta1+beta2) ) +(deltat)*exp( -(beta11+beta21) )
!              end do
!          end do
!      end do


      call scram(tempreal,tempcomp)
      call xyzfft(tempcomp,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
      call norm(tempcomp)


      call c0derbw(tempcomp,wrk1,1)


      do k=1,nxh
          do j=1,nz
              bcscltop(j,k)= cmplx(0.0,0.0) 
              bcsclbot(j,k)= cmplx(0.0,0.0) ! wrk1(1,j,k) 
          end do
      end do

end



subroutine polyrhs(scalar,wrkc,wrk1,cnl,bcbot,bctop)
!
!   if theta.le..99 crank-nickolson,
!   works on both odd and even parts simultaneously
!
      use grid_size
      complex scalar(nyp,nz,nxh),wrkc(nyp,nz,nxh),wrk1(nyp,nz,nxh)
      complex cnl(nyp,nz,nxh)
      complex bcbot(nz,nxh),bctop(nz,nxh)

      ! Common variables
      real pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
      real re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
      real zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z
      common/data1/ pi,dt,theta,wavz,wavx,c,yc
      common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
      common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z
!
!**********************************************************************
!                                                                     *
!     this subroutine evaluates the rhs of the poisson eqn. for the   *
!     particular phi.                                                 *
!                                                                     *
!     see p. 139 jfm v.177 (1987) kmm                                 *
!                                                                     *
!     rhs=-(1.-thta)(phi"(n)-(kx*kx+kz*kz)*phi(n))/thta               *
!         -re*phi(n)/(dt*thta)                                        *
!         -(3fn-fnm1)*re/(2*thta)                                     *
!                                                                     *
!**********************************************************************
!
!    first reset bc's to 0 since they may have been altered in 
!    pressure sol'n.


      do k=1,nxh
          do j=1,nz
            bctop(j,k)=(0.0,0.0)  !these arrays are defined as bcscltop
            bcbot(j,k)=(0.0,0.0)  ! and bclsclbot in main 
          end do
      end do

      do k=1,nxh
          do j=1,nz
            do i=1,nyp
                wrkc(i,j,k)= -(re*diffpoly)/(dt*theta)*scalar(i,j,k)-(re*diffpoly)/(2.*theta)*(cnl(i,j,k))
            end do
          end do
      end do    

      if(theta.gt..99) then
       do k=1,nxh
           do j=1,nz
             do i=1,nyp
               scalar(i,j,k) = wrkc(i,j,k)
             end do
           end do
       end do
      end if



      do k=1,nxh
          do j=1,nz
            w2 = wavz(j)**2 +wavx(k)**2
            do i=1,nyp
                wrkc(i,j,k) = wrkc(i,j,k) +(1.-theta)*w2*scalar(i,j,k)/theta
            end do
          end do
      end do

!  first deriv.
      call cderiv(scalar,wrk1)

!  second deriv.
      call cderiv(wrk1,scalar)

      do k=1,nxh
          do j=1,nz
              do i=1,nyp
                scalar(i,j,k)=wrkc(i,j,k) -(1.-theta)*scalar(i,j,k)/theta
              end do
          end do
      end do
end



subroutine polynl(c11nl,c12nl,c13nl,c22nl,c23nl,c33nl,    &
           c11n,c12n,c13n,c22n,c23n,c33n,c11nm1,c12nm1,   &
           c13nm1,c22nm1,c23nm1,c33nm1,str11n,str12n,     &
           str13n,str22n,str23n,str33n,str11nm1,str12nm1, &
           str13nm1,str22nm1,str23nm1,str33nm1)

      use grid_size
!        nonlinear terms for conformation tensor
      complex c11nl(nyp,nz,nxh),c12nl(nyp,nz,nxh),c13nl(nyp,nz,nxh)
      complex c22nl(nyp,nz,nxh),c23nl(nyp,nz,nxh)
      complex c33nl(nyp,nz,nxh) 

      complex c11n(nyp,nz,nxh),c12n(nyp,nz,nxh),c13n(nyp,nz,nxh)
      complex c22n(nyp,nz,nxh),c23n(nyp,nz,nxh),c33n(nyp,nz,nxh) ! nonlinear terms for conformation
      complex c11nm1(nyp,nz,nxh),c12nm1(nyp,nz,nxh),c13nm1(nyp,nz,nxh)
      complex c22nm1(nyp,nz,nxh),c23nm1(nyp,nz,nxh),c33nm1(nyp,nz,nxh) ! nonlinear terms for conformation
      complex str11n(nyp,nz,nxh),str12n(nyp,nz,nxh),str13n(nyp,nz,nxh)
      complex str22n(nyp,nz,nxh),str23n(nyp,nz,nxh),str33n(nyp,nz,nxh) ! nonlinear terms for conformation
      complex str11nm1(nyp,nz,nxh),str12nm1(nyp,nz,nxh),str13nm1(nyp,nz,nxh)
      complex str22nm1(nyp,nz,nxh),str23nm1(nyp,nz,nxh),str33nm1(nyp,nz,nxh) ! nonlinear terms for conformation
      common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
      common/polymer1/ zlmax,tpoly,qbeta
      common/polymer3/ ipolyflag
        

      tpolyinv = 1./tpoly
      do k=1,nxh
          do j=1,nz
              do i=1,nyp
                  c11nl(i,j,k)= 3.*(c11n(i,j,k)-tpolyinv*str11n(i,j,k))-(c11nm1(i,j,k)-tpolyinv*str11nm1(i,j,k))  ! note factor of 3 for ab 2nd order 
                  c12nl(i,j,k)= 3.*(c12n(i,j,k)-tpolyinv*str12n(i,j,k))-(c12nm1(i,j,k)-tpolyinv*str12nm1(i,j,k))
                  c13nl(i,j,k)= 3.*(c13n(i,j,k)-tpolyinv*str13n(i,j,k))-(c13nm1(i,j,k)-tpolyinv*str13nm1(i,j,k))
                  c22nl(i,j,k)= 3.*(c22n(i,j,k)-tpolyinv*str22n(i,j,k))-(c22nm1(i,j,k)-tpolyinv*str22nm1(i,j,k))
                  c23nl(i,j,k)= 3.*(c23n(i,j,k)-tpolyinv*str23n(i,j,k))-(c23nm1(i,j,k)-tpolyinv*str23nm1(i,j,k))
                  c33nl(i,j,k)= 3.*(c33n(i,j,k)-tpolyinv*str33n(i,j,k))-(c33nm1(i,j,k)-tpolyinv*str33nm1(i,j,k))
              end do
          end do
      end do

end


subroutine sclrhs(scalar,wrkc,wrk1,scn,scnm1,bcbot,bctop)

!   if theta.le..99 crank-nickolson,
!   works on both odd and even parts simultaneously
!
     use grid_size
     complex scalar(nyp,nz,nxh),wrkc(nyp,nz,nxh),wrk1(nyp,nz,nxh)
     complex scn(nyp,nz,nxh),scnm1(nyp,nz,nxh)
     complex bcbot(nz,nxh),bctop(nz,nxh)
     common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
     common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
     common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z
!
!**********************************************************************
!                                                                     *
!     this subroutine evaluates the rhs of the poisson eqn. for the   *
!     particular phi.                                                 *
!                                                                     *
!     see p. 139 jfm v.177 (1987) kmm                                 *
!                                                                     *
!     rhs=-(1.-thta)(phi"(n)-(kx*kx+kz*kz)*phi(n))/thta               *
!         -re*phi(n)/(dt*thta)                                        *
!         -(3fn-fnm1)*re/(2*thta)                                     *
!                                                                     *
!**********************************************************************
!
!    first reset bc's to 0 since they may have been altered in 
!    pressure sol'n.


        temp = 0.0 ! -1.0
  
     do k=1,nxh
         do j=1,nz
             bctop(j,k)=(0.0,0.0)  !these arrays are defined as bcscltop
             bcbot(j,k)=(0.0,0.0)  ! and bclsclbot in main 
         end do
     end do

     bctop(1,1) = cmplx(temp,0.0)
     bcbot(1,1) = cmplx(0.0,0.0)

     do k=1,nxh
         do j=1,nz
              do i=1,nyp
                   wrkc(i,j,k)= -(re*diff)/(dt*theta)*scalar(i,j,k)-(re*diff)/(2.*theta)*(3.*scn(i,j,k) -scnm1(i,j,k))
              end do
         end do
     end do


     if (theta .gt. 0.99) then
         do k=1,nxh
             do j=1,nz
                 do i=1,nyp
                     scalar(i,j,k) = wrkc(i,j,k)
                     return
                 end do
             end do
         end do
     end if



     do k=1,nxh
         do j=1,nz
            w2 = wavz(j)**2 +wavx(k)**2
            do i=1,nyp
               wrkc(i,j,k) = wrkc(i,j,k) +(1.-theta)*w2*scalar(i,j,k)/theta
            end do
         end do
     end do

!  first deriv.
     call cderiv(scalar,wrk1)

!  second deriv.
     call cderiv(wrk1,scalar)


     do k=1,nxh
         do j=1,nz
             do i=1,nyp
               scalar(i,j,k)=wrkc(i,j,k) -(1.-theta)*scalar(i,j,k)/theta
             end do
         end do
     end do
end


subroutine bforce1(a,b)
!****************************************************************
!  calculate a force array to be used to initiate a vortex pair.
!
!
!  the total force specified is
!
!                f_total = amp * f(y) * f(x) 
!
!          modified from original jds code by rah on 7/26/96
!****************************************************************

    use grid_size
    real a(nyp,nz,nx), b(nyp,nz,nx)
    common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
    common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z
    common/polymer2/ ampbdy1,ampbdy2,c11amp,tfrac
    common/data22/ forbeta
    real bdyfx, radx, radius, xi, yi, zi
    common/vortRing/ bdyfx, radx, radius, xi, yi, zi

      rny = float(ny)
      rnz = float(nz)
      rnx = float(nx)
      pi =  4.* atan(1.)
      twopi = 2.*pi

! Specify the force application parameters

!   Loading these variables from setup file - Ryan 7/19/22
!       xi = xl/2.      !centroid of volume over which force is applied
!       zi = zl/2. 
!       yi = 0.0 !-.9 !-.8
!
!       radx = 0.25     !radii of core region for force application
!       radius = 0.25

        yi = yi - 1.0 ! Need to adjust y inside this routine only - Ryan 7/19/22


! zero out force array

      do i =1,nyp
        do j =1,nz 
          do k =1,nx 
            a(i,j,k) =  0.0
            b(i,j,k) =  0.0
          end do
        end do
      end do

! force array tapering in y-direction
      zk2 = pi !wavenumber in y = 2.*pi/yhl

      do i=1,nyp
        arg = pi*(float(i-1)/rny)
        y =  cos(arg) ! full channel !0.5*(cos(arg)-1.) half-channel
        ysq = (y - yi)**2

        fykol = sin(zk2*y) 

! force array tapering in r-direction (f1 & f2)

        do j=1,nz
            z = zl*(float(j-1)/rnz)
            zsq = (z - zi)**2

            radsq = ysq + zsq
            zr  = sqrt(radsq)
            argrad = forbeta*(radius - zr)
            fr = 0.5*(1.0 + tanh(argrad))


          do k=1,nx
            x = xl*(float(k-1)/rnx)
            xsq = (x - xi)**2
            xx  = sqrt(xsq)
            argx = forbeta*(radx - xx)
            fx = 0.5*(1.0 + tanh(argx))   !force profile in x-direction 


! total force distribution combining tanh profiles

            a(i,j,k) = ampbdy1*fr*fx
            b(i,j,k) = ampbdy2*fykol
          end do
        end do
      end do

end

subroutine c11init(a)
!**************************************************************** !
!  calculate a force array to be used to initiate a vortex pair.
!
!
!  the total force specified is
!
!                f_total = amp * f(y) * f(x) 
!
!          modified from original jds code by rah on 7/26/96
!**************************************************************** !

    use grid_size
      real a(nyp,nz,nx)
      common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
      common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z
      common/polymer2/ ampbdy1,ampbdy2,c11amp,tfrac

! some stuff

      rny = float(ny)
      rnz = float(nz)
      rnx = float(nx)
      pi =  4.* atan(1.)
      twopi = 2.*pi

! specify the force application parameters

       xi = xl/2.      !centroid of volume over which force is applied
       zi = zl/2. 
       yi = 0.0 !-.9 !-.8

       radx = 0.25     !radii of core region for force application
       radz = 0.25
       rady = 0.25

       betax = 10.0        !tanh decay factor for x
       betaz = 10.0
       betay = 50.0        !tanh decay factor for y

              
! zero out force array

      do i =1,nyp
        do j =1,nz 
          do k =1,nx 
            a(i,j,k) =  0.0
          end do
        end do
      end do

! force array tapering in y-direction

      do i=1,nyp
        arg = pi*(float(i-1)/rny)
        y =  cos(arg) ! full channel !0.5*(cos(arg)-1.) half-channel
        ysq = (y - yi)**2
        yy = sqrt(ysq)
        argy = betay*(rady - yy)
        fy = 0.5*(1.0 + tanh(argy))  !force profile in fs-normal direction

        ! force array tapering in r-direction (f1 & f2)

        do j=1,nz
            z = zl*(float(j-1)/rnz)
            zp = z - zi
            zsq = (z - zi)**2
            zz  = sqrt(zsq)
            argz = betaz*(radz - zz)
            fz = 0.5*(1.0 + tanh(argz))   !force profile in z-direction 

          do k=1,nx
            x = xl*(float(k-1)/rnx)
            xp = x - xi
            xsq = (x - xi)**2
            xx  = sqrt(xsq)
            argx = betax*(radx - xx)
            fx = 0.5*(1.0 + tanh(argx))   !force profile in x-direction 


            ! total force distribution combining tanh profiles

            a(i,j,k) = c11z + c11amp*fx*fz  ! c11z + c11amp*fy*fx*fz
          end do
        end do
      end do

end subroutine
