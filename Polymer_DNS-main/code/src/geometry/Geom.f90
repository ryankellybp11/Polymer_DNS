program geom
    use grid_size

    implicit none

    integer perturb
    integer kwall, kmax, kbuff
    integer geomtype, flow_select
    integer bfhead            ! Index # for the Buffer zone head face.
    integer bftail            ! Index # for the Buffer zone tail face.
    integer bfwidth, bw2, bw2p1
    integer i, k, j, n, skip
    integer kgeom, imax, imin, jmax, jmin
    integer kmin2(2)
    real  ygeom, xmax, xmin, zmax, zmin, geomstart, ywidth, ymin, ymax
    real xi, zj, L, rad, xcenter, ycenter, zcenter
    real rsq
  
    integer imatrix(nyp,mzp,mxp)
    real yl, zl, xl, ywall, ytop
    real delxm, delzm, pi
    real ybuff, xlocation, ylocation, zlocation
  
    real geomheight, geomwidth
  
    ! Grid spacing in y direction
    integer grid_spacing
    real yfrac, a_spacing, r_spacing
      
    !******************************************************************************************
  
    pi = 2.0*acos(0.0)

    open(110, file = '../../../setup/dns.config', status = 'old')
 
    do skip = 1,12
        read(110,*) 
    end do
    read(110,*) bftail
    read(110,*) yl
    read(110,*) zl
    read(110,*) xl
    read(110,*) 
    read(110,*)
    read(110,*) 
    read(110,*) kwall
    do skip = 1,14
        read(110,*) 
    end do
    read(110,*) geomtype
    read(110,*)
    read(110,*)
    read(110,*) flow_select
    do skip = 1,43
        read(110,*) 
    end do
    read(110,*) xcenter, ycenter, zcenter
    read(110,*) L, rad
    close(110)
  

    ! Set y coordinate
    do k = 1, nyp
        ycoord(k) = (1. - cos(float(k-1) * pi / float(ny))) * (yl / 2.0)
    end do
  
    delxm = xl /float(mx)
    delzm = zl /float(mz)
    
    bfhead  = 2
    bfwidth = (bftail - bfhead)
    bw2 = bfwidth * 2
    bw2p1 = bw2 + 1
  
    ywall = ycoord(kwall)
    ytop = yl - ywall
  
   

    !------------------------------------------------------------!
    !                           Wall                             !
    !------------------------------------------------------------!
    do k = 1,nyp
       do j = 1,mz
          do i = 1,mx
             if (k .lt. kwall) then 
                imatrix(k,j,i) = 0
             else if (k .eq. kwall) then
                imatrix(k,j,i) = 1
             else 
                imatrix(k,j,i) = 0
             end if
          end do
       end do
    end do


    !------------------------------------------------------------!
    !                       Buffer Region                        !
    !------------------------------------------------------------!
    do k = 1,nyp
       do j = 1,mz
          do i = bfhead,bftail
             if (imatrix(k,j,i) .eq. 0) then
                imatrix(k,j,i) = 6
             end if
          end do
       end do
    end do


    !------------------------------------------------------------!
    !                         3D Brick                           !
    !------------------------------------------------------------!
    if(geomtype.eq.3) then
      xmin = XL*0.1
      xmax = xmin + 0.1*XL
      zmin = ZL*0.2
      zmax = ZL*0.8
      ygeom = 0.1*YL
      kbuff = kwall
      ybuff = ycoord(kbuff)
      ywidth = 0.25*YL
      ymin = ygeom
      ymax = ygeom + ywidth
      do i = bftail+1,mx
           xlocation = (float(i-1))*delxm
           if ((xlocation .le. xmax) .and. (xlocation .ge. xmin)) then
             do j = 1,mz
               zlocation = (float(j-1))*delzm
               if ((zlocation.le.zmax).and.(zlocation.ge.zmin)) then
                   do k = kbuff, nyp
                     ylocation = ycoord(k) - ybuff
                     if((ylocation .le. ymax) .and. (ylocation .ge. ymin)) then
                         imatrix(k,j,i) = 4
                     end if
                   end do
               end if
             end do
           end if
      end do
    end if
  
  
    !------------------------------------------------------------!
    !                         2D Brick                           !
    !------------------------------------------------------------!
    if(geomtype.eq.2) then
        geomstart = 0.4*XL
        geomwidth = 0.1*YL
        geomheight = 0.1*YL
        imin = floor(geomstart / delxm) + 1  
        imax =  floor((geomstart+geomwidth) / delxm) + 1
        kmin2 = [floor(0.1*float(nyp)), floor(0.7*float(nyp))]
        do n = 1,2
            kmax = kmin2(n) + floor(0.1*float(nyp))
            do j = 1,mz
             do i = imin,imax
               do k = kmin2(n), kmax
                 imatrix(k,j,i) = 4
               end do ! k
             end do ! i 
            end do ! j
        end do ! n
    end if


    ! This was my original plan to generate a vortex ring, but it turned out
    ! better by coding it into dns.f90 with the flow_select switch. However,
    ! this could be a useful reference for an additional type of forcefield
    ! later on, so I'll leave it commented out - Ryan 7/11/22
    
    !------------------------------------------------------------!
    !                       Vortex Ring                          !
    !------------------------------------------------------------!
!    if (geomtype .eq. 10) then
!        do i = 1,mx
!          xi = float(i-1)*delxm
!          do j = 1,mz
!            zj = float(j-1)*delzm
!            do k = 1,nyp
!              if (abs(xi - xcenter) .lt. L/2.) then 
!                 rsq = (zj - zcenter)**2 + (ycoord(k) - ycenter)**2
!                 if (rsq .lt. rad**2) then
!                     imatrix(k,j,i) = 100
!                 end if
!              end if
!            end do
!          end do
!        end do
!    end if
 
    
    !------------------------------------------------------------!
    !                      Suction Region                        !
    !                   (Based on flow type)                     !
    !------------------------------------------------------------!
    if (flow_select .eq. 4) then
         k = 0.703125*ny
    end if

    do i = 1,mx
       if (i .gt. (bfwidth/2 + 1)) then
          do j = 1,mz
             if (flow_select .eq. 4) then
                 imatrix(k,j,i) = 7
             else
                 imatrix(k,j,i) = 0
             end if
          end do
       end if
    end do


    
    open(111, file = 'geometry', status = 'replace', form = 'unformatted')
    write(111) imatrix
    close(111)
    
end program geom
