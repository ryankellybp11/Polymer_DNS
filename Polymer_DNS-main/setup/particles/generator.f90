program generator

    implicit none

    real dt,yl,zl,xl,delay
    real xp, yp, zp, up, vp, wp
    real xc,yc,zc,L,rad,pi
    real r, theta,phi
    integer npart,j,tp,skip,new_part
    integer i,k,npx,npy,npz
    character(len = 126) label
    real vyc, vzc

    open(110, file='../dns.config', status = 'old')
    do skip = 1,5
        read(110,*)
    end do
    read(110,*) dt
    do skip = 1,7
        read(110,*)
    end do
    read(110,*) yl
    read(110,*) zl
    read(110,*) xl
    do skip = 1,31
        read(110,*)
    end do
    read(110,*) npart
    read(110,*) new_part
    do skip = 1,33
        read(110,*)
    end do
    read(110,*) xc,yc,zc
    read(110,*) L, rad
    close(110) 

    open(111, file='../vort.config',status = 'old')
    do skip = 1,15
        read(111,*) 
    end do
    read(111,*) vyc, vzc
    close(111)
 
    pi = acos(-1.0)

    if (new_part .eq. 11) then

        ! Non-random, uniform distribution (regular grid)
        ! Assumes npart = 2^16, not sure the best way to generalize
        npx = npart**(3./8)
        npy = npart**(5./16)
        npz = npart**(5./16)

        up = 0
        vp = 0
        wp = 0
        tp = 1
 
        do k = 1,npz
          do j = 1,npy
            do i = 1,npx
              xp = (i-1)*(xl/npx)
              yp = (j-1)*(yl/npy)
              zp = (k-1)*(zl/npz)

              ! Write these to init_particle.dat file
                  if (j .eq. 1 .and. i .eq. 1 .and. k .eq. 1) then
                      call system('rm particles.dat')
                      open(22, file = 'particles.dat')
                      label = '!      x       |      y       |      z       |     vpx      |     vpy      |      vpz     |  t0 (timestep)'
                      write(22,"(a126)") label 
                  else
                      open(22, file = 'particles.dat', position = "append")
                  end if
                  write(22,"(e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,i5)") xp,yp,zp,up,vp,wp,tp
                  close(22) 
            end do
          end do
        end do
                    
    else if (new_part .eq. 12) then
        ! Uniform distrubtion on a given y-plane assumes npart is a 5th power,
        ! e.g., 4^5 = 1024
        npx = npart**(3./5.)
        npz = npart**(2./5.)
  
        up = 0
        vp = 0
        wp = 0
        tp = 1

        print *,'Please enter the y-coordinate of the particle plane: '
        read *,yp 

        if (yp .ne. 1.0) then
            npx = npx/2
            do k = 1,npz
                do i = 1,npx
                    xp = (i-1)*(xl/npx)
                    zp = (k-1)*(zl/npz)
    
                    ! Write these to init_particle.dat file
                    if (i .eq. 1 .and. k .eq. 1) then
                        call system('rm particles.dat')
                        open(22, file = 'particles.dat')
                        label = '!      x       |      y       |      z       |     vpx      |     vpy      |      vpz     |  t0 (timestep)'
                        write(22,"(a126)") label 
                    else
                        open(22, file = 'particles.dat', position = "append")
                    end if
                    write(22,"(e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,i5)") xp,yp,zp,up,vp,wp,tp
                    close(22) 
                end do
            end do
            yp = yl - yp
        end if
        do k = 1,npz
            do i = 1,npx
                xp = (i-1)*(xl/npx)
                zp = (k-1)*(zl/npz)
    
                open(22, file = 'particles.dat', position = "append")
                write(22,"(e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,i5)") xp,yp,zp,up,vp,wp,tp
                close(22) 
            end do
        end do
    else
        call random_seed()
        do j = 1,npart
            ! Generate random numbers for each component
            call random_number(xp)
            call random_number(yp)
            call random_number(zp)
            if (new_part .eq. 1) then ! Uniform distribution
       
                ! Scale these according to domain size
                xp = xl*xp
                yp = yl*yp
                zp = zl*zp
    
    
            else if (new_part .eq. 2) then ! Sphere
    
                r = rad*xp**(1/3.)
                theta = yp*2*pi
                phi = acos(2*zp - 1)
    
                xp = xc + 4*L + r*cos(theta)*sin(phi)
                yp = yc + r*sin(theta)*sin(phi)
                zp = zc + r*cos(phi)
            
            else if (new_part .eq. 3) then ! Circle
    
    !            theta = 2*pi*float(j)/float(npart)
                theta = pi/2.0
                xp = xc + 4*L
                yp = yc + 0.25*rad*sin(theta)
                zp = zc + 0.25*rad*cos(theta)
            
            else if (new_part .eq. 4) then ! Thin y-z sheet
    
                xp = xc
                yp = yl*yp
                zp = zl*zp
    
            else if (new_part .eq. 5) then ! Near-wall
        
                if (j .le. floor(float(npart)/2.0)) then
                    yp = yl/10.0*yp
                else
                    yp = yl*(1.0 - yp/10.0)
                end if
    
                xp = xl*xp
                zp = zl*zp
    
            else if (new_part .eq. 6) then ! cylinder
                
               call random_number(theta) 
               theta = theta*2*pi
               xp = xc + 4*L + xp*2.0
               yp = yc + 0.25*rad*sqrt(yp)*sin(theta)
               zp = zc + 0.25*rad*sqrt(zp)*cos(theta)
        
            else if (new_part .eq. 7) then ! horizontal line
                
                xp = xl/2.0
                yp = vyc - mod(j,npart/10)*0.05
                zp = (j*10/npart)*zl/10.0
            else if (new_part .eq. 10) then ! Steady release
    
                ! All particles release from the same spot (or spots)
                xp = 1.0
                yp = yl/8.0
                zp = zl/2.0
    
            end if
            
            ! Particle Velocities - As of now, just released from rest
            up = 0.0
            vp = 0.0
            wp = 0.0
        
        ! Particle release time - testing regular particle release
            delay = 10 ! Arbitrary delay for testing - number of timesteps between particle release
            if (new_part .eq. 10) then
               tp = (j-1)*dt*delay 
            else
                tp = 1
            end if
     
        ! Write these to init_particle.dat file
            if (j .eq. 1) then
                call system('rm particles.dat')
                open(22, file = 'particles.dat')
                label = '!      x       |      y       |      z       |     vpx      |     vpy      |      vpz     |  t0 (timestep)'
                write(22,"(a126)") label 
            else
                open(22, file = 'particles.dat', position = "append")
            end if
            write(22,"(e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6,1x,i5)") xp,yp,zp,up,vp,wp,tp
            close(22) 
        
        end do
    end if

end program
