module derivs
implicit none
! This module contains derivative-related subroutines used in main dns code
contains

    subroutine gradv(u,v,w,u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw)
    !---------------------------------------------------------------------!
    ! This subroutine computes velocity gradient terms as well as the     !
    ! Laplacian of the velocity field. These are used for calculation of  !
    ! the swirl criterion as well as for particle motion integration.     !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: u,v,w,Lu,Lv,Lw
    complex, dimension(nyp,nz,nxh) :: u11,u12,u13
    complex, dimension(nyp,nz,nxh) :: u21,u22,u23
    complex, dimension(nyp,nz,nxh) :: u31,u32,u33
    
    ! Calculation variables
    complex, dimension(nyp,nz,nxh) :: wrk1,wrk2,wrk3
    complex :: im
    integer :: i,j,k
    real    :: wn2
    
    ! Solver variables
    real :: wavx(nxh),wavz(nz) 
    real :: c(nyp)
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/      wavx,wavz,c
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    im = (0.0,1.0)
    
    ! Calculate y-derivatives
    call cderiv(u,u12) ! du/dy
    call cderiv(v,u22) ! dv/dy
    call cderiv(w,u32) ! dw/dy
    
    call cderiv(u12,wrk1) ! d^2(u)/dy^2
    call cderiv(u22,wrk2) ! d^2(v)/dy^2
    call cderiv(u32,wrk3) ! d^2(w)/dy^2
    
    ! Calculate x- and z-derivatives spectrally
    do k = 1,nxh
        do j = 1,nz
            wn2 = wavx(k)**2 + wavz(j)**2
            do i = 1,nyp
                u11(i,j,k) = im*wavx(k)*u(i,j,k) ! du/dx
                u21(i,j,k) = im*wavx(k)*v(i,j,k) ! dv/dx
                u31(i,j,k) = im*wavx(k)*w(i,j,k) ! dw/dx
    
                u13(i,j,k) = im*wavz(j)*u(i,j,k) ! du/dz
                u23(i,j,k) = im*wavz(j)*v(i,j,k) ! dv/dz
                u33(i,j,k) = im*wavz(j)*w(i,j,k) ! dw/dz
    
                ! Laplacian terms
                Lu(i,j,k) = wrk1(i,j,k) - wn2*u(i,j,k)
                Lv(i,j,k) = wrk2(i,j,k) - wn2*v(i,j,k)
                Lw(i,j,k) = wrk3(i,j,k) - wn2*w(i,j,k)
    
            end do
        end do
    end do
    
    end subroutine gradv
    
    !---------------------------------------------------------------------!
#IFDEF SCALAR    
    subroutine gradscl(scalar,sclx,scly,sclz,wrkc)
    !---------------------------------------------------------------------!
    ! This subroutine computes gradient of the scalar and stores them in  !
    ! sclx, scly, and sclz                                                !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: scalar,sclx,scly,sclz,wrkc
    
    ! Calculation variables
    complex :: im
    integer :: i,j,k
    
    ! Solver variables
    real :: wavx(nxh),wavz(nz) 
    real :: c(nyp)
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/      wavx,wavz,c
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    im = (0.0,1.0)
   
    call cderiv(scalar,scly)

    !$omp parallel do
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                sclx(i,j,k) = im*wavx(k)*scalar(i,j,k)
                sclz(i,j,k) = im*wavz(j)*scalar(i,j,k)
            end do
        end do
    end do
    !$omp end parallel do 
    
    end subroutine gradscl
#ENDIF    
    !---------------------------------------------------------------------!
#IFDEF POLYMER   
    subroutine derivscji(c11,c12,c13,c21,c22,c23,c31,c32,c33,                   &
                         dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313, &
                         dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323, &
                         dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333)  
    !---------------------------------------------------------------------!
    ! This subroutine computes derivatives of the conformation tensor     !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    use omp_lib
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: c11,c12,c13,c21,c22,c23,c31,c32,c33

    complex, dimension(nyp,nz,nxh) :: dc111,dc112,dc113,dc121,dc122,dc123,dc131,dc132,dc133
    complex, dimension(nyp,nz,nxh) :: dc211,dc212,dc213,dc221,dc222,dc223,dc231,dc232,dc233
    complex, dimension(nyp,nz,nxh) :: dc311,dc312,dc313,dc321,dc322,dc323,dc331,dc332,dc333
    
    ! Calculation variables
    complex :: im
    integer :: i,j,k
    
    ! Solver variables
    real :: wavx(nxh),wavz(nz) 
    real :: c(nyp)
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/      wavx,wavz,c
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    im = (0.0,1.0)
  
    ! Calculate y-derivatives
    call cderiv(c11,dc112)
    call cderiv(c12,dc122)
    call cderiv(c13,dc132)
    call cderiv(c21,dc212)
    call cderiv(c22,dc222)
    call cderiv(c23,dc232)
    call cderiv(c31,dc312)
    call cderiv(c32,dc322)
    call cderiv(c33,dc332)

    ! Calculate x- and z- derivatives
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                dc111(i,j,k) = im*wavx(k)*c11(i,j,k) 
                dc113(i,j,k) = im*wavz(j)*c11(i,j,k)

                dc121(i,j,k) = im*wavx(k)*c12(i,j,k) 
                dc123(i,j,k) = im*wavz(j)*c12(i,j,k)

                dc131(i,j,k) = im*wavx(k)*c13(i,j,k) 
                dc133(i,j,k) = im*wavz(j)*c13(i,j,k)

                dc211(i,j,k) = im*wavx(k)*c21(i,j,k) 
                dc213(i,j,k) = im*wavz(j)*c21(i,j,k)

                dc221(i,j,k) = im*wavx(k)*c22(i,j,k) 
                dc223(i,j,k) = im*wavz(j)*c22(i,j,k)

                dc231(i,j,k) = im*wavx(k)*c23(i,j,k) 
                dc233(i,j,k) = im*wavz(j)*c23(i,j,k)

                dc311(i,j,k) = im*wavx(k)*c31(i,j,k) 
                dc313(i,j,k) = im*wavz(j)*c31(i,j,k)

                dc321(i,j,k) = im*wavx(k)*c32(i,j,k) 
                dc323(i,j,k) = im*wavz(j)*c32(i,j,k)

                dc331(i,j,k) = im*wavx(k)*c33(i,j,k) 
                dc333(i,j,k) = im*wavz(j)*c33(i,j,k)
            end do
        end do
    end do

    end subroutine derivscji
#ENDIF    
    !---------------------------------------------------------------------!
    
    subroutine cderiv(f,df)
    !---------------------------------------------------------------------!
    !  This subroutine computes the wall-normal derivative (df) of f      !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: f,df
    
    ! Calculation variables
    integer :: i,j,k
    
    ! Solver variables
    real :: gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/solver/     gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    do i = 1,2 ! Dummy loop to fool tau. Without this, tau-instrumented code doesn't work
    end do     ! Editor's note: idk what that means, but I'm not messing with it
    
    ! Top wall
    do k = 1,nxh
        do j = 1,nz
            df(nyp,j,k) = 0.0
            df(ny,j,k)  = 2.0*float(ny)*f(nyp,j,k)*dyde
        end do
    end do
    
    ! Interior
    do j = 1,nz
        do i = nym,1,-1
            do k = 1,nxh
                df(i,j,k) = df(i+2,j,k) + 2.0*float(i)*f(i+1,j,k)*dyde
            end do
        end do
    end do
    
    ! Bottom wall
    do k = 1,nxh
        do j = 1,nz
            df(1,j,k) = 0.5*df(1,j,k) 
        end do
    end do
    
    end subroutine cderiv   
    
    !---------------------------------------------------------------------!
    
    subroutine c0derbw(wrkc,wrk1,iy)
    !---------------------------------------------------------------------!
    !  This subroutine evaluates the function at the y = -1 boundary. It  !
    !  assumes that wrkc is spectral in both z and y. wrk1 is the normal  !
    !  derivative at y = -1 and is physical in y.                         !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: wrkc,wrk1
    integer :: iy
    
    ! Calculation variables
    integer :: i,j,k
    real    :: sgn
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    do k = 1,nxh
        do j = 1,nz
            wrk1(iy,j,k) = 0.0
        end do
    end do
    
    sgn = -1.0
    
    do i = 1,nyp
        sgn = -sgn
        do k = 1,nxh
            do j = 1,nz
                wrk1(iy,j,k) = wrk1(iy,j,k) + sgn*wrkc(i,j,k)
            end do
        end do
    end do
    
    end subroutine c0derbw
    
    !---------------------------------------------------------------------!
    
    subroutine c1derbw(wrkc,wrk1,iy)
    !---------------------------------------------------------------------!
    !  This subroutine evaluates the normal derivative at the y = -1      !
    !  boundary. It assumes that wrkc is spectral in both z and y. wrk1   !
    !  is the normal derivative at y = -1 and is physical in y.           !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: wrkc,wrk1
    integer :: iy
    
    ! Calculation variables
    integer :: i,j,k
    real    :: sgn,rp
    
    ! Solver variables
    real :: gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/solver/ gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    do k = 1,nxh
        do j = 1,nz 
            wrk1(iy,j,k) = 0.0
        end do
    end do
    
    sgn = 1.0
    
    do i = 1,nyp
        sgn = -sgn
        rp = float(i-1)
        do k = 1,nxh
            do j = 1,nz 
                wrk1(iy,j,k) = wrk1(iy,j,k) + rp*rp*sgn*wrkc(i,j,k)*dyde
            end do
        end do
    end do
    
    end subroutine c1derbw
    
    
    !---------------------------------------------------------------------!
    
    subroutine c2derbw(wrkc,wrk1,iy)
    !---------------------------------------------------------------------!
    !  This subroutine evaluates the second derivative at the y = -1      !
    !  boundary. It assumes that wrkc is spectral in both z and y. wrk1   !
    !  is the normal derivative at y = -1 and is physical in y.           !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: wrkc,wrk1
    integer :: iy
    
    ! Calculation variables
    integer :: i,j,k
    real    :: sgn,rp2
    
    ! Solver variables
    real :: gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/solver/ gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    do k = 1,nxh
        do j = 1,nz 
            wrk1(iy,j,k) = 0.0
        end do
    end do
    
    sgn = 1.0
    
    do i = 1,nyp
        sgn = -sgn
        rp2 = float(i-1)*float(i-1)
        do k = 1,nxh
            do j = 1,nz 
                wrk1(iy,j,k) = wrk1(iy,j,k) + (rp2 - 1.0)*rp2*sgn*wrkc(i,j,k)*dyde*dyde/3.0
            end do
        end do
    end do
    
    end subroutine c2derbw
    
    !---------------------------------------------------------------------!
    
    subroutine c0dertw(wrkc,wrk1,iy)
    !---------------------------------------------------------------------!
    !  This subroutine evaluates the function at the y = +1 boundary. It  !
    !  assumes that wrkc is spectral in both z and y. wrk1 is the normal  !
    !  derivative at y = -1 and is physical in y.                         !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: wrkc,wrk1
    integer :: iy
    
    ! Calculation variables
    integer :: i,j,k
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    do k = 1,nxh
        do j = 1,nz
            wrk1(iy,j,k) = 0.0
        end do
    end do
    
    do i = 1,nyp
        do k = 1,nxh
            do j = 1,nz
                wrk1(iy,j,k) = wrk1(iy,j,k) + wrkc(i,j,k)
            end do
        end do
    end do
    
    end subroutine c0dertw
    
    !---------------------------------------------------------------------!
    
    subroutine c1dertw(wrkc,wrk1,iy)
    !---------------------------------------------------------------------!
    !  This subroutine evaluates the normal derivative at the y = +1      !
    !  boundary. It assumes that wrkc is spectral in both z and y. wrk1   !
    !  is the normal derivative at y = -1 and is physical in y.           !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: wrkc,wrk1
    integer :: iy
    
    ! Calculation variables
    integer :: i,j,k
    real    :: rp
    
    ! Solver variables
    real :: gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/solver/ gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    do k = 1,nxh
        do j = 1,nz 
            wrk1(iy,j,k) = 0.0
        end do
    end do
    
    do i = 1,nyp
        rp = float(i-1)
        do k = 1,nxh
            do j = 1,nz 
                wrk1(iy,j,k) = wrk1(iy,j,k) + rp*rp*wrkc(i,j,k)*dyde
            end do
        end do
    end do
    
    end subroutine c1dertw
    
    !---------------------------------------------------------------------!
    
    subroutine c2dertw(wrkc,wrk1,iy)
    !---------------------------------------------------------------------!
    !  This subroutine evaluates the second derivative at the y = +1      !
    !  boundary. It assumes that wrkc is spectral in both z and y. wrk1   !
    !  is the normal derivative at y = -1 and is physical in y.           !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: wrkc,wrk1
    integer :: iy
    
    ! Calculation variables
    integer :: i,j,k
    real    :: rp2
    
    ! Solver variables
    real :: gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/solver/ gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    do j = 1,nz 
        do k = 1,nxh
            wrk1(iy,j,k) = 0.0
        end do
    end do
    
    do i = 1,nyp
        rp2 = float(i-1)*float(i-1)
        do j = 1,nz 
            do k = 1,nxh
                wrk1(iy,j,k) = wrk1(iy,j,k) + (rp2 - 1.0)*rp2*wrkc(i,j,k)*dyde*dyde/3.0
            end do
        end do
    end do
    
    end subroutine c2dertw
    
    !---------------------------------------------------------------------!
    
    subroutine cderiv1d(f,df)
    !---------------------------------------------------------------------!
    !  This subroutine computes the wall-normal derivative (df) of f in 1D!
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    real, dimension(nyp) :: f,df
    
    ! Calculation variables
    integer :: i
    
    ! Solver variables
    real :: gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/solver/     gain,ugain,theta,alpha,beta,dyde
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    ! Top wall
    df(nyp) = 0.0
    df(ny) = 2.0*float(ny)*f(nyp)*dyde
    
    ! Interior
    do i = nym,1,-1
        df(i) = df(i+2) + 2.0*float(i)*f(i+1)*dyde
    end do
    
    ! Bottom wall
    df(1) = 0.5*df(1)
    
    end subroutine cderiv1d


end module derivs
