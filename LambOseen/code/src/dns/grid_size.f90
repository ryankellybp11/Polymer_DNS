module grid_size
    integer, parameter :: nx_ = 2, ny_ = 256, nz_ = 256, bftail_ = 0, npart = 0, qn = 65536 
    integer, parameter :: ny = ny_
    integer, parameter :: nz = nz_
    integer, parameter :: nx = nx_ 
    integer, parameter :: nyp = ny + 1
    integer, parameter :: nym = ny - 1 
    integer, parameter :: nxh = nx/2
    integer, parameter :: nyh = ny/2
    integer, parameter :: nzh = nz/2
    integer, parameter :: nyhp = nyh + 1
    integer, parameter :: nyhm = nyh - 1 
 
    integer, parameter :: mx = 3*nx/2
    integer, parameter :: mz = 3*nz/2
    integer, parameter :: mz2 = 10
    integer, parameter :: mx2 = 10     

    integer, parameter :: mx4 = mx/4
    integer, parameter :: mx4p = mx4 + 1
    integer, parameter :: mzp = mz + 1
    integer, parameter :: mxp = mx + 1

    integer, parameter :: nyp2 = 2*(ny + 1)

    integer, parameter :: nxp = nx + 1
    integer, parameter :: nxp2 = nx + 2
    integer, parameter :: nzhp = nzh + 1
    integer, parameter :: mxp2 = mx + 2 
    integer, parameter :: mzp2 = mz + 2
    integer, parameter :: mzm1 = mz - 1
    integer, parameter :: mzm2 = mz - 2
    integer, parameter :: nwrk = 4*mzp2*mxp2 
    real, parameter :: rmz=1.0/mz

    real, dimension(nyp)   :: ycoord,seght

    ! Particle variables
    real, dimension(npart) :: xpart,ypart,zpart,upart,vpart,wpart,swirl_part
    
    integer :: particle_flag,CD_switch
	real    :: ratio,ap,C_mu,gravity

    ! Logger 
    integer :: cadence

    ! Lamb Oseen variable
    real    :: vortR

end module grid_size
