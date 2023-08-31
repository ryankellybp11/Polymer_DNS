module grid_size
    integer, parameter :: nx_ = 256, ny_ = 128, nz_ = 128, bftail_ = 0, npart = 1024 
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
 
    integer, parameter :: nmax = 8*nz*ny*nx        !nmax should be at least 4 times the product of the two largest grid dimensions.
    integer, parameter :: nx32 = (3*nx)/2
    integer, parameter :: nz32 = (3*nz)/2
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

    real :: ycoord(nyp)
end module grid_size
