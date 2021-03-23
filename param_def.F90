module param_def
	integer, save :: ix0, jx0, kx0
	integer, save :: ib,jb,kb
	integer, dimension(:,:,:), allocatable :: ijkb2rank

	type :: geometry
		real(KIND(0.d0)) :: xmin,xmax,ymin,ymax,zmin,zmax
		real(KIND(0.d0)), dimension(:), allocatable :: x,y,z
		real(KIND(0.d0)), dimension(:), allocatable :: x00,y00,z00
		integer :: nx ,ny ,nz     ! number of grid in each MPI thread
		integer :: nxg,nyg,nzg    ! number of grid in each MPI thread with margin
		integer :: ix0,jx0,kx0    ! number of MPI thread
		integer :: ix00,jx00,kx00 ! number of 
		integer :: margin,mtype,npe
		integer :: ix_ununi
		logical :: ununiform_flag
		real(KIND(0.d0)) :: dx00
		real(KIND(0.d0)), dimension(:,:,:,:), allocatable :: qq
		integer, dimension(:), allocatable :: xstt,xend,xwdt
		integer, dimension(:), allocatable :: ystt,yend,ywdt
		integer, dimension(:), allocatable :: zstt,zend,zwdt
	end type
	
	type(geometry), save :: orgl, upgd
		
	real(KIND(0.d0)), parameter :: rstar = 6.9598947d+10
	logical, save :: ununiform_flag 
	integer, save :: ix_ununi
	real(KIND(0.d0)), save :: dx00

	real(KIND(0.d0)), dimension(:), allocatable, save :: dx0,dx1,dy0,dy1,dz0,dz1
	integer, dimension(:), allocatable, save :: iloc,jloc,kloc
	! MPI parameter
	integer, save :: myrank,npe
	integer, save :: imin,imax,jmin,jmax,kmin,kmax
	integer, save :: merr

end module param_def