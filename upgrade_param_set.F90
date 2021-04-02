subroutine upgrade_param_set
	use param_def, only: geometry, orgl, upgd, rstar &
		& , ununiform_flag, ix_ununi, dx00 &
		& , myrank &
		& , ix0, jx0, kx0, ib, jb, kb, rstar, merr
	implicit none
	include 'mpif.h'

	integer :: i,j,k
	real(KIND(0.d0)) :: dx_unif,dy_unif,dz_unif
	integer :: nxx
	real(KIND(0.d0)) :: xrange,xrange0,xrange1,dx11,fdx
  integer :: idf = 11
  integer :: ix00,jx00,kx00
  namelist /upgrade/ix00,jx00,kx00,ununiform_flag,ix_ununi,dx00

!==================================================================
  ! specify calculation domain in upgraded data
!	upgd%xmax = 0.94d0*rstar; upgd%xmin = 0.71d0*rstar
!	upgd%xmax = rstar + 0.7d8 ; upgd%xmin = rstar - 5.444d8
	upgd%xmax = rstar + 0.7d8 ; upgd%xmin = 0.71d0*rstar
	upgd%ymax = orgl%ymax; upgd%ymin = orgl%ymin
	upgd%zmax = orgl%zmax; upgd%zmin = orgl%zmin
!==================================================================

	upgd%margin = orgl%margin

  if(myrank == 0) then
	  open(idf,file='./namelist.in',status='old',action='read')
	  read(idf,nml=upgrade)
	  close(idf)
  endif

	call MPI_BCAST(ix00,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(jx00,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(kx00,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)
  call MPI_BCAST(ununiform_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,merr)
  call MPI_BCAST(ix_ununi,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)
  call MPI_BCAST(dx00,1,MPI_REAL8,0,MPI_COMM_WORLD,merr)

	upgd%ix00 = ix00 + 2*upgd%margin
	upgd%jx00 = jx00 + 2*upgd%margin
	upgd%kx00 = kx00 + 2*upgd%margin
  upgd%ununiform_flag = ununiform_flag
  upgd%ix_ununi = ix_ununi
  upgd%dx00 = dx00

	upgd%mtype = orgl%mtype
	allocate(upgd%x00(upgd%ix00))
	allocate(upgd%y00(upgd%jx00))
	allocate(upgd%z00(upgd%kx00))

	if(ununiform_flag) then
!-----
    ! inhomogeneous spacing
    !
    ! homogeneous at top mpi thread
    ! This inhomogeneous setting width of the domain and number of grid
    ! Then the grid spacing increases downward linearly
    ! see an explanation in Powerpoint file provided by H. Hotta

    xrange = upgd%xmax - upgd%xmin             ! x domain size
    xrange0 = dx00*dble(ix_ununi)              ! x domain size in top MPI thread, in which grid is uniform
    xrange1 = xrange - xrange0                 ! x domain size in the other MPI thread than top
    nxx = upgd%ix00 - 2*upgd%margin - ix_ununi ! number of grid for inhomogeneous grid

    fdx = 2.d0*(xrange1 - dx00*dble(nxx))/dble((nxx - 4)*nxx) ! increase rate of the grid

    !-------------------------------------
    ! grid for uniform area
    ! set uniform area
    ! set the geometry at the top margin

    ! set the geometry around the top boundary
    ! the geometry is defined downward from the top
    upgd%x00(upgd%ix00 - upgd%margin) = upgd%xmax - 0.5d0*dx00

    do i = upgd%ix00 - upgd%margin + 1, upgd%ix00
      upgd%x00(i) = upgd%x00(i - 1) + dx00
    enddo

    do i = upgd%ix00 - upgd%margin - 1, upgd%ix00 - upgd%margin - ix_ununi + 1, -1
      upgd%x00(i) = upgd%x00(i + 1) - dx00
    enddo

    ! uniform area finish
    !-------------------------------------

    ! first 2 uniform grids in ununiform area
    do i = upgd%ix00 - upgd%margin - ix_ununi, upgd%ix00 - upgd%margin - ix_ununi - 1, -1
      upgd%x00(i) = upgd%x00(i + 1) - dx00
    enddo

    ! ordinary ununiform gridding
    dx11 = dx00
    do i = upgd%ix00 - upgd%margin - ix_ununi - 2, 5, -1
      upgd%x00(i) = upgd%x00(i + 1) - dx11
      dx11 = dx11 + fdx
    enddo

    ! final 4 uniform grid including margin
    do i = 4, 1, -1
      upgd%x00(i) = upgd%x00(i + 1) - dx11
    enddo

	else
		! Define uniform grid

  	! define global geometry (Samples are for uniform)
  	dx_unif = (upgd%xmax - upgd%xmin)/real(upgd%ix00 - 2*upgd%margin)
  	upgd%x00(1) = upgd%xmin + (0.5d0 - dble(upgd%margin))*dx_unif
  	do i = 2, upgd%ix00
      upgd%x00(i) = upgd%x00(i - 1) + dx_unif
    enddo
	endif

  dy_unif = (upgd%ymax - upgd%ymin)/real(upgd%jx00 - 2*upgd%margin)
  upgd%y00(1) = upgd%ymin + (0.5d0 - dble(upgd%margin))*dy_unif
  do j = 2, upgd%jx00
    upgd%y00(j) = upgd%y00(j - 1) + dy_unif
  enddo

  dz_unif = (upgd%zmax - upgd%zmin)/real(upgd%kx00 - 2*upgd%margin)
  upgd%z00(1) = upgd%zmin + (0.5d0 - dble(upgd%margin))*dz_unif
  do k = 2, upgd%kx00
    upgd%z00(k) = upgd%z00(k - 1) + dz_unif
 enddo
 
	allocate(upgd%xstt(0:ix0-1))
	allocate(upgd%xend(0:ix0-1))
	allocate(upgd%xwdt(0:ix0-1))
	allocate(upgd%ystt(0:jx0-1))
	allocate(upgd%yend(0:jx0-1))
	allocate(upgd%ywdt(0:jx0-1))
	allocate(upgd%zstt(0:kx0-1))
	allocate(upgd%zend(0:kx0-1))
	allocate(upgd%zwdt(0:kx0-1))

  upgd%xstt(0) = 1
  upgd%xend(ix0 - 1) = upgd%ix00

  upgd%ystt(0) = 1
  upgd%yend(jx0 - 1) = upgd%jx00

  upgd%zstt(0) = 1
  upgd%zend(kx0 - 1) = upgd%kx00

  call specify_location(ix0,orgl%ix00,upgd%ix00,orgl%margin,upgd%margin,orgl%nx,orgl%x00,upgd%x00,upgd%xstt,upgd%xend,upgd%xwdt)
  call specify_location(jx0,orgl%jx00,upgd%jx00,orgl%margin,upgd%margin,orgl%ny,orgl%y00,upgd%y00,upgd%ystt,upgd%yend,upgd%ywdt)
  call specify_location(kx0,orgl%kx00,upgd%kx00,orgl%margin,upgd%margin,orgl%nz,orgl%z00,upgd%z00,upgd%zstt,upgd%zend,upgd%zwdt)

  upgd%nxg = upgd%xwdt(ib)
  upgd%nyg = upgd%ywdt(jb)
  upgd%nzg = upgd%zwdt(kb)

  upgd%nx = upgd%nxg - 2*upgd%margin
  upgd%ny = upgd%nyg - 2*upgd%margin
  upgd%nz = upgd%nzg - 2*upgd%margin

  allocate(upgd%x(upgd%nxg))
  allocate(upgd%y(upgd%nyg))
  allocate(upgd%z(upgd%nzg))

  upgd%x = upgd%x00(upgd%xstt(ib):upgd%xend(ib))
  upgd%y = upgd%y00(upgd%ystt(jb):upgd%yend(jb))
  upgd%z = upgd%z00(upgd%zstt(kb):upgd%zend(kb))

	return
end subroutine upgrade_param_set
