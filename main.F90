!===========================================================
program main
	use param_def, only: geometry, orgl, upgd &
		& ,myrank,npe &
		& ,ix0,jx0,kx0, ib, jb, kb &
		& ,dx0,dx1,dy0,dy1,dz0,dz1,iloc,jloc,kloc &
		& ,imin,imax,jmin,jmax,kmin,kmax &
		& ,ijkb2rank, merr, rstar
	use change, only: change_judge
	implicit none
	include 'mpif.h'

	integer :: i,j,k,m,ibt,jbt,kbt

	integer :: nd
	real(KIND(0.d0)), parameter :: pi = 3.141592653589d0
	character*16 :: R2D2_ver
	character*4  :: caseid,caseid_out
	character*8  :: c_nd_read,c_nd
	integer :: idf = 11
	integer, parameter :: ndims = 3
	integer, dimension(ndims) :: dims
	logical, dimension(ndims) :: periodic
	logical :: reorder
	integer, dimension(:,:), allocatable :: xyz
	integer :: mpi_comm_world_cart
	integer :: np
	integer :: xdcheck,ydcheck,zdcheck
	integer :: mmx
	character*11 :: geometry_char
	integer, dimension(4) :: gsize,ssize,start
	integer(Kind=MPI_OFFSET_KIND), parameter :: disp = 0
	integer :: fh
	integer, dimension(mpi_status_size) :: mstatus

	integer :: global_orgl,global_upgd
	logical :: exist
	real(KIND(0.d0)) :: time = 0.d0
	integer :: itmp
	real(KIND(0.d0)) :: rtmp
	logical :: ltmp
	character*20 :: server_name
	character*11 :: change_char
	
	namelist /info/caseid,caseid_out,c_nd_read,ix0,jx0,kx0
!===========================================================
	! initialize MPI
	call MPI_INIT(merr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, npe, merr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, merr)
!-----------------------------------------------------------
	! read basic parameters
	if(myrank == 0) then
		open(idf,file='./namelist.in',status='old',action='read')
		read(idf,nml=info)
		close(idf)
	endif

	call MPI_BCAST(caseid    ,4,MPI_CHARACTER,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(caseid_out,4,MPI_CHARACTER,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(c_nd_read,8,MPI_CHARACTER,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(ix0,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(jx0,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(kx0,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)

	! check caseid and caseid_out
	if(caseid == caseid_out) then
		if(myrank == 0) then
			write(*,*) "The origin and the destination are the same. stop..."
		endif
		call mpi_finalize(merr)		
		stop
	endif
	
	! check destination directory
	inquire(file="../run/"//caseid_out,exist=exist)
	if(exist) then
		if(myrank == 0) then
			write(*,*) "The destination directory "//caseid_out//" already exists. stop..."
		endif
		call mpi_finalize(merr)
		stop
	endif

	if(npe /= ix0*jx0*kx0) then
		if(myrank == 0) then
			write(*,*) "Required and prepared MPI is different"
			write(*,*) "Required: ",ix0*jx0*kx0
			write(*,*) "Prepared: ",npe
		endif
		call mpi_finalize(merr)
		stop
	endif

	! copy original directory
	if(myrank == 0) then
		call system('rsync -a --exclude="data" ../run/'//caseid//'/ ../run/'//caseid_out)
		call system('mkdir -p ../run/'//caseid_out//'/data/param')
		call system('mkdir -p ../run/'//caseid_out//'/data/qq')
		call system('mkdir -p ../run/'//caseid_out//'/data/slice')
		call system('mkdir -p ../run/'//caseid_out//'/data/remap/qq')
		call system('mkdir -p ../run/'//caseid_out//'/data/remap/vl')
		call system('mkdir -p ../run/'//caseid_out//'/data/time/mhd')
		call system('mkdir -p ../run/'//caseid_out//'/data/time/tau')
		call system('mkdir -p ../run/'//caseid_out//'/data/tau')
	endif

!-----------------------------------------------------------

	! allocate 3D MPI thread
	dims(1) = ix0; dims(2) = jx0; dims(3) = kx0; periodic = .false.
	call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periodic, reorder, mpi_comm_world_cart, merr)
	allocate(xyz(0:npe-1,ndims))
	do np = 0,npe-1
		call MPI_CART_COORDS(mpi_comm_world_cart,np,ndims,xyz(np,:),merr)
	enddo

	ib = xyz(myrank,1)
	jb = xyz(myrank,2)
	kb = xyz(myrank,3)

	allocate(ijkb2rank(0:ix0-1,0:jx0-1,0:kx0-1))
	do np = 0,npe-1
		ijkb2rank(xyz(np,1),xyz(np,2),xyz(np,3)) = np
	enddo
	
!-----------------------------------------------------------
	if(myrank == 0) then
		open(idf,file='../run/'//caseid//'/data/param/nd.dac',form='formatted')
		read(idf,*) nd
		close(idf)
	endif	

	call MPI_BCAST(nd,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)

	if(myrank == 0) then
		if (c_nd_read == 'end_step') then
			if (nd == (nd/2)*2) then
				c_nd = 'e'
			else
				c_nd = 'o'
			endif
		else
			read(c_nd_read,*) nd
		  	write(c_nd,'(i8.8)') nd
		endif

		write(*,'(A)') "### Data info ###"
		write(*,'(A)') "Read information from namelist.in"
		write(*,'(A,A)') "Input  caseid is: ",caseid
		write(*,'(A,A)') "Output caseid is: ",caseid_out
		write(*,'(A,I0)') "Read time step: ",nd
	endif	
	call MPI_BCAST(c_nd,8,MPI_CHARACTER,0,MPI_COMM_WORLD,merr)

!-----------------------------------------------------------
	! read basic parameters of original data
	if(myrank == 0) then
		open(idf,file='../run/'//caseid//'/data/param/params.dac',form='formatted')
		read(idf,'(A16)') R2D2_ver
		read(idf,*) xdcheck
		read(idf,*) ydcheck
		read(idf,*) zdcheck
		read(idf,*) orgl%margin
		read(idf,*) orgl%nx
		read(idf,*) orgl%ny
		read(idf,*) orgl%nz
		read(idf,*) orgl%npe
		read(idf,*) orgl%ix0
		read(idf,*) orgl%jx0
		read(idf,*) orgl%kx0
		read(idf,*) orgl%mtype
		read(idf,*) orgl%xmax
		read(idf,*) orgl%xmin
		read(idf,*) orgl%ymax
		read(idf,*) orgl%ymin
		read(idf,*) orgl%zmax
		read(idf,*) orgl%zmin
		read(idf,*) rtmp ! dtout
		read(idf,*) rtmp ! tend
		read(idf,*) itmp ! swap
		read(idf,*) itmp ! ix_e
		read(idf,*) itmp ! jx_e
		read(idf,*) itmp ! ixr
		read(idf,*) itmp ! jxr
		read(idf,*) itmp ! m_in
		read(idf,*) itmp ! m_tu
		read(idf,*) rtmp ! dtout_tau
		read(idf,*) rtmp ! ifac
		read(idf,*) rtmp ! omfac
		read(idf,*) rtmp ! om00
		read(idf,*) itmp ! jc
		read(idf,*) itmp ! kc
		read(idf,*) rtmp ! rstar
		read(idf,*) rtmp ! lstar
		read(idf,*) server_name
		read(idf,*) ltmp ! rte_multiray_flag
		read(idf,*) rtmp ! potential_alpha
		read(idf,*) orgl%ununiform_flag ! ununiform_flag
		read(idf,*) orgl%ix_ununi
		read(idf,*) orgl%dx00
		read(idf,*) geometry_char

		close(idf)

		orgl%ix00 = orgl%nx*orgl%ix0 + 2*orgl%margin
		orgl%jx00 = orgl%ny*orgl%jx0 + 2*orgl%margin
		orgl%kx00 = orgl%nz*orgl%kx0 + 2*orgl%margin
	endif

	call MPI_BCAST(orgl%margin,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(orgl%mtype ,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(orgl%ix00  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(orgl%jx00  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(orgl%kx00  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,merr)
		
	allocate(orgl%x00(orgl%ix00))
	allocate(orgl%y00(orgl%jx00))
	allocate(orgl%z00(orgl%kx00))

	if(myrank == 0) then
		open(idf,file='../run/'//caseid//'/data/param/back.dac',form='unformatted')
		read(idf) orgl%x00,orgl%y00,orgl%z00
		close(idf)
	endif

	call MPI_BCAST(orgl%x00,orgl%ix00,MPI_REAL8,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(orgl%y00,orgl%jx00,MPI_REAL8,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(orgl%z00,orgl%kx00,MPI_REAL8,0,MPI_COMM_WORLD,merr)

	orgl%nx = (orgl%ix00 - 2*orgl%margin)/ix0
	orgl%ny = (orgl%jx00 - 2*orgl%margin)/jx0
	orgl%nz = (orgl%kx00 - 2*orgl%margin)/kx0

	orgl%nxg = orgl%nx + 2*orgl%margin
	orgl%nyg = orgl%ny + 2*orgl%margin
	orgl%nzg = orgl%nz + 2*orgl%margin

	allocate(orgl%x(orgl%nxg))
	allocate(orgl%y(orgl%nyg))
	allocate(orgl%z(orgl%nzg))

	call MPI_BCAST(orgl%xmin,1,MPI_REAL8,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(orgl%xmax,1,MPI_REAL8,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(orgl%ymin,1,MPI_REAL8,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(orgl%ymax,1,MPI_REAL8,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(orgl%zmin,1,MPI_REAL8,0,MPI_COMM_WORLD,merr)
	call MPI_BCAST(orgl%zmax,1,MPI_REAL8,0,MPI_COMM_WORLD,merr)

	allocate(orgl%xstt(0:ix0-1))
	allocate(orgl%xend(0:ix0-1))
	allocate(orgl%xwdt(0:ix0-1))
	allocate(orgl%ystt(0:jx0-1))
	allocate(orgl%yend(0:jx0-1))
	allocate(orgl%ywdt(0:jx0-1))
	allocate(orgl%zstt(0:kx0-1))
	allocate(orgl%zend(0:kx0-1))
	allocate(orgl%zwdt(0:kx0-1))

	do ibt = 0,ix0 - 1
		orgl%xstt(ibt) = orgl%nx*ibt + 1
		orgl%xend(ibt) = orgl%xstt(ibt) + orgl%nxg - 1
		orgl%xwdt(ibt) = orgl%nxg
	enddo

	do jbt = 0,jx0 - 1
		orgl%ystt(jbt) = orgl%ny*jbt + 1
		orgl%yend(jbt) = orgl%ystt(jbt) + orgl%nyg - 1
		orgl%ywdt(jbt) = orgl%nyg
	enddo

	do kbt = 0,kx0 - 1
		orgl%zstt(kbt) = orgl%nz*kbt + 1
		orgl%zend(kbt) = orgl%zstt(kbt) + orgl%nzg - 1
		orgl%zwdt(kbt) = orgl%nzg
	enddo

	orgl%x = orgl%x00(orgl%xstt(ib):orgl%xend(ib))
	orgl%y = orgl%y00(orgl%ystt(jb):orgl%yend(jb))
	orgl%z = orgl%z00(orgl%zstt(kb):orgl%zend(kb))

!-----------------------------------------------------------
	! geometry for upgraded data

	call upgrade_param_set

!-----------------------------------------------------------

	allocate(iloc(upgd%nxg))
	allocate(jloc(upgd%nyg))
	allocate(kloc(upgd%nzg))

	allocate(dx0(upgd%nxg))
	allocate(dx1(upgd%nxg))
	allocate(dy0(upgd%nyg))
	allocate(dy1(upgd%nyg))
	allocate(dz0(upgd%nzg))
	allocate(dz1(upgd%nzg))

	allocate(upgd%qq(upgd%nxg,upgd%nyg,upgd%nzg,upgd%mtype))

	call id_loc(orgl%x,upgd%x,orgl%nxg,upgd%nxg,imin,imax,iloc,dx0,dx1)
	call id_loc(orgl%y,upgd%y,orgl%nyg,upgd%nyg,jmin,jmax,jloc,dy0,dy1)	
	call id_loc(orgl%z,upgd%z,orgl%nzg,upgd%nzg,kmin,kmax,kloc,dz0,dz1)

!-----------------------------------------------------------

	gsize(1) = orgl%ix00
	gsize(2) = orgl%jx00
	gsize(3) = orgl%kx00
	gsize(4) = orgl%mtype

	ssize(1) = orgl%nxg
	ssize(2) = orgl%nyg
	ssize(3) = orgl%nzg
	ssize(4) = orgl%mtype

	start(1) = orgl%nx*ib
	start(2) = orgl%ny*jb
	start(3) = orgl%nz*kb
	start(4) = 0

	call MPI_TYPE_CREATE_SUBARRAY(4,gsize,ssize,start, &
		& MPI_ORDER_FORTRAN,MPI_REAL8,global_orgl,merr)
	call MPI_TYPE_COMMIT(global_orgl,merr)

	allocate(orgl%qq(orgl%nxg,orgl%nyg,orgl%nzg,orgl%mtype))
	call mpi_file_open(mpi_comm_world, '../run/'//caseid//'/data/qq/qq.dac.'//trim(c_nd), &
    & mpi_mode_rdonly, mpi_info_null, fh, merr)
	call mpi_file_set_view(fh, disp, mpi_real8, global_orgl, "native", mpi_info_null, merr)
	mmx = orgl%nxg*orgl%nyg*orgl%nzg*orgl%mtype
	!call mpi_file_read(fh, orgl%qq, mmx, mpi_real8, mstatus, merr)
        call mpi_file_read_all(fh, orgl%qq, mmx, mpi_real8, mstatus, merr)
	call mpi_file_close(fh, merr)

	upgd%qq = 0.d0

        !$omp parallel do private(i,j,k,m)
	do k = kmin,kmax
	do j = jmin,jmax
	do i = imin,imax
	do m = 1,orgl%mtype
           upgd%qq(i,j,k,m) = ( &
          & + orgl%qq(iloc(i)+0,jloc(j)+0,kloc(k)+0,m)*dz1(k)*dy1(j)*dx1(i) &
          & + orgl%qq(iloc(i)+0,jloc(j)+0,kloc(k)+1,m)*dz0(k)*dy1(j)*dx1(i) &
          & + orgl%qq(iloc(i)+0,jloc(j)+1,kloc(k)+0,m)*dz1(k)*dy0(j)*dx1(i) &
          & + orgl%qq(iloc(i)+0,jloc(j)+1,kloc(k)+1,m)*dz0(k)*dy0(j)*dx1(i) &
          & + orgl%qq(iloc(i)+1,jloc(j)+0,kloc(k)+0,m)*dz1(k)*dy1(j)*dx0(i) &
          & + orgl%qq(iloc(i)+1,jloc(j)+0,kloc(k)+1,m)*dz0(k)*dy1(j)*dx0(i) &
          & + orgl%qq(iloc(i)+1,jloc(j)+1,kloc(k)+0,m)*dz1(k)*dy0(j)*dx0(i) &
          & + orgl%qq(iloc(i)+1,jloc(j)+1,kloc(k)+1,m)*dz0(k)*dy0(j)*dx0(i) &
          & )/(dx0(i) + dx1(i))/(dy0(j) + dy1(j))/(dz0(k) + dz1(k))
	enddo
	enddo
	enddo
	enddo

	call bc_mpi(upgd%qq,upgd%nxg,upgd%nyg,upgd%nzg,upgd%mtype,upgd%margin,ib,jb,kb,ix0,jx0,kx0,ijkb2rank,merr)

	gsize(1) = upgd%ix00
	gsize(2) = upgd%jx00
	gsize(3) = upgd%kx00
	gsize(4) = upgd%mtype

	ssize(1) = upgd%nxg
	ssize(2) = upgd%nyg
	ssize(3) = upgd%nzg
	ssize(4) = upgd%mtype

	start(1) = upgd%xstt(ib) - 1
	start(2) = upgd%ystt(jb) - 1
	start(3) = upgd%zstt(kb) - 1
	start(4) = 0

	call MPI_TYPE_CREATE_SUBARRAY(4,gsize,ssize,start, &
		& MPI_ORDER_FORTRAN,MPI_REAL8,global_upgd,merr)
	call MPI_TYPE_COMMIT(global_upgd,merr)

	call mpi_file_open(mpi_comm_world, '../run/'//caseid_out//'/data/qq/qq.dac.e', &
    & ior(mpi_mode_create, mpi_mode_wronly), mpi_info_null, fh, merr)
	call mpi_file_set_view(fh, disp, mpi_real8, global_upgd, "native", mpi_info_null, merr)
	mmx = upgd%nxg*upgd%nyg*upgd%nzg*upgd%mtype
	!call mpi_file_write(fh, upgd%qq, mmx, mpi_real8, mstatus, merr)
        call mpi_file_write_all(fh, upgd%qq, mmx, mpi_real8, mstatus, merr)
	call mpi_file_close(fh, merr)

	if(myrank == 0) then
		open(idf,file='../run/'//caseid_out//'/data/time/mhd/t.dac.e',form='unformatted',access='stream')
		write(idf) time
		close(idf)

		open(idf,file='../run/'//caseid_out//'/data/param/nd.dac',form='formatted')
		write(idf,'(i8,i8)') 0,0
		close(idf)

		open(idf,file='../run/'//caseid_out//'/data/cont_log.txt',form='formatted')
		write(idf,'(A)') 'The initial condition of this run is upgraded data from the other directory.'
		write(idf,'(A)') 'This file describes the original and the upgraded data'
		write(idf,'(A)') '### Original data ###'
		write(idf,*)
		write(idf,'(A)') 'Server: '//trim(server_name)
		write(idf,'(A)') 'Datadir: ../run/'//caseid//'/data'
		write(idf,'(A,i6)') 'Output step: ',nd
		call change_judge(orgl%xmin,upgd%xmin,change_char)
		write(idf,'(A,SP,ES10.3,A,SS,F5.3,A)') 'xmin = rsun',orgl%xmin-rstar,' or ',orgl%xmin/rstar,'rsun '//trim(change_char)
		call change_judge(orgl%xmax,upgd%xmax,change_char)
		write(idf,'(A,SP,ES10.3,A,SS,F5.3,A)') 'xmax = rsun',orgl%xmax-rstar,' or ',orgl%xmax/rstar,'rsun '//trim(change_char)
		if (trim(geometry_char) == 'Spherical' .or. trim(geometry_char) == 'YinYang') then
			call change_judge(orgl%ymin,upgd%ymin,change_char)
			write(idf,'(A,F6.2,A)') 'ymin = ',orgl%ymin/pi*180.d0,' [rad] '//trim(change_char)
			call change_judge(orgl%ymax,upgd%ymax,change_char)
			write(idf,'(A,F6.2,A)') 'ymax = ',orgl%ymax/pi*180.d0,' [rad] '//trim(change_char)
			call change_judge(orgl%zmin,upgd%zmin,change_char)
			write(idf,'(A,F6.2,A)') 'zmin = ',orgl%zmin/pi*180.d0,' [rad] '//trim(change_char)
			call change_judge(orgl%zmax,upgd%zmax,change_char)
			write(idf,'(A,F6.2,A)') 'zmax = ',orgl%zmax/pi*180.d0,' [rad] '//trim(change_char)
		else
			call change_judge(orgl%ymin,upgd%ymin,change_char)
			write(idf,'(A,ES10.3,A)') 'ymin = ',orgl%ymin,' [cm] '//trim(change_char)
			call change_judge(orgl%ymax,upgd%ymax,change_char)
			write(idf,'(A,ES10.3,A)') 'ymax = ',orgl%ymax,' [cm] '//trim(change_char)
			call change_judge(orgl%zmin,upgd%zmin,change_char)
			write(idf,'(A,ES10.3,A)') 'zmin = ',orgl%zmin,' [cm] '//trim(change_char)
			call change_judge(orgl%zmax,upgd%zmax,change_char)
			write(idf,'(A,ES10.3,A)') 'zmax = ',orgl%zmax,' [cm] '//trim(change_char)
		endif

		call change_judge(orgl%ix00,upgd%ix00,change_char)
		write(idf,'(A,I0,A)') 'nx0*ix0 = ',orgl%ix00 - 2*orgl%margin,trim(change_char)
		call change_judge(orgl%jx00,upgd%jx00,change_char)
		write(idf,'(A,I0,A)') 'ny0*jx0 = ',orgl%jx00 - 2*orgl%margin,trim(change_char)
		call change_judge(orgl%kx00,upgd%kx00,change_char)
		write(idf,'(A,I0,A)') 'nz0*kx0 = ',orgl%kx00 - 2*orgl%margin,trim(change_char)

		call change_judge(orgl%ununiform_flag,upgd%ununiform_flag,change_char)
		if(orgl%ununiform_flag) then
			write(idf,'(A)') 'ununiform_flag = .true.'//trim(change_char)
		else
			write(idf,'(A)') 'ununiform_flag = .false.'//trim(change_char)
		endif
		call change_judge(orgl%ix_ununi,upgd%ix_ununi,change_char)
		write(idf,'(A,I0,A)') 'ix_ununi = ',orgl%ix_ununi,trim(change_char)
		call change_judge(orgl%dx00,upgd%dx00,change_char)
		write(idf,'(A,ES10.3,A)') 'dx00 = ',orgl%dx00,trim(change_char)
		write(idf,*)
		write(idf,'(A)') '### Upgraded data ###'
		write(idf,*)
		call change_judge(orgl%xmin,upgd%xmin,change_char)
		write(idf,'(A,SP,ES10.3,A,SS,F5.3,A)') 'xmin = rsun',upgd%xmin-rstar,' or ',upgd%xmin/rstar,'rsun '//trim(change_char)
		call change_judge(orgl%xmax,upgd%xmax,change_char)
		write(idf,'(A,SP,ES10.3,A,SS,F5.3,A)') 'xmax = rsun',upgd%xmax-rstar,' or ',upgd%xmax/rstar,'rsun '//trim(change_char)
		if (trim(geometry_char) == 'Spherical' .or. trim(geometry_char) == 'YinYang') then
			call change_judge(orgl%ymin,upgd%ymin,change_char)
			write(idf,'(A,F6.2,A)') 'ymin = ',upgd%ymin/pi*180.d0,' [rad] '//trim(change_char)
			call change_judge(orgl%ymax,upgd%ymax,change_char)
			write(idf,'(A,F6.2,A)') 'ymax = ',upgd%ymax/pi*180.d0,' [rad] '//trim(change_char)
			call change_judge(orgl%zmin,upgd%zmin,change_char)
			write(idf,'(A,F6.2,A)') 'zmin = ',upgd%zmin/pi*180.d0,' [rad] '//trim(change_char)
			call change_judge(orgl%zmax,upgd%zmax,change_char)
			write(idf,'(A,F6.2,A)') 'zmax = ',upgd%zmax/pi*180.d0,' [rad] '//trim(change_char)
		else
			call change_judge(orgl%ymin,upgd%ymin,change_char)
			write(idf,'(A,ES10.3,A)') 'ymin = ',upgd%ymin,' [cm] '//trim(change_char)
			call change_judge(orgl%ymax,upgd%ymax,change_char)
			write(idf,'(A,ES10.3,A)') 'ymax = ',upgd%ymax,' [cm] '//trim(change_char)
			call change_judge(orgl%zmin,upgd%zmin,change_char)
			write(idf,'(A,ES10.3,A)') 'zmin = ',upgd%zmin,' [cm] '//trim(change_char)
			call change_judge(orgl%zmax,upgd%zmax,change_char)
			write(idf,'(A,ES10.3,A)') 'zmax = ',upgd%zmax,' [cm] '//trim(change_char)
		endif

		call change_judge(orgl%ix00,upgd%ix00,change_char)
		write(idf,'(A,I0,A)') 'nx0*ix0 = ',upgd%ix00 - 2*upgd%margin,trim(change_char)
		call change_judge(orgl%jx00,upgd%jx00,change_char)
		write(idf,'(A,I0,A)') 'ny0*jx0 = ',upgd%jx00 - 2*upgd%margin,trim(change_char)
		call change_judge(orgl%kx00,upgd%kx00,change_char)
		write(idf,'(A,I0,A)') 'nz0*kx0 = ',upgd%kx00 - 2*upgd%margin,trim(change_char)

		call change_judge(orgl%ununiform_flag,upgd%ununiform_flag,change_char)
		if(upgd%ununiform_flag) then
			write(idf,'(A)') 'ununiform_flag = .true.'//trim(change_char)
		else
			write(idf,'(A)') 'ununiform_flag = .false.'//trim(change_char)
		endif
		call change_judge(orgl%ix_ununi,upgd%ix_ununi,change_char)
		write(idf,'(A,I0,A)') 'ix_ununi = ',upgd%ix_ununi,trim(change_char)
		call change_judge(orgl%dx00,upgd%dx00,change_char)
		write(idf,'(A,ES10.3,A)') 'dx00 = ',upgd%dx00,trim(change_char)
		
		close(idf)
	endif
	
	call mpi_finalize(merr)
	stop
end program main
