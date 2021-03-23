subroutine bc_mpi(qq,nxg,nyg,nzg,mtype,margin,ib,jb,kb,ix0,jx0,kx0,ijkb2rank,merr)
	implicit none
	include 'mpif.h'

	integer :: mmx,nerank
	integer, intent(in) :: nxg,nyg,nzg,mtype,margin,ib,jb,kb,ix0,jx0,kx0,merr
	integer, dimension(0:ix0-1,0:jx0-1,0:kx0-1), intent(in) :: ijkb2rank
	real(KIND(0.d0)), dimension(nxg   ,nyg   ,nzg   ,mtype), intent(inout) :: qq
	real(KIND(0.d0)), dimension(margin,nyg   ,nzg   ,mtype) :: bufsndx_up,bufsndx_dw,bufrcvx_up,bufrcvx_dw
	real(KIND(0.d0)), dimension(nxg   ,margin,nzg   ,mtype) :: bufsndy_up,bufsndy_dw,bufrcvy_up,bufrcvy_dw
	real(KIND(0.d0)), dimension(nxg   ,nyg   ,margin,mtype) :: bufsndz_up,bufsndz_dw,bufrcvz_up,bufrcvz_dw

	integer :: ireqsndx_up,ireqsndx_dw,ireqrcvx_up,ireqrcvx_dw
	integer :: ireqsndy_up,ireqsndy_dw,ireqrcvy_up,ireqrcvy_dw
	integer :: ireqsndz_up,ireqsndz_dw,ireqrcvz_up,ireqrcvz_dw
	integer :: mstatus(mpi_status_size)

	! x direction
	mmx = margin*nyg*nzg*mtype
	if(ib /= ix0 - 1) then
		bufsndx_up = qq(nxg-2*margin+1:nxg-margin,:,:,:)
		nerank = ijkb2rank(ib+1,jb,kb)
		call mpi_isend(bufsndx_up,mmx,MPI_REAL8,nerank,10,mpi_comm_world,ireqsndx_up,merr)
		call mpi_irecv(bufrcvx_up,mmx,MPI_REAL8,nerank,20,mpi_comm_world,ireqrcvx_up,merr)
	endif

	if(ib /= 0)	 then
		bufsndx_dw = qq(margin+1:2*margin,:,:,:)
		nerank = ijkb2rank(ib-1,jb,kb)
		call mpi_isend(bufsndx_dw,mmx,MPI_REAL8,nerank,20,mpi_comm_world,ireqsndx_dw,merr)
		call mpi_irecv(bufrcvx_dw,mmx,MPI_REAL8,nerank,10,mpi_comm_world,ireqrcvx_dw,merr)
	endif

	if(ib /= ix0 - 1) then
		call mpi_wait(ireqsndx_up,mstatus,merr)		
		call mpi_wait(ireqrcvx_up,mstatus,merr)
		qq(nxg-margin+1:nxg,:,:,:) = bufrcvx_up
	endif

	if(ib /= 0) then
		call mpi_wait(ireqsndx_dw,mstatus,merr)		
		call mpi_wait(ireqrcvx_dw,mstatus,merr)
		qq(1:margin,:,:,:) = bufrcvx_dw
	endif

	!---------------------------------------------------------------------------------------
	! y direction
	mmx = nxg*margin*nzg*mtype	
	if(jb /= jx0 - 1) then
		bufsndy_up = qq(:,nyg-2*margin+1:nyg-margin,:,:)
		nerank = ijkb2rank(ib,jb+1,kb)
		call mpi_isend(bufsndy_up,mmx,MPI_REAL8,nerank,30,mpi_comm_world,ireqsndy_up,merr)
		call mpi_irecv(bufrcvy_up,mmx,MPI_REAL8,nerank,40,mpi_comm_world,ireqrcvy_up,merr)
	endif

	if(jb /= 0)	 then
		bufsndy_dw = qq(:,margin+1:2*margin,:,:)
		nerank = ijkb2rank(ib,jb-1,kb)
		call mpi_isend(bufsndy_dw,mmx,MPI_REAL8,nerank,40,mpi_comm_world,ireqsndy_dw,merr)
		call mpi_irecv(bufrcvy_dw,mmx,MPI_REAL8,nerank,30,mpi_comm_world,ireqrcvy_dw,merr)
	endif

	if(jb /= jx0 - 1) then
		call mpi_wait(ireqsndy_up,mstatus,merr)
		call mpi_wait(ireqrcvy_up,mstatus,merr)
		qq(:,nyg-margin+1:nyg,:,:) = bufrcvy_up
	endif

	if(jb /= 0) then
		call mpi_wait(ireqsndy_dw,mstatus,merr)		
		call mpi_wait(ireqrcvy_dw,mstatus,merr)
		qq(:,1:margin,:,:) = bufrcvy_dw
	endif

	!---------------------------------------------------------------------------------------
	! z direction
	mmx = nxg*nyg*margin*mtype
	if(kb /= kx0 - 1) then
		bufsndz_up = qq(:,:,nzg-2*margin+1:nzg-margin,:)
		nerank = ijkb2rank(ib,jb,kb+1)
		call mpi_isend(bufsndz_up,mmx,MPI_REAL8,nerank,50,mpi_comm_world,ireqsndz_up,merr)
		call mpi_irecv(bufrcvz_up,mmx,MPI_REAL8,nerank,60,mpi_comm_world,ireqrcvz_up,merr)
	endif

	if(kb /= 0)	 then
		bufsndz_dw = qq(:,:,margin+1:2*margin,:)
		nerank = ijkb2rank(ib,jb,kb-1)
		call mpi_isend(bufsndz_dw,mmx,MPI_REAL8,nerank,60,mpi_comm_world,ireqsndz_dw,merr)
		call mpi_irecv(bufrcvz_dw,mmx,MPI_REAL8,nerank,50,mpi_comm_world,ireqrcvz_dw,merr)
	endif

	if(kb /= kx0 - 1) then
		call mpi_wait(ireqsndz_up,mstatus,merr)
		call mpi_wait(ireqrcvz_up,mstatus,merr)
		qq(:,:,nzg-margin+1:nzg,:) = bufrcvz_up
	endif

	if(kb /= 0) then
		call mpi_wait(ireqsndz_dw,mstatus,merr)
		call mpi_wait(ireqrcvz_dw,mstatus,merr)
		qq(:,:,1:margin,:) = bufrcvz_dw
	endif

return
end subroutine bc_mpi