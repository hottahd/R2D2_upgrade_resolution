subroutine specify_location(ix0,ix00_orgl,ix00_upgd,margin_orgl,margin_upgd,nx,x00_orgl,x00_upgd,xstt,xend,xwdt)
  implicit none

	integer, intent(in) :: ix0, nx
	integer, intent(in) :: ix00_orgl,ix00_upgd
	integer, intent(in) :: margin_orgl,margin_upgd
	real(KIND(0.d0)), dimension(ix00_orgl), intent(in) :: x00_orgl
	real(KIND(0.d0)), dimension(ix00_upgd), intent(in) :: x00_upgd
	integer, dimension(0:ix0-1),  intent(out) :: xstt,xend,xwdt
	real(KIND(0.d0)), parameter :: rstar = 6.9598947d+10

  integer :: i,ibt
	real(KIND(0.d0)) :: xbot,xtop

	do ibt = 0,ix0 - 1
    xbot = 0.5d0*(x00_orgl(margin_orgl + nx* ibt   ) + x00_orgl(margin_orgl + nx* ibt + 1))
    xtop = 0.5d0*(x00_orgl(margin_orgl + nx*(ibt+1)) + x00_orgl(margin_orgl + nx*(ibt+1) + 1))
		
    if(ibt /= 0) then
    	do i = 1,ix00_upgd
        if(x00_upgd(i) > xbot) then
          xstt(ibt) = i - margin_upgd
          goto 1000
        endif
      enddo
    endif
1000 continue

    if(ibt /= ix0 - 1) then
      do i = ix00_upgd,1,-1
        if(x00_upgd(i) <= xtop) then
          xend(ibt) = i + margin_upgd
          goto 2000
        endif
      enddo
    endif
2000 continue

		xwdt(ibt) = xend(ibt) - xstt(ibt) + 1
	enddo

return
end subroutine specify_location