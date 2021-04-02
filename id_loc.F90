subroutine id_loc(x_orgl,x_upgd,nxg_orgl,nxg_upgd,imin,imax,iloc,dx0,dx1,margin)
  implicit none
  
  integer :: i,ii
  integer, intent(in) :: nxg_orgl, nxg_upgd, margin
  real(KIND(0.d0)), dimension(nxg_orgl), intent(in) :: x_orgl
  real(KIND(0.d0)), dimension(nxg_upgd), intent(in) :: x_upgd

  integer, intent(out) :: imin,imax
  integer, dimension(nxg_upgd), intent(out) :: iloc
  real(8), dimension(nxg_upgd), intent(out) :: dx0,dx1

  ! serch minmum
  do i = 1+margin,nxg_upgd-margin
     if(x_upgd(i) > x_orgl(1)) then
        imin = i
        goto 2000
     endif
  enddo
  2000 continue

  ! serch maximum
  do i = nxg_upgd-margin,1+margin,-1
     if(x_upgd(i) < x_orgl(nxg_orgl)) then
        imax = i
        goto 3000
     endif
  enddo
3000 continue

  do i = 1+margin,nxg_upgd-margin
     do ii = nxg_orgl-1,1,-1
        if(x_upgd(i) >= x_orgl(ii)) then
           iloc(i) = ii
           dx0(i) = x_upgd(i   ) - x_orgl(ii)
           dx1(i) = x_orgl(ii+1) - x_upgd(i )
           goto 1000
        endif
     enddo
1000 continue
  enddo
  
  return
end subroutine id_loc
