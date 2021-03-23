module change
  private
  interface change_judge
    module procedure change_judge_d, change_judge_i,change_judge_l
  end interface
  public :: change_judge

contains
 
!--------------------------------------------------------
 subroutine change_judge_d(a,b,change)
		implicit none

		real(KIND(0.d0)), intent(in) :: a,b
		character*11, intent(inout) :: change

		if (a == b) then
			change = ' (unchange)'
		else
			change = ' (change)  '
		endif

		return
	end subroutine change_judge_d
!--------------------------------------------------------
 subroutine change_judge_i(a,b,change)
		implicit none

		integer, intent(in) :: a,b
		character*11, intent(inout) :: change

		if (a == b) then
			change = ' (unchange)'
		else
			change = ' (change)  '
		endif

		return
	end subroutine change_judge_i

!--------------------------------------------------------
 subroutine change_judge_l(a,b,change)
		implicit none

		logical, intent(in) :: a,b
		character*11, intent(inout) :: change

		if (a .eqv. b) then
			change = ' (unchange)'
		else
			change = ' (change)  '
		endif

		return
	end subroutine change_judge_l

end module change