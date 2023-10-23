module particle_snapshot
  use amr_parameters, only: dp, i8b
  use pm_commons, only: dp

  implicit none

  integer(i8b), allocatable, dimension(:) :: particle_ids
  integer :: npart_subset
  integer :: selected_particle_type = -999  ! Dummy value to prevent any error

contains

  pure function binary_search(id, ids_list, np)
    ! Uses a binary search to find `id` in `ids_list`, assuming the latter is sorted in increasing order
    integer(i8b), intent(in) :: id
    integer(i8b), dimension(np), intent(in) :: ids_list
    integer, intent(in) :: np

    integer :: binary_search

    integer(i8b) :: v
    integer :: L, M, R
    L = 1; R = np

    ! Early break if the id is not within min/max
    if ((id < ids_list(L)) .or. (id > ids_list(R))) then
      binary_search = 0
      return
    end if

    do while (L <= R)
       M = (L + R) / 2
       v = ids_list(M)
       if (v < id) then
          L = M + 1
       else if (v > id) then
          R = M - 1
       else
          binary_search = M
          return
       end if
    end do

    binary_search = 0

 end function binary_search

 recursive subroutine quicksort(a, first, last)
    integer(i8b), intent(inout), dimension(:) :: a
    integer(i8b) :: x, t
    integer, intent(in) :: first, last
    integer :: i, j

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
       do while (a(i) < x)
          i=i+1
       end do
       do while (x < a(j))
          j=j-1
       end do
       if (i >= j) exit
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    end do
    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
 end subroutine quicksort

  subroutine read_particle_ids(particle_source_file, n, ids, particle_type)
    character(len=*), intent(in) :: particle_source_file
    integer, intent(out) :: n
    integer(i8b), intent(out), allocatable, dimension(:) :: ids
    integer, intent(out) :: particle_type

    integer :: read_unit, ios
    integer(i8b) :: idpart

    open(newunit=read_unit, file=trim(particle_source_file), form='formatted', iostat=ios, status="old")
    if (ios /= 0) then
      print*, "ERROR: Could not read file `", trim(particle_source_file), "`"
      call clean_stop
    end if

    read(read_unit, *, iostat=ios) particle_type
    if (particle_type < -5 .or. particle_type > 5) then
      print*, "ERROR: the first line of", trim(particle_source_file), "should be a valid particle family. Got", particle_type
      call clean_stop
    end if
    n = 0
    do
      read(read_unit, *, iostat=ios) idpart
      if (ios /= 0) exit
      n = n + 1
    end do

    allocate(ids(1:n))
    rewind(read_unit)
    n = 0
    read(read_unit, *) particle_type ! Pass first line
    do
      read(read_unit, *, iostat=ios) idpart
      if (ios /= 0) exit
      n = n + 1
      ids(n) = idpart
    end do

    close(read_unit)

    call quicksort(ids, 1, n)

  end subroutine read_particle_ids


end module particle_snapshot
