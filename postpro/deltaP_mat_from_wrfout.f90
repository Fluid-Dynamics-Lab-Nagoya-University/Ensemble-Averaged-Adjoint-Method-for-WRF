program deltaP_mat_from_wrfout
  use netcdf
  implicit none

  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: NCASE = 100
  integer, parameter :: T_INDEX_1B = 7   ! 1-based (Python t_idx=6)

  character(len=1024) :: NL_ROOT, REF_ROOT, OUT_DIR
  character(len=4096) :: SG_LIST, IG_LIST
  character(len=64)   :: NACTR_FIELD
  integer :: cx0, cy0  ! 0-based in input
  integer :: cx,  cy   ! 1-based internal (we,sn)

  integer :: nsg, nig
  character(len=32), allocatable :: sg_vals(:), ig_vals(:)

  real(sp), allocatable :: dP_center(:,:), dP_g3(:,:)

  character(len=1024) :: out_center, out_g3

  ! defaults for center (same as your Python)
  cx0 = 34
  cy0 = 16

  call getenv_req('NL_ROOT',  NL_ROOT)
  call getenv_req('REF_ROOT', REF_ROOT)
  call getenv_req('OUT_DIR',  OUT_DIR)
  call getenv_req('SG_LIST',  SG_LIST)
  call getenv_req('IG_LIST',  IG_LIST)
  call getenv_allow_empty('NACTR_FIELD', NACTR_FIELD)
  call getenv_int_default('CENTER_X', cx0)
  call getenv_int_default('CENTER_Y', cy0)

  cx = cx0 + 1
  cy = cy0 + 1

  call split_list(SG_LIST, sg_vals, nsg)
  call split_list(IG_LIST, ig_vals, nig)

  if (nsg <= 0 .or. nig <= 0) then
    write(*,*) '[ERR] empty SG_LIST or IG_LIST.'
    stop 2
  end if

  allocate(dP_center(nsg, nig))
  allocate(dP_g3(nsg, nig))
  dP_center = 0.0_sp
  dP_g3     = 0.0_sp

  call compute_all(nsg, nig, sg_vals, ig_vals, cx, cy, NACTR_FIELD, dP_center, dP_g3)

  out_center = trim(OUT_DIR)//'/deltaP_mat_absG_i0M_NOV0.02_ngEQsg'//trim(NACTR_FIELD)//'_sg.npy'
  out_g3     = trim(OUT_DIR)//'/deltaP_mat_absG_i0M_NOV0.02_ngEQsg'//trim(NACTR_FIELD)//'_sg_G3x3.npy'

  call save_npy_f32_2d(out_center, dP_center, .true.)
  call save_npy_f32_2d(out_g3,     dP_g3,     .true.)

  write(*,*) '[OK] wrote:'
  write(*,*) '  ', trim(out_center)
  write(*,*) '  ', trim(out_g3)

contains

  subroutine compute_all(nsg, nig, sg_vals, ig_vals, cx, cy, nactr_field, dP_center, dP_g3)
    integer, intent(in) :: nsg, nig, cx, cy
    character(len=32), intent(in) :: sg_vals(nsg), ig_vals(nig)
    character(len=*),  intent(in) :: nactr_field
    real(sp), intent(inout) :: dP_center(nsg, nig), dP_g3(nsg, nig)

    integer :: isg, iig, nl_ncase, ref_ncase
    character(len=1024) :: ref_dir, nl_dir
    real(dp) :: pref_c, pctrl_c, pref_g, pctrl_g
    real(dp) :: rv_c, rv_g

    do isg = 1, nsg
      ref_dir = trim(REF_ROOT)//'/woinput_absG_i0M_NOV0.02_sg_'//trim(sg_vals(isg))
      ! When sg=0, all ensemble members are identical; only member 1 exists.
      if (trim(sg_vals(isg)) == '0') then
        nl_ncase  = 1
        ref_ncase = 1
      else
        nl_ncase  = NCASE
        ref_ncase = NCASE
      end if
      do iig = 1, nig
        nl_dir  = trim(NL_ROOT)//'/absG_i0M_NOV0.02_ng'//trim(sg_vals(isg))// &
                  '_ig'//trim(ig_vals(iig))//trim(nactr_field)//'_sg'//trim(sg_vals(isg))

        pref_c  = mean_weighted_over_members(ref_dir, 'AD', cx, cy, 0, ref_ncase)
        pctrl_c = mean_weighted_over_members(nl_dir,  'NL', cx, cy, 0, nl_ncase)
        rv_c    = pref_c - pctrl_c

        pref_g  = mean_weighted_over_members(ref_dir, 'AD', cx, cy, 1, ref_ncase)
        pctrl_g = mean_weighted_over_members(nl_dir,  'NL', cx, cy, 1, nl_ncase)
        rv_g    = pref_g - pctrl_g

        dP_center(isg, iig) = real(rv_c, sp)
        dP_g3(isg, iig)     = real(rv_g, sp)

        write(*,'(A,1X,A,1X,A,1X,A,F10.4,1X,A,F10.4)') &
          '[DONE]', 'sg='//trim(sg_vals(isg)), 'ig='//trim(ig_vals(iig)), &
          'RV_center=', real(rv_c,sp), 'RV_G3x3=', real(rv_g,sp)
      end do
    end do
  end subroutine compute_all

  function mean_weighted_over_members(dir, tag, cx, cy, kernel_mode, ncase) result(meanv)
    character(len=*), intent(in) :: dir, tag
    integer, intent(in) :: cx, cy, kernel_mode, ncase
    real(dp) :: meanv
    integer :: i
    character(len=1024) :: fpath
    real(dp) :: sumv

    sumv = 0.0_dp
    !$omp parallel do private(i,fpath) reduction(+:sumv) schedule(dynamic,1)
    do i = 1, ncase
      fpath = trim(dir)//'/wrfout_d01_2018-07-05_120000_woinput_'//trim(tag)//trim(itoa(i))
      sumv  = sumv + weighted_value_one_file(fpath, cx, cy, kernel_mode)
    end do
    !$omp end parallel do
    meanv = sumv / real(ncase, dp)

  end function mean_weighted_over_members

  function weighted_value_one_file(fpath, cx, cy, kernel_mode) result(val)
    character(len=*), intent(in) :: fpath
    integer, intent(in) :: cx, cy, kernel_mode
    real(dp) :: val

    logical :: ex
    integer :: ncid, varid, e
    integer :: nd, dimids(3)
    integer :: time_pos, we_pos, sn_pos
    character(len=NF90_MAX_NAME) :: dname
    integer :: dlen
    integer :: we_len, sn_len, t_len

    integer :: start3(3), count3(3)
    real(sp) :: raw(3,3)         ! as read (order depends on dim order)
    real(dp) :: win(3,3)         ! mapped to (we,sn) order

    inquire(file=trim(fpath), exist=ex)
    if (.not. ex) then
      write(*,*) '[ERR] missing file: ', trim(fpath)
      stop 10
    end if

    e = nf90_open(trim(fpath), NF90_NOWRITE, ncid); call chk(e, 'open '//trim(fpath))
    e = nf90_inq_varid(ncid, 'RAINNC', varid);       call chk(e, 'inq_varid RAINNC')
    e = nf90_inquire_variable(ncid, varid, ndims=nd, dimids=dimids); call chk(e, 'inq_var RAINNC')
    if (nd /= 3) then
      write(*,*) '[ERR] RAINNC ndims != 3 in ', trim(fpath)
      stop 11
    end if

    call find_dim_positions(ncid, dimids, time_pos, we_pos, sn_pos, t_len, we_len, sn_len)

    if (T_INDEX_1B < 1 .or. T_INDEX_1B > t_len) then
      write(*,*) '[ERR] Time index out of range. T=', T_INDEX_1B, ' TimeLen=', t_len
      stop 12
    end if

    ! boundary check for 3x3 when kernel_mode==1
    
    ! IMPORTANT: follow Python indexing:
    ! center_x -> first spatial axis, center_y -> second spatial axis
    if (kernel_mode == 1) then
      if (cx < 2) then
        write(*,*) '[ERR] cx too small for 3x3'
        stop 13
      end if
      if (cy < 2) then
        write(*,*) '[ERR] cy too small for 3x3'
        stop 13
      end if
      if (cx > sn_len-1) then
        write(*,*) '[ERR] cx too large for 3x3'
        stop 13
      end if
      if (cy > we_len-1) then
        write(*,*) '[ERR] cy too large for 3x3'
        stop 13
      end if
    else
      if (cx < 1 .or. cx > sn_len) then
        write(*,*) '[ERR] cx out of bounds'
        stop 14
      end if
      if (cy < 1 .or. cy > we_len) then
        write(*,*) '[ERR] cy out of bounds'
        stop 14
      end if
    end if


    start3 = 1
    count3 = 1
    start3(time_pos) = T_INDEX_1B
    count3(time_pos) = 1
    start3(sn_pos) = cx - 1
    count3(sn_pos) = 3
    start3(we_pos) = cy - 1
    count3(we_pos) = 3

    e = nf90_get_var(ncid, varid, raw, start=start3, count=count3); call chk(e, 'get_var RAINNC 3x3')
    call chk(nf90_close(ncid), 'close')

    ! Map raw(3,3) into win(we,sn) = win(i_we, j_sn)
    call map_raw_to_win(raw, time_pos, we_pos, sn_pos, win)

    if (kernel_mode == 0) then
      val = win(2,2)
    else
      val = sum(win * kernel_g3()) / 4.0_dp
    end if
  end function weighted_value_one_file

  subroutine find_dim_positions(ncid, dimids, time_pos, we_pos, sn_pos, t_len, we_len, sn_len)
    integer, intent(in) :: ncid
    integer, intent(in) :: dimids(3)
    integer, intent(out) :: time_pos, we_pos, sn_pos
    integer, intent(out) :: t_len, we_len, sn_len
    character(len=NF90_MAX_NAME) :: nm
    integer :: lenv
    integer :: p

    time_pos = 0; we_pos = 0; sn_pos = 0
    t_len = -1; we_len = -1; sn_len = -1

    do p = 1, 3
      call dim_name_len(ncid, dimids(p), nm, lenv)
      if (trim(nm) == 'Time') then
        time_pos = p
        t_len = lenv
      else if (trim(nm) == 'west_east' .or. trim(nm) == 'west_east_stag') then
        we_pos = p
        we_len = lenv
      else if (trim(nm) == 'south_north' .or. trim(nm) == 'south_north_stag') then
        sn_pos = p
        sn_len = lenv
      end if
    end do

    if (time_pos == 0 .or. we_pos == 0 .or. sn_pos == 0) then
      write(*,*) '[ERR] cannot identify dims (Time/we/sn).'
      write(*,*) '       pos: time,we,sn=', time_pos, we_pos, sn_pos
      stop 15
    end if
  end subroutine find_dim_positions

  subroutine map_raw_to_win(raw, time_pos, we_pos, sn_pos, win)
    real(sp), intent(in) :: raw(3,3)
    integer, intent(in) :: time_pos, we_pos, sn_pos
    real(dp), intent(out) :: win(3,3)
    ! raw is 2D container for the two non-time dims in their order.
    ! We need win(i_we, j_sn) in (we,sn) order.
    !
    ! There are only two possible non-time orders: (we,sn) or (sn,we).
    ! If we_pos < sn_pos among the 3 dims and time_pos is the remaining one, nf90_get_var will
    ! fill raw with first index = dim with smaller position among non-time dims.
    ! The simplest robust approach: detect whether the first non-time dim is we or sn.

    integer :: first_non_time_pos, second_non_time_pos
    logical :: first_is_we
    integer :: i, j

    ! Determine first and second non-time positions
    if (time_pos == 1) then
      first_non_time_pos = 2
      second_non_time_pos = 3
    else if (time_pos == 2) then
      first_non_time_pos = 1
      second_non_time_pos = 3
    else
      first_non_time_pos = 1
      second_non_time_pos = 2
    end if

    first_is_we = (first_non_time_pos == we_pos)

    if (first_is_we) then
      ! raw(i,j) = raw(we,sn)
      do j = 1, 3
        do i = 1, 3
          win(i,j) = real(raw(i,j), dp)
        end do
      end do
    else
      ! raw(i,j) = raw(sn,we) -> transpose into (we,sn)
      do j = 1, 3
        do i = 1, 3
          win(i,j) = real(raw(j,i), dp)
        end do
      end do
    end if
  end subroutine map_raw_to_win

  pure function kernel_g3() result(k)
    real(dp) :: k(3,3)
    k(1,1)=0.25_dp; k(2,1)=0.5_dp;  k(3,1)=0.25_dp
    k(1,2)=0.5_dp;  k(2,2)=1.0_dp;  k(3,2)=0.5_dp
    k(1,3)=0.25_dp; k(2,3)=0.5_dp;  k(3,3)=0.25_dp
  end function kernel_g3

  subroutine dim_name_len(ncid, dimid, name, lenout)
    integer, intent(in) :: ncid, dimid
    character(len=*), intent(out) :: name
    integer, intent(out) :: lenout
    integer :: e
    e = nf90_inquire_dimension(ncid, dimid, name=name, len=lenout)
    call chk(e, 'inq_dim')
  end subroutine dim_name_len

  subroutine chk(e, ctx)
    integer, intent(in) :: e
    character(len=*), intent(in) :: ctx
    if (e /= nf90_noerr) then
      write(*,*) '[NETCDF ERR] ', trim(ctx), ': ', trim(nf90_strerror(e))
      stop 99
    end if
  end subroutine chk

  subroutine getenv_req(key, val)
    character(len=*), intent(in) :: key
    character(len=*), intent(out):: val
    integer :: st
    call get_environment_variable(key, val, status=st, trim_name=.true.)
    if (st /= 0 .or. len_trim(val) == 0) then
      write(*,*) '[ERR] missing env: ', trim(key)
      stop 3
    end if
  end subroutine getenv_req

  subroutine getenv_allow_empty(key, val)
    character(len=*), intent(in) :: key
    character(len=*), intent(out):: val
    integer :: st
    call get_environment_variable(key, val, status=st, trim_name=.true.)
    if (st /= 0) val = ''
  end subroutine getenv_allow_empty

  subroutine getenv_int_default(key, x)
    character(len=*), intent(in) :: key
    integer, intent(inout) :: x
    character(len=64) :: s
    integer :: st, iosr, tmp
    call get_environment_variable(key, s, status=st, trim_name=.true.)
    if (st == 0 .and. len_trim(s) > 0) then
      read(s, *, iostat=iosr) tmp
      if (iosr == 0) x = tmp
    end if
  end subroutine getenv_int_default

  pure function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=32) :: s
    write(s,'(I0)') i
  end function itoa

  subroutine split_list(str, items, n)
    character(len=*), intent(in) :: str
    character(len=32), allocatable, intent(out) :: items(:)
    integer, intent(out) :: n
    integer :: i, L, cnt, p, q
    character(len=4096) :: s

    s = trim(str)
    L = len_trim(s)

    cnt = 0
    i = 1
    do while (i <= L)
      do while (i <= L .and. s(i:i) == ' ')
        i = i + 1
      end do
      if (i > L) exit
      cnt = cnt + 1
      do while (i <= L .and. s(i:i) /= ' ')
        i = i + 1
      end do
    end do

    n = cnt
    allocate(items(n))

    i = 1
    cnt = 0
    do while (i <= L)
      do while (i <= L .and. s(i:i) == ' ')
        i = i + 1
      end do
      if (i > L) exit
      p = i
      do while (i <= L .and. s(i:i) /= ' ')
        i = i + 1
      end do
      q = i - 1
      cnt = cnt + 1
      items(cnt) = adjustl(s(p:q))
    end do
  end subroutine split_list

  ! -----------------------------
  ! NPY writer: float32, 2D, v1.0
  ! -----------------------------
  subroutine save_npy_f32_2d(path, A, fortran_order)
    character(len=*), intent(in) :: path
    real(sp), intent(in) :: A(:,:)
    logical, intent(in) :: fortran_order

    integer :: u, iosw
    integer :: n1, n2
    character(len=6) :: magic
    integer(kind=1) :: vmaj, vmin
    integer(kind=2) :: hlen
    character(len=512) :: header
    character(len=5) :: fflag
    integer :: hlen_i, pad, i

    n1 = size(A,1)
    n2 = size(A,2)

    if (fortran_order) then
      fflag = 'True'
    else
      fflag = 'False'
    end if

    write(header,'(A,I0,A,I0,A)') "{'descr': '<f4', 'fortran_order': "//fflag//", 'shape': (", n1, ', ', n2, '), }'
    hlen_i = len_trim(header)

    ! pad so that (10 + header_len) is multiple of 16, and header ends with '\n'
    pad = mod(16 - mod(10 + hlen_i + 1, 16), 16)
    do i = hlen_i+1, hlen_i+pad
      header(i:i) = ' '
    end do
    header(hlen_i+pad+1:hlen_i+pad+1) = new_line('a')
    hlen_i = hlen_i + pad + 1
    hlen = int(hlen_i, kind=2)

    magic = char(147)//'NUMPY'
    vmaj = 1
    vmin = 0

    open(newunit=u, file=trim(path), access='stream', form='unformatted', &
         status='replace', action='write', iostat=iosw)
    if (iosw /= 0) then
      write(*,*) '[ERR] cannot open for write: ', trim(path)
      stop 20
    end if

    write(u) magic
    write(u) vmaj, vmin
    write(u) hlen
    write(u) header(1:hlen_i)
    write(u) A

    close(u)
  end subroutine save_npy_f32_2d

end program deltaP_mat_from_wrfout

