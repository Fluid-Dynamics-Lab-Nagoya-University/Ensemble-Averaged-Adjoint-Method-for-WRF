program rand_pert_gen
  implicit none

  integer i,j,k
  integer e_we,e_sn,e_vert
  real amp
  real, allocatable :: q(:,:,:)

  read(*,*) e_we, e_sn, e_vert
  read(*,*) amp
  print*, e_we, e_sn, e_vert, amp

  allocate(q   (e_we,e_sn,e_vert))

  call random_number(q)

  print*, shape(q)
  open(1,file='aux_field_tl',status='unknown',form='formatted')

  write(1,*) e_we, e_sn, e_vert
  write(1,*) (q-0.5)*amp
  close(1)

  j=1
  k=1
  do i=1,e_we
    print*,(q(i,j,k)-0.5)*amp
  enddo
  deallocate(q)
end
