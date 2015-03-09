program main
  use omp_lib
  implicit none

  double precision, parameter :: pi = 3.141592653589793d+0
  double precision, parameter :: ydistance = 10 !meters Distance from
                                                 !viewpoint to view
                                                 !window
  double precision, parameter :: window_max = 10 !meters
  double precision, parameter :: radius = 6 !meters Size of sphere to
                                                 !view
  
  double precision, parameter, dimension(3) :: clocation =  (/0, 12, 0/)  !meters
  double precision, parameter, dimension(3) :: llocation =  (/4, 4, -1/)  !meters

  
  double precision, dimension(:,:), allocatable :: G !!  View Window
  integer Nrays !! Number of rays to generate
  double precision, dimension(3) :: V=0 !!View Array
  double precision, dimension(3) :: W=0 !!Window Array 
  double precision, dimension(3) :: C=0 !!Sphere Location Array
  double precision, dimension(3) :: Isec=0 !!Intersection of view array and
                                                 !sphere
  double precision, dimension(3) :: L=0 !!Light source array
  double precision, dimension(3) :: N=0 !!Unit norml vecotr at I
  double precision, dimension(3) :: S=0 !!Direction of light source at I
  double precision b !!brightness at specific point
  integer n_grid !!nXn grid size for G
  double precision theta !! polar angle
  double precision phi !! azimuthal angle
  double precision t, t1, t2 !scalar
  double precision start, finish
  real*8 rando !random number
  integer i,j,k
  character(10) prog_name, arg1_str, arg2_str 
 
  !Get command line arguments.
  call getarg(0, prog_name)
  call getarg(1, arg1_str)
  call getarg(2, arg2_str)

  read (arg1_str, *) Nrays
  read (arg2_str, *) n_grid
  
  allocate(G(n_grid,n_grid))
  G=0
  
  !set y distance for window
  W(2) = ydistance
  C=clocation
  L=llocation
  start = omp_get_wtime();
!$OMP PARALLEL DO &
!$OMP PRIVATE(theta, phi, V, t, t1, t2, Isec, N, S, i, j, b) &
!$OMP FIRSTPRIVATE(W) &
!$OMP REDUCTION(+:G)

  do k=1, Nrays
    call set_random_angle(theta)
    call set_random_angle(phi)

    call set_cart(V, theta, phi)
    call inter_view_wind(W,V)


    !multiplier = W(2)/V(2)
    !write(*,*) W(2)
    !W = multiplier*V
    !write(*,*) W(1), W(2), W(3)
    
    if(abs(W(1)) .lt. window_max .and. abs(W(3)) .lt. window_max) then
      t1 = calc_scalar_1(V,C)
      t2 = calc_scalar_2_sqr(t1, C, radius)
      !write(*,*) t1, t2
      if(t2 .ge. 0) then
        t=t1-sqrt(t2)
        Isec = t*V
        call normal_vector(Isec, C, N)
        call normal_vector(L, Isec, S)
        
        !write(*,*) S
        b = brightness(S,N)
        call get_i_j(W, window_max, n_grid, i,j)
        !write(*,*)i,j
        G(i,j) = G(i,j) + b

        !write(*,*) i, ':', j, '=', G(i,j) 
      end if
    end if
  end do

!$OMP END PARALLEL DO
  finish = omp_get_wtime();
  write(*,*) finish-start
  call write_G(G, n_grid, Nrays)

  deallocate(G)

  contains

  subroutine set_cart(field, theta, phi)
    double precision, intent(inout) :: field(*)
    double precision, intent(in) :: theta
    double precision, intent(in) :: phi

    field(1) = sin(theta) * cos(phi)
    field(2) = sin(theta) * sin(phi)
    field(3) = cos(theta)
  end subroutine set_cart

  subroutine set_random_angle(angle)
    double precision, intent(inout) :: angle
    double precision rando

    call random_number(rando)
    angle = rando * 2.*pi   
  end subroutine set_random_angle

  subroutine inter_view_wind(window,view)
    double precision, intent(inout) :: window(3)
    double precision, intent(in) :: view(3)
   
    double precision multiplier

    multiplier = window(2)/view(2)
    window = multiplier*view
  end subroutine inter_view_wind

  double precision function calc_scalar_1(view, sphere)
    double precision, intent(in) :: view(3)
    double precision, intent(in) :: sphere(3)

    double precision v_dot_s
    
    v_dot_s = dot_product(view, sphere)
    calc_scalar_1 = v_dot_s
  end function calc_scalar_1 
  
  double precision function calc_scalar_2_sqr(scalar_1, sphere, radius)
    double precision, intent(in) :: scalar_1
    double precision, intent(in) :: sphere(3)
    double precision, intent(in) :: radius

    double precision s_dot_s
    double precision scalar_1_sqr

    s_dot_s = dot_product(sphere, sphere)
    scalar_1_sqr = scalar_1**2
    calc_scalar_2_sqr = scalar_1_sqr + radius**2 - s_dot_s
 
  end function calc_scalar_2_sqr
  
  subroutine normal_vector(array1, array2, result_array)
    double precision, intent(in) :: array1(3)
    double precision, intent(in) :: array2(3)
    double precision, intent(inout) :: result_array(3)
    
    double precision  difference(3)
    
    difference = array1 - array2
    result_array = difference / abs(difference)

  end subroutine normal_vector
  
  double precision function brightness(S, N)
    double precision, intent(in) :: S(3)
    double precision, intent(in) :: N(3)

    double precision s_dot_n
    
    s_dot_n = dot_product(S,N)

    if(s_dot_n .lt. 0) then
      brightness=0
    else
      brightness = s_dot_n
    end if
  end function brightness
 
  subroutine get_i_j(window, window_max, grid_size, i, j)
    double precision, intent(in) :: window(*)
    double precision, intent(in) :: window_max
    integer, intent(in) :: grid_size
    integer, intent(out) :: i, j

    double precision step
    double precision half_grid

    step = (2*window_max)/grid_size
    half_grid = grid_size/2

    if(window(1)<0) then
      i = 1+half_grid+int(window(1)/step)
    else
      i = half_grid+int(window(1)/step)
    end if

    if(window(3)<0) then
      j = 1+half_grid+int(window(3)/step)
    else
      j = half_grid+int(window(3)/step)
    end if 
  end subroutine get_i_j

  subroutine write_G(G_array, n, num_rays)
    integer, intent(in) :: n
    double precision, dimension(n,n), intent(in) :: G_array
    integer, intent(in) :: num_rays

    character(20) :: outfile, rays, n_size
    integer i
    write(rays, *) num_rays
    write(n_size,*) n

    outfile = "G_"//trim(adjustl(rays))//"_"//trim(adjustl(n_size))//"X"//trim(adjustl(n_size))//".out"
    open(99, file=outfile)

    do i=1, n
      write(99,*) G_array(i,:)
    end do
    
    close(99)
  end subroutine write_G
  end program main
