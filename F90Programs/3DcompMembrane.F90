program composite_3D
  implicit none

  integer :: n, i, info
  real(8), allocatable :: E1(:), E2(:), E3(:)
  real(8), allocatable :: G12(:), G13(:), G23(:)
  real(8), allocatable :: nu12(:), nu13(:), nu23(:)
  real(8), allocatable :: theta(:), t(:), V(:)

  real(8) :: A(3,3), Ainv(3,3)
  real(8) :: Q11, Q22, Q12, Q66, nu21, denom
  real(8) :: Qb11, Qb22, Qb12, Qb66, Qb16, Qb26
  real(8) :: c, s, h

  ! Final properties
  real(8) :: Ex, Ey, Ez
  real(8) :: Gxy, Gxz, Gyz
  real(8) :: nuxy, nuxz, nuyz
  real(8) :: Ez_inv, Gxz_inv, Gyz_inv

OPEN(UNIT=20, FILE='LayerData.txt', STATUS='OLD', ACTION='READ')
Read(20,*) n
     print *, "Layers: ", n

  allocate(E1(n),E2(n),E3(n),G12(n),G13(n),G23(n))
  allocate(nu12(n),nu13(n),nu23(n),theta(n),t(n),V(n))

  print *, "Enter per layer:"
  print *, "(E1 E2 E3 G12 G13 G23 nu12 nu13 nu23 theta(deg) thickness)"

  do i = 1,n
    print *, "Layer ", i
    read (20,*) E1(i),E2(i),E3(i),G12(i),G13(i),G23(i), &
             nu12(i),nu13(i),nu23(i),theta(i),t(i)
  end do

  ! Thickness
  h = sum(t)
  V = t / h

  A = 0.0d0
  Ez_inv = 0.0d0
  Gxz_inv = 0.0d0
  Gyz_inv = 0.0d0
  nuxz = 0.0d0
  nuyz = 0.0d0

  do i = 1,n

    ! Angle
    c = cos(theta(i)*acos(-1.0d0)/180.0d0)
    s = sin(theta(i)*acos(-1.0d0)/180.0d0)

    ! Plane stress Q
    nu21 = nu12(i)*E2(i)/E1(i)
    denom = 1.0d0 - nu12(i)*nu21

    Q11 = E1(i)/denom
    Q22 = E2(i)/denom
    Q12 = nu12(i)*E2(i)/denom
    Q66 = G12(i)

    ! Transform
    Qb11 = Q11*c**4 + Q22*s**4 + 2d0*(Q12+2d0*Q66)*s**2*c**2
    Qb22 = Q11*s**4 + Q22*c**4 + 2d0*(Q12+2d0*Q66)*s**2*c**2
    Qb12 = (Q11+Q22-4d0*Q66)*s**2*c**2 + Q12*(s**4+c**4)
    Qb66 = (Q11+Q22-2d0*Q12-2d0*Q66)*s**2*c**2 + Q66*(s**4+c**4)

    Qb16 = (Q11-Q12-2d0*Q66)*c**3*s - (Q22-Q12-2d0*Q66)*c*s**3
    Qb26 = (Q11-Q12-2d0*Q66)*c*s**3 - (Q22-Q12-2d0*Q66)*c**3*s

    ! Assemble A
    A(1,1)=A(1,1)+Qb11*t(i)
    A(2,2)=A(2,2)+Qb22*t(i)
    A(1,2)=A(1,2)+Qb12*t(i)
    A(2,1)=A(1,2)

    A(3,3)=A(3,3)+Qb66*t(i)

    A(1,3)=A(1,3)+Qb16*t(i)
    A(3,1)=A(1,3)

    A(2,3)=A(2,3)+Qb26*t(i)
    A(3,2)=A(2,3)

    ! --- Transverse ---
    Ez_inv  = Ez_inv  + V(i)/E3(i)
    Gxz_inv = Gxz_inv + V(i)/G13(i)
    Gyz_inv = Gyz_inv + V(i)/G23(i)

    nuxz = nuxz + V(i)*nu13(i)
    nuyz = nuyz + V(i)*nu23(i)

  end do
A(1,3) = 0.0d0
A(3,1) = 0.0d0
A(2,3) = 0.0d0
A(3,2) = 0.0d0
  call invert3x3(A,Ainv)

  ! In-plane
  Ex   = 1.0d0/(Ainv(1,1)*h)
  Ey   = 1.0d0/(Ainv(2,2)*h)
  Gxy  = 1.0d0/(Ainv(3,3)*h)
  nuxy = A(1,2)/A(2,2)

  ! Transverse
  Ez  = 1.0d0/Ez_inv
  Gxz = 1.0d0/Gxz_inv
  Gyz = 1.0d0/Gyz_inv

  print *, "-----------------------------"
  print *, "3D Equivalent Properties"
  print *, "Ex  =",Ex
  print *, "Ey  =",Ey
  print *, "Ez  =",Ez
  print *, "Gxy =",Gxy
  print *, "Gxz =",Gxz
  print *, "Gyz =",Gyz
  print *, "nu_xy =",nuxy
  print *, "nu_xz =",nuxz
  print *, "nu_yz =",nuyz
  print *, "-----------------------------"

  print *, "-------------------------------------"
OPEN(UNIT=10, FILE='Compliance.txt', status='UNKNOWN')
  write(10,fmt='(A23)') 'Youngs Modulus(6,6) = \'
  write(10,fmt='(6e13.5,A2)') Ex,1/Ey,1/Ez,0.,0.,0.,' \'
  write(10,fmt='(6e13.5,A2)') 1/Ex,Ey,1/Ez,0.,0.,0.,' \'
  write(10,fmt='(6e13.5,A2)') 1/Ex,1/Ey,Ez,0.,0.,0.,' \'
  write(10,fmt='(6e13.5,A2)') 0.,0.,0.,Gxy,0.,0.,' \'
  write(10,fmt='(6e13.5,A2)') 0.,0.,0.,0.,Gxz,0.,' \'
  write(10,fmt='(6e13.5)') 0.,0.,0.,0.,0.,Gyz


contains

  subroutine invert3x3(A,Ainv)
    real(8) :: A(3,3),Ainv(3,3),det

    det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
        - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
        + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

    Ainv(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
    Ainv(1,2)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))/det
    Ainv(1,3)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))/det

    Ainv(2,1)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))/det
    Ainv(2,2)=(A(1,1)*A(3,3)-A(1,3)*A(3,1))/det
    Ainv(2,3)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))/det

    Ainv(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
    Ainv(3,2)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))/det
    Ainv(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
  end subroutine
subroutine invert6(A, Ainv, info)
  implicit none
  integer, parameter :: n = 6
  real(8), intent(in)  :: A(n,n)
  real(8), intent(out) :: Ainv(n,n)
  integer, intent(out) :: info

  real(8) :: aug(n,2*n)
  real(8) :: temp(2*n)
  real(8) :: factor, pivot
  integer :: i, j, k, maxrow
  real(8) :: maxval

  info = 0

  ! Form augmented matrix [A | I]
  aug(:,1:n) = A
  aug(:,n+1:2*n) = 0.0d0
  do i = 1, n
     aug(i,n+i) = 1.0d0
  end do

  ! Gauss-Jordan elimination
  do i = 1, n

     ! Partial pivoting
     maxrow = i
     maxval = abs(aug(i,i))
     do k = i+1, n
        if (abs(aug(k,i)) > maxval) then
           maxval = abs(aug(k,i))
           maxrow = k
        end if
     end do

     ! Check for singular matrix
     if (maxval == 0.0d0) then
        info = 1
        return
     end if

     ! Swap rows if needed
     if (maxrow /= i) then
        temp = aug(i,:)
        aug(i,:) = aug(maxrow,:)
        aug(maxrow,:) = temp
     end if

     ! Normalize pivot row
     pivot = aug(i,i)
     aug(i,:) = aug(i,:) / pivot

     ! Eliminate other rows
     do j = 1, n
        if (j /= i) then
           factor = aug(j,i)
           aug(j,:) = aug(j,:) - factor * aug(i,:)
        end if
     end do

  end do

  ! Extract inverse matrix
  Ainv = aug(:,n+1:2*n)

end subroutine invert6
end program composite_3D
