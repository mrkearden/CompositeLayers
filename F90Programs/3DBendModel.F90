program composite_3D
  implicit none

  integer :: n, i
  real(8), allocatable :: E1(:), E2(:), E3(:)
  real(8), allocatable :: G12(:), G13(:), G23(:)
  real(8), allocatable :: nu12(:), nu13(:), nu23(:)
  real(8), allocatable :: theta(:), t(:), V(:), z(:)

  real(8) :: A(3,3), Ainv(3,3), D(3,3)
  real(8) :: Q11, Q22, Q12, Q66, nu21, denom
  real(8) :: Qb11, Qb22, Qb12, Qb66, Qb16, Qb26
  real(8) :: c, s, h

  ! Final properties
  real(8) :: Ex_mem, Ex_bend, Ey_mem, Ey_bend
  real(8) :: Gxy_mem
  real(8) :: nuxy
  ! Final properties
  real(8) :: Ex, Ey, Ez
  real(8) :: Gxy, Gxz, Gyz
  real(8) :: nuxz, nuyz
  real(8) :: Ez_inv, Gxz_inv, Gyz_inv

OPEN(UNIT=20, FILE='LayerData.txt', STATUS='OLD', ACTION='READ')
Read(20,*) n
     print *, "Layers: ", n

  allocate(E1(n),E2(n),E3(n),G12(n),G13(n),G23(n))
  allocate(nu12(n),nu13(n),nu23(n),theta(n),t(n),V(n))
  allocate(z(n))
  print *, "Enter per layer:"
  print *, "(E1 E2 E3 G12 G13 G23 nu12 nu13 nu23 theta(deg) thickness)"

  do i = 1,n
    print *, "Layer ", i
    read (20,*) E1(i),E2(i),E3(i),G12(i),G13(i),G23(i), &
             nu12(i),nu13(i),nu23(i),theta(i),t(i)
  end do

  ! Thickness + volume fractions
  h = sum(t)
  V = t / h

  ! z-coordinates (critical for D)
  z(0) = -h/2.0d0
  do i = 1,n
     z(i) = z(i-1) + t(i)
  end do
print *, "Z complete"
  A = 0.0d0
  D = 0.0d0

  do i = 1,n

    c = cos(theta(i)*acos(-1.0d0)/180.0d0)
    s = sin(theta(i)*acos(-1.0d0)/180.0d0)

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

    ! --- A matrix ---
    A(1,1)=A(1,1)+Qb11*t(i)
    A(2,2)=A(2,2)+Qb22*t(i)
    A(1,2)=A(1,2)+Qb12*t(i)
    A(2,1)=A(1,2)

    A(3,3)=A(3,3)+Qb66*t(i)

    A(1,3)=A(1,3)+Qb16*t(i)
    A(3,1)=A(1,3)

    A(2,3)=A(2,3)+Qb26*t(i)
    A(3,2)=A(2,3)

    ! --- D matrix (NEW) ---
    D(1,1)=D(1,1)+Qb11*(z(i)**3 - z(i-1)**3)/3.0d0
    D(2,2)=D(2,2)+Qb22*(z(i)**3 - z(i-1)**3)/3.0d0
    D(1,2)=D(1,2)+Qb12*(z(i)**3 - z(i-1)**3)/3.0d0
    D(2,1)=D(1,2)

    D(3,3)=D(3,3)+Qb66*(z(i)**3 - z(i-1)**3)/3.0d0

    D(1,3)=D(1,3)+Qb16*(z(i)**3 - z(i-1)**3)/3.0d0
    D(3,1)=D(1,3)

    D(2,3)=D(2,3)+Qb26*(z(i)**3 - z(i-1)**3)/3.0d0
    D(3,2)=D(2,3)

  end do

  call invert3x3(A,Ainv)

  ! --- Membrane ---
  Ex_mem = 1.0d0/(Ainv(1,1)*h)
  Ey_mem = 1.0d0/(Ainv(2,2)*h)
  Gxy_mem = 1.0d0/(Ainv(3,3)*h)
  nuxy = -Ainv(1,2)/Ainv(1,1)

  ! --- Bending (CRITICAL FIX) ---
  Ex_bend = 12.0d0*D(1,1)/(h**3)
  Ey_bend = 12.0d0*D(2,2)/(h**3)

  print *, "-----------------------------"
  print *, "MEMBRANE:"
  print *, "Ex =",Ex_mem
  print *, "Ey =",Ey_mem
  print *, "Gxy=",Gxy_mem
  print *, "nu_xy=",nuxy

  print *, "-----------------------------"
  print *, "BENDING:"
  print *, "Ex_bend =",Ex_bend
  print *, "Ey_bend =",Ey_bend
  print *, "-----------------------------"

  ! In-plane
  Ex   = Ex_bend
  Ey   = Ey_bend
  Gxy  = 1.0d0/(Ainv(3,3)*h)
  nuxy = A(1,2)/A(2,2)

  ! Transverse
  Ez  = E3(1)
  Gxz = G13(1)
  Gyz = G23(1)
  print *, "-----------------------------"
  print *, "3D Equivalent Properties"
  print *, "Ex  =",Ex
  print *, "Ey  =",Ey
  print *, "Ez  =",E3(1)
  print *, "Gxy =",G13(1)
  print *, "Gxz =",G23(1)
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

open(unit=11,file="abaqus_section.inp", status='UNKNOWN')

write(11,*) "*MATERIAL, NAME=PLYMAT"
write(11,*) "*ELASTIC, TYPE=ENGINEERING CONSTANTS"
write(11,'(9(E12.5,1X))') E1(1),E2(1),E3(1),nu12(1),nu13(1),nu23(1), &
                         G12(1),G13(1),G23(1)

write(11,*) "*SHELL SECTION, COMPOSITE, ELSET=BEAM"

do i=1,n
   write(11,'(F10.5,1X,A,1X,F10.5)') t(i),"PLYMAT",theta(i)
end do

close(11)

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

end program composite_3D
