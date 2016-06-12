program main
! {{ MAIN }}
 use vector_module 
 implicit none
 type (vector_array)      ::  r,v,a
 real                     ::  Uvdw, Ukin,virial,dt=0.001,dens=1.0,press
 integer                  ::  natoms = 10
 real                     ::  vr(3,3),rv(3,3),cell_0(6)
 integer                  ::  i
 cell_0(1)=10.0
 cell_0(2)=10.0
 cell_0(3)=10.0
 cell_0(4)=90.0
 cell_0(5)=90.0
 cell_0(6)=90.0
 call cell(rv,vr,cell_0) 
 allocate(r%x(natoms),r%y(natoms),r%z(natoms))
 allocate(v%x(natoms),v%y(natoms),v%z(natoms))
 allocate(a%x(natoms),a%y(natoms),a%z(natoms))
 !
 do i =1,natoms
  r%x(i) = rand()*cell_0(1); r%y(i) = rand()*cell_0(2); r%z(i) = rand()*cell_0(3)
  v%x(i) = rand()*200-1; v%y(i) = rand()*200-1; v%z(i) = rand()*200-1
  a%x(i) = 0.0; a%y(i) = 0.0; a%z(i) = 0.0
 end do
 !
 do i=1,1000
  call verlet(r,v,a,natoms,dt,dens,Ukin,Uvdw,press)
  write(6,*)i*dt,Ukin,Uvdw,press
 end do
 !
 deallocate(r%x,r%y,r%z)
 deallocate(v%x,v%y,v%z)
 deallocate(a%x,a%y,a%z)
 contains
!
 subroutine verlet(r,v,a,natoms,dt,dens,Ukin,Uvdw,press)
  implicit none
  integer             :: natoms,iatom,i
  type(vector_array)  :: r,v,a
  real                :: dt !,dt2
  real                :: dens,press,Ukin,Uvdw
!
  !dt2=0.5*dt
  do iatom = 1,natoms
   v%x(iatom) = v%x(iatom) + a%x(iatom)*0.5*dt
   v%y(iatom) = v%y(iatom) + a%y(iatom)*0.5*dt
   v%z(iatom) = v%z(iatom) + a%z(iatom)*0.5*dt
!
   r%x(iatom) = r%x(iatom) + v%x(iatom)*dt
   r%y(iatom) = r%y(iatom) + v%y(iatom)*dt
   r%z(iatom) = r%z(iatom) + v%z(iatom)*dt
  end do
!
  call BoundCond(r,natoms,vr,rv)
  call forces(r,a,natoms,Uvdw,virial)
!
  Ukin=0.0
  do iatom=1,natoms
   v%x(iatom) = v%x(iatom) + a%x(iatom)*0.5*dt
   v%y(iatom) = v%y(iatom) + a%y(iatom)*0.5*dt
   v%z(iatom) = v%z(iatom) + a%z(iatom)*0.5*dt
   Ukin = Ukin + v%x(iatom)**2+v%y(iatom)**2+v%z(iatom)**2
  end do
  Ukin  = 0.5*Ukin/natoms
  press = dens*(2.0*Ukin+virial)/3.0
  !
  !do i=1,natoms
  ! write(6,'(4(f14.7,2x))')r%x(i),r%y(i),r%z(i),Uvdw
  ! write(6,'(3(f14.7,2x))')v%x(i),v%y(i),v%z(i)
  ! write(6,'(3(f14.7,2x))')a%x(i),a%y(i),a%z(i)
  !end do
  return
 end subroutine verlet
!
 subroutine BoundCond(r,natoms,vr,rv)
  implicit none
  integer             :: natoms,iatom 
  type (vector_array) :: r
  type (vector)       :: r1,rc
  real                :: vr(3,3),rv(3,3)
  do iatom = 1,natoms
   r1%x = r%x(iatom) ; r1%y = r%y(iatom) ; r1%z = r%z(iatom)
! 
   rc%x = mod(vr(1,1)*r1%x + vr(1,2)*r1%y + vr(1,3)*r1%z,1.0)
   rc%y = mod(vr(2,1)*r1%x + vr(2,2)*r1%y + vr(2,3)*r1%z,1.0)
   rc%z = mod(vr(3,1)*r1%x + vr(3,2)*r1%y + vr(3,3)*r1%z,1.0) 
!
   r%x(iatom) = rv(1,1)*rc%x + rv(1,2)*rc%y + rv(1,3)*rc%z
   r%y(iatom) = rv(2,1)*rc%x + rv(2,2)*rc%y + rv(2,3)*rc%z
   r%z(iatom) = rv(3,1)*rc%x + rv(3,2)*rc%y + rv(3,3)*rc%z
  end do
  return
 end subroutine BoundCond
!
!============================================================================
 subroutine forces(r,a,natoms,Uvdw,virial)
 implicit none
 integer              :: natoms
 type (vector_array)  :: r,a
 real                 :: Uvdw,virial
! local
 real,parameter       :: epsilon=1.0,sigma=1.0,sigma2=sigma*sigma
 real                 :: dist,dist2,fr2,fr6,fpr,fxi,fyi,fzi
 real,parameter       :: Rcut = 10.0,Rcut2=Rcut*Rcut
 integer              :: iatom,jatom
 type (vector)        :: r1,r2,dr
!
 do iatom = 1,natoms
  a%x(iatom) = 0.0
  a%y(iatom) = 0.0
  a%z(iatom) = 0.0
 end do
 do iatom =1,natoms-1
  do jatom=iatom+1,natoms
   ! minimum image criterion:
   r1%x = r%x(iatom)
   r1%y = r%y(iatom)
   r1%z = r%z(iatom)
   r2%x = r%x(jatom)
   r2%y = r%y(jatom)
   r2%z = r%z(jatom)
   call make_distances(.true.,r2,r1,dist,dr)
   !
   dist2 = dist*dist
   if (dist2 < Rcut2)then
! Potential:
    fr2 = sigma2 / dist2
    fr6 = fr2*fr2*fr2
    fpr = 48.0 * epsilon*fr6*(fr6-0.5)/(dist2)  ! f/r
! Newton:   
    a%x(iatom) = a%x(iatom) + fpr*dr%x
    a%y(iatom) = a%y(iatom) + fpr*dr%y
    a%z(iatom) = a%z(iatom) + fpr*dr%z
    a%x(jatom) = a%x(jatom) - fpr*dr%x
    a%y(jatom) = a%y(jatom) - fpr*dr%y
    a%z(jatom) = a%z(jatom) - fpr*dr%z
!
    Uvdw   = Uvdw   + 4.0*epsilon*fr6*(fr6-1.0)
    virial = virial + fpr*dist2
   end if
   !write(6,'(3(f10.5,2x))') a%x(iatom),a%y(iatom),a%z(iatom)
   !write(6,'(10(f10.5,2x))')r1%x,r1%y,r1%z,r2%x,r2%y,r2%z,dr%x,dr%y,dr%z,dist
  end do
 end do
 Uvdw = Uvdw/real(natoms)
 virial=virial/real(natoms)
 return
 end subroutine forces
!
 subroutine make_distances(flag,r2,r1,dist,dr)
 implicit none
 type (vector)  :: r1,r2,dr
 real           :: image(3,27),r3(3)
 REAL           :: dist                      ! matriz de distancias N X N
 REAL           :: d_image(1:27)             !,image(3,27)        ! array de distancias
 REAL           :: phi = 1000.0
 INTEGER        :: k,l,m,n,o,i,j             ! variables mudas
 REAL           :: atom(3),ouratom(3)        ! coordenadas preparadas
 logical        :: flag
 ouratom(1) = vr(1,1)*r1%x + vr(1,2)*r1%y + vr(1,3)*r1%z
 ouratom(2) = vr(2,1)*r1%x + vr(2,2)*r1%y + vr(2,3)*r1%z
 ouratom(3) = vr(3,1)*r1%x + vr(3,2)*r1%y + vr(3,3)*r1%z
 k=0
 do l=-1,1
  do m=-1,1
     do n=-1,1
        k = k + 1
        atom(1) = l + vr(1,1)*r2%x + vr(1,2)*r2%y + vr(1,3)*r2%z
        atom(2) = m + vr(2,1)*r2%x + vr(2,2)*r2%y + vr(2,3)*r2%z
        atom(3) = n + vr(3,1)*r2%x + vr(3,2)*r2%y + vr(3,3)*r2%z
        d_image(k) = distance(atom,ouratom,dr)
        forall ( i=1:3)
         image(i,k) = atom(i)
        end forall
    enddo
  enddo
 enddo
 if(flag)then
  k=1
  do l=1,27
   if(d_image(l)<=phi)then
     phi=d_image(l) ! seleccionamos el parametro menor
     k=l            ! y el contador correspondiente.
   endif
  enddo
  atom(1)=image(1,k) 
  atom(2)=image(2,k)
  atom(3)=image(3,k)
  dist = distance(atom,ouratom,dr)
 else
  dist = MINVAL(d_image)
  dr%x = 0.0 
  dr%y = 0.0
  dr%z = 0.0
 end if
 return
end subroutine
!
REAL FUNCTION distance(atom,ouratom,dr)
 IMPLICIT NONE
 type (vector)  :: dr
 INTEGER        :: j
 REAL  :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
 FORALL ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
 END FORALL
 dr%x = dist(1)
 dr%y = dist(2)
 dr%z = dist(3)
 distance = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
END FUNCTION
!
SUBROUTINE cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: pi,DEGTORAD
 pi = ACOS(-1.0)
 DEGTORAD=pi/180.0
 IF(cell_0(4) == 90.0) THEN
   cosa = 0.0
 ELSE
   ALP=cell_0(4)*degtorad
   COSA=cos(ALP)
 ENDIF
 IF(cell_0(5) == 90.0) THEN
   cosb = 0.0
 ELSE
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 ENDIF
 IF(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 ELSE
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 ENDIF
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3)) 
 call inverse(rv,vr,3)
 !print*,'Cell:'
 !WRITE(*,'(6F14.7)')( cell_0(j), j=1,6 )
 !print*,'Box:'
 !DO i=1,3
 !   WRITE(*,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
 !ENDDO
 !WRITE(*,*)'----------------------------------------'
 !WRITE(*,*)'bOX:'
 !DO i=1,3
 !   WRITE(*,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
 !ENDDO
 RETURN
END SUBROUTINE cell
!
SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out) :: cell_0(6)
  real,intent(in)  :: rv(3,3)
  integer  :: i,j
  real     :: temp(6)
  REAL :: radtodeg,PI
  PI=ACOS(-1.0)
  radtodeg=180.0/PI
!
  do i = 1,3
    temp(i) = 0.0
    do j = 1,3
      temp(i) = temp(i) + rv(j,i)**2
    enddo
    temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
    temp(3+i) = 0.0
  enddo
  do j = 1,3
    temp(4) = temp(4) + rv(j,2)*rv(j,3)
    temp(5) = temp(5) + rv(j,1)*rv(j,3)
    temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
!  Avoid round off errors for 90.0 and 120.0 degrees
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
!
  return
end subroutine uncell
!
subroutine inverse(a,c,n)
 implicit none 
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
 L=0.0
 U=0.0
 b=0.0
! step 1: forward elimination
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
 do i=1,n
  L(i,i) = 1.0
 end do
! U matrix is the upper triangular part of A
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
!
! Step 3: compute columns of the inverse matrix C
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 return
END SUBROUTINE inverse
!
real function volume(rv)
  implicit none
  real, intent(in)  :: rv(3,3)
  real       :: r1x
  real       :: r1y
  real       :: r1z
  real       :: r2x
  real       :: r2y
  real       :: r2z
  real       :: r3x
  real       :: r3y
  real       :: r3z
  real       :: vol
!
  r1x = rv(1,1)
  r1y = rv(2,1)
  r1z = rv(3,1)
  r2x = rv(1,2)
  r2y = rv(2,2)
  r2z = rv(3,2)
  r3x = rv(1,3)
  r3y = rv(2,3)
  r3z = rv(3,3)
  vol = r1x*(r2y*r3z - r2z*r3y) + r1y*(r3x*r2z - r3z*r2x) + r1z*(r2x*r3y - r2y*r3x)
  volume = abs(vol)
  RETURN
end function 
end program
