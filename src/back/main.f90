program main
! {{ MAIN }}
 use vector_module 
 implicit none
 integer                  ::  natoms = 10
 type (vector_array)      ::  r,v,f
 real                     ::  Uvdw, virial
 allocate(r%x(natoms),r%y(natoms),r%z(natoms))
 allocate(v%x(natoms),v%y(natoms),v%z(natoms))
 allocate(f%x(natoms),f%y(natoms),f%z(natoms))
 r = 0.0
 v = 0.0
 f = 0.0
 call forces(r,f,natom,Udvdw,virial)
end program main
!============================================================================
subroutine forces(r,f,natom,Uvdw,virial)
 implicit none
 type (vector_array), intent(in) :: r
 type (vector_array),intent(out) :: f
 real                            :: Uvdw = 0.0,virial = 0.0
 integer,intent(in)              :: natom
! local
 real,parameter                  :: (epsilon=1.0,sigma=1.0,sigma2=sigma*sigma)
 real,parameter                  :: (Rcut = 5.0,Rcut2=Rcut*Rcut)
 integer                         :: iatom,jatom
 type (vector)                   :: d,r1,r2
!
 do iatom = 1,natom
  f%x = 0.0; f%y = 0.0; f%z = 0.0
 end do
 !do iatom =1,natom-1
 ! do jatom=iatom+1,natom
 !  r1 = rcryst%x(
 !  ! minimum image criterion
 !  if(abs(d%x)>0.5*
 ! end do
 end do
 return
end subroutine forces
!
subroutine make_distances(flag,cell_0,r2,r1,rv,r3,dist)
! {{ prepara para el calculo de distancias en una celda triclinica }}
 implicit none
 type (vector),intent(in)  :: r1,r2
 type (vector),intent(out) :: dist
 REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
 REAL,    intent(out) :: dist,r3(1:3)                     ! matriz de distancias N X N
 REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
 REAL                 :: distance,rcm(3),phi
 INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
 REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
 LOGICAL              :: flag
! {{ calculamos la matriz de distancias }}
  k=0
  do l=-1,1
   do m=-1,1
      do n=-1,1
         k = k + 1
         ouratom(1) = r1(1)
         ouratom(2) = r1(2)
         ouratom(3) = r1(3)
         atom(1) = r2(1) + l
         atom(2) = r2(2) + m
         atom(3) = r2(3) + n
         d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)
           image(i,k) = atom(i)
         end forall
     enddo
   enddo
  enddo
  dist=MINVAL(d_image)
  if(flag)then
   phi=1000.0
   k=1
   do l=1,27
    if(d_image(l)<=phi)then
 !     PRINT*,d_image(l),( image(m,l), m=1,3 )
      phi=d_image(l) ! seleccionamos el parametro menor
      k=l            ! y el contador correspondiente.
    endif
   enddo
   forall ( l=1:3)
     r3(l)=image(l,k)
   end forall
  else
   r3(1:3)=0.0
  end if
 return
end subroutine

