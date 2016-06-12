module vector_module
! ...
 implicit none
 type  :: vector
  sequence
  real :: x
  real :: y
  real :: z
 end type vector
! 
 type  :: vector_array
  sequence
  real,allocatable :: x(:) :: length
  real,allocatable :: y(:) :: length
  real,allocatable :: z(:) :: length
 end type vector_array
 contains
! 
   TYPE (vector) FUNCTION vector_add (v1,v2)
    IMPLICIT NONE
    TYPE (vector), INTENT(IN) :: v1
    TYPE (vector), INTENT(IN) :: v2
    vector_add%x = v1%x + v2%x
    vector_add%y = v1%y + v2%y
    vector_add%z = v1%z + v2%z
   END FUNCTION vector_add
!
   TYPE (vector) FUNCTION vector_sub (v1,v2)
    IMPLICIT NONE
    TYPE (vector), INTENT(IN) :: v1
    TYPE (vector), INTENT(IN) :: v2
    vector_sub%x = v1%x - v2%x
    vector_sub%y = v1%y - v2%y
    vector_sub%z = v1%z - v2%z
   END FUNCTION vector_sub
!
   TYPE (vector) FUNCTION cross(a,b) 
    IMPLICIT NONE
    TYPE (vector), INTENT (in) :: a, b
    cross%x = a%y * b%z - a%z * b%y
    cross%y = a%z * b%x - a%x * b%z
    cross%z = a%x * b%y - a%y * b%x
   END FUNCTION cross
!
   REAL FUNCTION absvec(a)
    TYPE (vector), INTENT (in) :: a
    absvec = sqrt(a%x**2 + a%y**2 + a%z**2)
   END FUNCTION absvec
end module vector_module
