subroutine GULP(xfrac,n,label,cell_0,CORE_SHELL)
 implicit none
 integer, intent(in)              :: n
 real, intent(in)                 :: xfrac(0:3,1:n),cell_0(6)
 CHARACTER (LEN=2), intent(in)    :: label(1:n,2)
 LOGICAL, intent(in)              :: CORE_SHELL
 integer :: k
!
 OPEN(1,FILE='gin')
 WRITE(1,'(A)')'opti conp phon freq noden'
 WRITE(1,'(A)')'title'
 WRITE(1,'(A,1X,I3)')' Al/Si substitutions',n
 WRITE(1,'(A)')'end'
 WRITE(1,'(A)')'cell'
 WRITE(1,'(6(F14.7))')(cell_0(k), k=1,6)
 WRITE(1,'(A)')'frac'
 write_: do k=1,n
  WRITE(1,'(A,1X,A,1X,F14.7,F14.7,F14.7)')LABEL(k,1),'core',xfrac(1,k),xfrac(2,k),xfrac(3,k)
  IF(CORE_SHELL.and.label(k,2)=='O') THEN
   WRITE(1,'(A,1X,A,1X,F14.7,F14.7,F14.7)')LABEL(K,1),'shel',xfrac(1,k),xfrac(2,k),xfrac(3,k)
  ENDIF
 enddo write_
 close(1)
 RETURN
end subroutine
!
SUBROUTINE escritura_cif(xcryst,n_atoms,label,cell_0,rv)
 IMPLICIT NONE
 INTEGER           :: n_atoms
 REAL              :: xcryst(0:3,n_atoms)
 REAL              :: volume,cell_0(6),rv(1:3,1:3)
 INTEGER           :: I,U
 CHARACTER (LEN=4) :: label(1:n_atoms,1:2)
 U=1000
 OPEN(U,FILE="P1.cif")
 WRITE(U,'(A)')'data_subtitutions'
 WRITE(U,'(A)')'_audit_creation_method    iGOR'
 WRITE(U,'(A)')"_audit_author_name 'Sponge Bob'"
 WRITE(U,'(A,F14.7)')'_cell_length_a',cell_0(1)
 WRITE(U,'(A,F14.7)')'_cell_length_b',cell_0(2)
 WRITE(U,'(A,F14.7)')'_cell_length_c',cell_0(3)
 WRITE(U,'(A,F14.7)')'_cell_angle_alpha',cell_0(4)
 WRITE(U,'(A,F14.7)')'_cell_angle_beta',cell_0(5)
 WRITE(U,'(A,F14.7)')'_cell_angle_gamma',cell_0(6)
 WRITE(U,'(A,F14.7)')'_cell_volume',volume(rv)
 WRITE(U,'(A)')'_symmetry_cell_setting cubic'
 WRITE(U,'(A)')"_symmetry_space_group_name_Hall 'P 1'"
 WRITE(U,'(A)')"_symmetry_space_group_name_H-M 'P 1'"
 WRITE(U,'(A)')'_symmetry_Int_Tables_number 1'
 WRITE(U,'(A)')"_symmetry_equiv_pos_as_xyz 'x,y,z'"
 WRITE(U,'(A)')'loop_'
 WRITE(U,'(A)')'_atom_site_label'
 WRITE(U,'(A)')'_atom_site_fract_x'
 WRITE(U,'(A)')'_atom_site_fract_y'
 WRITE(U,'(A)')'_atom_site_fract_z' 
 atoms_: DO I=1,n_atoms
   WRITE(U,'(A3,3x,3(f12.8,2x))')label(i,1),xcryst(1,i),xcryst(2,i),xcryst(3,i)
 ENDDO atoms_
 CLOSE(U)
 RETURN
END SUBROUTINE escritura_cif
