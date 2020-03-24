! Programa para leer un archivo pdb y copiar los carbonos alfa en un fichero nuevo
      PROGRAM ATOM_CA.f90
      IMPLICIT NONE
      INTEGER :: i,n,icont,FILE,LEN
	  CHARACTER :: LINES

      OPEN(UNIT=10,FILE='Ingrese nombre del archivo pdb: ',STATUS=OLD,ACTION='READ',POSITION='REWIND',IOSTAT=ierror)
      OPEN(UNIT=20,FILE='carbonos_alfa.txt',STATUS=NEW,ACTION='WRITE')

      icont=0

      DO i = 1, n
	     

100   FORMAT(A)

!-----SUBROUTINES

      INTEGER FUNCTION READSTRING(FILE,LINES,LEN)
	  INTEGER :: FILE,LEN=a
	  CHARACTER :: LINES
	  
	  DO WHILE (a <= 80)
	     READ(10,100) LINES
	  ENDDO