      INTEGER FUNCTION READSTRING(FILE,LINES,LEN)
	  INTEGER :: FILE,LEN=a
	  CHARACTER :: LINES
	  
	  OPEN(UNIT=10,FILE='Ingrese nombre del archivo pdb: ',STATUS=OLD,ACTION='READ',POSITION='REWIND',IOSTAT=ierror)

	  DO WHILE (a <= 80)
	     READ(10,100) LINES
	  ENDDO
	  
100   FORMAT(A)
