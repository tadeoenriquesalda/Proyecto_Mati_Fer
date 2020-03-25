	PROGRAM CA.f
	IMPLICIT NONE
		CHARACTER :: string,substring1,substring2
		
		string = 'ATOM      2  CA  LEU A   4      -6.696  22.003  26.447  1.00 48.82           C'
		substring1 = 'ATOM'
		substring2 = 'CA'
	  
		IF(INDEX(string, substring1) = 1)THEN
			IF(INDEX(string, substring2) = 14)THEN
				PRINT(*,*)string
			ELSE
				PRINT (*,*)'No CA found'
		ELSE
			PRINT (*,*)'No CA found'
			ENDIF
		ENDIF
		
	END PROGRAM CA.f