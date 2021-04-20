      IMPLICIT NONE

      integer readstring
      integer i,slen,flen,long,j,k,m
      integer icont
      character*1000 card,card1,card2

      icont=0

      open(50,file='tempXX')
      open(31,file='dimXX')
      open(40,file='coordXX')
      do while (readstring(5,card,slen).ge.0 )
         long=slen
         if(card(1:4).eq.'ATOM') then
          if(card(14:15).eq.'CA'.or. card(14:15).eq.'O ') then
           if(card(17:17).eq.' '.or. card(17:17).eq.'A' ) then
            if(card(22:22).eq.'A') then
            card1=card(1:15)
            card2=card(31:54)
            write(50,111) card(1:80)
            write(40,113) card1,card2
            icont=icont+1
            endif
           endif
          endif
         endif
      enddo
      write(31,112) 'Cantidad de átomos: ',icont
      write(50,112) 'Cantidad de átomos: ',icont

      close(50)
      close(31)
      close(40)
111   FORMAT(A80)
112   FORMAT(/,A,1x,I6)
113   FORMAT(A15,1x,A24)

      stop
      end


c------ SUBROUTINES

      integer function readstring(file,card,flen)
      integer       file, flen
      character*1000 card
      if(file.gt.200)STOP'ERROR: file number too large'
      read(file,'(a)',err=100,end=100)card
      flen=1000
      do while(card(flen:flen).eq.' ')
         flen=flen-1
      enddo
      readstring=flen
      return
 100  readstring=-1

      return
      end