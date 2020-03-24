      IMPLICIT NONE

      integer readstring
      integer i,slen,flen,long,j,k,m
      integer icont
      character*1000 card

      icont=0

      open(50,file='tempXX')
      open(31,file='dimXX')
      do while (readstring(5,card,slen).ge.0 )
         long=slen
         if(card(1:4).eq.'ATOM') then
          if(card(14:15).eq.'CA') then
           if(card(17:17).eq.' '.or. card(17:17).eq.'A' ) then
           write(50,111) card(1:80)
           icont=icont+1
           endif
          endif
         endif
      enddo
      write(31,*) icont

      close(50)
      close(31)
111   FORMAT(A80)
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

