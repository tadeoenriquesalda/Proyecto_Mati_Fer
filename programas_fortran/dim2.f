

      IMPLICIT NONE

      integer readstring
      integer i,icont,slen,flen,long,j,k,m
      character*1000 card

      open(50,file='dim2')
   

      icont=0
      do while (readstring(5,card,slen).ge.0 )
         long=slen
         if(card(1:5).ne.'     ') then
               icont=icont+1
         endif
      enddo

c      write(50,222) '      parameter(icont2=',icont,')'
      write(50,221) icont
      close(50)

222   FORMAT(A23,I6,A1) 
221   FORMAT(I6)
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
