      IMPLICIT NONE

      integer readstring
      integer i,slen,flen,long,j,k,m,p,q
      integer icont,icont2,icontmax
      character*1000 card
      parameter(icontmax=10000)
      character*1000 coordx,cardX,cardY
      character*1000 cardZ
      double precision cmx(1),cmy(1),cmz(1)
      double precision x(icontmax),y(icontmax),z(icontmax)
      character*1000 cardA

C*************************************************************************
C Limpio el pdb y dejo sólo los átomos con todos sus parámetros
C*************************************************************************

      icont=0

      open(50,file='tempXX')
      do while (readstring(5,card,slen).ge.0 )
       long=slen
        if(card(1:4).eq.'ATOM') then
         write(50,111) card(1:80)
         icont=icont+1
        endif
      enddo
      write(50,112) 'Cantidad de átomos: ',icont

      close(50)
111   FORMAT(A80)
112   FORMAT(/,A,1x,I6)

C*************************************************************************
C Calculo el centro de masa
C*************************************************************************
      open(50,file='tempXX')
      open(52,file='tempY')

      icont2=0
C      coordx=ichar(cardB(31:38))
      do p=1,icont
       read(52,*) cardA(1:30),cardX(31:38) !,cardY(39:46),cardZ(47:54)
       x(p)=dble(ichar(cardX(31:38)))
C       y(p)=cardY
C       z(p)=cardZ
       do q=1,icont
        cmx(q)=cmx(q)+x(p)/icont
C        cmy(q)=cmy(q)+y(p)/icont
C        cmz(q)=cmy(q)+z(p)/icont

        icont2=icont2+1

       enddo
      enddo

      write(50,113) 'Centro de masa: ',cmx(1),cmy(1),cmz(1)

113   FORMAT(/,A,1x,I6,1x,I6,1x,I6)
114   FORMAT(A30,F8.3,F8.3,F8.3)
      close(50)

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