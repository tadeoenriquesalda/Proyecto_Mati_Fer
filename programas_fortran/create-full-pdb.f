      IMPLICIT NONE

      integer readstring
      integer i,slen,flen,long,j,k,m
      integer icont,icont2
      character*1000 card


      icont=0
      open(48,file='dat_scwrl')
      open(49,file='dat_pdbs')
      do while (readstring(48,card,slen).ge.0 )
         long=slen
         if(card(1:4).eq.'ATOM') then
           write(49,109) card(1:80)
           icont=icont+1
         endif
      enddo
      close(48)
      close(49)

      icont2=0
      open(48,file='archivo.pdb-original')
      open(49,file='dat_bf')
      do while (readstring(48,card,slen).ge.0 )
         long=slen
         if(card(1:4).eq.'ATOM') then
c|ATOM      2  CA  CYS A  10      15.972  30.146  25.450  1.00 35.29 
c|ATOM   1806  CA  ASN D 124       8.178  49.693  67.736  1.00 45.54
           if(card(17:17).eq.' '.or. card(17:17).eq.'A' ) then
             if(card(14:15).eq.'N ') then
               write(49,109) card(1:80)
               icont2=icont2+1
             elseif(card(14:15).eq.'CA') then
               write(49,109) card(1:80)
               icont2=icont2+1
             elseif(card(14:15).eq.'C ') then
               write(49,109) card(1:80)
               icont2=icont2+1
             elseif(card(14:15).eq.'O ') then
               write(49,109) card(1:80)
               icont2=icont2+1
             elseif(card(14:15).eq.'CB') then
               write(49,109) card(1:80)
               icont2=icont2+1
             elseif((card(14:15).ne.'N '.or.card(14:15).ne.'CA'
     &.or. card(14:15).ne.'C '.or.card(14:15).eq.'O '
     &.or.card(14:15).ne.'CB')
     &.and.card(14:14).ne.'H'.and.card(13:13).ne.'H') then
               write(49,109) card(1:80)
               icont2=icont2+1
             endif
           endif
         endif
      enddo
      close(48)
      close(49)

      call calculo(icont,icont2)

109   FORMAT(A80)
      stop
      end

C Fin del programa principal
C**************************************
C Subroutine calculo

      subroutine calculo(icont,icont2)
      integer ifail,ii
      integer i,icont,slen,flen,long,j,k,m
      parameter(icontmax=10000)
c### Coord de Main Chain ###
      character*7 priAT(icontmax)
      integer nlin(icontmax)
      character*3 AAtom(icontmax)
      character*1 AAocu(icontmax)
      character*3 AAid(icontmax)
      character*1 CHid(icontmax)
      integer npdb(icontmax)
      double precision r(icontmax,3)
      double precision bfactor(icontmax),xx(icontmax)
      character*2 AA2(icontmax),AA3(icontmax)

      character*7 priATO(icontmax)
      integer nlinO(icontmax)
      character*3 AAtomO(icontmax)
      character*1 AAocuO(icontmax)
      character*3 AAidO(icontmax)
      character*1 CHidO(icontmax)
      integer npdbO(icontmax)
      double precision rO(icontmax,3)
      double precision bfactorO(icontmax),xxO(icontmax)
      character*2 AA2O(icontmax),AA3O(icontmax)


      open(53,file='dat_pdbs')
      do i=1,icont
        read(53,110) priAT(i),nlin(i),AAtom(i),AAocu(i),AAid(i),
     &CHid(i),npdb(i),r(i,1),r(i,2),r(i,3),xx(i),AA2(i),AA3(i)
      enddo
      close(53)

      open(54,file='dat_bf')
        do i=1,icont2
        read(54,112) priATO(i),nlinO(i),AAtomO(i),AAocuO(i),AAidO(i),
     &CHidO(i),npdbO(i),rO(i,1),rO(i,2),rO(i,3),xxO(i),bfactorO(i),
     &AA2O(i),AA3O(i)
        enddo
      close(54)

      if ( icont .eq. icont2 ) then
        do i=1,icont
           do j=1,icont
             if ( AAtomO(j).eq.AAtom(i).and.
     &AAidO(j).eq.AAid(i).and.
     &npdbO(j).eq.npdb(i) ) then
             bfactor(i)=bfactorO(j)
             endif
           enddo
        enddo
        open(55,file='full_pdbs')      
        do i=1,icont
          write(55,112) priAT(i),nlin(i),AAtom(i),AAocu(i),AAid(i),
     &CHid(i),npdb(i),r(i,1),r(i,2),r(i,3),xx(i),bfactor(i),
     &AA2(i),AA3(i)
        enddo
          write(55,116) 'TER'
          write(55,116) 'END'
        close(55)
      elseif ( icont .ne. icont2 ) then
        do i=1,icont
           do j=1,icont2
             if ( AAtomO(j).eq.AAtom(i).and.
     &AAidO(j).eq.AAid(i).and.
     &npdbO(j).eq.npdb(i) ) then
               rO(j,1)=r(i,1)
               rO(j,2)=r(i,2)
               rO(j,3)=r(i,3)
               AA2O(j)=AA2(i)
               AA3O(j)=AA3(i)
             endif
           enddo
        enddo

       open(55,file='full_pdbs')
        do i=1,icont2
          write(55,112) priATO(i),nlinO(i),AAtomO(i),AAocuO(i),AAidO(i),
     &CHidO(i),npdbO(i),rO(i,1),rO(i,2),rO(i,3),xxO(i),bfactorO(i),
     &AA2O(i),AA3O(i) 
        enddo
          write(55,116) 'TER'
          write(55,116) 'END'
        close(55)
      endif

c|AAAAAABBBBB CCCCDEEE FGGGGH   IIIIIIIIJJJJJJJJKKKKKKKKLLLLLLMMMMMM          NNOO
c|         1         2         3         4         5         6         7         8
c|12345678901234567890123456789012345678901234567890123456789012345678901234567890
c|ATOM      2  CA  LEU A   1     107.125  21.792  60.131  0.00  0.00           C

c|Asi deberia llegar del pdb original
c|ATOM      2  CA  LEU A   1     107.125  21.792  60.131  0.00  0.00           C
c|ATOM   3118  CB  ALA   188       7.454  31.535  10.406  1.00 33.19

c|llega del SCWR4
c|ATOM   1652  N   GLY _ 214     -10.467  27.845  21.160  1.00                 N

c|Asi deberia quedar
c|ATOM      2  CA  LEU A   1     107.125  21.792  60.131  0.00  0.00           C

110   FORMAT(A7,I4,1X,A4,A1,A3,1x,A1,I4,4x,3(F8.3),F6.2,16x,2(A2))
112   FORMAT(A7,I4,1x,A4,A1,A3,1x,A1,I4,4x,3(F8.3),2(F6.2),10x,2(A2))
123   format(A6,I5,1x,A4,A1,A3,1x,A1,I4,4x,F8.3,F8.3,F8.3,1x,F5.2,F6.2)
115   FORMAT(F6.2)
116   FORMAT(A3)
      return
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

