C Programa para calcular los modos normales utilizando el metodo de ANM 

      IMPLICIT NONE

      integer readstring,icont,icontHET,icontR
      integer i,slen,flen,long,j,k,m
      character*1000 card
      character*8 saber

      open(48,file='archivo.pdb')
      open(49,file='temp-mainch')
      open(50,file='temp-grupoR')
c      open(51,file='temp-hetatm')

      icont=0
      icontR=0
      icontHET=0

      do while (readstring(48,card,slen).ge.0 )
         long=slen
         if(card(1:4).eq.'ATOM') then
c|ATOM      2  CA  CYS A  10      15.972  30.146  25.450  1.00 35.29 
c|ATOM   1806  CA  ASN D 124       8.178  49.693  67.736  1.00 45.54
           if(card(17:17).eq.' '.or. card(17:17).eq.'A' ) then
             if(card(14:15).eq.'N ') then
               write(49,110) card(8:68)
               icont=icont+1
             elseif(card(14:15).eq.'CA') then
               write(49,110) card(8:68)
               icont=icont+1
             elseif(card(14:15).eq.'C ') then
               write(49,110) card(8:68)
               icont=icont+1
             elseif(card(14:15).eq.'O ') then
               write(49,110) card(8:68)
               icont=icont+1
             elseif(card(14:15).eq.'CB') then
               write(49,110) card(8:68)
               icont=icont+1
             elseif((card(14:15).ne.'N '.or.card(14:15).ne.'CA'
     &.or. card(14:15).ne.'C '.or.card(14:15).eq.'O '
     &.or.card(14:15).ne.'CB')
     &.and.card(14:14).ne.'H'.and.card(13:13).ne.'H') then
                 write(50,110) card(8:68)
                 icontR=icontR+1
             endif
           endif
         endif
c         if(card(1:6).eq.'HETATM') then
c            if(card(18:20).ne.'HOH') then
c              if(card(14:14).ne.'H') then
c              if(card(13:13).ne.'H') then
cc|HETATM 1814  S   DMS A 201      24.610  19.854  47.132  1.00 50.73           S
cc|HETATM 4154  CAEA3MI D 128      -0.058   0.382 -30.131  0.25 17.22           C
c                write(51,110) card(8:68)
c                icontHET=icontHET+1
c              endif
c              endif
c            endif
c         endif
      enddo
      close(48)
      close(49)
      close(50)
c      close(51)

      call calculo(icont,icontHET,icontR)

c|1806  CA  ASN D 124       8.178  49.693  67.736  1.00 45.54
110   FORMAT(A59)
111   FORMAT(A80)
      stop
      end
C Fin del programa principal
C**************************************
C Subroutine calculo
      subroutine calculo(icont,icontHET,icontR)
      integer ifail,ii,icontmax
      integer i,icont,slen,flen,long,j,k,m,t
      parameter(icontmax=3000)
      double precision cutoff
      double precision w(icontmax,icontmax)
      double precision x(3*icontmax),col
      double precision SUMBFAC,SUMCORR,beta
      double precision bfacteor(icontmax)
      double precision bftprom,bfeprom
      double precision AA,AB,AC,CORSTD
      double precision s(icontmax,icontmax)
      double precision hess(3*icontmax,3*icontmax)
      double precision diag(3*icontmax,3*icontmax)
      double precision freq(3*icontmax),e(3*icontmax)
      double precision corr(3*icontmax,3*icontmax)
      character*1000 card
      integer icontMC
c### Coord de Main Chain ###
      integer nlin(icontmax)
      character*3 AAtom(icontmax)
      character*1 AAocu(icontmax)
      character*3 AAid(icontmax)
      character*1 CHid(icontmax)
      integer npdb(icontmax)
      double precision r(icontmax,3)
      double precision bfactor(icontmax),xx(icontmax)
c### Calculo de CM de la cadena lateral ###
      integer icontR
      integer h,kk(icontmax)
      integer nlinR(icontmax),tnlinR(icontmax)
      character*3 cardR(icontmax),tcardR(icontmax)
      character*1 ocuR(icontmax),tARocu(icontmax)
      character*3 ARid(icontmax),tARid(icontmax)
      character*1 Rchid(icontmax),tRchid(icontmax)
      integer Rnpdb(icontmax),tRnpdb(icontmax)
      double precision d(icontmax,3),cm1(icontmax,3)
      double precision xzx(icontmax),txzx(icontmax)
      double precision bfR(icontmax),tbfR(icontmax)
      integer icontGR,icontS1
c### Coordenadas del Ligando ###
      integer icontHET
      integer nlinHet(icontmax)
      character*3 AAtomHet(icontmax)
      character*1 AAocuHet(icontmax)
      character*3 ligHet(icontmax)
      character*1 CHidHet(icontmax)
      integer npdbHet(icontmax)
      double precision dHet(icontmax,3)
      double precision xzxHet(icontmax),bfHet(icontmax)
c### Definicion de contactos ###
c     RING
      integer icontCONT,z
      character*3 AAri1(icontmax),AAri2(icontmax)
      character*1 CHri1(icontmax),CHri2(icontmax)
      integer nring1(icontmax),nring2(icontmax)
      character*15 contac(icontmax)
      character*2 ATring1(icontmax),ATring2(icontmax)
c     LPC
      integer icontCONT2
      character*3 ATlpc1(icontmax),ATlpc2(icontmax)
      character*3 AAlpc1(icontmax),AAlpc2(icontmax)
      character*1 CHlpc1(icontmax)
      character*2 tipolpc(icontmax)
      integer nlpc1(icontmax)
      double precision dlpc(icontmax),asalpc(icontmax)
c     ATOMOS LIGANDOS
      integer icontCONT3
      character*3 ATliga1(icontmax),ATliga2(icontmax)
      character*3 AAliga1(icontmax),AAliga2(icontmax)
      character*1 CHliga1(icontmax),CHliga2(icontmax)
      integer npdbliga1(icontmax),npdbliga2(icontmax)
      character*3 Cov(icontmax)
      character*6 priAT
      character*3 termi,endi
c     Argumento de entrada
      character*13 cheqcorrelac
      character*10 cheqcolect

      read (5,*) cutoff

C###  Cadena Principal ###
      open(49,file='temp-mainch')
      do i=1,icont
         read(49,113) nlin(i),AAtom(i),AAocu(i),AAid(i),CHid(i),npdb(i),
     &r(i,1),r(i,2),r(i,3),xx(i),bfactor(i)
      enddo
      close(49)

C###  Calculo para el CM del grupo R ###
      open(50,file='temp-grupoR')
      do h=1,icontR
          read(50,113) nlinR(h),cardR(h),ocuR(h),ARid(h),Rchid(h),
     &Rnpdb(h),d(h,1),d(h,2),d(h,3),xzx(h),bfR(h)
      enddo
      close(50)

      ARid(icontR+1)='nad'
      Rnpdb(icontR+1)=9999
cDefino el numero de atomos por Grupo R 
      t=0
      do h=1,icontR
         if (Rnpdb(h).ne.Rnpdb(h-1).and.
     &Rnpdb(h-1).ne.Rnpdb(icontR+1)) then
           t=t+1
           kk(t)=0
         endif
         if (Rnpdb(h).ne.Rnpdb(h-1)) then
           kk(t)=kk(t)+1
         elseif (Rnpdb(h).eq.Rnpdb(h-1)) then
            kk(t)=kk(t)+1
         endif
      enddo
cDefino los CM de cada Grupo R
      t=0
      do h=1,icontR
         if (Rnpdb(h).ne.Rnpdb(h-1).and.
     &Rnpdb(h-1).ne.Rnpdb(icontR+1)) then
            t=t+1
            tnlinR(t)=nlinR(h)
            tcardR(t)=cardR(h)
            tARocu(t)=ocuR(h)
            tARid(t)=ARid(h)
            tRchid(t)=Rchid(h)
            tRnpdb(t)=Rnpdb(h)
            txzx(t)=xzx(h)
         endif
         if (Rnpdb(h).ne.Rnpdb(h-1)) then
            do j=1,3
               cm1(t,j)=cm1(t,j)+d(h,j)/dfloat(kk(t))
            enddo
            tbfR(t)=tbfR(t)+bfR(h)/dfloat(kk(t))
         elseif (Rnpdb(h).eq.Rnpdb(h-1)) then
            do j=1,3
               cm1(t,j)=cm1(t,j)+d(h,j)/dfloat(kk(t))
            enddo
            tbfR(t)=tbfR(t)+bfR(h)/dfloat(kk(t))
         endif
      enddo

      do i=1,t
        if (tcardR(i) .eq.'CG1'.or.tcardR(i).eq.'CG2') then
          tcardR(i)='CG '
        elseif (tcardR(i) .eq.'OG1'.or.tcardR(i).eq.'OG2') then
          tcardR(i)='OG '
        endif
      enddo

      OPEN(52, file='temp-grupoRNEW')
      do i=1,t
         write(52,113) tnlinR(i),tcardR(i),tARocu(i),tARid(i),tRchid(i),
     &tRnpdb(i),cm1(i,1),cm1(i,2),cm1(i,3),txzx(i),tbfR(i)
      enddo
      CLOSE(52)

      icontGR=t
      do i=1,icontGR
         nlin(icont+i)=tnlinR(i)
         AAtom(icont+i)=tcardR(i)
         AAocu(icont+i)=tARocu(i)
         AAid(icont+i)=tARid(i)
         CHid(icont+i)=tRchid(i)
         npdb(icont+i)=tRnpdb(i)
         r(icont+i,1)=cm1(i,1)
         r(icont+i,2)=cm1(i,2)
         r(icont+i,3)=cm1(i,3)
         xx(icont+i)=txzx(i)
         bfactor(icont+i)=tbfR(i)
      enddo

      icontS1=icont+icontGR

      icontMC=icont
      icont=icontMC+icontGR

C### ORDENA LAS LINEAS LUEGO DEL CM del grupo R ###
      CALL ORDEN(icont,nlin,AAtom,AAocu,AAid,CHid,npdb,
     &r,xx,bfactor)
      priAT='ATOM  '
      termi='TER'
      endi='END'
ctadeo check
      do i=1,icont
        if ( CHid(i).eq.' ') then
         CHid(i)='_'
        endif
      enddo
ctadeo check
      open(53,file='temp-coord')
      do i=1,icont
         write(53,123) priAT,nlin(i),AAtom(i),AAocu(i),AAid(i),CHid(i)
     &,npdb(i),r(i,1),r(i,2),r(i,3),xx(i),bfactor(i)
      enddo
         write(53,124) termi
         write(53,124) endi
      close(53)


C### Saber distancia entre atomos ###
      do i=1,icont
         do j=1,icont
             s(i,j)=dsqrt((r(j,1)-r(i,1))**2+
     $(r(j,2)-r(i,2))**2+(r(j,3)-r(i,3))**2)
         enddo
      enddo

C### CONTACTOS ############################
c                     O               "N "
c                     ||              "CA"
c  ...C-/-N---Calfa---C-/-N...        "C "
c             |                       "O "
c             Cbeta                   "CB"
c             |                       "CG" o "OG" o "SG"
c             R(CM)
c| nlin(i),AAtom(i),AAocu(i),AAid(i),CHid(i),npdb(i),r(i,1),r(i,2),r(i,3),xx(i),bfactor(i)
c| 7       N        A(oBoC)  PRO     A       11      14.368 30.118 27.192 1.00  25.63

      do i=1,icont
        do j=1,icont
          if (s(i,j).le.cutoff) then
            w(i,j)=0.001d0
            w(j,i)=0.001d0
          endif
        enddo
      enddo

      do i=1,icont
       do j=1,icont
        if (i.ne.j) then
        if (s(i,j).le.cutoff) then
          if (AAtom(i).eq.'N  ') then
            if (AAtom(j).eq.'CA '.and.npdb(i).eq.npdb(j)
     &.and.CHid(i).eq.CHid(j)) then
               w(i,j)=1.0d0
               w(j,i)=1.0d0
            endif
          elseif (AAtom(i).eq.'CA ') then
            if (AAtom(j).eq.'C  '.and.npdb(i).eq.npdb(j)
     &.and.CHid(i).eq.CHid(j).or.AAtom(j).eq.'CB '
     &.and.npdb(i).eq.npdb(j).and.CHid(i).eq.CHid(j)) then
               w(i,j)=1.0d0
               w(j,i)=1.0d0
            endif
          elseif (AAtom(i).eq.'C  ') then
            if (AAtom(j).eq.'O  '.and.npdb(i).eq.npdb(j)
     &.and.CHid(i).eq.CHid(j).or.AAtom(j).eq.'N  '
     &.and.npdb(i)+1.eq.npdb(j).and.CHid(i).eq.CHid(j)) then
               w(i,j)=1.0d0
               w(j,i)=1.0d0
            endif
          elseif (AAtom(i).eq.'CB ') then
            if (AAtom(j).eq.'CG '
     &.and.npdb(i).eq.npdb(j).and.CHid(i).eq.CHid(j)
     &.or.AAtom(j).eq.'OG '
     &.and.npdb(i).eq.npdb(j).and.CHid(i).eq.CHid(j)
     &.or.AAtom(j).eq.'SG '
     &.and.npdb(i).eq.npdb(j).and.CHid(i).eq.CHid(j)) then
               w(i,j)=1.0d0
               w(j,i)=1.0d0
            endif
          endif
        endif
        endif
       enddo
      enddo
c|ATOM     43  SG  CYS     5      13.025  27.990   5.327  1.00  5.75

c ### Defino las interacciones entre residuos del RING ###
c|A 11 PRO HBOND:MC_MC A 105 TYR O N
c|A 15 LYS HBOND:SC_MC A 52 SER NZ O
c|A 21 ARG HBOND:SC_SC A 23 SER NE OG
      icontCONT=0
      OPEN(42, file='temp-contactos')
      do i=1,icontmax
          read(42,*,end=101) CHri1(i),nring1(i),AAri1(i),contac(i),
     &CHri2(i),nring2(i),AAri2(i),ATring1(i),ATring2(i)
           icontCONT=icontCONT+1
      enddo
      CLOSE(42)
101   continue

c| nlin(i), AAtom(i)    AAid(i)   CHid(i)   npdb(i)
c| 7        N           PRO       A         11     
c|          ATring1(t)  AAri1(t)  CHri1(t)  nring1(t)  
c|          ATring2(t)  AAri2(t)  CHri2(t)  nring2(t)
      do i=1,icont
      do j=1,icont
       if (i.ne.j) then
       if (s(i,j).le.cutoff) then
        do t=1,icontCONT
        if (contac(t).eq.'VDW:MC_MC') then
              if ( ATring1(t).eq.AAtom(i).and.
     &AAri1(t).eq.AAid(i).and.
     &CHri1(t).eq.CHid(i).and.
     &nring1(t).eq.npdb(i).and.
     &ATring2(t).eq.AAtom(j).and.
     &AAri2(t).eq.AAid(j).and.
     &CHri2(t).eq.CHid(j).and.
     &nring2(t).eq.npdb(j) ) then
                w(i,j)=0.01d0
                w(j,i)=0.01d0
             endif
        elseif (contac(t).eq.'VDW:MC_SC') then
              if (ATring1(t).eq.AAtom(i).and.
     &AAri1(t).eq.AAid(i).and.
     &CHri1(t).eq.CHid(i).and.
     &nring1(t).eq.npdb(i).and.
     &AAri2(t).eq.AAid(j).and.
     &CHri2(t).eq.CHid(j).and.
     &nring2(t).eq.npdb(j)) then
                if (AAtom(j).eq.'CB '.and.
     &ATring2(t).eq.'CB ') then
                  w(i,j)=0.01d0
                  w(j,i)=0.01d0
                elseif ((AAtom(j).eq.'CG '.or.AAtom(j).eq.'OG '
     &.or.AAtom(j).eq.'SG ').and.ATring2(t).ne.'CB ') then
                  w(i,j)=0.01d0
                  w(j,i)=0.01d0
                endif
              endif
        elseif (contac(t).eq.'VDW:SC_MC') then
              if (ATring2(t).eq.AAtom(j).and.
     &AAri2(t).eq.AAid(j).and.
     &CHri2(t).eq.CHid(j).and.
     &nring2(t).eq.npdb(j).and.
     &AAri1(t).eq.AAid(i).and.
     &CHri1(t).eq.CHid(i).and.
     &nring1(t).eq.npdb(i)) then
                if (AAtom(i).eq.'CB '.and.
     &ATring1(t).eq.'CB ') then
                  w(i,j)=0.01d0
                  w(j,i)=0.01d0
                elseif ((AAtom(i).eq.'CG '.or.AAtom(i).eq.'OG '
     &.or.AAtom(i).eq.'SG ').and.ATring1(t).ne.'CB ') then
                  w(i,j)=0.01d0
                  w(j,i)=0.01d0
                endif
              endif
        elseif (contac(t).eq.'VDW:SC_SC') then
              if (AAri1(t).eq.AAid(i).and.
     &CHri1(t).eq.CHid(i).and.
     &nring1(t).eq.npdb(i).and.
     &AAri2(t).eq.AAid(j).and.
     &CHri2(t).eq.CHid(j).and.
     &nring2(t).eq.npdb(j))then
                if ((AAtom(i).eq.'CB '.and.AAtom(j).eq.'CB ').and.
     &(ATring1(t).eq.'CB '.and.ATring2(t).eq.'CB ')) then
                  w(i,j)=0.01d0
                  w(j,i)=0.01d0
                elseif ((AAtom(i).eq.'CB ').and.
     &(AAtom(j).eq.'CG '.or.AAtom(j).eq.'OG '.or.AAtom(j).eq.'SG ').and.
     &(ATring1(t).eq.'CB '.and.ATring2(t).ne.'CB ')) then
                  w(i,j)=0.01d0
                  w(j,i)=0.01d0
                elseif ((AAtom(j).eq.'CB ').and.
     &(AAtom(i).eq.'CG '.or.AAtom(i).eq.'OG '.or.AAtom(i).eq.'SG ').and.
     &(ATring1(t).ne.'CB '.and.ATring2(t).eq.'CB ')) then
                  w(i,j)=0.01d0
                  w(j,i)=0.01d0
                elseif (
     &(AAtom(i).eq.'CG '.or.AAtom(i).eq.'OG '.or.AAtom(i).eq.'SG ').and.
     &(AAtom(j).eq.'CG '.or.AAtom(j).eq.'OG '.or.AAtom(j).eq.'SG ').and.
     &(ATring1(t).ne.'CB '.or.ATring2(t).ne.'CB ')) then
                  w(i,j)=0.01d0
                  w(j,i)=0.01d0
                endif
              endif
        elseif (contac(t).eq.'HBOND:MC_MC') then
              if ( ATring1(t).eq.AAtom(i).and.
     &AAri1(t).eq.AAid(i).and.
     &CHri1(t).eq.CHid(i).and.
     &nring1(t).eq.npdb(i).and.
     &ATring2(t).eq.AAtom(j).and.
     &AAri2(t).eq.AAid(j).and.
     &CHri2(t).eq.CHid(j).and.
     &nring2(t).eq.npdb(j) ) then
                w(i,j)=0.1d0
                w(j,i)=0.1d0
              endif
        elseif (contac(t).eq.'HBOND:MC_SC') then
              if (ATring1(t).eq.AAtom(i).and.
     &AAri1(t).eq.AAid(i).and.
     &CHri1(t).eq.CHid(i).and.
     &nring1(t).eq.npdb(i).and.
     &AAri2(t).eq.AAid(j).and.
     &CHri2(t).eq.CHid(j).and.
     &nring2(t).eq.npdb(j).and.
     &(AAtom(j).eq.'CG '.or.AAtom(j).eq.'OG '.or.AAtom(j).eq.'SG ')
     &) then
                w(i,j)=0.1d0
                w(j,i)=0.1d0
              endif
        elseif (contac(t).eq.'HBOND:SC_MC') then
              if (ATring2(t).eq.AAtom(j).and.
     &AAri2(t).eq.AAid(j).and.
     &CHri2(t).eq.CHid(j).and.
     &nring2(t).eq.npdb(j).and.
     &AAri1(t).eq.AAid(i).and.
     &CHri1(t).eq.CHid(i).and.
     &nring1(t).eq.npdb(i).and.
     &(AAtom(i).eq.'CG '.or.AAtom(i).eq.'OG '.or.AAtom(i).eq.'SG ')
     &) then
                w(i,j)=0.1d0
                w(j,i)=0.1d0
              endif
        elseif (contac(t).eq.'HBOND:SC_SC') then
              if (AAri1(t).eq.AAid(i).and.
     &CHri1(t).eq.CHid(i).and.
     &nring1(t).eq.npdb(i).and.
     &AAri2(t).eq.AAid(j).and.
     &CHri2(t).eq.CHid(j).and.
     &nring2(t).eq.npdb(j).and.
     &(AAtom(i).eq.'CG '.or.AAtom(i).eq.'OG '.or.AAtom(i).eq.'SG ').and.
     &(AAtom(j).eq.'CG '.or.AAtom(j).eq.'OG '.or.AAtom(j).eq.'SG ')
     &) then
                w(i,j)=0.1d0
                w(j,i)=0.1d0
              endif
        elseif (contac(t).eq.'SSBOND:SC_SC') then
              if (AAri1(t).eq.AAid(i).and.
     &CHri1(t).eq.CHid(i).and.
     &nring1(t).eq.npdb(i).and.
     &AAri2(t).eq.AAid(j).and.
     &CHri2(t).eq.CHid(j).and.
     &nring2(t).eq.npdb(j).and.
     &(AAtom(i).eq.'CG '.or.AAtom(i).eq.'OG '.or.AAtom(i).eq.'SG ').and.
     &(AAtom(j).eq.'CG '.or.AAtom(j).eq.'OG '.or.AAtom(j).eq.'SG ')
     &) then
                w(i,j)=1.0d0
                w(j,i)=1.0d0
              endif
        endif
        enddo
       endif
       endif
      enddo
      enddo

c#### Comento porque saco LIGANDO #####
cc ### Defino las interacciones entre residuos del LPC ###
cc|OE2 GLU A 54 Hb D75 O2 2.9 9.4
cc|NZ LYS C 15 Hb D75 O2 3.9 11.8
cc|O SER D 52 Hb D75 O3 5.5 0.2
cc Los N u O del grupo R tienen que cambiarse a CG u OG en el archivo temp-LPC
cc|OG GLU A 54 Hb D75 O2 2.9 9.4    ! Grupo R !
cc|CG LYS C 15 Hb D75 O2 3.9 11.8   ! Grupo R !
cc|O SER D 52 Hb D75 O3 5.5 0.2
c      icontCONT2=0
c      OPEN(42, file='temp-LPC')
c      do i=1,icontmax
c          read(42,*,end=102) ATlpc1(i),AAlpc1(i),CHlpc1(i),nlpc1(i)
c     &,tipolpc(i),AAlpc2(i),ATlpc2(i),dlpc(i),asalpc(i)
c           icontCONT2=icontCONT2+1
c      enddo
c      CLOSE(42)
c102   continue
cc| nlin(i), AAtom(i)    AAid(i)   CHid(i)   npdb(i)
cc| 7        N           PRO       A         11     
cc|          ATlpc1(t)  AAlpc1(t)  CHlpc1(t)  nlpc1(t)  
cc|          ATlpc2(t)  AAlpc2(t) 
c      do i=1,icont
c       do j=1,icont
c        if (i.ne.j) then
c        if (s(i,j).le.cutoff) then
c         do t=1,icontCONT2
c          if (tipolpc(t).eq.'Hb') then
c           if (AAlpc2(t).eq.AAid(j).and.ATlpc2(t).eq.AAtom(j)) then
c            if (ATlpc1(t).eq.AAtom(i).and.AAlpc1(t).eq.AAid(i).and.
c     &CHlpc1(t).eq.CHid(i).and.nlpc1(t).eq.npdb(i)) then
c             if (dlpc(t).le.3.5d0) then
c                w(i,j)=0.1d0
c                w(j,i)=0.1d0
c             endif
c            endif
c           endif
c          endif
c         enddo
c        endif
c        endif
c       enddo
c      enddo

cc ### Defino los contactos entre Atomos del Ligando ###
cc|O10 D75 A 1126 Cov C12 D75 A 1126
cc| FORMAT(A3,1x,A3,1x,A1,1x,I4,1x,A3,1x,A3,1x,A3,1x,A1,1x,I4)
cc|            ATliga1(i)              AAliga1(i) CHliga1(i) npdbliga1(i)
cc|            OAB                     3MI        A          128         Cov 
cc|            ATliga2(i)              AAliga2(i) CHliga2(i) npdbliga2(i)  
cc|            CAM                     3MI        A          128
cc| 4155       CAF         A           3MI        D          128         0.308     1.090     -28.998   0.25      17.47
cc| nlin(i),   AAtom(i),   AAocu(i),   AAid(i),   CHid(i),   npdb(i),    r(i,1),   r(i,2),   r(i,3),   xx(i),    bfactor(i)
cc| nlinHet(h),AAtomHet(h),AAocuHet(h),ligHet(h), CHidHet(h),npdbHet(h), dHet(h,1),dHet(h,2),dHet(h,3),xzxHet(h),bfHet(h)
c      icontCONT3=0
c      OPEN(42, file='temp-AtomLig')
c      do i=1,100
c       read(42,*,end=103) ATliga1(i),AAliga1(i),CHliga1(i),npdbliga1(i),
c     &Cov(i),ATliga2(i),AAliga2(i),CHliga2(i),npdbliga2(i)
c           icontCONT3=icontCONT3+1
c      enddo
c      CLOSE(42)
c103   continue

c      do i=1,icont
c        do j=1,icont
c          if (i.ne.j) then
c            if (s(i,j).le.cutoff) then
c              do t=1,icontCONT3
c                if (AAliga1(t).eq.AAid(i).and.ATliga1(t).eq.AAtom(i)
c     &.and.CHliga1(t).eq.CHid(i)) then
c                  if (AAliga2(t).eq.AAid(j).and.ATliga2(t).eq.AAtom(j)
c     &.and.CHliga2(t).eq.CHid(j)) then
c                    w(i,j)=1.0d0
c                    w(j,i)=1.0d0
c                  endif
c                endif
c              enddo
c            endif
c          endif
c        enddo
c      enddo
c##################################################
      OPEN(54, file='temp-chequeoCONT')
c      OPEN(74, file='temp-chequeoCONT-MainChain-Ring')
      do i=1,icont
        do j=1,icont
          if (i.ne.j) then
          if (s(i,j).le.cutoff) then
            write(54,114) AAtom(i),AAid(i),CHid(i),npdb(i),
     &AAtom(j),AAid(j),CHid(j),npdb(j),w(i,j),s(i,j)
c            if (w(i,j).ne.0.001d0) then
c                write(74,114) AAtom(i),AAid(i),CHid(i),npdb(i),
c     &AAtom(j),AAid(j),CHid(j),npdb(j),w(i,j),s(i,j)
c            endif
          endif
          endif
        enddo
      enddo
      CLOSE(54)
c      CLOSE(74)

c|   7  N   PRO A  11      14.368  30.118  27.192  1.00 25.63
c|4075  CB  ASN D 124       6.233  49.107  67.149  1.00 38.22
113   format(I4,2x,A3,A1,A3,1x,A1,I4,4x,F8.3,F8.3,F8.3,1x,F5.2,F6.2)
114   format(2(A3,1x,A3,1x,A1,I4,1x),F5.3,1x,F5.2)
123   format(A6,I5,1x,A4,A1,A3,1x,A1,I4,4x,F8.3,F8.3,F8.3,1x,F5.2,F6.2)
c123   FORMAT(A6,I5,1x,A4,A1,A3,1x,A1,I4,A1,3x,F8.3,F8.3,F8.3,F6.2,F6.2)
c|ATOM     30  N   VAL A  14      -0.400 -10.074 -15.073  1.00 13.27
124   format(A3)

      OPEN(53,file='dim-anmh')
      write(53,*) 'icontTOTAL;icontTOTAL*3;icontMC;icontGR;icontHET'
         write(53,*) icont,icont*3,icontMC,icontGR,icontHET
      CLOSE(53)


C### COMIENZA CALCULO DE MODOS NORMALES ###

****************************************************
* Calculo de los elementos triangulares superiores *
*  de la matriz Hess 3Nx3N.- Aclaracion: i.NE.j    *
****************************************************      
      do i=1,icont*3
         do j=1,icont*3
            hess(i,j)=0.0d0
         enddo
      enddo
      do i=1,icont
        do j=1,icont
          if(i.ne.j.and.i.lt.j) then
            if(s(i,j).le.cutoff) then
              do m=0,2
                do k=0,2
                  If(k.eq.m) then
                  hess(i*3-2+m,j*3-2+k)=w(i,j)*(r(j,k+1)-r(i,m+1))**2
     $/s(i,j)**2
                             else
                  hess(i*3-2+m,j*3-2+k)=-w(i,j)*((r(j,m+1)-r(i,m+1))
     $*(r(j,k+1)-r(i,k+1)))/s(i,j)**2
                  endif
                enddo
              enddo
            endif
          endif
        enddo
      enddo

*************************************************
* Calculo de los elementos fuera de la diagonal,* 
*    dentro de las submatrices con i=j.         *
*************************************************
      do i=1,icont
        do j=1,icont
          if(i.ne.j) then
            if(s(i,j).le.cutoff) then            
              do m=0,2
                do k=0,2
                  if(k.ne.m) then
                    hess(i*3-2+m,i*3-2+k)=hess(i*3-2+m,i*3-2+k)
     $+w(i,j)*((r(j,m+1)-r(i,m+1))*(r(j,k+1)-r(i,k+1)))/s(i,j)**2
                  endif
                enddo
              enddo
            endif            
          endif
        enddo
      enddo

**************************************************
* Calculo de los elementos diagonales dentro     *
*         de las submatrices i=j.-               *
**************************************************
      do i=1,icont
        do j=1,icont
          if(i.ne.j) then
            if(s(i,j).le.cutoff) then
              do k=0,2
                  hess(i*3-2+k,i*3-2+k)=hess(i*3-2+k,i*3-2+k)
     $+w(i,j)*((r(j,k+1)-r(i,k+1))*(r(j,k+1)-r(i,k+1)))/s(i,j)**2
              enddo
            endif
          endif
        enddo
      enddo

************************************************
* Calculo de los elementos transpuestos fuera * 
*   de la diagonal de las submatrices i.ne.j   *
************************************************

      do i=1,icont
        do j=1,icont
          if(i.ne.j.and.i.gt.j) then
            if(s(i,j).le.cutoff) then
              do m=0,2
                do k=0,2
                  if(k.eq.m) then
                    hess(i*3-2+m,j*3-2+k)=-hess(j*3-2+k,i*3-2+m)
                        else
**********************************************************
* Calculo de los elem. transp. de las submatrices i.ne.j *
*     fuera de la diagonal de esa submatriz              *
**********************************************************  
                    hess(i*3-2+m,j*3-2+k)=hess(j*3-2+k,i*3-2+m)
                  endif
                enddo
              enddo
            endif
          endif
        enddo
      enddo
**********************     
**********************
 
      ifail=0
      do i=1,icont*3
        do j=1,icont*3
           diag(i,j)=0.0d0
        enddo
      enddo

c      open(10,file='hessiano')
c      do i=1,icont*3
c        write(10,333) (hess(i,j),j=1,icont*3)
c      enddo
c      close(10)
 
      call f02abf(hess,icontmax*3,icont*3,freq,diag,icontmax*3,e,ifail)

      open(10,file='modos')
      do i=1,icont*3
        write(10,333) (diag(i,j),j=7,icont*3)
      enddo
      close(10)


*********************
*   CORRELACIONES   *
*********************

      do i=1,icont*3
         do j=1,icont*3
             corr(i,j)=0.0d0
         enddo
      enddo

      do i=1,icont*3
         do j=1,icont*3
           do k=7,icont*3
               corr(i,j)=corr(i,j)+(diag(i,k)*diag(j,k))/freq(k)
           enddo
        enddo
      enddo

      open(11,file='cheq-correlac')
      read(11,*) cheqcorrelac
       if (cheqcorrelac.eq.'yes_correlac') then
         open(10,file='correlac')
         do i=1,icont*3
          write(10,339) (corr(i,j),j=1,icont*3)
         enddo 
         close(10)
       endif
      close(11)

      open(11,file='cheq-colect')
      read(11,*) cheqcolect
       if (cheqcolect.eq.'yes_colect') then
        open(10,file='colect')
         do j=7,3*icont     
           do i=1,3*icont
             x(i)=diag(i,j)
           enddo
           call colec(x,icont,col)
           write(10,*) j,col,freq(j)
         enddo
        close(10)
       endif
      close(11)

      SUMCORR=0.0d0
      SUMBFAC=0.0d0
      ii=1
      do i=1,icont*3,3
         bfacteor(ii)=(corr(i,i)+corr(i+1,i+1)+corr(i+2,i+2))
         ii=ii+1
      enddo 

      do i=1,icont
          SUMCORR = SUMCORR+bfacteor(i)
      enddo

      do i=1,icont
          SUMBFAC = SUMBFAC+bfactor(i) 
      enddo

      open(29,file='bfacteorico')
      do i=1,icont
          write(29,778) i,bfacteor(i)
      enddo
      close(29)
778   FORMAT(I5,1x,ES17.9E3)     
*******************************************************
*  CORRELACION ENTRE bfactor teorico y experimental   *
*******************************************************
*     AA es la suma de los (Xi-Xprom)(Yi-Yprom)
*     AB es la suma de los (Xi-Xprom) al cuadrado
*     AC es la suma de los (Yi-Yprom) al cuadrado

      bftprom=0.0d0
      bfeprom=0.0d0

      bftprom=SUMCORR/dfloat(icont)
      bfeprom=SUMBFAC/dfloat(icont)

      AA=0.0d0
      AB=0.0d0
      AC=0.0d0

      do i=1,icont
           AA = AA + ((bfactor(i)-bfeprom)*(bfacteor(i)-bftprom))
           AB = AB + (bfactor(i)-bfeprom)**2
           AC = AC + (bfacteor(i)-bftprom)**2
      enddo

      CORSTD= AA/(dsqrt(AB)*dsqrt(AC))   
   
      beta=SUMBFAC/SUMCORR
      if (beta.eq.0.0d0) then
        beta=1.0d0
        CORSTD=0.0d0
      endif

      open(10,file='freq')
c        write(10,*)'beta=',SUMBFAC/SUMCORR
        write(10,*)'beta=',beta
      do i=1,icont*3
        write(10,*) freq(i)
      enddo
      close(10)

      open(10,file='bfactor')
      do i=1,icont
c          write(10,444)i,bfacteor(i)*SUMBFAC/SUMCORR
c     $,bfactor(i),bfacteor(i)*CORSTD,SUMBFAC/SUMCORR,CORSTD
          write(10,444)i,bfacteor(i)*beta
     $,bfactor(i),bfacteor(i)*CORSTD,beta,CORSTD
      enddo
      close(10)
c678   FORMAT(ES12.4E3,1x,F5.3)

cc## Si quiero que el pdb random tenga un bfactor teorico escalado en la columna
cc      open(53,file='temp-coord2')
cc      do i=1,icont
cc         write(53,123) priAT,nlin(i),AAtom(i),AAocu(i),AAid(i),CHid(i)
cc     &,npdb(i),r(i,1),r(i,2),r(i,3),xx(i),bfacteor(i)*SUMBFAC/SUMCORR
cc      enddo
cc         write(53,124) termi
cc         write(53,124) endi
cc      close(53)

c      open(10,file='bfactorcorr')
c      write(10,335) CORSTD
c      close(10)


222   FORMAT(3(1x,F7.3),1x,F6.2,1x,F6.2) 
333   FORMAT(10000(1x,ES12.4E3))
c333   FORMAT(10000(1x,F9.6)) 
c|  -1.0485116050807597
c| -1.0485E+000
339   FORMAT(10000(1x,ES12.4E3))
444   FORMAT(I4,2(1x,F7.2),1x,F8.2,2(1x,ES12.4E3))
445   FORMAT(I5,1x,F8.2)
555   FORMAT(3(1x,I3))
776   FORMAT(2(1x,I4))
775   FORMAT(I6)
774   FORMAT(I2)
c778   FORMAT(I4,1x,I4,A)
c|  27  125   27  125  2.09
335   FORMAT(10000(1x,F14.11))
      return
      end

****************************************
* Subrutina de calculo de colectividad *
****************************************

      SUBROUTINE colec(x,n,col)

      IMPLICIT NONE
      INTEGER n,i,ii
      DOUBLE PRECISION x(3*n), alfa, col, peq
      DOUBLE PRECISION r2(n) 
      PARAMETER (peq = 1.e-30)


      alfa=0.0d0

      ii=1
      do i=1,n*3,3
        r2(ii)=x(i)**2+x(i+1)**2+x(i+2)**2
        ii=ii+1
      enddo
    
      do i=1,3*n
          alfa=alfa+x(i)*x(i)
      enddo


      alfa=1.0d0/alfa
      col=0.0d0
      do i=1,n
          col=col+alfa*r2(i)*log(alfa*r2(i))
      enddo

      col=exp(-col)/dfloat(n)

      RETURN
      END


C      do i=1,icont*3,3
C         write(17,333) corr(i,i)-
C     $((corr(i,i)+corr(i+1,i+1)+corr(i+2,i+2)))/3
C         write(18,333) corr(i+1,i+1)-((corr(i,i)+corr(i+1,i+1)+
C     $corr(i+2,i+2)))/3
C         write(19,333) corr(i+2,i+2)-((corr(i,i)+corr(i+1,i+1)+
C     $corr(i+2,i+2)))/3
C         write(20,333) (corr(i,i)+corr(i+1,i+1)+corr(i+2,i+2))
C      enddo 

c      close(50,STATUS='DELETE')
     
c222   FORMAT(3(1x,F7.3),1x,F5.2) 
c333   FORMAT(1000(1x,F7.3)) 

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

*---------------------------------------------------
      SUBROUTINE ORDEN(NELEM,AG1,AG2,AG11,AG3,AG4,AG5,
     &AG6,AG9,AG0)
c AG1=nlin(i)  AG2=AAtom(i)  AG11=AAocu(i)  AG3=AAid(i)  AG4=CHid(i)  AG5=npdb(i) 
c AG6=r1(i,1)  AG7=r2(i,2)  AG8=r3(i,3)  AG9=xx  AG0=bfactor(i)
* --------------------------------------------------
* SUBROUTINE ORDEN(NELEM,ARREG)
* ORDENACION POR BURBUJA ("buble sort") de un arreglo
* unidimensional, de menor a mayor.
*
* NELEM = Número de elementos del arreglo
* ARREG = Arreglo unidimensional a ordenar
* --------------------------------------------------
      IMPLICIT NONE
      integer icontmax
      parameter(icontmax=3000)
      INTEGER NELEM
      INTEGER AG1(icontmax)
      CHARACTER*3 AG2(icontmax)
      CHARACTER*1 AG11(icontmax)
      CHARACTER*3 AG3(icontmax)
      CHARACTER*1 AG4(icontmax)
      INTEGER AG5(icontmax)
      DOUBLE PRECISION AG6(icontmax,3)
      DOUBLE PRECISION AG9(icontmax),AG0(icontmax)
* --------------------------------------------------
      INTEGER I,J
      INTEGER AUX
      CHARACTER*4 AUX2
      CHARACTER*1 AUX11
      CHARACTER*3 AUX3
      CHARACTER*1 AUX4
      INTEGER AUX5
      DOUBLE PRECISION AUX6,AUX7,AUX8
      DOUBLE PRECISION AUX9,AUX0
* --------------------------------------------------
      IF (NELEM.LT.2) RETURN
      DO I=1,NELEM-1
       DO J=1,NELEM-I
        IF (AG1(J).GT.AG1(J+1)) THEN
          AUX = AG1(J)
          AG1(J) = AG1(J+1)
          AG1(J+1) = AUX
C-----------------------
          AUX2 = AG2(J)
          AG2(J) = AG2(J+1)
          AG2(J+1) = AUX2
C-----------------------
          AUX11 = AG11(J)
          AG11(J) = AG11(J+1)
          AG11(J+1) = AUX11
C-----------------------
          AUX3 = AG3(J)
          AG3(J) = AG3(J+1)
          AG3(J+1) = AUX3
C-----------------------
          AUX4 = AG4(J)
          AG4(J) = AG4(J+1)
          AG4(J+1) = AUX4
C-----------------------
          AUX5 = AG5(J)
          AG5(J) = AG5(J+1)
          AG5(J+1) = AUX5
C-----------------------
          AUX6 = AG6(J,1)
          AG6(J,1) = AG6(J+1,1)
          AG6(J+1,1) = AUX6
C-----------------------
          AUX7 = AG6(J,2)
          AG6(J,2) = AG6(J+1,2)
          AG6(J+1,2) = AUX7
C-----------------------
          AUX8 = AG6(J,3)
          AG6(J,3) = AG6(J+1,3)
          AG6(J+1,3) = AUX8
C-----------------------
          AUX9 = AG9(J)
          AG9(J) = AG9(J+1)
          AG9(J+1) = AUX9
C-----------------------
          AUX0 = AG0(J)
          AG0(J) = AG0(J+1)
          AG0(J+1) = AUX0
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END




c-------------------------------------------------------      
      SUBROUTINE f02abf (A, IA, N, R, V, IV, E, IFAIL)
C
C  diagonalizes a real, symmetric matrix
C  A(IA,N) : matrix to be diagonalized,
C    IA      its leading (first) dimension
C    N     : order of the matrix (.LE. IA!!!)
C             (only the NxN part of the matrix is actually needed)
C  R(N)    : on output, contains the eigenvalues (in ascending order)
C  V(IV,N) : on output, contains the eigen vectors (in columns)
C    IV      its leading dimension
C   E(N)   : working array
C  IFAIL   : indicator for errors: 0 on output if everything worked...
C
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C
C     EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIX MATRIX
C     1ST AUGUST 1971
      INTEGER P01AAF, ISAVE, IFAIL, N, IA, IV
C$P 1
      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION TOL, XXXX, A(IA,N), R(N), V(IV,N), E(N), X02ADF,
     * X02AAF
      DATA SRNAME /8H F02ABF /
      ISAVE = IFAIL
      IFAIL = 1
      TOL = X02ADF(1.0D0)
      CALL F01AJF(N, TOL, A, IA, R, E, V, IV)
      DO 1 I=1,N
C     PRINT 100,IA,IV,R(I),E(I),(A(I,J),V(I,J),J=1,N)
  100 FORMAT(/,' AJF :IA,IV,R,E:',2I2,2(D15.8),/,' A',3D15.8,/,' D',
     X3D15.8)
   1  CONTINUE
      TOL = X02AAF(1.0D0)
      CALL F02AMF(N, TOL, R, E, V, IV, IFAIL)
      DO 2 I=1,N
C     PRINT 200,R(I),E(I),IV,IFAIL,(V(I,J),J=1,3)
  200 FORMAT(/,' AMF :R,E,IV,IFAIL=',2D15.8,2I3,/,' V',3D15.8)
   2  CONTINUE
      IF (IFAIL.NE.0) IFAIL = P01AAF(ISAVE,IFAIL,SRNAME)
      RETURN
      END

      SUBROUTINE F01AJF(N, ATOL, A, IA, D, E, Z, IZ)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C
C     TRED2
C     THIS SUBROUTINE REDUCES THE GIVEN LOWER TRIANGLE OF A
C     SYMMETRIC MATRIX, A, STORED IN THE ARRAY A(N,N), TO
C     TRIDIAGONAL FORM USING HOUSEHOLDERS REDUCTION. THE DIAGONAL
C     OF THE RESULT IS STORED IN THE ARRAY D(N) AND THE
C     SUB-DIAGONAL IN THE LAST N - 1 STORES OF THE ARRAY E(N)
C     (WITH THE ADDITIONAL ELEMENT E(1) = 0). THE TRANSFORMATION
C     MATRICES ARE ACCUMULATED IN THE ARRAY Z(N,N). THE ARRAY
C     A IS LEFT UNALTERED UNLESS THE ACTUAL PARAMETERS
C     CORRESPONDING TO A AND Z ARE IDENTICAL.
C     1ST AUGUST 1971
C
      INTEGER I, IA, II, IZ, J1, J, K, L, N
      DOUBLE PRECISION ATOL, F, G, H, HH, A(IA,N), D(N), E(N), Z(IZ,N)
      DO 40 I=1,N
         DO 20 J=1,I
            Z(I,J) = A(I,J)
   20    CONTINUE
   40 CONTINUE
      IF (N.EQ.1) GO TO 280
      DO 260 II=2,N
         I = N - II + 2
         L = I - 2
         F = Z(I,I-1)
         G = 0.0D0
         IF (L.EQ.0) GO TO 80
         DO 60 K=1,L
            G = G + Z(I,K)*Z(I,K)
   60    CONTINUE
   80    H = G + F*F
C     IF G IS TOO SMALL FOR ORTHOGONALITY TO BE
C     GUARANTEED THE TRANSFORMATION IS SKIPPED
         IF (G.GT.ATOL) GO TO 100
         E(I) = F
         H = 0.0D0
         GO TO 240
  100    L = L + 1
         G = DSQRT(H)
         IF (F.GE.0.0D0) G = -G
         E(I) = G
         H = H - F*G
         Z(I,I-1) = F - G
         F = 0.0D0
         DO 180 J=1,L
            Z(J,I) = Z(I,J)/H
            G = 0.0D0
C     FORM ELEMENT OF A*U
            DO 120 K=1,J
               G = G + Z(J,K)*Z(I,K)
  120       CONTINUE
            J1 = J + 1
            IF (J1.GT.L) GO TO 160
            DO 140 K=J1,L
               G = G + Z(K,J)*Z(I,K)
  140       CONTINUE
C     FORM ELEMENT OF P
  160       E(J) = G/H
            F = F + G*Z(J,I)
  180    CONTINUE
C     FORM K
         HH = F/(H+H)
C     FORM REDUCED A
         DO 220 J=1,L
            F = Z(I,J)
            G = E(J) - HH*F
            E(J) = G
            DO 200 K=1,J
               Z(J,K) = Z(J,K) - F*E(K) - G*Z(I,K)
  200       CONTINUE
  220    CONTINUE
  240    D(I) = H
  260 CONTINUE
  280 E(1) = 0.0D0
      D(1) = 0.0D0
C     ACCUMULATION OF TRANSFORMATION MATRICES
      DO 400 I=1,N
         L = I - 1
         IF (D(I).EQ.0.0D0) GO TO 360
         DO 340 J=1,L
            G = 0.0D0
            DO 300 K=1,L
               G = G + Z(I,K)*Z(K,J)
  300       CONTINUE
            DO 320 K=1,L
               Z(K,J) = Z(K,J) - G*Z(K,I)
  320       CONTINUE
  340    CONTINUE
  360    D(I) = Z(I,I)
         Z(I,I) = 1.0D0
         IF (L.EQ.0) GO TO 400
         DO 380 J=1,L
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  380    CONTINUE
  400 CONTINUE
      RETURN
      END


      SUBROUTINE F02AMF(N, ACHEPS, D, E, Z, IZ, IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C
C     TQL2
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS OF A
C     TRIDIAGONAL MATRIX, T, GIVEN WITH ITS DIAGONAL ELEMENTS IN
C     THE ARRAY D(N) AND ITS SUB-DIAGONAL ELEMENTS IN THE LAST N
C     - 1 STORES OF THE ARRAY E(N), USING QL TRANSFORMATIONS. THE
C     EIGENVALUES ARE OVERWRITTEN ON THE DIAGONAL ELEMENTS IN THE
C     ARRAY D IN ASCENDING ORDER. THE EIGENVECTORS ARE FORMED IN
C     THE ARRAY Z(N,N), OVERWRITING THE ACCUMULATED
C     TRANSFORMATIONS AS SUPPLIED BY THE SUBROUTINE F01AJF. THE
C     SUBROUTINE WILL FAIL IF ANY ONE EIGENVALUE TAKES MORE THAN 30
C     ITERATIONS.
C     1ST APRIL 1972
C
      INTEGER P01AAF, ISAVE, IFAIL, N, I, L, J, M, I1, M1, II, K, IZ
C$P 1
      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION B, F, H, ACHEPS, G, P, R, C, S, D(N), E(N), Z(IZ,
     *N)
      DATA SRNAME /8H F02AMF /
      ISAVE = IFAIL
      IF (N.EQ.1) GO TO 40
      DO 20 I=2,N
         E(I-1) = E(I)
   20 CONTINUE
   40 E(N) = 0.0D0
      B = 0.0D0
      F = 0.0D0
      DO 300 L=1,N
         J = 0
         H = ACHEPS*(DABS(D(L))+DABS(E(L)))
         IF (B.LT.H) B = H
C     LOOK FOR SMALL SUB-DIAG ELEMENT
         DO 60 M=L,N
            IF (DABS(E(M)).LE.B) GO TO 80
   60    CONTINUE
   80    IF (M.EQ.L) GO TO 280
  100    IF (J.EQ.30) GO TO 400
         J = J + 1
C     FORM SHIFT
         G = D(L)
         H = D(L+1) - G
         IF (DABS(H).GE.DABS(E(L))) GO TO 120
         P = H*0.5D0/E(L)
         R = DSQRT(P*P+1.0D0)
         H = P + R
         IF (P.LT.0.0D0) H = P - R
         D(L) = E(L)/H
         GO TO 140
  120    P = 2.0D0*E(L)/H
         R = DSQRT(P*P+1.0D0)
         D(L) = E(L)*P/(1.0D0+R)
  140    H = G - D(L)
         I1 = L + 1
         IF (I1.GT.N) GO TO 180
         DO 160 I=I1,N
            D(I) = D(I) - H
  160    CONTINUE
  180    F = F + H
C     QL TRANSFORMATION
         P = D(M)
         C = 1.0D0
         S = 0.0D0
         M1 = M - 1
         DO 260 II=L,M1
            I = M1 - II + L
            G = C*E(I)
            H = C*P
            IF (DABS(P).LT.DABS(E(I))) GO TO 200
            C = E(I)/P
            R = DSQRT(C*C+1.0D0)
            E(I+1) = S*P*R
            S = C/R
            C = 1.0D0/R
            GO TO 220
  200       C = P/E(I)
            R = DSQRT(C*C+1.0D0)
            E(I+1) = S*E(I)*R
            S = 1.0D0/R
            C = C/R
  220       P = C*D(I) - S*G
            D(I+1) = H + S*(C*G+S*D(I))
C     FORM VECTOR
            DO 240 K=1,N
               H = Z(K,I+1)
               Z(K,I+1) = S*Z(K,I) + C*H
               Z(K,I) = C*Z(K,I) - S*H
  240       CONTINUE
  260    CONTINUE
         E(L) = S*P
         D(L) = C*P
         IF (DABS(E(L)).GT.B) GO TO 100
  280    D(L) = D(L) + F
  300 CONTINUE
C     ORDER EIGENVALUES AND EIGENVECTORS
      DO 380 I=1,N
         K = I
         P = D(I)
         I1 = I + 1
         IF (I1.GT.N) GO TO 340
         DO 320 J=I1,N
            IF (D(J).GE.P) GO TO 320
            K = J
            P = D(J)
  320    CONTINUE
  340    IF (K.EQ.I) GO TO 380
         D(K) = D(I)
         D(I) = P
         DO 360 J=1,N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  360    CONTINUE
  380 CONTINUE
      IFAIL = 0
      RETURN
  400 IFAIL = P01AAF(ISAVE,1,SRNAME)
      RETURN
      END

      INTEGER FUNCTION P01AAF(IFAIL, ERROR, SRNAME)
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971
C     MARK 3 REVISED
C     MARK 4A REVISED, IER-45
C     MARK 4.5 REVISED
C     MARK 7 REVISED (DEC 1978)
C     RETURNS THE VALUE OF ERROR OR TERMINATES THE PROGRAM.
      INTEGER ERROR, IFAIL, NOUT
C$P 1
      DOUBLE PRECISION SRNAME
C     TEST IF NO ERROR DETECTED
      IF (ERROR.EQ.0) GO TO 20
C     DETERMINE OUTPUT UNIT FOR MESSAGE
      CALL X04AAF (0,NOUT)
C     TEST FOR SOFT FAILURE
      IF (MOD(IFAIL,10).EQ.1) GO TO 10
C     HARD FAILURE
      WRITE (NOUT,99999) SRNAME, ERROR
C     STOPPING MECHANISM MAY ALSO DIFFER
      STOP
C     SOFT FAIL
C     TEST IF ERROR MESSAGES SUPPRESSED
   10 IF (MOD(IFAIL/10,10).EQ.0) GO TO 20
      WRITE (NOUT,99999) SRNAME, ERROR
   20 P01AAF = ERROR
      RETURN
99999 FORMAT (1H0, 38HERROR DETECTED BY NAG LIBRARY ROUTINE , A8,
     * 11H - IFAIL = , I5//)
      END


      DOUBLEPRECISION FUNCTION X02AAF(X)
      DOUBLEPRECISION X,Z
C     NAG COPYRIGHT 1975
C     IBM DOUBLE PRECISION
C     MARK 4.5 RELEASE
C     * MACHEPS *
C     RETURNS THE VALUE MACHEPS WHERE MACHEPS IS THE SMALLEST POSITIVE
C     NUMBER SUCH THAT 1.0 + EPS > 1.0
C     THE X PARAMETER IS NOT USED
C     FOR IBM RISC 6000
      X02AAF = 2.0D0**(-52.0D0)
C     SET IN HEX FOR ACCURACY
C      DATA Z/Z/
C      X02AAF=Z
      RETURN
      END


      DOUBLEPRECISION FUNCTION X02ADF(X)
      DOUBLEPRECISION X,Z
C     NAG COPYRIGHT 1975
C     MARK 4.5 RELEASE
C     IBM DOUBLE PRECISION
C     * MINTOEPS *
C     RETURNS THE RATIO OF THE SMALLEST POSITIVE REAL FLOATING-
C     POINT NUMBER REPRESENTABLE ON THE COMPUTER TO MACHEPS
C     FOR IBM RISC 6000
      X02ADF = 2.0D0**(-448.0)
C     FOR IBM 370
C     X02ADF = 16.0D0**(-52)
C     SET IN HEX FOR ACCURACY
C      DATA Z/Z0D000000/
C      X02ADF=Z
      RETURN
      END


      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     *** NOTE ***
C     THIS ROUTINE ASSUMES THAT THE VALUE OF NERR1 IS SAVED
C     BETWEEN CALLS.  IN SOME IMPLEMENTATIONS IT MAY BE
C     NECESSARY TO STORE NERR1 IN A LABELLED COMMON
C     BLOCK /AX04AA/ TO ACHIEVE THIS.
C
C     .. SCALAR ARGUMENTS ..
      INTEGER I, NERR
C     ..
C     .. LOCAL SCALARS ..
      INTEGER NERR1
C     ..
      DATA NERR1 /6/
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
