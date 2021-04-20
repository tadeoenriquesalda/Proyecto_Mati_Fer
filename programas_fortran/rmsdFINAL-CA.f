C Programa para calcular el promedio termalizado 
C de el RMSD  


      IMPLICIT NONE      
 
      INTEGER i,j
      INTEGER i2,j2,ii,kk2,k1,k2,k,kk
      INTEGER icont,icont2,icontmax
      parameter (icontmax=10000)
      DOUBLE PRECISION cm1(3),cm2(3) 
      DOUBLE PRECISION r1(icontmax,3) 
      DOUBLE PRECISION r2(icontmax,3) 
      double precision r1ini(icontmax,3)
      double precision r2ini(icontmax,3)
      double precision xxyx,xxyy,xxyz,xyyx,xyyy,xzyx,xzyy,xzyz,xyyz 
      DOUBLE PRECISION dist,dist1,dist2,distdiff
      DOUBLE PRECISION c(4,4),d(4) 
      DOUBLE PRECISION rmsd,x,x2,rmsdp 
      DOUBLE PRECISION xrot,yrot,zrot 
      DOUBLE PRECISION rot(3,3),q(4),v(4,4) 
      DOUBLE PRECISION work1(4),work2(4)
      DOUBLE PRECISION xx,bfactor(icontmax),xx2,bfactor2(icontmax)
      character*30 cardA(icontmax),cardB(icontmax)
      character*14 cardAX(icontmax),cardBX(icontmax)
C****LA Kb ESTA EN ANGSTROMS CUADRADO

c      icont=0
      OPEN(21, FILE='dimAA')
      read(21,*) icont
      CLOSE(21)

      OPEN(21, FILE='tempAA')
      do i=1,icont
      READ(21,111) cardA(i),r1ini(i,1),r1ini(i,2)
     $,r1ini(i,3),xx,bfactor(i),cardAX(i)
      enddo
      close(21)

      OPEN(21, FILE='dimBB')
      read(21,*) icont2
      CLOSE(21)

      OPEN(21, FILE='tempBB')
c      icont2=0
      do i=1,icont2
      READ(21,111) cardB(i),r2ini(i,1),r2ini(i,2)
     $,r2ini(i,3),xx2,bfactor2(i),cardBX(i) 
      enddo
      CLOSE(21)

      if(icont.ne.icont2) write(99,*)"distinto nº de residuos"

C*****************************************************
C Inicio del bucle para calcular todas las combinaciones 
C de a pares
C*****************************************************
c      OPEN (22,FILE='rmsd.out')

C*****************************************************
C Calculo del CM de los atomos de la estructura 1
C*****************************************************
      do j=1,3
        cm1(j)=0.0d0
      enddo
      do i=1,icont
        do j=1,3
           cm1(j)=cm1(j)+r1ini(i,j)/dfloat(icont)
        enddo
      enddo
C*****************************************************
C Traslacion de los atomos de la estructura 1 
C*****************************************************
      open(13,file='pdb.modSL')
      do j=1,icont
         do kk=1,3
            r1(j,kk)=r1ini(j,kk)-cm1(kk)
         enddo
        write(13,111) cardA(j),(r1(j,kk),kk=1,3),xx,bfactor(j),cardAX(j)
      enddo 
      CLOSE(13)
C*****************************************************
C Calculo del CM de los atomos de la estructura 2
C*****************************************************
      do j=1,3
        cm2(j)=0.0d0
      enddo
      do i=1,icont
        do j=1,3
           cm2(j)=cm2(j)+r2ini(i,j)/dfloat(icont)
        enddo
      enddo
C*****************************************************
C Traslacion de los atomos de la estructura 2 
C*****************************************************
      do j=1,icont
        do kk=1,3
          r2(j,kk)=r2ini(j,kk)-cm2(kk)
        enddo
      enddo 
C***************************************************
C Rotacion 
C*****************************************************
      xxyx=0.0d0
      xxyy=0.0d0
      xxyz=0.0d0
      xyyx=0.0d0
      xyyy=0.0d0
      xyyz=0.0d0
      xzyx=0.0d0
      xzyy=0.0d0
      xzyz=0.0d0
      do j = 1,icont
         xxyx = xxyx + r1(j,1)*r2(j,1)
         xxyy = xxyy + r1(j,2)*r2(j,1)
         xxyz = xxyz + r1(j,3)*r2(j,1)
         xyyx = xyyx + r1(j,1)*r2(j,2)
         xyyy = xyyy + r1(j,2)*r2(j,2)
         xyyz = xyyz + r1(j,3)*r2(j,2)
         xzyx = xzyx + r1(j,1)*r2(j,3)
         xzyy = xzyy + r1(j,2)*r2(j,3)
         xzyz = xzyz + r1(j,3)*r2(j,3)
      end do
      c(1,1) = xxyx + xyyy + xzyz
      c(1,2) = xzyy - xyyz
      c(2,2) = xxyx - xyyy - xzyz
      c(1,3) = xxyz - xzyx
      c(2,3) = xxyy + xyyx
      c(3,3) = xyyy - xzyz - xxyx
      c(1,4) = xyyx - xxyy
      c(2,4) = xzyx + xxyz
      c(3,4) = xyyz + xzyy
      c(4,4) = xzyz - xxyx - xyyy
c
c     diagonalize the quadratic form matrix
c
      call jacobi (4,4,c,d,v,work1,work2)
c
c     extract the desired quaternion
c
      q(1) = v(1,4)
      q(2) = v(2,4)
      q(3) = v(3,4)
      q(4) = v(4,4)
c
c     assemble rotation matrix that superimposes the molecules
c
      rot(1,1) = q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
      rot(2,1) = 2.0d0 * (q(2) * q(3) - q(1) * q(4))
      rot(3,1) = 2.0d0 * (q(2) * q(4) + q(1) * q(3))
      rot(1,2) = 2.0d0 * (q(3) * q(2) + q(1) * q(4))
      rot(2,2) = q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
      rot(3,2) = 2.0d0 * (q(3) * q(4) - q(1) * q(2))
      rot(1,3) = 2.0d0 * (q(4) * q(2) - q(1) * q(3))
      rot(2,3) = 2.0d0 * (q(4) * q(3) + q(1) * q(2))
      rot(3,3) = q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2
c
c     rotate second molecule to best fit with first molecule
c
      open(14,file='pdb.modCL')
      do j = 1,icont 
         xrot = r2(j,1)*rot(1,1) + r2(j,2)*rot(1,2) + r2(j,3)*rot(1,3)
         yrot = r2(j,1)*rot(2,1) + r2(j,2)*rot(2,2) + r2(j,3)*rot(2,3)
         zrot = r2(j,1)*rot(3,1) + r2(j,2)*rot(3,2) + r2(j,3)*rot(3,3)
         r2(j,1) = xrot
         r2(j,2) = yrot
         r2(j,3) = zrot
      write(14,111) cardB(j),(r2(j,kk),kk=1,3),xx2,bfactor2(j),cardBX(j)
      enddo
      CLOSE(14)
C***************************************************
C Calculo de la diferencia de posiciones en las 2 estructuras 
C***************************************************
      open(50,file='rmsd') 
      open(51,file='vecdif')
c *RMSD total      rmsd=0.0d0
      rmsd=0.0d0
      do j=1,icont
         dist=0.0d0
c *RMSD posicional         rmsdp=0.0d0
         do kk=1,3
            x=r1(j,kk)-r2(j,kk)
            x2=x*x
            dist=dist+x2
         write(51,224) x
         enddo
****Si quiero RMSD posisional uso esta parte***********
c         rmsdp=dsqrt(dist)
c         write(50,124) j,rmsdp
c       enddo
*******************************************************

***Si quiero RMSD total uso esta parte*****************
         rmsd=rmsd+dist/dfloat(icont)
      enddo
      rmsd=dsqrt(rmsd)
      write(50,776) rmsd
******************************************************
      close(50)
      close(51)

776   FORMAT(F7.4)
111   FORMAT(A30,3(F8.3),F6.2,F6.2,A14)
112   FORMAT(3(F8.3))
113   FORMAT(A30)
114   FORMAT(F6.2,F6.2)
123   FORMAT(I4,1X,F18.10)
222   FORMAT(3(F8.3))
333   FORMAT(10000(1x,F7.4))
555   FORMAT(10000(1x,F9.4))
777   FORMAT(F18.10)
888   FORMAT(F18.10,1X,I3,1X,I3)
889   FORMAT(F18.10,1X,F18.10,1X,F18.10,100(1X,I3))
444   FORMAT(1x,F8.4)
666   FORMAT(10000(1x,F16.10))
124   FORMAT(10000(1x,I3,1x,F11.8))
224   FORMAT(ES12.4E3)
      STOP
      END

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine jacobi  --  jacobi matrix diagonalization  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "jacobi" performs a matrix diagonalization of a real
c     symmetric matrix by the method of Jacobi rotations
c
c     variables and parameters:
c
c     n     logical dimension of the matrix to be diagonalized
c     np    physical dimension of the matrix storage area
c     a     input with the matrix to be diagonalized; only
c              the upper triangle and diagonal are required
c     d     returned with the eigenvalues in ascending order
c     v     returned with the eigenvectors of the matrix
c     b     temporary work vector
c     z     temporary work vector
c
c
      subroutine jacobi (n,np,a,d,v,b,z)
      implicit none
      integer i,j,k
      integer n,np,ip,iq
      integer nrot,maxrot
      real*8 sm,tresh,s,c,t
      real*8 theta,tau,h,g,p
      real*8 d(np),b(np),z(np)
      real*8 a(np,np),v(np,np)
c
c
c     setup and initialization
c
      maxrot = 100
      nrot = 0
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0d0
         end do
         v(ip,ip) = 1.0d0
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0d0
      end do
c
c     perform the jacobi rotations
c
      do i = 1, maxrot
         sm = 0.0d0
         do ip = 1, n-1
            do iq = ip+1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm .eq. 0.0d0)  goto 10
         if (i .lt. 4) then
            tresh = 0.2d0*sm / n**2
         else
            tresh = 0.0d0
         end if
         do ip = 1, n-1
            do iq = ip+1, n
               g = 100.0d0 * abs(a(ip,iq))
               if (i.gt.4 .and. abs(d(ip))+g.eq.abs(d(ip))
     &                    .and. abs(d(iq))+g.eq.abs(d(iq))) then
                  a(ip,iq) = 0.0d0
               else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h)+g .eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5d0*h / a(ip,iq)
                     t = 1.0d0 / (abs(theta)+sqrt(1.0d0+theta**2))
                     if (theta .lt. 0.0d0)  t = -t
                  end if
                  c = 1.0d0 / sqrt(1.0d0+t**2)
                  s = t * c
                  tau = s / (1.0d0+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0d0
                  do j = 1, ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = ip+1, iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = iq+1, n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h+g*tau)
                     a(iq,j) = h + s*(g-h*tau)
                  end do
                  do j = 1, n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h+g*tau)
                     v(j,iq) = h + s*(g-h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
         end do
      end do
c
c     print warning if not converged
c
   10 continue
      if (nrot .eq. maxrot) then
         write (99,20)
   20    format (/,' JACOBI  --  Matrix Diagonalization not Converged')
      end if
c
c     sort the eigenvalues and vectors
c
      do i = 1, n-1
         k = i
         p = d(i)
         do j = i+1, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            end if
         end do
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
            end do
         end if
      end do
      return
      end

