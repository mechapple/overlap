      program GCCM_GenerateCrystal
      implicit none

!     ####################################################################
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!     ####################################################################
!
!     Written by: Matthew P. Kroonblawd
!                 Department of Chemistry
!                 University of Missouri - Columbia
!                 May 2016
!
!     This program is the second part of the Generalized Crystal-Cutting
!     Method (GCCM) that can be used to generate crystals that are oriented
!     with an arbitrary direction S parallel to the normal of a face of
!     the simulation cell. This program takes cell edge vectors found
!     by the program GCCM_SingleCrystal.f90 that are copied into the file
!     chosen.solution.txt and the original unit cell and performs the
!     necessary rotations and cuts to populate a new simulation cell.
!     The direction x1 x x2 || S is aligned to be exactly parallel to  
!     the +z axis.
!
!     ####################################################################

      double precision :: a, b, c, alpha, beta, gamma, volunit
      double precision :: lx1, lx2, lx3, theta, phi, psi
      double precision :: av(3), bv(3), cv(3)
      double precision, allocatable :: ro(:,:,:)
      double precision, allocatable :: mass(:,:)
      double precision, allocatable :: molmass(:)
      double precision :: xcm, ycm, zcm, unitmass
      double precision :: x1(3), x2(3), x3(3)
      double precision :: x1p(3), x2p(3), x3p(3)
      double precision :: volcell
      double precision :: theta1, theta2, theta3
      double precision :: cost1, sint1
      double precision :: cost2, sint2
      double precision :: cost3, sint3
      double precision :: oldvol, newvol
      double precision :: hinv(3,3), Rx(3,3), Ry(3,3), Rz(3,3), Rall(3,3), Rtemp(3,3)
      double precision, allocatable :: r(:,:,:), rp(:,:,:), rbox(:,:,:)
      double precision, allocatable :: rcmo(:,:), rcm(:,:), rcmp(:,:)
      double precision :: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
      double precision :: ri(3), rj(3), rijmin(3)
      double precision :: rijmag
      double precision :: cellcom(3), geocenter(3)
      double precision :: totalmass
      double precision, parameter :: pi=4.0d0*datan(1.0d0)
      double precision, parameter :: zero=1.0d-12

      integer :: molperunit
      integer, allocatable :: atpermol(:)
      integer :: i, j, k, l, index
      integer :: m, mmin, mmax, allmin, allmax
      integer :: n, nmin, nmax
      integer :: p, pmin, pmax
      integer :: x3scale, type_atom, mol_id
      integer :: mnp(3,3)
      integer :: nmolEST, nmol, natom, iter, newnmol, newnatom
      integer, allocatable :: moltype(:), duplicates(:)

      character*2, allocatable :: el(:,:)
      character*1 :: char

!     Read control file

      open(30,action='read',file='control.txt')

      read(30,*)
      read(30,*) char
      read(30,*) char
      read(30,*) char, mmin, mmax
      read(30,*) char, nmin, nmax
      read(30,*) char, pmin, pmax
      read(30,*) char, x3scale
      read(30,*) char
      read(30,*)
      read(30,*)
      read(30,*) char
      read(30,*) char
      read(30,*) char
      read(30,*)
      read(30,*)
      read(30,*) char, a            !get unit cell lattice parameters
      read(30,*) char, b
      read(30,*) char, c
      read(30,*) char, alpha
      read(30,*) char, beta
      read(30,*) char, gamma
      read(30,*) char, molperunit

      allocate (atpermol(molperunit))

      read(30,*) char, (atpermol(i),i=1,molperunit)
      close(30)

      allocate (ro(molperunit,maxval(atpermol,DIM=1),3))
      allocate (el(molperunit,maxval(atpermol,DIM=1)))

!     Compute lattice vectors from lattice parameters
!
!     Unit cell H-matrix is H = (av, bv, cv) 
!     where av, bv, and cv are column vectors
!     Unit cell volume is volunit = Tr[H]

      alpha=alpha*pi/180.0d0
      beta=beta*pi/180.0d0
      gamma=gamma*pi/180.0d0

      av(1)=a               !x component
      av(2)=0.0d0           !y component
      av(3)=0.0d0           !z component

      bv(1)=b*dcos(gamma)   !same format as above
      bv(2)=b*dsin(gamma)
      bv(3)=0.0d0

      cv(1)=c*dcos(beta)
      cv(2)=c*(dcos(alpha)-dcos(beta)*dcos(gamma))/dsin(gamma)
      cv(3)=c*dsqrt(1.0d0-(dcos(alpha))**2-(dcos(beta))**2-(dcos(gamma))**2 &
                    +2.0d0*dcos(alpha)*dcos(beta)*dcos(gamma))/dsin(gamma)

      volunit=av(1)*bv(2)*cv(3)

!     Optional dump to get rotation matrices.
!     Make sure to uncomment all lines dealing with file 50. 

      open(50,file='rotation_matrices.txt')
      write(50,*) '#Unitcell H-Matrix'
      write(50,5000) av(1), bv(1), cv(1)
      write(50,5000) av(2), bv(2), cv(2)
      write(50,5000) av(3), bv(3), cv(3)
5000  format(3F24.10)

!     Get basis molecules/atoms for original unit cell.
!     This code assumes that the unit cell is oriented
!     such that a || x, b is in the xy plane, and c is
!     in the +z direction.

      open(30,file='unitcell.xyz')
      read(30,*) 
      read(30,*)

      do i=1,molperunit
       do j=1,atpermol(i)
        read(30,*) el(i,j), (ro(i,j,k),k=1,3)
       enddo
      enddo

      close(30)

!     Assign atom masses and get molecular mass.
!     This assumes that the original unit cell was molecule-ordered.
!     Currently only coded for H, C, N, and O.

      allocate (mass(molperunit,maxval(atpermol,DIM=1)))
      allocate (molmass(molperunit))

      molmass=0.0d0

      do i=1,molperunit
       do j=1,atpermol(i)
        if(el(i,j).eq.'h'.or.el(i,j).eq.'H') mass(i,j)=1.0079d0
        if(el(i,j).eq.'c'.or.el(i,j).eq.'C') mass(i,j)=12.011d0
        if(el(i,j).eq.'n'.or.el(i,j).eq.'N') mass(i,j)=14.007d0
        if(el(i,j).eq.'o'.or.el(i,j).eq.'O') mass(i,j)=15.999d0
        molmass(i)=molmass(i)+mass(i,j)
       enddo
      enddo

!     Place origin at unitcell center-of-mass

      xcm=0.0d0
      ycm=0.0d0
      zcm=0.0d0

      unitmass=0.0d0

      do i=1,molperunit
       unitmass=unitmass+molmass(i)
      enddo

      do i=1,molperunit
       do j=1,atpermol(i)
        xcm=xcm+mass(i,j)*ro(i,j,1)
        ycm=ycm+mass(i,j)*ro(i,j,2)
        zcm=zcm+mass(i,j)*ro(i,j,3)
       enddo
      enddo

      xcm=xcm/unitmass
      ycm=ycm/unitmass
      zcm=zcm/unitmass

      do i=1,molperunit
       do j=1,atpermol(i)
        ro(i,j,1)=ro(i,j,1)-xcm !??
        ro(i,j,2)=ro(i,j,2)-ycm
        ro(i,j,3)=ro(i,j,3)-zcm
       enddo
      enddo

!     Get x1, x2, and x3

      open(40,file='chosen.solution.txt')
      read(40,*)
      read(40,*) (mnp(1,i),i=1,3),(mnp(2,i),i=1,3),(mnp(3,i),i=1,3)

      do j=1,3
       x1(j)=dfloat(mnp(1,1))*av(j)+dfloat(mnp(1,2))*bv(j)+dfloat(mnp(1,3))*cv(j)
       x2(j)=dfloat(mnp(2,1))*av(j)+dfloat(mnp(2,2))*bv(j)+dfloat(mnp(2,3))*cv(j)
       x3(j)=dfloat(mnp(3,1))*av(j)+dfloat(mnp(3,2))*bv(j)+dfloat(mnp(3,3))*cv(j)
      enddo

      write(*,9092)
      write(*,9093) x1(1), x2(1), x3(1)
      write(*,9093) x1(2), x2(2), x3(2)
      write(*,9093) x1(3), x2(3), x3(3)

!     Determine rotations to align x1 -> x, x2 in xy plane, x3 in +z

      !rotate about z by theta1

      theta1=-datan(x1(2)/x1(1))

      cost1=dcos(theta1)
      sint1=dsin(theta1)

      x1p(1)=x1(1)*cost1-x1(2)*sint1
      x1p(2)=x1(1)*sint1+x1(2)*cost1
      x1p(3)=x1(3)

      if(x1p(1).lt.0.0d0) then
       theta1=theta1+pi
       cost1=dcos(theta1)
       sint1=dsin(theta1)
       x1p(1)=x1(1)*cost1-x1(2)*sint1
       x1p(2)=x1(1)*sint1+x1(2)*cost1
       x1p(3)=x1(3)
      endif

      x2p(1)=x2(1)*cost1-x2(2)*sint1
      x2p(2)=x2(1)*sint1+x2(2)*cost1
      x2p(3)=x2(3)

      x3p(1)=x3(1)*cost1-x3(2)*sint1
      x3p(2)=x3(1)*sint1+x3(2)*cost1
      x3p(3)=x3(3)

      do i=1,3
       x1(i)=x1p(i)
       x2(i)=x2p(i)
       x3(i)=x3p(i)
      enddo

      write(50,*) '#Rotation Rz(theta1) about z'
      Rz(1,:) = (/  cost1, -sint1, 0.0d0 /)
      Rz(2,:) = (/  sint1,  cost1, 0.0d0 /)
      Rz(3,:) = (/  0.0d0,  0.0d0, 1.0d0 /)
      
      do i=1,3
         write(50,5000) ( Rz(i,j), j=1,3 )
      enddo

      !now rotate about y by theta2

      theta2=datan(x1(3)/x1(1))

      cost2=dcos(theta2)
      sint2=dsin(theta2)

      x1p(1)=x1(1)*cost2+x1(3)*sint2
      x1p(2)=x1(2)
      x1p(3)=-x1(1)*sint2+x1(3)*cost2

      x2p(1)=x2(1)*cost2+x2(3)*sint2
      x2p(2)=x2(2)
      x2p(3)=-x2(1)*sint2+x2(3)*cost2

      x3p(1)=x3(1)*cost2+x3(3)*sint2
      x3p(2)=x3(2)
      x3p(3)=-x3(1)*sint2+x3(3)*cost2

      do i=1,3
       x1(i)=x1p(i)
       x2(i)=x2p(i)
       x3(i)=x3p(i)
      enddo

      write(50,*) '#Rotation Ry(theta2) about y'
      Ry(1,:) = (/  cost2,  0.0d0, sint2 /)
      Ry(2,:) = (/  0.0d0,  1.0d0, 0.0d0 /)
      Ry(3,:) = (/ -sint2,  0.0d0, cost2 /)
      
      do i=1,3
         write(50,5000) ( Ry(i,j), j=1,3 )
      enddo

      !now rotate about x by theta3

      theta3=-datan(x2(3)/x2(2))

      cost3=dcos(theta3)
      sint3=dsin(theta3)

      x3p(1)=x3(1)
      x3p(2)=x3(2)*cost3-x3(3)*sint3
      x3p(3)=x3(2)*sint3+x3(3)*cost3

      if(x3p(3).lt.0.0d0) then
       theta3=theta3+pi
       cost3=dcos(theta3)
       sint3=dsin(theta3)
       x3p(1)=x3(1)
       x3p(2)=x3(2)*cost3-x3(3)*sint3
       x3p(3)=x3(2)*sint3+x3(3)*cost3
      endif

      x1p(1)=x1(1)
      x1p(2)=x1(2)*cost3-x1(3)*sint3
      x1p(3)=x1(2)*sint3+x1(3)*cost3

      x2p(1)=x2(1)
      x2p(2)=x2(2)*cost3-x2(3)*sint3
      x2p(3)=x2(2)*sint3+x2(3)*cost3

      do i=1,3
       x1(i)=x1p(i)
       x2(i)=x2p(i)
       x3(i)=x3p(i)
      enddo

      write(50,*) '#Rotation Rx(theta3) about x'
      Rx(1,:) = (/ 1.0d0,  0.0d0,  0.0d0 /)
      Rx(2,:) = (/ 0.0d0,  cost3, -sint3 /)
      Rx(3,:) = (/ 0.0d0,  sint3,  cost3 /)
      
      do i=1,3
         write(50,5000) ( Rx(i,j), j=1,3 )
      enddo

!     Now that we have the rotations, we generate the crystal, rotate,
!     and place molecules into the new cell.

      Rtemp = matmul(Ry,Rz)
      Rall = matmul(Rx,Rtemp)
      write(50,*) '#Net Rotation matrix'
      
! 	  Rall*chosen_h_matrix = newcell_h-matrix (tested)

      do i=1,3
         write(50,5000) ( Rall(i,j), j=1,3 )
      enddo
      
      write(*,*) '#Net Rotation matrix'
      
      do i=1,3
         write(*,5000) ( Rall(i,j), j=1,3 )
      enddo
      
      close(50)
      !Get inverse H matrix for new cell

      lx1=x1(1)
      lx2=dsqrt(x2(1)**2+x2(2)**2)
      lx3=dsqrt(x3(1)**2+x3(2)**2+x3(3)**2)

      theta=dacos(x2(1)/lx2)
      phi=dacos(x3(1)/lx3)
      psi=dacos((x2(1)*x3(1)+x2(2)*x3(2))/(lx2*lx3))

      volcell=x1(1)*x2(2)*x3(3)/(lx1*lx2*lx3)

!     write out new cell parameters and h matrix

      open(90,file='newcell.parameters.txt')
      write(90,9090)
9090  format('#|x1| |x2| |x3| theta phi psi')
      write(90,9091) lx1, lx2, lx3, theta*180.0d0/pi, phi*180.0d0/pi, psi*180.0d0/pi
9091  format(6F16.8)

      write(90,9092)
9092  format('#H-Matrix')   
      write(90,9093) x1(1), x2(1), x3(1)
      write(90,9093) x1(2), x2(2), x3(2)
      write(90,9093) x1(3), x2(3), x3(3)
9093  format(3F16.8)

      

      hinv(1,1)=1.0d0/lx1     !define inverse h-matrix
      hinv(2,1)=0.0d0
      hinv(3,1)=0.0d0

      hinv(1,2)=-dcos(theta)/(lx1*dsin(theta))
      hinv(2,2)=1.0d0/(lx2*dsin(theta))
      hinv(3,2)=0.0d0

      hinv(1,3)=(dcos(psi)*dcos(theta)-dcos(phi))/(lx1*volcell*dsin(theta))
      hinv(2,3)=(dcos(phi)*dcos(theta)-dcos(psi))/(lx2*volcell*dsin(theta))
      hinv(3,3)=dsin(theta)/(lx3*volcell)

      write(90,9094)
9094  format('#invH-Matrix')   
      write(90,9095) hinv(1,1), hinv(1,2), hinv(1,3)
      write(90,9095) hinv(2,1), hinv(2,2), hinv(2,3)
      write(90,9095) hinv(3,1), hinv(3,2), hinv(3,3)
9095  format(3F16.8)

	   close(90)
	
      !Populate new cell

      allocate (r(molperunit,maxval(atpermol,DIM=1),3))
      allocate (rp(molperunit,maxval(atpermol,DIM=1),3))
      allocate (rcmo(molperunit,3),rcm(molperunit,3),rcmp(molperunit,3))

      oldvol=av(1)*bv(2)*cv(3)  !estimate new number of molecules
      newvol=x1(1)*x2(2)*x3(3)

      nmolEST=4*molperunit*nint(newvol/oldvol)  !somewhat arbitrary upper bound on number of mols in new cell

      write(*,100) oldvol, newvol, molperunit*nint(newvol/oldvol)
100   format(2F10.4,I8) 

      allocate (rbox(nmolEST,maxval(atpermol,DIM=1),3))
      allocate (moltype(nmolEST))   !holds info on each mol type for printing atpermol quantities later

      mmin=2*mmin*x3scale
      mmax=2*mmax*x3scale

      nmin=2*nmin*x3scale
      nmax=2*nmax*x3scale

      pmin=2*pmin*x3scale
      pmax=2*pmax*x3scale

      nmol=0

!     Get molecular centers-of-mass

      rcmo=0.0d0

      do i=1,molperunit
       do j=1,atpermol(i)
        do k=1,3
         rcmo(i,k)=rcmo(i,k)+ro(i,j,k)*mass(i,j)
        enddo
       enddo
       do k=1,3
        rcmo(i,k)=rcmo(i,k)/molmass(i)
       enddo
      enddo

      do m=mmin,mmax
       do n=nmin,nmax
        do p=pmin,pmax

         !get translated molecular center-of-mass coordinates

         do i=1,molperunit
          do j=1,3
           rcm(i,j)=dfloat(m)*av(j)+dfloat(n)*bv(j)+dfloat(p)*cv(j)+rcmo(i,j)
          enddo
         enddo

         !Rotate the translated unit cell by theta1, theta2, and theta3

         do i=1,molperunit     !rotate by theta1 about z
          rcmp(i,1)=rcm(i,1)*cost1-rcm(i,2)*sint1
          rcmp(i,2)=rcm(i,1)*sint1+rcm(i,2)*cost1
          rcmp(i,3)=rcm(i,3)
         enddo

         do i=1,molperunit
          do j=1,3
           rcm(i,j)=rcmp(i,j)
          enddo
         enddo

         do i=1,molperunit     !rotate by theta2 about y
          rcmp(i,1)=rcm(i,1)*cost2+rcm(i,3)*sint2
          rcmp(i,2)=rcm(i,2)
          rcmp(i,3)=-rcm(i,1)*sint2+rcm(i,3)*cost2
         enddo

         do i=1,molperunit
          do j=1,3
           rcm(i,j)=rcmp(i,j)
          enddo
         enddo

         do i=1,molperunit    !rotate by theta3 about x
          rcmp(i,1)=rcm(i,1)
          rcmp(i,2)=rcm(i,2)*cost3-rcm(i,3)*sint3
          rcmp(i,3)=rcm(i,2)*sint3+rcm(i,3)*cost3
         enddo

         do i=1,molperunit
          do j=1,3
           rcm(i,j)=rcmp(i,j)
          enddo
         enddo

         !Transform to fractional coordinates

         do i=1,molperunit
          rcmp(i,1)=rcm(i,1)*hinv(1,1)+rcm(i,2)*hinv(1,2)+rcm(i,3)*hinv(1,3)
          rcmp(i,2)=rcm(i,1)*hinv(2,1)+rcm(i,2)*hinv(2,2)+rcm(i,3)*hinv(2,3)
          rcmp(i,3)=rcm(i,1)*hinv(3,1)+rcm(i,2)*hinv(3,2)+rcm(i,3)*hinv(3,3)
         enddo

         !Find mols that are in the new cell
         !Note that we take the domain to be [0,1.1] so we get
         !molecules that lie on or near the edges of the cell.
         !We'll remove any periodically equivalent molecules later.

         do i=1,molperunit
          if(rcmp(i,1).ge.0.0d0 .and. rcmp(i,1).le.1.01d0)then
           if(rcmp(i,2).ge.0.0d0 .and. rcmp(i,2).le.1.01d0)then
            if(rcmp(i,3).ge.0.0d0 .and. rcmp(i,3).le.1.01d0)then
             nmol=nmol+1
             moltype(nmol)=i

             !Now that we have identified a molecule that belongs
             !in the cell, we need its translated-rotated atomic
             !coordinates.

             !translate the coordinates of molecule i 

             do j=1,atpermol(i)
              do k=1,3
               r(i,j,k)=dfloat(m)*av(k)+dfloat(n)*bv(k)+dfloat(p)*cv(k)+ro(i,j,k)
              enddo
             enddo

             !Rotate the atomic coordinates of molecule i by theta1, theta2, and theta3
             
             do j=1,atpermol(i) !rotate by theta1 about z
              rp(i,j,1)=r(i,j,1)*cost1-r(i,j,2)*sint1
              rp(i,j,2)=r(i,j,1)*sint1+r(i,j,2)*cost1
              rp(i,j,3)=r(i,j,3)
             enddo
                  
             do j=1,atpermol(i)
              do k=1,3
               r(i,j,k)=rp(i,j,k)
              enddo
             enddo
         
             do j=1,atpermol(i) !rotate by theta2 about y
              rp(i,j,1)=r(i,j,1)*cost2+r(i,j,3)*sint2
              rp(i,j,2)=r(i,j,2)
              rp(i,j,3)=-r(i,j,1)*sint2+r(i,j,3)*cost2
             enddo
                  
             do j=1,atpermol(i)
              do k=1,3
               r(i,j,k)=rp(i,j,k)
              enddo
             enddo
         
             do j=1,atpermol(i) !rotate by theta3 about x
              rp(i,j,1)=r(i,j,1)
              rp(i,j,2)=r(i,j,2)*cost3-r(i,j,3)*sint3
              rp(i,j,3)=r(i,j,2)*sint3+r(i,j,3)*cost3
             enddo
         
             do j=1,atpermol(i)
              do k=1,3
               r(i,j,k)=rp(i,j,k)
              enddo
             enddo
         
             !Now put the translated-rotated molecule into an array 
             !that holds all of the molecules in the cell

             do j=1,atpermol(i)
              do k=1,3
               rbox(nmol,j,k)=r(i,j,k)
              enddo
             enddo

            endif
           endif
          endif
         enddo

        !m,n,p loops end here

        enddo
       enddo
      enddo

!     Remove any molecules that overlap with periodic images.

      allocate (duplicates(nmol*nmol))
      duplicates=0
      index=1

      do i=1,nmol-1
        ri(1)=rbox(i,1,1)
        ri(2)=rbox(i,1,2)
        ri(3)=rbox(i,1,3)
        do j=i+1,nmol
          rj(1)=rbox(j,1,1)
          rj(2)=rbox(j,1,2)
          rj(3)=rbox(j,1,3)
          call minimage(ri,rj,x1,x2,x3,rijmin)
          rijmag=dsqrt(rijmin(1)**2+rijmin(2)**2+rijmin(3)**2)
          if(rijmag.lt.1.0d0)then  !setting the distance threshold to 1.0 A is arbitrary 
            duplicates(index)=j    !and can be changed if problems are encountered. 
            index=index+1
          endif
        enddo
      enddo

!     Get new nmol after removing overlapping periodic images

      newnmol=0
      newnatom=0
      cellcom=0.0d0
      totalmass=0.0d0

      do i=1,nmol
        if(ANY(duplicates.eq.i))then
          !found a duplicate, don't count it
        else
          newnmol=newnmol+1
          newnatom=newnatom+atpermol(moltype(i))

          !get molecule contribution to COM

          do j=1,atpermol(moltype(i))
            do k=1,3
              cellcom(k)=cellcom(k)+mass(moltype(i),j)*rbox(i,j,k)
            enddo
            totalmass=totalmass+mass(moltype(i),j)
          enddo

        endif
      enddo

      do k=1,3
        cellcom(k)=cellcom(k)/totalmass
      enddo

!     Get center of cell, which is (0.5,0.5,0.5) in fractional coords.
!     Note that the matrix (x1,x2,x3) gives this transformation and
!     should be an upper-triangular matrix.

      geocenter(1)=0.5d0*x1(1)+0.5d0*x2(1)+0.5d0*x3(1)
      geocenter(2)=            0.5d0*x2(2)+0.5d0*x3(2)
      geocenter(3)=                        0.5d0*x3(3)

!     Translate molecules so the system com is at the cell center

      do i=1,nmol
        if(ANY(duplicates.eq.i))then
          !found a duplicate, skip it
        else

          do j=1,atpermol(moltype(i))
            do k=1,3
              rbox(i,j,k)=rbox(i,j,k)-cellcom(k)+geocenter(k)
              
            enddo
          enddo

        endif
      enddo

!     Check that new number of molecules makes sense compared to
!     the unit cell and issue the user a warning if it does not.

      volcell=x1(1)*x2(2)*x3(3)

      if(newnmol.ne.nint(volcell/volunit)*molperunit)then
        write(*,1090)
        write(*,1091) newnmol
        write(*,1092) nint(volcell/volunit)*molperunit
      endif

1090  format('Incorrect number of molecules in cell!')
1091  format('There are ', I12, ' molecules in the cell,')
1092  format('but there should only be ', I12)

!     Print new box to files

      write(*,1010) newnmol
1010  format(' Number of molecules in new simulation cell: ', I12)

      open(94,file='newcell.xyz')
      write(94,9401) newnatom
      write(94,9402)
9401  format(I16)
9402  format('Atoms')

!     LAMMPS style dump with pbc box

	  xy=x2(1)
      xz=x3(1)
      yz=x3(2)

!	  major bug fixed here ! (Apr 17, 2017)
!	  do while(abs(xy).gt.(x1(1)*0.5))
!		xy = xy - x1(1)*sign(1.0d+00,xy)
!      end do	
!
!	  do while(abs(xz).gt.(x1(1)*0.5))
!		xz = xz - x1(1)*sign(1.0d+00,xz)
!      end do	
!
!	  do while(abs(yz).gt.(x2(2)*0.5))
!		yz = yz - x2(2)*sign(1.0d+00,yz)
!      end do	


      xlo=min(0.0d0,xy,xz,(xy+xz))
      xhi=x1(1)+max(0.0d0,xy,xz,(xy+xz))

      ylo=min(0.0d0,yz)
      yhi=x2(2)+max(0.0d0,yz)

      zlo=0.0d0
      zhi=x3(3)
      
      open(95,file='newcell.lammpstrj')
	  open(96,file='newcell.lmp')

	  write(96,9601)
	  write(96,9602)
	  write(96,9603) newnatom
	  write(96,9604)
	  write(96,9605) 4
	  write(96,9606)
	  write(96,9607)
	  write(96,9608) 0.0d0,x1(1)
	  write(96,9609) 0.0d0,x2(2)
	  write(96,9610) 0.0d0,x3(3)
	  write(96,9611) xy,xz,yz
	  write(96,9612)
	  write(96,9613)
	  write(96,9614)
	  write(96,9615)
	  write(96,9616)
	  write(96,9617)
	  write(96,9618)
	  write(96,9619)
	  write(96,9620)
	  write(96,9621)
	  

9601  format('LAMMPS Description')
9602  format('')
9603  format(I6, ' atoms')
9604  format('')
9605  format(I6, ' atom types')
9606  format('')
9607  format('# Cell: monoclinic')
9608  format(2F16.8, '    xlo xhi')
9609  format(2F16.8, '    ylo yhi')
9610  format(2F16.8, '    zlo zhi')
9611  format(3F16.8, '    xy xz yz')
9612  format('')
9613  format('Masses')
9614  format('')
9615  format('1 12.011 # C')
9616  format('2 14.007 # N')
9617  format('3 15.999 # O')
9618  format('4 1.0079 # H')
9619  format('')
9620  format('Atoms')
9621  format('')

      write(95,9501)
      write(95,9502)
      write(95,9503)
      write(95,9504) newnatom
      write(95,9505)
      write(95,9506) xlo, xhi, xy
      write(95,9506) ylo, yhi, xz
      write(95,9506) zlo, zhi, yz
      write(95,9507)

9501  format('ITEM: TIMESTEP')
9502  format('0')
9503  format('ITEM: NUMBER OF ATOMS')
9504  format(I16)
9505  format('ITEM: BOX BOUNDS xy xz yz pp pp pp')
9506  format(3F16.8)
9507  format('ITEM: ATOMS id type xu yu zu')

      iter=1

      do i=1,nmol
        if(ANY(duplicates.eq.i))then
          !found a duplicate molecule, don't write it out
        else
          do j=1,atpermol(moltype(i))
            write(94,9403) el(moltype(i),j), (rbox(i,j,k),k=1,3)
9403        format(A6,3ES24.16)
			
			if(el(moltype(i),j).eq.'c'.or.el(moltype(i),j).eq.'C') type_atom=1
			if(el(moltype(i),j).eq.'n'.or.el(moltype(i),j).eq.'N') type_atom=2
			if(el(moltype(i),j).eq.'o'.or.el(moltype(i),j).eq.'O') type_atom=3
			if(el(moltype(i),j).eq.'h'.or.el(moltype(i),j).eq.'H') type_atom=4
			mol_id = iter
			
            write(95,9508) iter, type_atom, (rbox(i,j,k),k=1,3)
            write(96,9622) iter, type_atom, (rbox(i,j,k),k=1,3)
9508        format(I6,I6,3F16.8)
9622        format(I6,I6,3F16.8)
            iter=iter+1
          enddo
        endif
      enddo

      close(94)
      close(95)
      close(96)

      deallocate (atpermol,ro,el,mass,molmass,r,rp,rcmo,rcm,rcmp,rbox,moltype,duplicates)

      end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine minimage(ri,rj,x1,x2,x3,rijmin)
      implicit none

!     This subroutine is a modified F90 version of the c++ LAMMPS  
!     routine closest_image in domain.cpp. License information for  
!     the original version is immediately below. 
!
!     /* ----------------------------------------------------------------------
!     LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
!     http://lammps.sandia.gov, Sandia National Laboratories
!     Steve Plimpton, sjplimp@sandia.gov
!
!     Copyright (2003) Sandia Corporation.  Under the terms of Contract
!     DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
!     certain rights in this software.  This software is distributed under
!     the GNU General Public License.
!
!     See the README file in the top-level LAMMPS directory.
!     ------------------------------------------------------------------------- */
!
!     /* ----------------------------------------------------------------------
!     Contributing author (triclinic) : Pieter in 't Veld (SNL)
!     ------------------------------------------------------------------------- */

      double precision, intent(in) :: ri(3), rj(3)
      double precision, intent(in out) :: rijmin(3)
      double precision, intent(in) :: x1(3), x2(3), x3(3)
      double precision :: dx, dy, dz

      dx = rj(1) - ri(1)
      dy = rj(2) - ri(2)
      dz = rj(3) - ri(3)

      !First handle x3

      if(dz.lt.0.0d0)then
        do while(dz.lt.0.0d0)
          dx = dx + x3(1)
          dy = dy + x3(2)
          dz = dz + x3(3)
        enddo
        if(dz.gt.0.5d0*x3(3))then
          dx = dx - x3(1)
          dy = dy - x3(2)
          dz = dz - x3(3)
        endif
      else
        do while(dz.gt.0.0d0)
          dx = dx - x3(1)
          dy = dy - x3(2)
          dz = dz - x3(3)
        enddo
        if(dz.lt.-0.5d0*x3(3))then
          dx = dx + x3(1)
          dy = dy + x3(2)
          dz = dz + x3(3)
        endif
      endif

      !Now x2

      if(dy.lt.0.0d0)then
        do while(dy.lt.0.0d0)
          dx = dx + x2(1)
          dy = dy + x2(2)
        enddo
        if(dy.gt.0.5d0*x2(2))then
          dx = dx - x2(1)
          dy = dy - x2(2)
        endif
      else
        do while(dy.gt.0.0d0)
          dx = dx - x2(1)
          dy = dy - x2(2)
        enddo
        if(dy.lt.-0.5d0*x2(2))then
          dx = dx + x2(1)
          dy = dy + x2(2)
        endif
      endif

      !And last x1
      if(dx.lt.0.0d0)then
        do while(dx.lt.0.0d0)
          dx = dx + x1(1)
        enddo
        if(dx.gt.0.5d0*x1(1))then
          dx = dx - x1(1)
        endif
      else
        do while(dx.gt.0.0d0)
          dx = dx - x1(1)
        enddo
        if(dx.lt.-0.5d0*x1(1))then
          dx = dx + x1(1)
        endif
      endif

      rijmin(1) = dx
      rijmin(2) = dy
      rijmin(3) = dz

      end subroutine minimage















