      program GCCM_SingleCrystal
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
!     This program is the first part of the single-crystal version of
!     the Generalized Crystal-Cutting Method (GCCM) and generates crystals
!     that are oriented with an arbitrary direction S parallel to the
!     normal of a face of the simulation cell. It searches possible
!     new lattice vectors x1, x2, and x3, that define the new cell
!     and have form x = ma + nb + pc, where (a,b,c) are the crystal
!     lattice vectors and (m,n,p) are integers. Chosing the x's to be
!     of this form guarantees that the new cell can generate the
!     original crystal through simple translations. Periodicity, even
!     with a complicated set of basis atoms, is exactly preserved.
!
!     The general proceedure is as follows:
!
!     1) First choose x1 and x2 such that x1 x x2 is parallel to S
!        within a user-defined tolerance. This ensures that the
!        normal of a face of the new simulation cell is parallel to
!        the direction S.
!
!     2) The parameter angledev12 controls the maximum deviation
!        of the angle between x1 and x2 away from 90 degrees.
!        That is, setting angledev12 = 10.0 will allow the angle
!        x1 and x2 to vary between 80.0 and 100.0. It is possible
!        that no set of (x1,x2) exist in the searched domain.
!
!     3) Having hopefully found possible pairs of (x1,x2), we then
!        search for the best x3. In this case, we explicitly look 
!        for the x3 that is most closely aligned to S, which is
!        equivalent to finding the that x3 which is closest to
!        being orthogonal to x1 and x2.
!
!     4) The solutions (x1,x2,x3) are then printed to a file with
!        geometric details on each solution. 
!
!     5) Pick the best solution for your problem and copy it to the
!        file chosen.solution.txt to generate a new simulation cell
!        using the second GCCM program, GCCM_GenerateCrystal.f90.
!
!     ####################################################################

      double precision :: S(3), Sm(3)
      double precision :: a, b, c, alpha, beta, gamma, volume
      double precision :: tol, angledev12, angledevS3
      double precision :: av(3), bv(3), cv(3)
      double precision :: Rav(3), Rbv(3), Rcv(3), STP
      double precision :: dpcount
      double precision :: dprod0, costheta, area, length, dprod1
      double precision :: lx1, lx2, lx3, theta, phi, psi
      double precision :: V(3), V_orig(3), V_orign(3)
      double precision :: x1(3), x2(3), x3(3), x1_orig(3), x2_orig(3), x3_orig(3), x1n(3), x2n(3)
      double precision :: oldmax, oldmin, min_angle12, min_area12
      double precision :: Maxangledev12, Aveangledev12
      double precision :: MinDevMinVol
      double precision :: MinVol
      double precision :: statistics(2,7)
      double precision, parameter :: pi=4.0d0*datan(1.0d0)

      integer :: mmin, mmax
      integer :: nmin, nmax
      integer :: pmin, pmax
      integer :: maxpairs, x3scale, max3
      integer :: mcount, ncount, pcount
      integer :: m, mm, n, nn, p, pp
      integer :: i, j, k, l
      integer :: scount, scount3, gen_count, count_X12, count_X3
      integer, allocatable :: mnp(:,:,:), mnp3(:,:,:), X_12(:,:,:), X_3(:,:,:)
      integer :: mostOrtho, leastVol

      character*1 :: char

!     Read control file

      open(30,action='read',file='control.txt')

      read(30,*)
      read(30,*) char, tol
      read(30,*) char, angledev12, angledevS3
      read(30,*) char, mmin, mmax
      read(30,*) char, nmin, nmax
      read(30,*) char, pmin, pmax
      read(30,*) char, x3scale
      read(30,*) char, maxpairs, max3
      read(30,*)
      read(30,*)
      read(30,*) char, Sm(1)
      read(30,*) char, Sm(2)
      read(30,*) char, Sm(3)
      read(30,*)
      read(30,*)
      read(30,*) char, a
      read(30,*) char, b
      read(30,*) char, c
      read(30,*) char, alpha
      read(30,*) char, beta
      read(30,*) char, gamma

      close(30)

!     Prepare for solution search

      alpha=alpha*pi/180.0d0
      beta=beta*pi/180.0d0
      gamma=gamma*pi/180.0d0

      av(1)=a           !x component
      av(2)=0.0d0       !y component
      av(3)=0.0d0       !z component

      bv(1)=b*dcos(gamma)   !same format as above
      bv(2)=b*dsin(gamma)
      bv(3)=0.0d0

      cv(1)=c*dcos(beta)
      cv(2)=c*(dcos(alpha)-dcos(beta)*dcos(gamma))/dsin(gamma)
      cv(3)=c*dsqrt(1.0d0-(dcos(alpha))**2-(dcos(beta))**2-(dcos(gamma))**2+2.0d0*dcos(alpha)*dcos(beta)*dcos(gamma))/dsin(gamma)

      ! get reciprocal lattice vectors
            
      call cross(bv,cv,V)
      !call dot(av,V,STP)
      STP = av(1)*bv(2)*cv(3)
      Rav(1) = V(1)/STP
      Rav(2) = V(2)/STP
      Rav(3) = V(3)/STP
      
      call cross(cv,av,V)
      Rbv(1) = V(1)/STP
      Rbv(2) = V(2)/STP
      Rbv(3) = V(3)/STP
      
      call cross(av,bv,V)
      Rcv(1) = V(1)/STP
      Rcv(2) = V(2)/STP
      Rcv(3) = V(3)/STP
      
      ! express direction normal in terms of miller indices and reciprocal lattice vectors
      
      do i=1,3
          S(i)=Sm(1)*Rav(i)+Sm(2)*Rbv(i)+Sm(3)*Rcv(i)
      enddo

      call norm(S)  !normalize S

      write(*,800) (Rav(i),i=1,3)
800   format(3F10.6)
      write(*,800) (Rbv(i),i=1,3)
      write(*,800) (Rcv(i),i=1,3)

	  write(*,801) (Sm(i),i=1,3)
801   format(3F6.2)

	  write(*,802) (S(i),i=1,3)
802   format(3F10.6)

      write(*,900) (av(i),i=1,3)
900   format(3F10.6)
      write(*,900) (bv(i),i=1,3)
      write(*,900) (cv(i),i=1,3)

      allocate (mnp(maxpairs,3,3))
      allocate (mnp3(max3,1,3))

      allocate (X_12(100,2,3))
      allocate (X_3(100,1,3))

      mcount=0  !get total number of (x1,x2) solutions searched
      ncount=0
      pcount=0

      do m=mmin,mmax
       mcount=mcount+1
      enddo
      do n=nmin,nmax
       ncount=ncount+1
      enddo
      do p=pmin,pmax
       pcount=pcount+1
      enddo

      dpcount=dfloat(mcount*mcount)*dfloat(ncount*ncount)*dfloat(pcount*pcount)

      angledev12=dcos((90.0d0-angledev12)*pi/180.0d0)  !take cos of 90-angledev12 for numerically cheap
                                                   !comparisions to the dot product x1.x2

      angledevS3=dcos((90.0d0-angledevS3)*pi/180.0d0)     !comparisions to the dot product x3.S
      
!     Search (x1,x2) pairs

      write(*,1000) dpcount
1000  format('Searching over ', ES9.2, ' (x1,x2) pairs.')

      scount=0  !solution counting variable
      min_area12=1.0d5

      open(100,file='pairs.txt')


      mloop: do m=mmin,mmax
       nloop: do n=nmin,nmax
        ploop: do p=pmin,pmax

!        form x1

         do i=1,3
          x1(i)=dfloat(m)*av(i)+dfloat(n)*bv(i)+dfloat(p)*cv(i)
          x1_orig(i)=dfloat(m)*av(i)+dfloat(n)*bv(i)+dfloat(p)*cv(i)
         enddo

         call norm(x1)

         mmloop: do mm=mmin,mmax
          nnloop: do nn=nmin,nmax
           pploop: do pp=pmin,pmax

            do i=1,3
              x2(i)=dfloat(mm)*av(i)+dfloat(nn)*bv(i)+dfloat(pp)*cv(i)
              x2_orig(i)=dfloat(mm)*av(i)+dfloat(nn)*bv(i)+dfloat(pp)*cv(i)
            enddo

            call norm(x2)

!           if x1 and x2 are nearly parallel, skip that solution as it
!           will be unstable and likely result in a highly skewed box.
!           we also

            call dot(x1,x2,costheta)

            if (dabs(costheta).le.angledev12)then

!            compute cross product, V = x1 x x2
!            if x1 || x2, then cross product is the 0 vector

             call cross(x1,x2,V)
             call norm(V)
             call dot(V,S,dprod0)
             
             call cross(x1_orig,x2_orig,V_orig)

             if(dprod0.ge.(1.0d0-tol))then   !This if expression ensures
              scount=scount+1               !that we have a right-handed
              mnp(scount,1,1)=m             !lattice coordinate system
              mnp(scount,1,2)=n             !wrt the direction S
              mnp(scount,1,3)=p
              mnp(scount,2,1)=mm            !store the solutions for finding x3
              mnp(scount,2,2)=nn
              mnp(scount,2,3)=pp
              mnp(scount,3,1)=0             !set initial x3 solutions to zero vector
              mnp(scount,3,2)=0
              mnp(scount,3,3)=0
              area=dsqrt(V_orig(1)*V_orig(1)+V_orig(2)*V_orig(2)+V_orig(3)*V_orig(3))

              if(area.le.min_area12)then
               min_area12=area
              endif

              write(100,901) (mnp(scount,1,i),i=1,3), (mnp(scount,2,i),i=1,3), &
                             (x1(i),i=1,3), (x2(i),i=1,3), (V(i),i=1,3), dprod0, costheta, area
901           format(6I4,12F10.2)


              if(scount.ge.maxpairs)then
                write(*,1001) maxpairs
1001            format('Minimum number of (x1,x2) pairs found: ', I16)
                write(*,1002)
1002            format('Increase maxpairs to consider more pairs.')
                exit mloop
              endif

             endif

            endif

           end do pploop
          end do nnloop
         end do mmloop

        end do ploop
       end do nloop
      end do mloop

      close(100)

      write(*,1003) scount
1003  format('Number of (x1,x2) pairs found: ', I16)

      if(scount.eq.0)then
        write(*,1005)
1005    format('No (x1,x2) found. Increase possible m, n, and p.')
        stop
      endif
      
      
      gen_count=0
 
      do i=1,scount
       do j=1,3
        x1(j)=dfloat(mnp(i,1,1))*av(j)+dfloat(mnp(i,1,2))*bv(j)+dfloat(mnp(i,1,3))*cv(j)
        x2(j)=dfloat(mnp(i,2,1))*av(j)+dfloat(mnp(i,2,2))*bv(j)+dfloat(mnp(i,2,3))*cv(j)
       enddo
	   call cross(x1,x2,V_orig)
	   area=dsqrt(V_orig(1)*V_orig(1)+V_orig(2)*V_orig(2)+V_orig(3)*V_orig(3))
	   call dot(x1,x2,costheta)
	   
	   if(area.le.(min_area12*1.2))then
	    gen_count=gen_count+1
	    X_12(gen_count,1,1)=mnp(i,1,1)
	    X_12(gen_count,1,2)=mnp(i,1,2)
	    X_12(gen_count,1,3)=mnp(i,1,3)
	    
	    X_12(gen_count,2,1)=mnp(i,2,1)
	    X_12(gen_count,2,2)=mnp(i,2,2)
	    X_12(gen_count,2,3)=mnp(i,2,3)
	    
	    call norm(V_orig)
	    
	    call dot(V_orig,S,dprod0)
	    
	    !write(101,905) (mnp(i,1,k),k=1,3), (mnp(i,2,k),k=1,3), (S(k),k=1,3), (V_orig(k),k=1,3), dprod0, costheta, area
905     format(6I8,6F8.4,3F10.4)
		!write(101,907) (x1(k),k=1,3), (x2(k),k=1,3) 
907     format(6F10.4)

	   endif
		
      enddo
      
      count_X12 = gen_count;
      
!     Now that we found (x1,x2) pairs, find best x3.
!     Note that when defining the best x3 to be that which is 
!     most closely aligned with S, each (x1,x2) pair has the
!     same 'best x3'.

      write(*,1004)
1004  format('Determining the best x3.')

      oldmax=-1.0d3  !set initial oldmax to some large negative value
      oldmin=1.0d5   !set initial oldmin to some large positive value

      scount3=0

      mloop2: do m=mmin*x3scale,mmax*x3scale      !loop over possible x3
       nloop2: do n=nmin*x3scale,nmax*x3scale
        ploop2: do p=pmin*x3scale,pmax*x3scale

         do j=1,3
          x3(j)=dfloat(m)*av(j)+dfloat(n)*bv(j)+dfloat(p)*cv(j)
          x3_orig(j)=dfloat(m)*av(j)+dfloat(n)*bv(j)+dfloat(p)*cv(j)
         enddo

         call norm(x3)
         call dot(x3,S,dprod0)
         length = dsqrt(x3_orig(1)*x3_orig(1)+x3_orig(2)*x3_orig(2)+x3_orig(3)*x3_orig(3))


         if((1.0-dabs(dprod0)).le.angledevS3)then
          scount3=scount3+1     !Pick out a set of x3 that are within angledevS3 of S
          mnp3(scount3,1,1)=m                 
          mnp3(scount3,1,2)=n                 
          mnp3(scount3,1,3)=p                
          if(length.le.oldmin)then
           oldmin=length
          endif
          
          if(scount3.ge.max3)then
                write(*,903) max3
903             format('Minimum number of (x3) pairs found: ', I16)
                write(*,904)
904             format('Increase max3 to consider more x3.')
                exit mloop2
          endif
            
         endif
        
        end do ploop2
       end do nloop2
      end do mloop2

      
      gen_count=0
      
      do i=1,scount3
       do j=1,3
        x3(j)=dfloat(mnp3(i,1,1))*av(j)+dfloat(mnp3(i,1,2))*bv(j)+dfloat(mnp3(i,1,3))*cv(j)
       enddo
       
       length = dsqrt(x3(1)*x3(1)+x3(2)*x3(2)+x3(3)*x3(3))
        
       if(length.le.(oldmin*1.2))then
         gen_count=gen_count+1
         X_3(gen_count,1,1)=mnp3(i,1,1)
	     X_3(gen_count,1,2)=mnp3(i,1,2)
 	     X_3(gen_count,1,3)=mnp3(i,1,3)
 	     
 	     !write(101,906) (mnp3(i,1,k),k=1,3), (x3(k),k=1,3), length
906     format(3I8,4F10.4)
 	     
       endif
      enddo

      open(101,file='select_pairs.txt')

      count_X3 = gen_count;
      
      write(*,902) min_area12, oldmin, scount3, count_X12, count_X3
902   format('Min area x1_x2: ', F16.6, ' Min length x3: ', F8.4, 3I16)

      do i=1,count_X12
       do j=1,3
        x1(j)=dfloat(X_12(i,1,1))*av(j)+dfloat(X_12(i,1,2))*bv(j)+dfloat(X_12(i,1,3))*cv(j)
        x2(j)=dfloat(X_12(i,2,1))*av(j)+dfloat(X_12(i,2,2))*bv(j)+dfloat(X_12(i,2,3))*cv(j)
       enddo
       
       call cross(x1,x2,V_orig)
       
       call norm1(x1,x1n)
       call norm1(x2,x2n)
       call dot(x1n,x2n,costheta)
       
       if(dabs(costheta).le.angledev12) then
	       do l=1,count_X3
	        do j=1,3
	         x3(j)=dfloat(X_3(l,1,1))*av(j)+dfloat(X_3(l,1,2))*bv(j)+dfloat(X_3(l,1,3))*cv(j)
	        enddo
	       
	        call dot(x3,V_orig,dprod0)
	        if(dprod0.gt.0) then
	        
	         call norm1(V_orig,V_orign)
	         call dot(S,V_orign,dprod1)
	         
	         write(101,908) (X_12(i,1,k),k=1,3), (X_12(i,2,k),k=1,3), (X_3(l,1,k),k=1,3), costheta, abs(dprod0)/493.230d0, abs(dprod1)
908          format(9I8,3F10.4) 
             write(101,909) (x1(k),k=1,3), (x2(k),k=1,3), (x3(k),k=1,3)
909          format(9F10.4) 
	
	        endif
	       
	       enddo
	   endif 
      
      enddo
       
      
      close(101)

!     We found our cells (x1,x2,x3), now identify amoung those the cell   
!     that has the smallest volume, is right-handed and is not too skewed. 

      Maxangledev12=0.0d0
      MinDevMinVol=1.0d20
      MinVol=1.0d20

      do i=1,scount
       do k=1,scount3

!       do j=1,3
!        x1(j)=dfloat(mnp(i,1,1))*av(j)+dfloat(mnp(i,1,2))*bv(j)+dfloat(mnp(i,1,3))*cv(j)
!        x2(j)=dfloat(mnp(i,2,1))*av(j)+dfloat(mnp(i,2,2))*bv(j)+dfloat(mnp(i,2,3))*cv(j)
!        x3(j)=dfloat(mnp(i,3,1))*av(j)+dfloat(mnp(i,3,2))*bv(j)+dfloat(mnp(i,3,3))*cv(j)
!       enddo

!!      compute new cell edge vector lengths

!       lx1=0.0d0
!       lx2=0.0d0
!       lx3=0.0d0

!       do j=1,3
!        lx1=lx1+x1(j)*x1(j)
!        lx2=lx2+x2(j)*x2(j)
!        lx3=lx3+x3(j)*x3(j)
!       enddo

!       lx1=dsqrt(lx1)
!       lx2=dsqrt(lx2)
!       lx3=dsqrt(lx3)

!!      compute new cell volume

!       call cross(x1,x2,V)
!       call dot(x3,V,volume)
!       volume=dabs(volume)  !volume in A**3

!!      get angles between the edge vectors

!       call norm(x1)
!       call norm(x2)
!       call norm(x3)

!       call dot(x1,x2,theta)
!       call dot(x1,x3,phi)
!       call dot(x2,x3,psi)

!       theta=dacos(theta)*180.0d0/pi
!       phi=dacos(phi)*180.0d0/pi
!       psi=dacos(psi)*180.0d0/pi

!!      Find the closest to orthorhombic cell

!       call cross(x1,x2,V)
!       call dot(x3,V,AveAngleDev)
!       AveAngleDev=dabs(AveAngleDev)   !orthorhombic will yield 1

!       if(AveAngleDev.ge.MaxAngleDev)then
!        if(AveAngleDev.gt.MaxAngleDev)then
!         MaxAngleDev=AveAngleDev
!         MinDevMinVol=volume
!         mostOrtho=i !closest to orthogonal solution number
!         statistics(1,1)=lx1
!         statistics(1,2)=lx2
!         statistics(1,3)=lx3
!         statistics(1,4)=theta
!         statistics(1,5)=phi
!         statistics(1,6)=psi
!         statistics(1,7)=volume
!        elseif(AveAngleDev.eq.MaxAngleDev .and. volume.lt.MinDevMinVol)then
!         MinDevMinVol=volume
!         mostOrtho=i !closest to orthogonal solution number
!         statistics(1,1)=lx1
!         statistics(1,2)=lx2
!         statistics(1,3)=lx3
!         statistics(1,4)=theta
!         statistics(1,5)=phi
!         statistics(1,6)=psi
!         statistics(1,7)=volume
!        endif
!       endif

!!      Find the cell with smallest volume amoung those identified

!       if(volume.lt.MinVol)then
!        MinVol=volume
!        leastVol=i
!        statistics(2,1)=lx1
!        statistics(2,2)=lx2
!        statistics(2,3)=lx3
!        statistics(2,4)=theta
!        statistics(2,5)=phi
!        statistics(2,6)=psi
!        statistics(2,7)=volume
!       endif
       enddo
      enddo

      deallocate (mnp)
      deallocate (mnp3)

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine norm(v)

      double precision, intent(in out) :: v(3)
      double precision :: mag

      mag=dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))

      v(1)=v(1)/mag
      v(2)=v(2)/mag
      v(3)=v(3)/mag

      end

      subroutine norm1(v1,v2)

      double precision, intent(in) :: v1(3)
      double precision, intent(in out) :: v2(3)
      double precision :: mag

      mag=dsqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))

      v2(1)=v1(1)/mag
      v2(2)=v1(2)/mag
      v2(3)=v1(3)/mag

      end

      subroutine dot(v1,v2,dprod0)

      double precision, intent(in) :: v1(3), v2(3)
      double precision, intent(in out) :: dprod0

      dprod0=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)

      end



      subroutine cross(v1,v2,v3)

      double precision, intent(in) :: v1(3), v2(3)
      double precision, intent(in out) :: v3(3)

      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)   !x component
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)   !y component
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)   !z component

      end
        












