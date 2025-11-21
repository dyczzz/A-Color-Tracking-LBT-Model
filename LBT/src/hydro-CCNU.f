      subroutine read_CCNU(dataFN_in)
       
      implicit none

      character (len=*) :: dataFN_in
      double precision CCNU_min_tau,CCNU_max_tau,CCNU_min_xy,CCNU_max_xy
      double precision CCNU_delta
      double precision CCNU_delta_xy !!added by WJ
      integer tau_nstep,x_nstep,y_nstep
      parameter(tau_nstep=100,x_nstep=65,y_nstep=65)
      double precision CCNU_temp(0:tau_nstep,0:x_nstep,0:y_nstep),
     &                 CCNU_vx(0:tau_nstep,0:x_nstep,0:y_nstep),
     &                 CCNU_vy(0:tau_nstep,0:x_nstep,0:y_nstep)
      common/CCNU_info/ CCNU_temp,CCNU_vx,CCNU_vy
      common/CCNU_grid/ CCNU_min_tau,CCNU_max_tau,CCNU_min_xy,
     &                  CCNU_max_xy,CCNU_delta,CCNU_delta_xy

      double precision dummy_float
      integer dummy_int,i,j,k

      CCNU_min_tau = 0.6d0
      CCNU_min_xy = -9.75d0
      CCNU_max_xy = 9.75d0
      CCNU_delta = 0.3d0 !!for PbPb5020 medium
!      CCNU_delta = 0.4d0  !!for XeXe5440 medium
      CCNU_delta_xy = 0.3d0


      open(unit=9,file=dataFN_in,status='old',form='formatted')

      i = 0
      do while(i.le.tau_nstep)
         j = 0
         do while(j.le.x_nstep)
            k = 0 
            do while(k.le.y_nstep)

               read(unit=9,fmt=*,err=101,end=102) dummy_float,
     &            dummy_float,dummy_float,dummy_float,CCNU_temp(i,j,k),
     &            CCNU_vx(i,j,k),CCNU_vy(i,j,k),dummy_float,dummy_float

!               write(unit=6,fmt=*) i*CCNU_delta+0.6,
!     &                                       j*CCNU_delta_xy-9.75,
!     &                                        k*CCNU_delta_xy-9.75
!               write(unit=6,fmt=*) i,j,k,i*(x_nstep+1)*(y_nstep+1)
!     &            +j*(x_nstep+1)+k
!               write(unit=6,fmt=*) "temperature: ", CCNU_temp(i,j,k)

               k = k+1
            end do
            j = j+1
         end do
         i = i+1
      end do 

      close(9)

      return
 

 101  continue
      write(6,*) 'ERROR reached in CCNU-hydro file'
      write(6,*) 'terminating ...'
      stop
      return

 102  continue

c      write(6,*) "CCNU hydro has been read in successfully."
c      write(6,*) CCNU_temp(i-1,x_nstep,y_nstep),
c     &           CCNU_vy(i-1,x_nstep,y_nstep)

      write(6,*) 'EOF reached in CCNU-hydro file'
      if (k.eq.0.and.j.eq.0) then
         CCNU_max_tau = (i-1)*CCNU_delta+CCNU_min_tau
         write(6,*) "Max tau in the hydro file: ",CCNU_max_tau
      else
         write(6,*) 'terminating ...'
         stop
      end if

      return

      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine hydroInfoCCNU(Ct,Cx,Cy,Cz,Ctemp,Cvx,Cvy,Cvz,Cflag)

      implicit none

      double precision CCNU_min_tau,CCNU_max_tau,CCNU_min_xy,CCNU_max_xy
      double precision CCNU_delta
      double precision CCNU_delta_xy
      integer tau_nstep,x_nstep,y_nstep
      parameter(tau_nstep=100,x_nstep=65,y_nstep=65)
      double precision CCNU_temp(0:tau_nstep,0:x_nstep,0:y_nstep),
     &                 CCNU_vx(0:tau_nstep,0:x_nstep,0:y_nstep),
     &                 CCNU_vy(0:tau_nstep,0:x_nstep,0:y_nstep)
      common/CCNU_info/ CCNU_temp,CCNU_vx,CCNU_vy
      common/CCNU_grid/ CCNU_min_tau,CCNU_max_tau,CCNU_min_xy,
     &                  CCNU_max_xy,CCNU_delta,CCNU_delta_xy

      double precision Ct,Cx,Cy,Cz,Ctemp,Cvx,Cvy,Cvz
      integer Cflag
      double precision Ctau,Cgamma
      integer i_tau,i_x,i_y

      Ctau =  sqrt(Ct**2-Cz**2)

      if(Ctau.lt.(CCNU_min_tau-1D-6).or.Ctau.gt.CCNU_max_tau) then
         Cflag = 1
         return
      endif

      if(Cx.lt.CCNU_min_xy.or.Cx.gt.CCNU_max_xy) then
         Cflag = 1
         return
      endif

      if(Cy.lt.CCNU_min_xy.or.Cy.gt.CCNU_max_xy) then
         Cflag = 1
         return
      endif

      i_tau = int((Ctau-CCNU_min_tau)/CCNU_delta+0.5d0)
!      i_x = int((Cx-CCNU_min_xy)/CCNU_delta+0.5d0)
!      i_y = int((Cy-CCNU_min_xy)/CCNU_delta+0.5d0)

      i_x = int((Cx-CCNU_min_xy)/CCNU_delta_xy+0.5d0)
      i_y = int((Cy-CCNU_min_xy)/CCNU_delta_xy+0.5d0)


      Ctemp = CCNU_temp(i_tau,i_x,i_y)    
      Cvx = CCNU_vx(i_tau,i_x,i_y)    
      Cvy = CCNU_vy(i_tau,i_x,i_y)

      Cvz = Cz/(Ct+1d-30)
      Cgamma = 1d0/(sqrt(1d0-Cvz*Cvz)+1D-30)
      Cvx = Cvx/Cgamma
      Cvy = Cvy/Cgamma
     
      Cflag = 0
      return

      end subroutine
