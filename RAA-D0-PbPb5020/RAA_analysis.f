
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     RAA_analysis.f version 1.0
c
c     This program will read data from HM.dat file and use those 
c     data to do some caculation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program CROSS_SECTION

	implicit none

	integer a,b,c,a8,b8,aa,bb,enumber1,enumber2,enumber3,enumber4
	integer SIZE1
	integer SIZE2
	parameter (SIZE1=2000000)
	parameter (SIZE2=2000000)
      double precision  id(SIZE1),KATT(SIZE1),p_px(SIZE1),p_py(SIZE1),p_pz(SIZE1),
     &                       p_p0(SIZE1),p_mass(SIZE1),p_rx(SIZE1),p_ry(SIZE1),
     &                       p_rz(SIZE1),p_r0(SIZE1),Thydro(SIZE1),c_vx(SIZE1),
     &                       c_vy(SIZE1),c_vz(SIZE1),p_ipT(SIZE1),p_wt(SIZE1),
     &                       col(SIZE1),acol(SIZE1),stat(SIZE1)

      double precision id8(SIZE2),KATT8(SIZE2),p_px8(SIZE2),p_py8(SIZE2),p_pz8(SIZE2),
     &                       p_p08(SIZE2),p_mass8(SIZE2),p_rx8(SIZE2),p_ry8(SIZE2),
     &                       p_rz8(SIZE2),p_r08(SIZE2),Thydro8(SIZE2),c_vx8(SIZE2),
     &                       c_vy8(SIZE2),c_vz8(SIZE2),p_ipT8(SIZE2),p_wt8(SIZE2)
	integer KATTT(SIZE1);
	double precision px(SIZE1),py(SIZE1),pz(SIZE1),p0(SIZE1),pmass(SIZE1),p5(SIZE1);
	double precision Vfrozenx(SIZE1),Vfrozeny(SIZE1),Vfrozenz(SIZE1),Vfrozen0(SIZE1);
	double precision vcfrozenx(SIZE1),vcfrozeny(SIZE1),vcfrozenz(SIZE1);
	double precision Tfrozen(SIZE1);
	double precision WT(SIZE1),dummyInt;
	double precision PT(SIZE1),PT8(SIZE2),PTT(SIZE1)
	double precision eta(SIZE1),eta8(SIZE2),etaa(SIZE1)
	integer i,j,k,l,m1
	integer i1,i2,i3
	double precision njsum,njsum8
	double precision ptsum(0:SIZE1),ptsum8(0:SIZE2)
	double precision njsum0
	double precision ptsum0(0:SIZE1)
	double precision ptmin
	double precision ptx
	double precision aaa, bbb, ccc

	integer N
	parameter (N=200)
	double precision countlength

	double precision weight(SIZE1)	
	double precision weight8(SIZE1)
	double precision pt_length
	integer ind_pt
	double precision weight_int
	integer ind_lt
	double precision weight_pt
	integer mmm,nnn
	integer status,status1
	double precision HF
	double precision weightsum8,weightsum

	double precision PTIMAX, PTIMIN
	double precision Nevent, Nevent8
        integer partonID,hadronID

        double precision PI

        double precision weightstart

        double precision ptbincms(29)
ccc        DIMENSION ptbincms(36)

        DATA ptbincms /
c     &  2d0, 3d0, 4d0, 5d0, 6d0, 8d0, 10d0, 12.5d0,
c     &  15d0, 20d0, 25d0, 30d0, 40d0, 60d0, 100d0, 120d0/
     &  0d0, 1d0, 2d0, 4d0, 6d0, 8d0, 10d0, 12d0,
     &  14d0, 16d0, 20d0, 24d0, 28d0, 32d0, 36d0, 40d0,
     &  50d0, 60d0, 70d0, 80d0, 100d0, 120d0, 140d0, 160d0,
     &  200d0, 240d0, 280d0, 320d0, 360d0/
	PTIMAX=2000d0
	PTIMIN=3d0
        Nevent=0d0
        Nevent8=0d0
        partonID=4
        hadronID=13
c        muonID=13

        weightsum8=0d0
        weightsum=0d0

        do i=0,SIZE1
           ptsum(i)=0d0
           ptsum8(i)=0d0
        enddo

	countlength = 600d0
	
	weight_int = 0.2d0
        weightstart=0.1d0

        ptmin=3d0	
        PI=3.141592653589793d0


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	OPEN(UNIT=40,FILE="../lpposi_Frag/hadrons-nega.dat")
	OPEN(UNIT=41,FILE="../lpposi_Frag/hadrons-posi.dat")
        OPEN(UNIT=42,FILE="../pp_lpposi_Frag/hadrons-posi.dat")      
 
c	OPEN(UNIT=40,FILE="/data/ycdang/particle/hadron/b/pt_hat_80-100/AA/hadrons-nega.1198210.dat")
c	OPEN(UNIT=41,FILE="/data/ycdang/particle/hadron/b/pt_hat_80-100/AA/hadrons-posi.1198210.dat")
c        OPEN(UNIT=42,FILE="/data/ycdang/particle/hadron/b/pt_hat_80-100/pp/hadrons-posi.1198210.dat")

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c For final hardon (after fragmentation)
c read input lines
	event: do i1=1,30000
		read(40,*) a,b
          	write(*,*) "a=",a,"b=",b
		do j=1,b
			read(40,*) KATT(j),p_px(j),p_py(j),p_pz(j),p_p0(j),p_mass(j),p_rx(j),p_ry(j),p_rz(j),p_r0(j)
c			read(40,*) id(j),KATT(j),p_px(j),p_py(j),p_pz(j),p_p0(j),p_mass(j),p_rx(j),p_ry(j),p_rz(j),p_r0(j),p_ipT(j),p_wt(j)
			PT(j) = sqrt(p_px(j)*p_px(j) + p_py(j)*p_py(j))
			eta(j) = (1.0 / 2.0) * dlog((p_p0(j) + p_pz(j)) / (p_p0(j) - p_pz(j)))
			
		enddo

	  	read(41,*) a8,b8
          	write(*,*) "a8=",a8,"b8=",b8
		do j=1,b8
			read(41,*) KATT8(j),p_px8(j),p_py8(j),p_pz8(j),p_p08(j),p_mass8(j),p_rx8(j),p_ry8(j),p_rz8(j),p_r08(j) 
c			read(41,*) id8(j),KATT8(j),p_px8(j),p_py8(j),p_pz8(j),p_p08(j),p_mass8(j),p_rx8(j),p_ry8(j),p_rz8(j),p_r08(j),p_ipT8(j),p_wt8(j)
			PT8(j) = sqrt(p_px8(j)*p_px8(j) + p_py8(j)*p_py8(j))
			eta8(j) = (1.0 / 2.0) * dlog((p_p08(j) + p_pz8(j)) / (p_p08(j) - p_pz8(j)))
		enddo

                read(42,*) aa,bb
                write(*,*) "aa=",aa,"bb=",bb
                do j=1,bb
                        read(42,*) KATTT(j),px(j),py(j),pz(j),p0(j),pmass(j),Vfrozenx(j),Vfrozeny(j),Vfrozenz(j),Vfrozen0(j) 
c                       read(42,*) dummyInt,KATTT(j),px(j),py(j),pz(j),p0(j),pmass(j),Vfrozenx(j),Vfrozeny(j),Vfrozenz(j),Vfrozen0(j),dummyInt,dummyInt 
                        PTT(j) = sqrt(px(j)*px(j) + py(j)*py(j))
                        etaa(j) = (1.0 / 2.0) * dlog((p0(j) + pz(j)) / (p0(j) - pz(j)))
                enddo

c count number(for P) for each bin

		do k=1,b
                        if((abs(KATT(k)).eq.321).or.(abs(KATT(k)).eq.321).or.(abs(KATT(k)).eq.211).or.(abs(KATT(k)).eq.2212)) then
                        if((eta(k).ge.-1).and.(eta(k).le.1)) then
			do l=1,28
				if( (PT(k).ge.ptbincms(l)).and.(PT(k).lt.ptbincms(l+1)) ) then
					njsum = njsum + 1.0
					ptsum(l) = ptsum(l) + 1.0
                                endif
			enddo
                        endif
                        endif
		enddo

		do k=1,b8
                        if((abs(KATT8(k)).eq.321).or.(abs(KATT8(k)).eq.321).or.(abs(KATT8(k)).eq.211).or.(abs(KATT8(k)).eq.2212)) then
                        if((eta8(k).ge.-1).and.(eta8(k).le.1)) then
			do l=1,28
				if( (PT8(k).ge.ptbincms(l)).and.(PT8(k).lt.ptbincms(l+1)) ) then

!				write(*,*) pt_length, weight_pt
				njsum8 = njsum8 + 1.0
				ptsum8(l) = ptsum8(l) + 1.0

                                endif
			enddo
                        endif
                        endif
		enddo

                do k=1,bb
                        if((abs(KATTT(k)).eq.321).or.(abs(KATTT(k)).eq.211).or.(abs(KATTT(k)).eq.321).or.(abs(KATTT(k)).eq.2212)) then
                        if((etaa(k).ge.-1).and.(etaa(k).le.1)) then
                        do l=1,28
                                if( (PTT(k).ge.ptbincms(l)).and.(PTT(k).lt.ptbincms(l+1)) ) then
                                        njsum0 = njsum0 + 1.0
                                        ptsum0(l) = ptsum0(l) + 1.0
                                endif
                        enddo
                        endif
                        endif
                enddo
                enumber1 = enumber2
                enumber2 = enumber3
                enumber3 = enumber4
                enumber4 = b8
                if ((enumber1 .eq. enumber2).and.(enumber2 .eq. enumber3).and.(enumber3 .eq. enumber4)) then
                    exit event
                endif
	end do event


	
        write(*,*) "Nevent is: ", Nevent
        write(*,*) "Nevent8 is: ", Nevent8
	
c normalization
	open(UNIT=22,FILE="pT-pp.dat")
        open(UNIT=23,FILE="pT-AA.dat")
        write(22,'(I5,1X,I1)') a8,0
        write(23,'(I5,1X,I1)') aa,0


	do m1=1,28
		ptx = ( ptbincms(m1)+ptbincms(m1+1) ) /2d0
                write(*,*) m1,ptsum(m1),ptsum8(m1)
ccc///		ptsum(m1)= ptsum(m1)/((countlength/(N*1.0)))*((PTIMAX-PTIMIN)/Nevent)
ccc///                ptsum(m1) = ptsum(m1) / (2d0*PI*ptx)
ccc///		ptsum8(m1)=ptsum8(m1)/((countlength/(N*1.0)))*((PTIMAX-PTIMIN)/Nevent8)
ccc///                ptsum8(m1) = ptsum8(m1) / (2d0*PI*ptx)
                aaa=ptsum(m1)/(ptbincms(m1+1)-ptbincms(m1))/(4d0*PI*ptx)
                bbb=ptsum8(m1)/(ptbincms(m1+1)-ptbincms(m1))/(4d0*PI*ptx)
                ccc=ptsum0(m1)/(ptbincms(m1+1)-ptbincms(m1))/(4d0*PI*ptx)
                write(22,'(E10.4,1X,E10.4)') ptx,ccc
                write(23,'(E10.4,1X,E10.4)') ptx,bbb-aaa
	enddo
	write(*,*) njsum
	write(*,*) njsum8

      write(6,*) "PROGRAM ENDS SUCCESSFULLY :)"

 7      format(6(1PE11.3))

      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
