!

!2007.9.20 the spectrum without Compton scattering  can be derived by 
! setting: Ntest=1.

	
C THIS CODE CALCULATE THE EMERGENT SPECTRUM FROM AN ADAF-SSD SYSTEM. 
C THE SOFT PHOTONS FROM SSD IS TAKEN INTO ACCOUNT AS SEED PHOTONS OF
C COMPTONIZATION; THE REPROCESSING OF X-RAY PHOTONS BY SSD IS TAKEN 
C INTO ACCOUNT 

c description of using this code:

c 1. the code should be run twice to get a spectrum, with different ranges of
c    frequencies. when focused on the Comptonization of BB photons,
c    you should set 'qbr=0'; and in the subroutine 'fls', appropriate
c    range of frequency should be set; when focused on the COmptonization
c    of the syn. photons, comment out 'qbr=0', and choose appropriate frequency 
c 2. Compared to the original code to calculate the spectrum of ADAF,
c    two changes here: 1. consider the seed photonn of BB from SSD; and 2. 
c    consider the reprocessing of ADAF emission by the SSD. The emission from
c    ADAF is absorbed and reflected. In this code, only the former case
c    is included. another seperate code has to be used to calculate the 
c   reflection spectrum. In this code, a parameter to describe the albedo is
c    used and set to 0.8, i.e., most of the emission from ADAF is aabsorbed.

c---------------2003, 6.9-----------------

	IMPLICIT REAL*8(A-H,O-Z)
        dimension r(0:3000),mach(3000),te(3000),ti(3000),
     $          s1(3000),s2(3000),tao(3000),be(3000),nnu(200),nu(200)
     $          ,nulu(200),f(200),rr(0:3000),sigma(0:1000),nulu2(200)
     $		,nulu3(200),ssum(200),sssum(200),ssssum(200),nulu4(200)   
     $		,nulu5(200),sssssum(200),tau(200),ssuumm(200),nulu6(200)
     $	,nulu7(200),ssssss(200),a1(30),b1(10),c1(30,10,1),z(1),su(200)
	dimension a2(100),b2(80),c2(100,80,1),z2(1),a3(100),b3(80),
     $		c3(100,80,1),z3(1)
	double precision mach,nu,nnu,nulu,mui,mue,m,mbsl4,kapa,mc2h
     $		,nulu2,nulu3,nulu4,nulu5,nulu6,nulu7,ssum

c *** Added y3 and y4 to a new common block. They are respectively y1 and y2 in
c the fls subroutine
	common /newnew/y3,y4
	common mui,mue,m,dotm0,rtran,beta,mc2h
	common /comp1/nu,ssum
	common /comp2/a1,b1,c1,a2,b2,c2,a3,b3,c3
  	data(nu(i),i=1,200)/200*0.d0/(ssum(i),i=1,200)/200*0.d0/

	open(33,file='romi-1.dat')
	open(34,file='romi-2.dat')
	open(37,file='romi-3.dat')

	do 41 i=1,30
41      read(33,*) a1(i)
        do 42 i=1,10
42      read(34,*) b1(i)

        do 43 j=1,10
        do 44 i=1,30
        read(37,*) c1(i,j,1)
44      continue
43      continue

	open(43,file='aomi-1.dat')
        open(44,file='aomi-2.dat')
        open(47,file='aomi-3.dat')

        do 51 i=1,100
51      read(43,*) a2(i)
        do 52 i=1,80
52      read(44,*) b2(i)

        do 53 j=1,80
        do 54 i=1,100
        read(47,*) c2(i,j,1)
54      continue
53      continue

	open(53,file='aomi2-1.dat')
        open(54,file='aomi2-2.dat')
        open(57,file='aomi2-3.dat')

        do 61 i=1,100
61      read(53,*) a3(i)
        do 62 i=1,80
62      read(54,*) b3(i)

        do 63 j=1,80
        do 64 i=1,100
        read(57,*) c3(i,j,1)
64      continue
63      continue

c ratio of gas to total pressure
c	print *,'beta:'
	read(*,*) beta
c	print *,'#beta',beta
c        beta=0.9d0	!---------------------------1--------------------
        ru=9.2d-14

c	print *,'Black hole mass (in Solar masses):'
	read(*,*) m
c	print *,'#m',m
c	m=1.d2	!---------------------------2-------------------

c	print *,'Distance (in pc):'
	read(*,*) distance
	distance=distance*3.09d18
c	distance=5.3d3*3.09d18	!-------------------(3)------------------
c !!!!!!!!!!!!!!

        mc2h=1.236d+20
	mui=1.24d0
        mue=1.13d0
c	nulu=0.d0
	n=0

c	print *,'alpha:'
	read(*,*) alfa
c	print *,'#alfa',alfa
c	alfa=0.3d0	!---------------------------4---------------------

c change below!!!!
	dotme=3.44d-15*m
! if there is no outflow, this value does not affect the spectrum.
c	print *,'Mdot_out (Eddington units):'
	read(*,*) dotm0
c	print *,'#dotm',dotm0
c       dotm0=dotme*1.d-1	!--------------------(5)-----------------------
	
c	print *,'R_out (units of R_S):'
	read(*,*) rtran	
	rtran=2.*rtran
c	print *,'#rtran',rtran
c	rtran=200.d0	!---------------------------(6)-----------------------

c *** Added: reads from stdinput the limits of second Comptonization (y1=y3 and y2=y4 defined
c inside subroutine fls).
	read(*,*) y3
	read(*,*) y4
c	print *,'#nui',y3
c	print *,'#nuf',y4
	
c *** Reads variable qbreset. 0 if qbr is untouched, 1 if qbr=0.
	read(*,*) qbreset
c	print *,'#qbreset',qbreset

	open(23,file='spectrum.dat')
!	open(23,file='s5.dat')
	open(24,file='ses5.dat')
        open(13,file='x.dat')
!        open(13,file='2-3.dat')
	open(14,file='hot.dat')
c 13 and 14 are two input files

! before plot the spectra, we need to see where the synchrotron peak is.
! Ntest is an integer to do this. ntest=1 the results is calculated without 
! considering the Compton heating. ntest=0 considering.
c	Ntest=0    

c *** Reads variable Ntest from STDIN
	read(*,*) Ntest  !-----------------------0-----------------------!
c	print *,'#ntest',Ntest
	nvmax=100 ! number of points in SED *** was 80!
	if (ntest .eq. 0) then
	   jmin=5
	else if (ntest .eq. 1) then
	   jmin=nvmax+1
	end if
c	print*,jmin

c *** Nlines below is the number of lines in the file created by the
c dynamics code. Previously the code had a fixed value for Nlines, but now
c I determine it self-consistently with the Perl script.
	read(*,*) Nlines
c	print *,'#nlines',Nlines

        do 11 i=1,Nlines	!-----------------------------7---------------------------
        n=n+1
c now input the data

        read(13,*) rr(i),mach(i),te(i),ti(i),s1(i),s2(i)
     $          ,tao(i),be(i),s1(i),s1(i),s1(i),s1(i)
        r(i)=10.d0**rr(i)*2.d0
	te(i)=10.d0**te(i)

        ti(i)=10.d0**ti(i)
	tao(i)=10.d0**tao(i)*1.0

c	if(r(i).lt.15.d0) then
c       te(i)=te(i)*1.44*(r(i)*2.)**(-0.2d0)
c       tao(i)=tao(i)*1.806*(r(i)*2.)**(-0.1738d0)
c       endif

	r0=r(1)
	do 12 j=1,i
	cs=sqrt(ru/beta*(ti(j)/mui+te(j)/mue))
        omigak=1./sqrt(r(j))/(r(j)-2.d0)
        h=cs/omigak/6.77d-12*m
c now h is in units of cm.g.s
        rho=tao(j)/h
        tau(i)=tau(i)+rho*(r0-r(j))*0.4/6.77d-12*m
	r0=r(j)
12	continue
c	print*,'tau(i)',tau(i)

c tao: the optical depth in the vertical direction
c tau: the optical depth in the radial direction

	aaa=sqrt((4.+1.33*alfa*alfa)/(8./3.))
	mach(i)=mach(i)*aaa
c-----------------------------
c note: tao is dimensionless
c n presents the number of the rings of the ADAF disk
c-----------------------------
11      continue

        r(0)=r(1)
c       r(n+1)=r(n)/2.
	    r(n+1)=2.d0


        ! biggest loop, goes through radial structure?
        ! ==============================================
        do 20 i=1,n
        print*,i

        cs=sqrt(1.d0/beta*ru*(ti(i)/mui+te(i)/mue))
        omigak=1./sqrt(r(i))/(r(i)-2.d0)
        h=cs/omigak

c now transfer them form c=G=M=1 units to cgs units

        hh=h/6.77d-12*m
        rhorho=tao(i)/hh

        taoes=rhorho*hh*0.4d0
        s=taoes+taoes*taoes

c	print*,cs,omigak,h,hh,rhorho,taoes,s
c	stop

c the first, the breamsstrahlung radiation

        setae=1.686d-10*te(i)
        setai=9.181d-14*ti(i)

        if(setae.lt.1.d0) then
        fseta=4.d0*sqrt(2.*setae/3.1416)/3.1416*(1.+1.781*setae**1.34d0)
     $          +1.73*setae**1.5d0*(1.+1.1*setae+setae*setae
     $          -1.25*setae**2.5d0)
        else
        fseta=9.d0*setae/2./3.1416*(log(1.123d0*setae+0.48)
     $          +1.5)+2.3d0*setae*
     $          (log(1.123d0*setae)+1.28)
        endif

        qbr=5.3d+25*rhorho*rhorho*fseta/mue/mue

c	print*,hh,rhorho,setae,setai
c	stop


c *** The statement "qbr=0" below must be commented out when focused on the 
c Comptonization of Syn. photons. To do this I introduced qbreset variable which is 
c read from stdinput. 
        if (qbreset .eq. 1.) then
		qbr=0.d0
	endif

c------------------------------------------------------------------
c the folowing calculate the luminosity at given frequency nu(j)  |
c emissioned from per accretion disk ring within r(i-1) --> r(i)  |
c------------------------------------------------------------------
c second big loop, goes through frequencies
c
        do 30 j=1,nvmax
! the value of dlog(nu) below must match the value of nvmax (number of
! steps in SED)
c *** previously was nnu(j)=10.d0**(10.1+0.15d0*j) 
        nnu(j)=10.d0**(9.0+0.12d0*j)	!----------------------(8)------------------
c now we consider the motion and the gravitational redshift
        if(mach(i)*cs.ge.1.) then
          redshift=9.d-1
c          print*,'mach*cs > 1!', mach(i)
        else
          redshift=sqrt(1.-2./r(i))*(1.-mach(i)*cs*mach(i)*cs)
        endif

	redshift=1.d0
        nu(j)=nnu(j)/redshift
c nu(j) present the 'unredened' frequency

c now calculate the black body radiation spectrum
c       bnu=1.47d-47*nu(j)*nu(j)*nu(j)/(dexp(4.8d-11*nu(j)/te(i))-1.d0)

	if(4.8d-11*nu(j)/te(i).lt.100.d0) then
        bnu=1.47d-47*nu(j)*nu(j)*nu(j)/(dexp(4.8d-11*nu(j)/te(i))-1.d0)
        else
	bnu=1.47d-47*nu(j)*nu(j)*nu(j)/((dexp(4.8d-11*nu(j)/te(i)/3.d0))
     $		**3.d0-1.d0)
        endif

c now calculate the bremstrahlung. gaunt factor
c       if(2.08d+10/nu(j)*te(i).lt.1.d0) then
c       gaunt=sqrt(3./3.1416d0*4.8d-11/nu(j)/te(i))
c       else
c       gaunt=4.8d-11/te(i)*sqrt(3.)/3.1416*
c    $          log(4./1.781*2.08d10*te(i)/nu(j))
c       endif

c now calculate the bremstrahlung. gaunt factor
 	g1=sqrt(3./3.1416d0*4.8d-11/nu(j)/te(i))
	g2=4.8d-11/te(i)*sqrt(3.)/3.1416*
     $          log(4./1.781*2.08d10*te(i)/nu(j))

        if(2.08d+10/nu(j)*te(i).lt.0.8d0) then
        gaunt=g1
        else
	   if(2.08d+10/nu(j)*te(i).lt.2.d0) then
           gaunt=sqrt(g1*g2)
	   else
	   gaunt=g2
	   endif
        endif

	if(4.8d-11*nu(j)/te(i).lt.50.d0) then
        qiabr=qbr*gaunt*dexp(-4.8d-11*nu(j)/te(i))
	else
	qiabr=qbr*gaunt*(dexp(-4.8d-11*nu(j)/te(i)/3.d0))**3.d0
	endif

c second, the synchrotron radiation

        xx=5.29d-12*nu(j)/setae/setae/sqrt(rhorho)/sqrt((1.-beta)/beta)
     $          /sqrt(ti(i)/mui+te(i)/mue)
        fxx=4.0505/xx**(1./6.d0)*(1.+0.4d0/xx**0.25d0+
     $          0.5316d0/sqrt(xx))*exp(-1.8899d0*xx**0.33333333333d0)

        fk2=mbsl4(2,1.d0/setae)

        qiasy=33.33482155d-6*rhorho*nu(j)*fxx/fk2/mue
c consider the radiation transfer
        qianu=qiabr+qiasy

        kapa=qianu/4.d0/3.1416/bnu
        taonux=sqrt(3.1416)/2.*kapa*hh
c here, taonux denotes taonu star
c fnu below denotes the luminosity per unit frequency radiated
c from per unit surface of the disk

        if(taonux.gt.1.d-6) then
        fnu=2.*3.1416d0/sqrt(3.d0)*bnu*
     $          (1.-dexp(-2.d0*sqrt(3.)*taonux))
        else
        fnu=2.*3.1416d0/sqrt(3.d0)*bnu*(2.*sqrt(3.)*taonux-
     $          6.*taonux*taonux+4.*sqrt(3.)*taonux**3.d0)
	endif


c	write(*,'(A,8e12.4)') 'fnu', nu(j),fk2,fxx,xx,qiabr,bnu,rhorho


c Transfer it to per ring the disk
c This breaks the independence between loop through the disk,
c since r[i] depends on r[i-1]
c
        sigma1=3.1416*((r(i-1)+r(i))**2.d0/4.d0-(r(i)+r(i+1))
     $          **2.d0/4.d0)/(6.77*6.77d-24)*m*m
        sigma2=3.1416*(r(i)*r(i)-r(i+1)*r(i+1))/(6.77*6.77d-24)*m*m
        f(j)=2.*fnu*sigma1
        sigma(i)=sigma1

c-------------------------------------------------------------------
c the following integration calculate the comptonization spectrum  |
c-------------------------------------------------------------------

      if(j.ge.jmin .and. j.le.(nvmax-10)) then
	  eps=1.d-2
	  fup=8.d+1		!-----------------------------(9)-----------------
	  flo=1.0d0*(1.+1.d-2)

c here fup, flo denotes the Lorentz factor of electrons

          call simps2(flo,fup,eps,kflag,n,sum,
     $		rhorho,hh,setae,setai,fk2,nu(j),r(i),tau(i))
	  if(sum.lt.0.) write(*,*) 'issue: sum<0; line 389 in spectrum.f'
c      pause 'sum<0!!!'
	  else
       sum=0.d0
      endif

c	print*,nulu(j),f(j),nu(j),sigma(i)

	nulu(j)=nulu(j)+(nu(j)*f(j)+nu(j)*sum*2.*sigma(i))
     $          *redshift**2.d0

	su(j)=su(j)+nu(j)*f(j)

	ssum(j)=sum*1.
30      continue
! end of big loop over frequencies
!

c---------------------------------------------------
c the following calculate the second componization |
c---------------------------------------------------
	if(jmin .ge. nvmax)  goto 20

!$OMP PARALLEL DO PRIVATE(jj,sum2) SHARED(sssum) 
!$OMP& REDUCTION(+:nulu2)
	do 34 jj=jmin,nvmax
c       call simp2(1.01d0,6.d+1,0.2d0,sum2,
c    $		rhorho,hh,setae,jj,nu,ssum)
	call simp2(flo,fup,0.2d0,sum2,
     $          rhorho,hh,setae,jj)

	nulu2(jj)=nulu2(jj)+nu(jj)*sum2*2.*sigma(i)
     $          *redshift**2.d0

 	sssum(jj)=sum2*1.
34	continue
!$ OMP END PARALLEL DO

	do 134 l=1,200
	ssum(l)=sssum(l)
134	continue
c---------------------------------------------------
c the following calculate the third comptonization |
c---------------------------------------------------

!$OMP PARALLEL DO PRIVATE(jj,sum3) SHARED(ssssum) 
!$OMP& REDUCTION(+:nulu3)
	do 35 jj=jmin+5,nvmax
c       call simp2(1.01d0,6.d+1,0.2d0,sum3,
c    $          rhorho,hh,setae,jj,nu,sssum)
	call simp2(flo,fup,0.2d0,sum3,
     $          rhorho,hh,setae,jj)

        nulu3(jj)=nulu3(jj)+nu(jj)*sum3*2.*sigma(i)
     $          *redshift**2.d0

 	ssssum(jj)=sum3*1.
35      continue
!$ OMP END PARALLEL DO

	do 135 l=1,200
        ssum(l)=ssssum(l)
135     continue
c---------------------------------------------------
c the following calculate the fourth comptonization |
c---------------------------------------------------

!$OMP PARALLEL DO PRIVATE(jj,sum4) SHARED(sssssum) 
!$OMP& REDUCTION(+:nulu4)
        do 36 jj=jmin+15,nvmax
c       call simp2(1.01d0,8.d+1,0.2d0,sum4,
c     $          rhorho,hh,setae,jj,nu,ssssum)
	call simp2(flo,fup,0.2d0,sum4,
     $          rhorho,hh,setae,jj)

        nulu4(jj)=nulu4(jj)+nu(jj)*sum4*2.*sigma(i)
     $          *redshift**2.d0
 	sssssum(jj)=sum4*1.
36      continue
!$ OMP END PARALLEL DO

c	print*,'nulu4'
	do 136 l=1,200
        ssum(l)=sssssum(l)
136     continue
c---------------------------------------------------
c the following calculate the fiveth comptonization |
c---------------------------------------------------

!$OMP PARALLEL DO PRIVATE(jj,sum5) SHARED(ssuumm) 
!$OMP& REDUCTION(+:nulu5)
        do 37 jj=jmin+15,nvmax
        call simp2(flo,fup,0.2d0,sum5,
     $          rhorho,hh,setae,jj)

        nulu5(jj)=nulu5(jj)+nu(jj)*sum5*2.*sigma(i)
     $          *redshift**2.d0
 	ssuumm(jj)=sum5*1.
37      continue
!$ OMP END PARALLEL DO

	do 137 l=1,200
        ssum(l)=ssuumm(l)
137     continue
c---------------------------------------------------
c the following calculate the sixth comptonization |
c---------------------------------------------------

!$OMP PARALLEL DO PRIVATE(jj,sum6) SHARED(ssssss) 
!$OMP& REDUCTION(+:nulu6)
        do 38 jj=jmin+20,nvmax
c       call simp2(1.01d0,8.d+1,0.2d0,sum6,
c    $          rhorho,hh,setae,jj,nu(jj),ssuumm(j))
	call simp2(flo,fup,0.2d0,sum6,
     $          rhorho,hh,setae,jj)

        nulu6(jj)=nulu6(jj)+nu(jj)*sum6*2.*sigma(i)
     $          *redshift**2.d0

 	ssssss(jj)=sum6*1.
38      continue
!$ OMP END PARALLEL DO

	do 138 l=1,200
        ssum(l)=ssssss(l)
138     continue

c---------------------------------------------------
c the following calculate the seventh comptonization |
c---------------------------------------------------

!$OMP PARALLEL DO PRIVATE(jj,sum7) REDUCTION(+:nulu7)
        do 39 jj=jmin+35,nvmax
        call simp2(flo,fup,0.2d0,sum7,
     $          rhorho,hh,setae,jj)

        nulu7(jj)=nulu7(jj)+nu(jj)*sum7*2.*sigma(i)
     $          *redshift**2.d0

39      continue
!$ OMP END PARALLEL DO


20      continue
!$ OMP END PARALLEL DO
! end of big loop over radial structure
!


	totlu=0.d0
        do 40 k=1,nvmax
	write(23,50) log10(nnu(k)),log10(nulu(k)+nulu2(k)
     $		+nulu3(k)+nulu4(k)+nulu5(k)+nulu6(k)+nulu7(k))
c     $			-log10(4.*3.1416*distance**2.)
	write(24,50) log10(nnu(k)),log10(su(k))
40      continue
50      format(2e15.6)

	Do k=2,nvmax
	totlu=totlu+(nulu(k)+nulu2(k)+nulu3(k)+nulu4(k)+nulu5(k)
     $		+nulu6(k)+nulu7(k))*log(nnu(k)/nnu(k-1))
	end do
	print*,'the total luminosity is: ',totlu

        END











	subroutine fls(x,y1,y2)
	implicit real*8(a-h,o-z)
	common /newnew/y3,y4
c x: the velocity; y: the frequency
c the following low and upper limits are very important to clculate the
c  flux of Comptonzization of, espessially, the second Comptonization
c       y1=4.d+13
c       y2=2.d+15

c *** Now y1 and y2 are taken from stdinput
	y1=y3	!----------------------------10-------------------
 	y2=y4

        return
        end











        function f(x,y,rhorho,hh,setae,setai,fk2,nu,radaf,tau)
	IMPLICIT REAL*8(A-H,O-Z)
	dimension a1(30),b1(10),c1(30,10,1),z(1)
	dimension a2(100),b2(80),c2(100,80,1),z2(1)
	dimension a3(100),b3(80),c3(100,80,1),z3(1)
        double precision nu,nunu,mui,mue,kapa,mc2h
     $  ,hea,aomi,m,aomi2
	common /comp2/a1,b1,c1,a2,b2,c2,a3,b3,c3
	common mui,mue,m,dotm0,rtran,beta,mc2h
        c=3.d+10
        te=setae/1.686d-10
        ti=setai/9.181d-14

        omip=y/mc2h
        omi=nu/mc2h
	
	gam=x
c gam is the Lorentz factor of electrons, be=v/c
	be=sqrt(1.-1./gam/gam)
c first, the scattering photo distribution P
c       aomi=4./3.*gam*gam*omip
c now use 2.12 to calculate aomi

c	if(omip.lt.5.d-4) then
c	  if(gam.gt.1.3) then
c	    aomi=4./3.*gam*gam*omip
c	  else
c	    aomi=(1.+4./3.*be*be-omip)*omip
c	  endif
c	else

c	 call lar2(100,80,1,omip,gam,a2,b2,c2,z2)
c         aomi=z2(1)

c        endif

 	aomi=4./3.*gam*gam*omip

c we found for the 1st scattering, we can use the above approximation

c now calculate aomi2 using eq. 2.15
	
 	if(omip.lt.2.d-6) then
          aomi2=14./5.*gam**4.d0*omip*omip*(1.-176/35.*gam*omip)
	aomi2=aomi2/3.5
         else
          call lar2(100,80,1,omip,gam,a3,b3,c3,z3)
          aomi2=z3(1)
         endif

c        aomi2=14./5.*gam**4.d0*omip*omip*(1.-176./35.*gam*omip)
c there is no good approximation to aomi2!!

        deomi=aomi2-aomi*aomi
        d1=sqrt(3.*deomi)
        d2=aomi
        if(d1.lt.d2) then
        dd=d1
        else
        dd=d2
        endif

        if(dd.gt.abs(omi-aomi)) then
        hea=1.d0
        else
        hea=0.d0
        endif

        p=1./2./dd*hea

c now we calculate the initial photon. first, synchrotron 
c  + bremsstrahlung (added 2004/09/14!!)

	nunu=omip*mc2h

        xx=5.29d-12*nunu/setae/setae/sqrt(rhorho)/sqrt((1.-beta)/beta)
     $          /sqrt(ti/mui+te/mue)
        fxx=4.0505/xx**(1./6.d0)*(1.+0.4d0/xx**0.25d0+
     $          0.5316d0/sqrt(xx))*exp(-1.8899d0*xx**0.33333333333d0)
        qiasy=33.33482155d-6*rhorho*nunu*fxx/fk2

c------------
	gg1=4.8d-11/te*sqrt(3./3.1416d0*2.08d10*te/nunu)
	gg2=4.8d-11/te*sqrt(3.)/3.1416*log(4./1.781*2.08d10*te/nunu)

	if(2.08d+10/nunu*te.lt.0.8d0) then
	gaunt=gg1
        else
	  if(2.08d+10/nunu*te.lt.3.d0) then
	  gaunt=sqrt(gg1*gg2)
	  else
          gaunt=gg2
	  endif
        endif

	if(setae.lt.1.) then
        fseta=4.*sqrt(2.*setae/3.1416)/3.1416*(1.+1.781*setae**1.34d0)
     $          +1.73*setae**1.5d0*(1.+1.1*setae+setae*setae
     $          -1.25*setae**2.5d0)
        else
        fseta=9.d0*setae/2./3.1416*(log(1.123*setae+0.48)+1.5)+2.3
     $          *setae*(log(1.123*setae)+1.28)
        endif

        qbr=5.3d+25*rhorho*rhorho*fseta/mue/mue
        qiabr=qbr*gaunt*exp(-4.8d-11*nunu/te)

        qianu=qiasy+qiabr

        taotao=4.8d-11*nunu/te
        if(taotao.lt.1.d-4) then
        bnu=1.47d-47*nunu*nunu*nunu/(taotao+taotao*taotao/2.d0)
        else
        bnu=1.47d-47*nunu*nunu*nunu/(dexp(4.8d-11*nunu/te)-1.d0)
        endif

        kapa=qianu/4.d0/3.1416/bnu
        taonux=sqrt(3.1416)/2.*kapa*hh
c here, taonux denotes taonu star
c fnu below denotes the luminosity per unit frequency radiated
c from per unit surface of the disk

c       if(taonux.gt.1.d-5) then
c       fnu=2.*3.1416d0/sqrt(3.d0)*bnu*
c    $          (1.-dexp(-2.d0*sqrt(3.)*taonux))
c       else
c       fnu=2.*3.1416d0/sqrt(3.d0)*bnu*(2.*sqrt(3.)*taonux-
c    $          6.*taonux*taonux+4.*sqrt(3.)*taonux**3.d0)
c       endif

	fnu1=2.*3.1416d0/sqrt(3.d0)*bnu*(1.-exp(-2.*sqrt(3.)*taonux))
	taoes=0.4*rhorho*hh
	fnu=2.*fnu1*(taoes+1.)

c the factor (taoes+1.) take into account the effect: the scatering
c between photons and electrons result in the duration of stay of
c photon in medium longer by a factor of taoes+1
c---------------------------------------------------------------
c now calculate the emission from the outside SSD, including
c the effect of reprocessing from the central ADAF
c----------------------------------------------------------------

c radaf: the location (scattering point) in ADAF
c rtran: the transition radius between ADAF and SSD

c fbb: the black body emission spectrum at nunu

c        sup=5.d0*rtran/6.77d-12*m
c        slow=rtran/6.77d-12*m
c	call simp(slow,sup,10,fbb,hh,radaf,tau,nunu)
c fbb: the spectrum of the black body seed photons
c	print*,fbb
c	pause
c----------------------------------------------------------------
c the following calculate the probility basing on the formula   |
c of Coppi & Blandford (2.3)--->2.6 & 2.4                                    |
c----------------------------------------------------------------
c use eq. 2.6 to calculate R(omi,gam):

c	pro1=3.*taoes/8./gam/omip*((1-2./gam/omip-2./(gam*omip)
c     $   /(gam*omip))*log(1+2*gam*omip)+0.5+4./gam/omip
c     $          -1./2./(1+2*gam*omip)/(1+2*gam*omip))
c	if(omip*2.*gam*(1+be).lt.0.1) then
c 	pro1=taoes*(1.-2.*gam*omip/3.*(3.+be*be))
c	endif
c 	if(1./(2.*gam*(1.+be)).lt.10.*omip) then
c         xl=2.*gam*(1.-be)*omip
c         xu=2.*gam*(1.+be)*omip
c         call simp3(xl,xu,1.d-1,sum3)
c         pro1=3.*taoes/32./gam/gam/be/omip/omip*sum3
c         endif

c	if(omip.lt.1.d-3.or.omip.gt.1.d1) then
c        pro1=3.*taoes/8./gam/omip*((1-2./gam/omip-2./(gam*omip)
c     $   /(gam*omip))*log(1+2*gam*omip)+0.5+4./gam/omip
c     $          -1./2./(1+2*gam*omip)/(1+2*gam*omip))

c        else

c        call lar(30,10,1,omip,gam,a1,b1,c1,z)
c        pro1=taoes*z(1)

c        endif

	if(omip.lt.1.d-3) then
c        pro1=3.*taoes/8./gam/omip*((1-2./gam/omip-2./(gam*omip)
c     $   /(gam*omip))*log(1+2*gam*omip)+0.5+4./gam/omip
c     $          -1./2./(1+2*gam*omip)/(1+2*gam*omip))
        pro1=taoes*(1.-2.*gam*omip/3.*(3.+be*be))

        else

        call lar(30,10,1,omip,gam,a1,b1,c1,z)
        pro1=taoes*z(1)

        endif

c for the 1st sccatering, the above approximation is ok

        if(pro1.lt.0.) pro1=0.d0

c	fgam=gam*gam/2.*(5.93d+9/te)**3.d0*exp(-gam*5.93d+9/te)
c now, the 'modified' relativistic Maxwell distribution, according
c to the formula in ApJ, 2000, 541. 234. Here, 'fk2' denotes the modified
c Besol function mbsl4(2,1./setae)

	fgam=gam*gam*(5.93d+9/te)/fk2*(1.-1./gam/gam)
     $		*exp(-gam*5.93d+9/te)

	finit=fnu+fbb*0.

c	if(nunu.lt.4.d14.and.nunu.gt.3.d14) then
c	print*,'ration,fbb=',log10(nunu),fnu/fbb,fbb
c	pause
c	endif

        f=1./mc2h*p*omi/omip*finit*pro1*fgam
c above `2': the volum of the accretion flow is 2*pi*rdr*2*H
c now `2` moved out of this subrotine

        RETURN

        END
	











	subroutine simp(a,b,n,s,hh,r,tau,nunu)
        implicit real*8(a-h,o-z)
	dimension repr(200),reprr(200),reph(200),emis(200),tau2(200)
	double precision nunu
        h=(b-a)/float(2*n)
        s=0.5*(f2(a,hh,r,tau,nunu)-f2(b,hh,r,tau,nunu))
        do 1 i=1,n
1       s=s+2.*f2(a+float(2*i-1)*h,hh,r,tau,nunu)
     $          +f2(a+2.*h*float(i),hh,r,tau,nunu)
        s=(b-a)*s/float(3*n)
        return
        end






        function f2(x,hh,radaf,tau,nunu)
c here small 'r' denotes the location in ADAF
        implicit real*8(a-h,o-z)
        dimension repr(200),reprr(200),reph(200),emis(200),tau2(200)
        double precision msun,nunu,mui,mue,m
	common mui,mue,m,dotm0,rtran,beta,mc2h

        msun=1.989d33
        g=6.67d-8
        pi=3.1416d0
        dotm=dotm0/3.44d-15/m*1.39d18*m*1.d6
        rg=g*m*1.d6*msun/9.d20

	r0=radaf*rg
c x: the location in SSD, in units of cm.g.s

        fvis=3./8./pi*g*m*1.d6*msun/x**3.d0*dotm*(1.-rtran*rg/x)

c now consider the reprocessing effect
c all the related quantites such as the location in ADAF, height
c and emissivity are started with 'rep'
c 'n' is the number of the integration zone of ADAF
c reph(i), repr(i): the hieght and radius of a point in ADAF

        repn=8 

        if(kkrep.eq.0) then

          do 103 i=1,int(repn)
          read(14,*) reprr(i),reph(i),emis(i),tau2(i)
c emis is in units of cm.g.s
          repr(i)=reprr(i)*rg
          reph(i)=reph(i)*rg
103       continue

	kkrep=1
        endif

	frep=0.d0
	do 50 i=1,int(repn-1)

        ff1=reph(i)/2.*pi/(reph(i)*reph(i)/4.+(x-repr(i))**2.d0)**1.5
        ff2=reph(i)*pi/(reph(i)*reph(i)/4.+repr(i)*repr(i)+x*x)**1.5
        ff3=reph(i)/2.*pi/(reph(i)*reph(i)/4.+(repr(i)+x)**2.d0)**1.5
        frep=frep+0.8d0*emis(i)*(ff1+ff2+ff3)
     $          *(repr(i)-repr(i+1))*repr(i)*reph(i)/exp(tau2(i)*1.0)
c frep: corresponds to fvis, the flux produce per unit area due to
c the reprocessing of the central X-ray source
50      continue

c	print*,'frep',ff3,emis(i),repr(i),i,frep

c	sigmarep=5.67d-5
c sigma*T**4 is the emergent flux from an isotropically emitting 
c surface (black body), but 'bnu' below is the specific intensity(I_nu).
        tbb=((frep*0.+fvis)/5.67d-5)**0.25d0

	if(frep.gt.fvis) then
c	print*, 'frep/fvis=',frep/fvis,x/rg
	endif

        fff1=hh/4.*pi/(hh*hh/4.+(x-r0)**2.d0)**1.5d0
        fff2=hh/2.*pi/(x*x+r0*r0+hh*hh/4.)**1.5d0
        fff3=hh/4.*pi/(hh*hh/4.+(x+r0)**2.d0)**1.5d0

	taotao=4.8d-11*nunu/tbb
        if(taotao.lt.1.d-4) then
        bnu=1.47d-47*nunu*nunu*nunu/(taotao+taotao*taotao/2.d0)
        else
        bnu=1.47d-47*nunu*nunu*nunu/(dexp(4.8d-11*nunu/tbb)-1.d0)
        endif

	f2=(fff1+fff2+fff3)*x*bnu/exp(tau)
c	print*,x/rg,r0/rg,tau,tbb
c	pause

        return
        end











	subroutine simps2(a,b,eps,kflag,n,sum,rhorho,
     $		hh,setae,setai,fk2,nu,radaf,tau)
	implicit real*8(a-h,o-z)
	double precision nu

	n2=1
	kflag=0
	c=abs(a)+abs(b)
	h2=0.5*(b-a)
	call simps1(a,eps,kflag,n1,ss1,rhorho,hh,setae,setai,fk2,
     $		nu,radaf,tau)
	call simps1(b,eps,kflag,n3,ss2,rhorho,hh,setae,setai,fk2,
     $		nu,radaf,tau)
	n=n1+n3+2
	tb1=h2*(ss1+ss2)
4	x=a-h2
	tb2=0.5*tb1
	do 5 j=1,n2
	x=x+2.*h2
	call simps1(x,eps,kflag,n1,s,rhorho,hh,setae,setai,fk2,
     $		nu,radaf,tau)
	n=n+n1+1
5	tb2=tb2+h2*s
	sum=0.3333333333*(4.*tb2-tb1)
	n2=n2+n2
	if(n2.lt.33) go to 6
c  111111111111
	if(abs(sum-sum0).le.eps*(abs(sum)+1.)) return
6	sum0=sum
	tb1=tb2
	h2=0.5*h2
	if(c+h2.ne.c) go to 4
	n=-n
	return
	end







	subroutine simps1(x,eps,kflag,n1,s,
     $		rhorho,hh,setae,setai,fk2,nu,radaf,tau)
	implicit real*8(a-h,o-z)
	double precision nu
	n1=1
	call fls(x,y1,y2)
	h1=0.5*(y2-y1)
	e=abs(y1)+abs(y2)
	t1=h1*(f(x,y1,rhorho,hh,setae,setai,fk2,nu,radaf,tau)+
     $	f(x,y2,rhorho,hh,setae,setai,fk2,nu,radaf,tau))
1	y=y1-h1
	t2=0.5*t1
	do 2 i=1,n1
	y=y+2 *h1
2	t2=t2+h1*f(x,y,rhorho,hh,setae,setai,fk2,nu,radaf,tau)
	s=0.3333333333*(4.*t2-t1)
	n1=n1+n1
	if(n1.lt.33) go to 3
c 222222222222222
	if(abs(s-s0).le.eps*(abs(s)+1.)) return
3	s0=s
	t1=t2
	h1=0.5*h1
	if(e+h1.ne.e) go to 1
	if(y1.eq.y2) return
	kflag=kflag+1
	return
	end





	function mbsl4(n,x)
        double precision mbsl4,x,mbsl3
        double precision y,p,b0,b1,a(7),b(7),c(7),d(7)
        data a/-0.57721566,0.4227842,0.23069756,0.0348859,
     $          0.00262698,0.0001075,0.0000074/
        data b/1.0,0.15443144,-0.67278579,-0.18156897,
     $          -0.01919402,-0.00110404,-0.00004686/
        data c/1.25331414,-0.07832358,0.02189568,-0.01062446,
     $          0.00587872,-0.0025154,0.00053208/
        data d/1.25331414,0.23498619,-0.0365562,0.01504268,
     $          -0.00780353,0.00325614,-0.00068245/
        if(n.lt.0) n=-n
        if(x.lt.0.) x=-x
        if(x+1.0.eq.1.) then
        mbsl4=1.d+35
        return
        endif
        if(n.ne.1) then
        if(x.le.2.) then
        y=x*x/4.
        p=a(7)
        do 10 i=6,1,-1
 10      p=p*y+a(i)
        p=p-mbsl3(0,x)*log(x/2.)
        else
        y=2.0/x
        p=c(7)
        do 20 i=6,1,-1
 20     p=p*y+c(i)
        p=p*exp(-x)/sqrt(x)
        endif
        endif
        if(n.eq.0) then
        mbsl4=p
        return
        endif
        b0=p
        if(x.le.2.0) then
        y=x*x/4.
        p=b(7)
        do 30 i=6,1,-1
 30     p=p*y+b(i)
        p=p/x+mbsl3(1,x)*log(x/2.)
        else
        y=2./x
        p=d(7)
        do 40 i=6,1,-1
 40     p=p*y+d(i)
	p=p*exp(-x)/sqrt(x)
        endif
        if(n.eq.1) then
        mbsl4=p
        return
        endif
        b1=p
        y=2./x
        do 50 i=1,n-1
        p=b0+i*y*b1
        b0=b1
        b1=p
 50     continue
        mbsl4=p
        return
        end


C THIS PROGRAM IS TO CALCULATE THE MODIFIED BASEL FUNCTION OF THE FIRST KIND
C IT IS CORRECT!!!!!!!!!!!!!!!
        function mbsl3(n,x)
        double precision mbsl3,x
        double precision t,y,p,b0,b1,q,a(7),b(7),c(9),d(9)
        data a/1.0, 3.5156229,3.0899424,1.2067492,
     $          0.2659732,0.0360768,0.0045813/
        data b/0.5,0.87890594,0.51498869,0.15084934,0.02658773,
     $          0.00301532,0.00032411/
        data c/0.39894228,0.01328592,0.00225319,-0.00157565,
     $          0.00916281,-0.02057706,0.02635537,
     $          -0.01647663,0.00392377/
        data d/0.39894228,-0.03988024,-0.00362018,0.00163801,
     $       -0.01031555,0.02282967,-0.02895312,0.01787654,-0.00420059/

        if(n.lt.0) n=-n
        t=abs(x)
        if (n.ne.1) then
         if(t.lt.3.75) then
           y=(x/3.75)*(x/3.75)
           p=a(7)
           do 10 i=6,1,-1
  10       p=p*y+a(i)
         else
           y=3.75/t
           p=c(9)
           do 20 i=8,1,-1
 20        p=p*y+c(i)
           p=p*exp(t)/sqrt(t)
         endif
        endif
        if(n.eq.0) then
         mbsl3=p
	return
        endif
        q=p
        if(t.lt.3.75) then
        y=(x/3.75)**2.d0
        p=b(7)
        do 30 i=6,1,-1
 30     p=p*y+b(i)
        p=p*t
        else
        y=3.75/t
        p=d(9)
        do 40 i=8,1,-1
 40     p=p*y+d(i)
        p=p*exp(t)/sqrt(t)
        endif
        if(x.lt.0.0) p=-p
        if(n.eq.1) then
        mbsl3=p
        return
        endif
        if(x+1.0.eq.1.0) then
        mbsl3=0.0
        return
        endif
        y=2.0/t
        t=0.0
        b1=1.0
        b0=0.0
        m=sqrt(40.0*n)
        m=m+n
        m=m+m
        do 50 i=m,1,-1
        p=b0+i*y*b1
        b0=b1
        b1=p
        if(abs(b1).gt.1.d+10) then
        t=t*1.d-10
        b0=b0*1.d-10
        b1=b1*1.d-10
        endif
        if(i.eq.n) t=b0
 50     continue
	p=t*q/b1
        if((x.lt.0.0).and.(mod(n,2).eq.1)) p=-p
        mbsl3=p
        return
        end













! This routine is called several times furing the Comptonization.
! It calls FCT several times.
! • the only variables that change between calls in the IC-calculation
!   are SUM and J
! 
	SUBROUTINE SIMP2(A,B,EPS,SUM,rhorho,hh,setae,j)
        implicit real*8(a-h,o-z)
c	double precision nu,ssum
        DIMENSION F(2,30),FM(2,30),E(2,30),KRTN(30),nu(200),ssum(200)
        double precision nu,ssum

        SUM=0.d0
        T=1.d0
        ABSA=1.0
        EST=1.d0
        DA=B-A
c       FA=FCT(A,rhorho,hh,setae,j,nu,ssum)
	FA=FCT(A,rhorho,hh,setae,j)
        FB=FCT(B,rhorho,hh,setae,j)
        FP=4.0*FCT(0.5*(A+B),rhorho,hh,setae,j)
        X=A
        L=0
 1      K=1
        L=L+1
        T=T*1.7
        DX=DA/3.0
        SX=DX/6.0
        FM1=4.0*FCT(X+0.5*DX,rhorho,hh,setae,j)
        F1=FCT(X+DX,rhorho,hh,setae,j)
        FM(1,L)=FP
	F(1,L)=FCT(X+2.0*DX,rhorho,hh,setae,j)
        FM(2,L)=4.0*FCT(X+2.5*DX,rhorho,hh,setae,j)
        F(2,L)=FB
        E1=SX*(FA+FM1+F1)
        E(1,L)=SX*(F1+FP+F(1,L))
        E(2,L)=SX*(F(1,L)+FM(2,L)+FB)
        S=E1+E(1,L)+E(2,L)
        ABSA=ABSA-ABS(EST)+ABS(E1)+ABS(E(1,L))+ABS(E(2,L))
        IF(EST.EQ.1.0) GOTO 5
        IF(T*ABS(EST-S).LE.EPS*ABSA) GOTO 2
        IF(L.LT.30) GOTO 5
2       SUM=SUM+S
3       L=L-1
        T=T/1.7
        K=KRTN(L)
        DX=DX*3.0
        IF(K.EQ.3) IF(L-1) 6,6,3
        EST=E(K,L)
        FP=FM(K,L)
        FA=FB
        FB=F(K,L)
        K=K+1
	X=X+DA
4       DA=DX
        KRTN(L)=K
        GOTO 1
5       EST=E1
        FP=FM1
        FB=F1
        GOTO 4
6       RETURN
        END






c	function fct(gam,rhorho,hh,setae,jj,nu,ssum)
 	function fct(gam,rhorho,hh,setae,jj)
        implicit real*8(a-h,o-z)
	dimension omip(200),nu(200),ssum(200)
	dimension a1(30),b1(10),c1(30,10,1),z(1)
	dimension a2(100),b2(80),c2(100,80,1),z2(1)
	dimension a3(100),b3(80),c3(100,80,1),z3(1)
c	common mui,mue,m,dotm0,rtran,beta,mc2h
	double precision nu,mc2h,mbsl4,ssum,s
	common /comp1/nu,ssum
	common /comp2/a1,b1,c1,a2,b2,c2,a3,b3,c3
	

c	print*,ssum(40),nu(140)

	c=3.d+10
	mc2h=1.236d+20

	te=setae/1.686d-10

	sig=0.d0

	do 91 i=1,140
	omip(i)=nu(i)/mc2h
	if(nu(i).lt.1.d0) goto 81
	omi=nu(jj)/mc2h
	be=sqrt(1.-1./gam/gam)

c 	if(gam*omip(i).lt.0.01d0) then
c	aomi=4./3.*gam*gam*omip(i)
c 	else
c now I use 2.12 in Coppi & Blandford to calculate aomi
c 	call simpsn(-1.d0,1.d0,1.d-4,6,n,s,omip(i),gam)
c 	aomi=s
c 	endif
c now use 2.12 to calculate aomi

	if(omip(i).lt.5.d-4) then
         if(gam.gt.1.1) then
            aomi=4./3.*gam*gam*omip(i)
          else
             aomi=(1.+4./3.*be*be-omip(i))*omip(i)
          endif
        else
          call lar2(100,80,1,omip(i),gam,a2,b2,c2,z2)
          aomi=z2(1)

        endif

c	aomi=4./3.*gam*gam*omip(i)
c the value of aomi is very important to determine the location 
c of the high-energy cut-off, should use the exact integration.

c 2006/07/30: now use eqs. 2.17 & 2.18 in Coppi & Blandford for aomi2:
c ------------------------------------------------------------------

c now calculate aomi2 using eq. 2.15

        if(omip(i).lt.2.d-6) then
         aomi2=14./5.*gam**4.d0*omip(i)*omip(i)*(1.-176/35.*gam*omip(i))
	aomi2=aomi2/4.5
        else
         call lar2(100,80,1,omip(i),gam,a3,b3,c3,z3)
         aomi2=z3(1)
        endif
c there is no good approximation to aomi2!!
c	aomi2=14./5.*gam**4.d0*omip(i)*omip(i)*(1.-176./35.*gam*omip(i))

c	if(gam*omip(i).gt.10.) then
c	romip=3./8./gam/omip(i)*log(4.*gam*omip(i))
c note the above formula is from eq. 2.5 in Coppi & Blandford, but
c c & sigma_T are omitted since there seems to be print error in eq. 2.18
c        q1=(6.*gam*omip(i)*(1.+be)*(2.*gam*gam+1)+6*gam*gam+3)
c     $          *log(2.*gam*omip(i)*(1.+be)+1.)
c        q2=(6.*gam*omip(i)*(1.-be)*(2.*gam*gam+1)+6*gam*gam+3)
c     $          *log(2.*gam*omip(i)*(1.-be)+1.)
c        q3=9./32.*(1-be*be)**2.d0*gam**3.*omip(i)-(58.*gam*gam+1)/
c     $          64./gam/omip(i)+7./32.*(2.*gam*gam*(1-be*be)+1)
c        aomi2=1./romip/64./be/gam/gam/omip(i)/omip(i)*(q1-q2)
c     $         +q3/romip
c	else
c	aomi2=14./5.*gam**4.d0*omip(i)*omip(i)*(1.-176./35.*gam*omip(i))
c	endif

c--------------------------------------------------------------------
        deomi=aomi2-aomi*aomi
        d1=sqrt(3.*deomi)
        d2=aomi
        if(d1.lt.d2) then
        dd=d1
        else
        dd=d2
        endif

        if(dd.gt.abs(omi-aomi)) then
        hea=1.d0
        else
        hea=0.d0
        endif

        p=1./2./dd*hea

	delomp=(nu(i+1)-nu(i-1))/mc2h/2.d0

c the following calculate R(omip,gam) according to eq. 2.3 
c	pro1=taoes*(1.-2.*gam*omip(i)/3.*(3.+be*be))
c        if(omip(i).gt.1./(2.*gam*(1.+be))/2000.) then
c changed above: from 10---> 500: 2006/07/30
c        xl=2.*gam*(1.-be)*omip(i)
c       xu=2.*gam*(1.+be)*omip(i)
c        call simp3(xl,xu,1.d-1,sum3)
c        pro1=3.*taoes/32./gam/gam/be/omip(i)/omip(i)*sum3
c        endif

	taoes=0.4*rhorho*hh
	if(omip(i).lt.1.d-3) then

c        pro1=3.*taoes/8./gam/omip(i)*((1-2./gam/omip(i)-2.
c     $		/(gam*omip(i))
c     $   /(gam*omip(i)))*log(1+2*gam*omip(i))+0.5+4./gam/omip(i)
c     $          -1./2./(1+2*gam*omip(i))/(1+2*gam*omip(i)))
	 pro1=taoes*(1.-2.*gam*omip(i)/3.*(3.+be*be))
        else

        call lar(30,10,1,omip(i),gam,a1,b1,c1,z)
        pro1=taoes*z(1)

        endif

        if(pro1.lt.0.) pro1=0.d0

	sig=sig+delomp*p*omi/omip(i)*ssum(i)*pro1

91	continue
81	taoes=0.4*rhorho*hh

        fgam=gam*gam*(5.93d+9/te)/mbsl4(2,5.93d+9/te)*(1.-1./gam/gam)
     $          *exp(-gam*5.93d+9/te)
c	fct=fgam*sig*pro1
	fct=fgam*sig
c	print*,'fct=',fct,jj,pro1
	RETURN

        END









	SUBROUTINE SIMP3(A,B,EPS,SUM)
        implicit real*8(a-h,o-z)
        DIMENSION F(2,30),FM(2,30),E(2,30),KRTN(30)

        SUM=0.d0
        T=1.d0
        ABSA=1.0
        EST=1.d0
        DA=B-A
        FA=FCT3(A)
        FB=FCT3(B)
        FP=4.0*FCT3(0.5*(A+B))
        X=A
        L=0
 1      K=1
        L=L+1
        T=T*1.7
        DX=DA/3.0
        SX=DX/6.0
        FM1=4.0*FCT3(X+0.5*DX)
        F1=FCT3(X+DX)
        FM(1,L)=FP
        F(1,L)=FCT3(X+2.0*DX)
        FM(2,L)=4.0*FCT3(X+2.5*DX)
        F(2,L)=FB
        E1=SX*(FA+FM1+F1)
        E(1,L)=SX*(F1+FP+F(1,L))
        E(2,L)=SX*(F(1,L)+FM(2,L)+FB)
        S=E1+E(1,L)+E(2,L)
        ABSA=ABSA-ABS(EST)+ABS(E1)+ABS(E(1,L))+ABS(E(2,L))
        IF(EST.EQ.1.0) GOTO 5
        IF(T*ABS(EST-S).LE.EPS*ABSA) GOTO 2
        IF(L.LT.30) GOTO 5
2       SUM=SUM+S
3       L=L-1
        T=T/1.7
        K=KRTN(L)
        DX=DX*3.0
        IF(K.EQ.3) IF(L-1) 6,6,3
        EST=E(K,L)
        FP=FM(K,L)
        FA=FB
        FB=F(K,L)
        K=K+1
        X=X+DA
4       DA=DX
	KRTN(L)=K
        GOTO 1
5       EST=E1
        FP=FM1
        FB=F1
        GOTO 4
6       RETURN
        END









	function fct3(x)
        implicit real*8(a-h,o-z)
	fct3=(1.-4./x-8./x/x)*log(1.+x)+0.5d0+8./x-1./2./(1.+x)/(1.+x)
	return
	end

	subroutine simpsn(a,b,eps,k,n,s,omip1,gam1)
        implicit real*8(a-h,o-z)
        n=1
        c=abs(a)+abs(b)
        h=0.5*(b-a)
        t1=h*(f3(a,omip1,gam1)+f3(b,omip1,gam1))
1       x=a-h
        t2=0.5*t1
        do 2 i=1,n
        x=x+2.*h
2       t2=t2+h*f3(x,omip1,gam1)
        s=0.33333333333d0*(4.*t2-t1)
        n=n+n
        if(abs(s-s0).le.(1.+abs(s))*eps) return
3       s0=s
        t1=t2
        h=h*0.5
        if(c+h.ne.c) goto 1
        n=-n
        return
        end




	function f3(xx,omip1,gam1)
        implicit real*8(a-h,o-z)
c	romip=3./8./gam1/omip1*log(4.*gam1*omip1)
	romip=3./8./gam1/omip1*((1-2./gam1/omip1-2./(gam1*omip1)
     $		/(gam1*omip1))*log(1+2*gam1*omip1)+0.5+4./gam1/omip1
     $		-1./2./(1+2*gam1*omip1)/(1+2*gam1*omip1))
	beta=1.-1./gam1/gam1
	x=gam1*omip1*(1-beta*xx)
	coef=3.*x/8./romip/gam1/omip1
	a=gam1*gam1*beta*omip1*(xx-beta)
	term1=2./(1.+2.*x)*(a+gam1*(x-2.)-(2.*gam1+a)/x-a*(6/x/x
     $		+3/x/x/x))
	term2=log(1.+2.*x)/x/x*(gam1-a+3.*a*(1./x+1./x/x))
	term3=1./2./x*(1.-1./(1+2*x)/(1+2*x))*(a*(3/x/x+1./x/x/x)
     $	   +(a+gam1)/x+2*gam1)+1./3*(1-1./(1+2*x)/(1+2*x)/(1+2*x))
     $     *(gam1+a*(1./x+1./x/x))-2*a/x/x/x
		
	f3=coef*(term1+term2+term3)/2.
	return
	end





	SUBROUTINE LAR(N,M,L,X,Y,A,B,C,Z)
	implicit real*8(a-h,o-z)
        DIMENSION A(N),B(M),C(N,M,L),Z(L),V(3),U(3)
        N1=N-2
        M1=M-2
        DO 1 I=1,N1
        IF(X.LE.A(I+1)) GO TO 2
1       CONTINUE
        I=N-2
2       DO 3 J=1,M1
        IF(Y.LE.B(J+1)) GOTO 4
3       CONTINUE
        J=M-2
4       IF(I.EQ.1) GOTO 5
        IF(X-A(I).GE.A(I+1)-X) GOTO 5
        I=I-1
5       IF(J.EQ.1) GOTO 6
        IF(Y-B(J).GE.B(J+1)-Y) GOTO 6
        J=J-1
6       A1=A(I)
        A2=A(I+1)
        A3=A(I+2)
        B1=B(J)
        B2=B(J+1)
        B3=B(J+2)
        U(1)=(X-A2)*(X-A3)/((A1-A2)*(A1-A3))
        U(2)=(X-A1)*(X-A3)/((A2-A1)*(A2-A3))
        U(3)=(X-A1)*(X-A2)/((A3-A1)*(A3-A2))
        V(1)=(Y-B2)*(Y-B3)/((B1-B2)*(B1-B3))
        V(2)=(Y-B1)*(Y-B3)/((B2-B1)*(B2-B3))
        V(3)=(Y-B1)*(Y-B2)/((B3-B1)*(B3-B2))
        DO 8 K=1,L
        W=0.0
        DO 7 II=1,3
        DO 7 JJ=1,3
        I1=I+II-1
        J1=J+JJ-1
7       W=W+U(II)*V(JJ)*C(I1,J1,K)
8       Z(K)=W
        RETURN
        END





	SUBROUTINE LAR2(N,M,L,X,Y,A,B,C,Z)
        implicit real*8(a-h,o-z)
        DIMENSION A(N),B(M),C(N,M,L),Z(L),V(3),U(3)
        N1=N-2
        M1=M-2
        DO 1 I=1,N1
        IF(X.LE.A(I+1)) GO TO 2
1       CONTINUE
        I=N-2
2       DO 3 J=1,M1
        IF(Y.LE.B(J+1)) GOTO 4
3       CONTINUE
        J=M-2
4       IF(I.EQ.1) GOTO 5
        IF(X-A(I).GE.A(I+1)-X) GOTO 5
        I=I-1
5       IF(J.EQ.1) GOTO 6
        IF(Y-B(J).GE.B(J+1)-Y) GOTO 6
        J=J-1
6       A1=A(I)
        A2=A(I+1)
        A3=A(I+2)
        B1=B(J)
        B2=B(J+1)
        B3=B(J+2)
        U(1)=(X-A2)*(X-A3)/((A1-A2)*(A1-A3))
        U(2)=(X-A1)*(X-A3)/((A2-A1)*(A2-A3))
        U(3)=(X-A1)*(X-A2)/((A3-A1)*(A3-A2))
        V(1)=(Y-B2)*(Y-B3)/((B1-B2)*(B1-B3))
        V(2)=(Y-B1)*(Y-B3)/((B2-B1)*(B2-B3))
        V(3)=(Y-B1)*(Y-B2)/((B3-B1)*(B3-B2))
        DO 8 K=1,L
        W=0.0
        DO 7 II=1,3
        DO 7 JJ=1,3
        I1=I+II-1
        J1=J+JJ-1
7       W=W+U(II)*V(JJ)*C(I1,J1,K)
8       Z(K)=W
        RETURN
        END

