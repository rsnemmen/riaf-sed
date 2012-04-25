c this code calculate the spectrum of a standard thin disk
c the code is based on the formula in the book "accretion
c power in astrophysics" by Frank et al. 1992
c in unit of cm.g.sec

     	IMPLICIT REAL*8(A-H,O-Z)
        dimension y(4),dy(4),w(16),nu(200),fnu(200)
	common dotm,h,k,mass,rs,solarmass,rin
        double precision mass,nu,k,h,nnu
	
	call cpu_time(start_time)

	open(20,file='ssd.dat')
	open(14,file='hottem.dat')

	h=6.626d-27
	c=2.9979d+10
	k=1.3807d-16
	pi=3.1416d0
	pc=3.086d+18
	solarmass=1.989d+33

        read(*,*) distance
        distance=distance*pc
c	distance=14.5d6*pc !--------------1-----------------!

        read(*,*) mass
        mass=mass*1.d6*solarmass
c	mass=1.2d8*solarmass !--------------2-----------------!

	g=6.672d-8
	rs=2.*g*mass/c/c
	dotme=1.39d18*mass/solarmass
c note in the above defination, the factor of 0.1 is included
c the unit is g/s

c***********************************
c *** I still did not figure this out. By comparing the output of
c this code with the SSD models in Lasota+96, Quataert+99, Nemmen+06,
c for the same mdot the SSD SEDs there are 10x weaker compared to the
c output of this code. By hand I introduced below a factor 0.1 for
c the correction of mdot. (9/01/2008)
c **********************************
c UPDATE Mar 8th 2012
c Still wondering about the 0.1 factor below
c **********************************
        read(*,*) dotm
        dotm=0.1d0*dotm*dotme
c        dotm=dotm*dotme
c	dotm=6.4d-4*dotme !--------------3-----------------!
	
        read(*,*) rin
        rin=rin*rs
c	rin=225.d0*rs  !--------------4-----------------!
c	rout=3.0d+4*rs  !--------------5-----------------!

        read(*,*) theta
	theta=theta/57.2958 !--------------6-----------------!
c*********************************
	do 10 i=1,50

c *** Changed the initial frequency from 12. to 11.
	nu(i)=10.d0**(11.+0.1d0*i)

	if(nu(i).lt.1.d16) rout=8.d+4*rs
	if(nu(i).gt.1.d16) rout=1.d+4*rs

	nnu=nu(i)
	s=0.d0
	call simpsn(rin,rout,1.d-3,8000,n,s,nnu)
	write(*,*) i,nu(i),n,s,exp(700.d0)

	fnu(i)=4.*pi*h*cos(theta)*nu(i)**3.d0/c/c*s 
c *** Changed the output below to nu*Lnu instead of nu*Fnu
c	write(20,*) log10(nu(i)),log10(nu(i)*fnu(i)/distance**2.)
	write(20,*) log10(nu(i)),log10(nu(i)*fnu(i)*4.*pi)

10	continue

	call cpu_time(finish_time)
	print*,'cpu times=',finish_time-start_time,' seconds'

	end
c=====================================================================c
	function f(x,nu)
	implicit real*8(a-h,o-z)
	dimension reprr(1000),repr(1000),rephh(1000),reph(1000),
     $		emis(1000),tau2(1000)
	common dotm,h,k,mass,rs,solarmass,rin
	double precision mass,nu,k,m8,msun

c	dotm26=dotm/1.d+26*(1.-rin/x)
 	dotm26=dotm/1.d+26

c now consider the reprocessing effect
c all the related quantites such as the location in ADAF, height
c and emissivity are started with 'rep'
c 'n' is the number of the integration zone of ADAF
c reph(i), repr(i): the hieght and radius of a point in ADAF

        repn=11

	msun=1.989d33
        g=6.67d-8
        pi=3.1416d0
c       rs=g*m*1.d6*msun/9.d20

        if(kkrep.eq.0) then

        do 103 i=1,repn
        read(14,*) reprr(i),rephh(i),emis(i),tau2(i)
c emis is in units of cm.g.s
        repr(i)=reprr(i)*rs/2.
        reph(i)=rephh(i)*rs/2.
103     continue

        kkrep=1
        endif

        frep=0.d0
        do 50 i=1,repn-1

c       ff1=reph(i)/2.*pi/(reph(i)*reph(i)/4.+(x-repr(i))**2.d0)**1.5

	ff1=reph(i)/2.*pi/(reph(i)*reph(i)/4.+
     $          (repr(i)*sin(22.5/57.))**2.d0+(x-repr(i))**2.d0)**1.5

        ff2=reph(i)*pi/(reph(i)*reph(i)/4.+repr(i)*repr(i)+x*x)**1.5
        ff3=reph(i)/2.*pi/(reph(i)*reph(i)/4.+(repr(i)+x)**2.d0)**1.5
        frep=frep+0.8d0*emis(i)*(ff1+ff2+ff3)
     $          *(repr(i)-repr(i+1))*repr(i)*reph(i)/exp(tau2(i)*1.0)
c frep: corresponds to fvis, the flux produce per unit area due to
c the reprocessing of the central X-ray source

50      continue

	fvis=3./8./pi*g*mass/x**3.d0*dotm26*1.d26*(1.-(3.*rs/x)**0.5d0)
c	fvis=3./8./pi*g*mass/x**3.d0*dotm26*1.d26

	m8=mass/solarmass/1.d+8
	r14=x/1.d+14
c	tem=2.2d+5*dotm26**0.25d0*m8**0.25d0*r14**(-0.75d0)

	tem=((frep*10.d-1+fvis)/5.67d-5)**0.25d0

c for the above formula, see the first page of Chapter 8 in <<accretion>> book

	if (h*nu/(k*tem).gt.700.) then
 	 f=0.d0
	else
         f=x/(dexp(h*nu/(k*tem))-1.d0)
	endif

	return
	end

c=============================================================c
	subroutine simpsn(a,b,eps,k,n,s,nu)
        implicit real*8(a-h,o-z)
	double precision nu
        n=1
        c=abs(a)+abs(b)
        h=0.5*(b-a)    
        t1=h*(f(a,nu)+f(b,nu))
1       x=a-h
        t2=0.5*t1
        do 2 i=1,n
        x=x+2.*h
2       t2=t2+h*f(x,nu)
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
