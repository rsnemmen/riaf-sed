
c-------------------------------------------------------
C THIS PROGRAM THE GLOBAL SOLUTIONS FOR TWO-TEMPERATURE
C ADAF USING SHOOTING METHOD. THE EQUATIONS ARE FROM MANMOTO ET AL.
C 1997 APJ, WITH THE FOLLOWING MODIFICATIONS:
C
C 1)THE MASS ACCRETION RATE IS THE FUNCTION OF RADIUS
C 2)IN THE ENERGY EQUATION OF ELECTRONS, THE VISOCUS DISSIPATION
C   IS ALSO INCLUDED, THE PARAMETER IS ``DELTA''


	IMPLICIT REAL*8(A-H,O-Z)
	dimension y(4),dy(4),w(16)
	dimension r(1000),ddhh(1000),sum(9000),sigma(1000)
	double precision m,mui,mue,mbsl4
	common alfa,beta,dotm0,rout,p0,gamai,ru,
     $          sl0,m,mui,mue,delta,rtran
	open(41,file='x.dat')
	open(42,file='hottem.dat')
	open(14,file='hot.dat')

	ru=9.2d-14
	mmm=0

c --------------------------------
c here, ru is k/M_mu
c---------------------------------
c	print *,'Adiabatic index gamma:'
	read(*,*) gamai
c        gamai=5./3.d0

c	print *,'Black hole mass (in Solar masses):'
	read(*,*) m
c	m=10.d0

	mui=1.24d0
	mue=1.13d0

c ratio of gas to total pressure
c	print *,'beta:'
	read(*,*) beta
c	beta=0.9d0

c	print *,'alpha:'
	read(*,*) alfa
c	alfa=3.d-1

c Fraction of turbulent dissipation that directly heats electrons
c	print *,'delta:'
	read(*,*) delta 
c	delta=5.d-1

	dotme=3.44d-15*m
c  dotme is the Eddington accretion rate

c	print *,'Mdot_out (Eddington units):'
	read(*,*) dotm0
	dotm0=dotme*dotm0
c dotm0 is the accretion rate at the outer boundary

c	print *,'R_out (units of R_S):'
	read(*,*) rout
c	rout=600.d0
c rout is the outer boundary

c	print *,'p_wind ("strength of wind"):'
	read(*,*) pp0
	p0=0.0d0
c	pp0=0.55d0
c  pp0 describes the strength of the outflow
	sint=0.d0
	
c	print *,' '
c	print *,'Boundary conditions:'
c	print *,'====================='
c	print *,'T_i (ion temperature):'
	read(*,*) ti
c	ti=5.d+9

c	print *,'T_e (electron temperature):'
	read(*,*) te
c	te=1.8d+9

	cs=sqrt(ru/beta*(ti/mui+te/mue))

c	print *,'v_R/c_s (radial velocity/sound speed):'
	read(*,*) vcs
c	vcs=3.d-1

c you need to adjust the boundary conditions ti, te, and vcs (the ratio
c between the radial velocity and local sound speed) to get a solution
c	print*,'sl0=? (eigenvalue of the problem)'
	read(*,*) sl0
c sl0 is the eignevalue of the problem

c Test if code is receiving its input correctly ***
c	print *,gamai,m,beta
c	print *,alfa,delta,dotm0
c	print *,rout,pp0,ti
c	print *,te,vcs,sl0

	I=2
c Added! ***
	y(1)=rout
c	y(1)=600.d0
	rtran=y(1)*1.001d0
	y(2)=-cs*vcs
	y(3)=te
	y(4)=ti

	sk=1./sqrt(y(1))/(y(1)-2.)
	sout=-alfa*y(1)*cs*cs/y(2)+sl0
	omiga=sout/y(1)/y(1)
	omik=sk

	print*,'# omiga/omigak=',omiga/omik,sout,100/y(1)/y(1)/sk
	if(omiga/omik.gt.1.) then
	   print*,'# omiga/omigak > 1!!!', omiga/omik
       	   print*,'# I must let vcs increase!'
	   stop
	endif

          k2=0
	inform=i+10
c------------------------
        tau=0.d0
	y0=rtran
c-------------------------

 222      if(k2.lt.160) hh=-1.d+1
	  if(k2.ge.160) hh=-1.d0
          if(k2.ge.500) hh=-1.d-1

	hh=-y(1)/100.d0
          
        call mers(4,inform,hh,1.d-4,y,dy,w,qn,qd,tau,sint)
c ***
	print *,'Breakpoint 1'

	cs=sqrt(ru/beta*(y(4)/mui+y(3)/mue))
	cs2=cs*cs
	omik=1./sqrt(y(1))/(y(1)-2.)
	slk=omik*y(1)*y(1)
	dlnok=-.5/y(1)-1./(y(1)-2.)
	ssll=-alfa*cs2*y(1)/y(2)+sl0
	omiga=ssll/y(1)/y(1)
c	rho=-dotm0*(y(1)/rout)**p0/4./3.1416/y(1)/cs*omik/y(2)
 	rho=-dotm0*exp(sint)/4./3.1416/y(1)/cs*omik/y(2)
	tau=tau+rho*(y0-y(1))*3.65d16/m
        y0=y(1)

	    ww=-alfa/y(2)*cs2
            w1=-ww/y(2)
            w2=ww/(y(3)/mue+y(4)/mui)
            w3=ww/y(1)-2.*omiga
            rdomiga=w1*dy(2)+w2*(dy(3)/mue+dy(4)/mui)+w3
            if(rdomiga.gt.0.) then
            print*,rdomiga,'# $$$$$$$$$$$$$$$'
            print*,'# FAILED!!!'
c           stop
            endif
	dh=cs/omik
        tao=dh*rho*3.65d16/m
        tao=tao/0.4d0
	dotmr=-4.*3.14*dh*rho*y(1)*y(2)/dotme

c now tao is in units of cgs
        poten=ssll*ssll/2./y(1)/y(1)-1./(y(1)-2.)
        Be=0.5*y(2)*y(2)+1.5/0.5*ru*y(3)+2.5*ru*y(4)+
     $          ssll*ssll/2./y(1)/y(1)-1./(y(1)-2.)

	
c the following calculates the advection factor

        a1=-1./y(2)
        a2=-1./2./(y(3)/mue+y(4)/mui)
c       a3=dlnok-1./y(1)
        a3=dlnok+(p0-1.)/y(1)
        rhodr=a1*dy(2)+a2*(1./mui*dy(4)+1./mue*dy(3))+a3
        setae=1.686d-10*y(3)
        setai=9.181d-14*y(4)
        qadv=rho*y(2)*ru*(dy(4)/mui/(gamai-1.)-y(4)/mui*rhodr)

	ate=1./setae*((3.*mbsl4(3,1.d0/setae)+mbsl4(1,1.d0/setae))/
     $       4./mbsl4(2,1.d0/setae)-1.d0)
	qadve=rho*y(2)*ru*(ate/mue*dy(3)-mue*y(3)*rhodr)
c qadve denotes the advection of electrons
	
        qie2=4./4.*37.85d0/m*rho*rho*(y(4)-y(3))*(sqrt(2./3.1416)+sqrt(
     $          setai+setae))/(setai+setae)**1.5d0

        qadv1=rho*y(2)*ru*dy(4)/mui/(gamai-1.)
        qadv2=rho*y(2)*ru*y(4)/mui*rhodr

        qvis=qadv+qie2
        advec=qadv/qvis
	advec2=(qadv+qadve)/qvis

        compy=4.*tao*1.686d-10*y(3)
	pgas=rho*6.d+5/m/m*1.38d-16/1.67d-24*(y(3)+y(4))
	bb=sqrt((1.-beta)/beta*pgas*8.*3.1416)

c--------------------------------------------------------------------------
c the following calculate the case of s is NOT equal to 0., i.e., outflow's case

	if(advec.lt.0.) p0=0.d0
 	if(advec.gt.0.) p0=advec*pp0

c the following calculate the integration for 's'
c sint denotes the integration for 's'

	if(k2.eq.0) p0=0.d0
	sint=sint+p0/y(1)*hh

c----------------------------------------------------------------------
	

	if(k2.eq.0) then
	seta=2.3536744d0
        ate=2.73821d0+0.12773*(setae-seta)-0.0901154*(setae-seta)**2.d0
     $      +0.0227991*(setae-seta)**3.d0-0.00219443*(setae-seta)**4.d0
     $      +7.07116d-5*(setae-seta)**5.d0

c	ate=1./setae*((3.*mbsl4(3,1.d0/setae)+mbsl4(1,1.d0/setae))/
c     $       4./mbsl4(2,1.d0/setae)-1.d0)

	qadve=rho*y(2)*ru*(ate/mue*dy(3)-mue*y(3)*rhodr)
	endif

          k2=k2+1

          if(abs(k2/10-k2/10.).lt.1.d-5.or.k2.eq.1) then

	  aaa=sqrt((4.+1.33*alfa*alfa)/(8./3.))
	  outflow=-rho*4.*3.14*y(1)*cs/omik*y(2)/dotme/0.05

c Modified the standard output ***
	    write(*,*) y(1),-y(2)/cs/aaa,log10(y(3))
     $	        ,log10(y(4)),advec,advec2,cs,dh,qn,tau,log10(slk)
     $		,log10(ssll),bb
c    $		log10(m*rho/1.01d-63*6.77d-12**3.d0/m**3.d0)
c the number density above is in cgs units 
 	    write(41,6) log10(y(1)/2.d0),-y(2),
     $		log10(y(3)),log10(y(4)),log10(slk),log10(ssll)
     $          ,log10(tao),bb,-8.d+8/(y(1)/y(2)/0.203*m)/bb/bb
c    $		,log10(tao),bb,log10(rho*6.1749d+5/m/m/1.67d-24/mue)
     $		,compy,advec,qn*y(1)*y(1)*dh*1.d+18
c    $		,log10(rho*6.1749d+5/m/m/1.67d-24/mue),advec,be
c    $          ,dotmr,advec,be

	write(42,*) y(1),dh,qn,tau

c the last quantity -8.d+8/(y(1)/y(2)/0.203*m)/bb/bb, is the coolong Lolentz
c factor of the electrons.
7	format(3e12.4)
6	format(12e12.4)
5	format(6e11.3)


         endif

	if(y(1).lt.1.6.and.mmm.eq.0) then
c	print*,y,rho

           y(4)=3./16./ru*y(2)*y(2)
           y(3)=y(4)
           y(2)=y(2)/4.d0
           inform=12
           mmm=1
           print*,y
           pause
        endif


c----------------------------------------------------------
c the follwoing code calculate the luminosity emitted by
c the hot disk
c----------------------------------------------------------

         r(k2)=y(1)
         ddhh(k2)=cs/omik
         rho=-dotm/4./3.1416/y(1)/cs*omik/y(2)
         sum(k2)=qn*8.87d-27*m*m*m

         if(y(1).ge.2.2d0) goto 222

         sumlum=0.d0

         do 11 i=2,k2-1

         sigma(i)=3.1416*((r(i-1)+r(i))**2.d0/4.d0-(r(i)+r(i+1))
     $          **2.d0/4.d0)
         sumlum=sumlum+sigma(i)*2.*ddhh(i)*sum(i)/2.d0*2.d0
c in the above formula, a factor '/2.d0' should be present. This is because
c  what we observe is only half of the total emitted radiation: no!!!!

         Eddlum=1.3d+38*1.d+6*m*5.59d-61/m/2.03d-1*m
 11      continue
         write(*,*) '# The total luminosity is:',sumlum/eddlum,
     $		sumlum/(5.59d-61/m/2.03d-1*m)
c *** added
     	 write(*,*) '# Run finished'


	  if(y(1).ge.2.2d0) goto 222
	
	end

	subroutine dery(n,y,dy,qn,qd,tau,sint)
	implicit real*8(a-h,o-z)
	double precision m,mui,mue,numax,mbsl4
	dimension y(4),dy(4)
	common alfa,beta,dotm0,rout,p0,gamai,ru,
     $          sl0,m,mui,mue,delta,rtran
c ------------------------------------
c y(1): r;  y(2): v; y(3):te, y(4):ti    |
c ------------------------------------	

	if(y(3).lt.0..or.y(4).lt.0.)  then
c	print*,y(3),y(4)
c	stop
	goto 232
	endif

	omik=1./sqrt(y(1))/(y(1)-2.)
	dlnok=-1./2.d0/y(1)-1./(y(1)-2.)
	cs2=ru/beta*(y(3)/mue+y(4)/mui)
	ssll=-alfa*cs2*y(1)/y(2)+sl0

	if(ssll.lt.0.) then
	print*,'# ssll < 0 '
	stop
	endif

	omiga=ssll/y(1)/y(1)

c	rho=-dotm/4./3.1416/y(1)/sqrt(cs2)*omik/y(2)
c	rho=-dotm0*(y(1)/rout)**p0/4./3.1416/y(1)/sqrt(cs2)*omik/y(2) 
 	rho=-dotm0*exp(sint)/4./3.1416/y(1)/sqrt(cs2)*omik/y(2)

        p=rho*cs2
	cs=sqrt(cs2)
	h=cs/omik

	ww=-alfa/y(2)*cs2
        w1=-ww/y(2)
        w2=ww/(y(3)/mue+y(4)/mui)
        w3=ww/y(1)-2.*omiga

	a1=-1./y(2)
	a2=-1./2./(y(3)/mue+y(4)/mui)
c	a3=dlnok-1./y(1)
	a3=dlnok+(p0-1.)/y(1)
	
	deltap=1.-delta
c note that in the following the viscous dissipation 
c heat electron is considered
	b1=rho*y(2)*ru*a1*y(4)/mui-alfa*p*w1*deltap
	b2=-alfa*p*w2*deltap/mue+rho*y(2)*ru*a2*y(4)/mui/mue
	b3=b2*mue/mui-1./(gamai-1.)*rho*y(2)*ru/mui

c--------------------------------
	setae=1.686d-10*y(3)
	setai=9.181d-14*y(4)

	if(1./setai.gt.7.d2) then
        qie=4./4.*37.85d0/m*rho*rho*(y(4)-y(3))*(sqrt(2./3.1416)+sqrt(
     $          setai+setae))/(setai+setae)**1.5d0
        else

        bel=1.d0*(setae+setai)/setae/setai
        qie=2.366*20.d0/m*rho*rho*(y(4)-y(3))/mbsl4(2,1.d0/setae)
     $          /mbsl4(2,1.d0/setai)*((2.*(setae+setai)**2.d0+1.)/
     $          (setae+setai)*mbsl4(1,bel)+2.*mbsl4(0,bel))
c ***
c	print *,'Breakpoint 4'
        endif

c---------------------------------
 	b4=qie+alfa*p*w3*deltap-rho*y(2)*ru*a3*y(4)/mui
	
	ppp=rho*y(2)*ru/mue
	c1=a1*y(3)-alfa*delta*p*w1/ppp
	
	ate=1./setae*((3.*mbsl4(3,1.d0/setae)+mbsl4(1,1.d0/setae))/
     $       4./mbsl4(2,1.d0/setae)-1.d0)

	c2=a2*y(3)/mue-ate-alfa*delta*p*w2/mue/ppp
	c3=a2*y(3)/mui-alfa*delta*p*w2/mui/ppp

c now we calculate c_4, i.e., the radiation rate.
c at first, we transfer the units of all the physical quantities 
c from c=G=M=1 units into cgs units.

	rhorho=rho*6.1749d+5/m/m
	hh=h/6.77d-12*m
	
c the first, the breamsstrahlung radiation
	
	if(setae.lt.1.) then
	fseta=4.*sqrt(2.*setae/3.1416)/3.1416*(1.+1.781*setae**1.34d0)
     $		+1.73*setae**1.5d0*(1.+1.1*setae+setae*setae
     $		-1.25*setae**2.5d0)
	else
	fseta=9.d0*setae/2./3.1416*(log(1.123*setae+0.48)+1.5)+2.3*setae*
     $		(log(1.123*setae)+1.28)
	endif

	qbr=5.3d+25*rhorho*rhorho*fseta/mue/mue

c	goto 232

c the other two radiation mechanisms, i.e., synchrotron and comptonization
c are all related with the frequency, so we must calculat them in a integration
c as follows. This integration is based on the fomular eq.(41) in Manmoto 
c et al. 1997, ApJ, 489, with the integartion variaty nu---> log10(nu), signed 
c as tt.

	numax=log10(6.d+10*y(3))

c in the above, '2' is added in 2005.09.14; in the below, 13.--->14 also
c the uppler limit of frequency in the integration is determined
c by the bremsstrahlung emission: h \nu = k T
c ** changed first argument of simp2 from 12d0 to 90, following Renyi's suggestion
c	print*,setae
	call simp2(9.d0,numax,1.d-2,sum,setae,setai,rhorho,hh,qbr)
c	print*,setae
c	pause

c now take into account the Compton cooling due to the soft photons
c from the outside SSD.
c currently the units are cm.g.s, except rtran and r

        rr=y(1)/6.77d-12*m
	sup=10.d0*rtran/6.77d-12*m
	slow=rtran/6.77d-12*m
c       call simp(slow,sup,10,ssd,hh,rr,tau,rhorho,setae)

c the above 'ssd' denotes the cooling rate due to the Comptonization of
c the  soft photon emitted by SSD from r_tr to 10.**r_tr
c the amplification factor has been included since it is a function
c of the blackbody temperature, i.e., R.
c -------- 2003, June 14-------------

	if(log10(y(1)/2.).gt.1.5.and.sum/ssd.lt.6.d-1) then
c	print*,'sum/ssd=',sum/ssd,ssd,sum
c	pause
	endif

c in the following, the factor '2' (in SSD) is because the disk is two-sided
        qrad=sum*1.d0+ssd*0.d0

c transfer qrad to the dissipation rate per volum
	qneg=qrad/2./hh
	qn=qneg

c then transfer it from cgs units to the c=G=M=1 units

	qneg=qn*8.87d-27*m*m*m

	c4=(qneg-qie)/rho/y(2)/ru*mue-a3*y(3)+alfa*delta*w3*p/ppp

	q1=cs2

	d1=y(2)+q1*a1
	d2=q1*a2/mue+1.5/beta*ru/mue

c	d2=q1*a2/mue+1.0/beta*ru/mue

	d3=d2*mue/mui
	d4=(omiga*omiga-omik*omik)*y(1)-q1*a3

	dd=c1*d2*b3+d1*b2*c3+c2*d3*b1-b1*d2*c3-c2*d1*b3-b2*d3*c1
	dd1=c4*d2*b3+d4*b2*c3+c2*d3*b4-b4*d2*c3-c2*d4*b3-b2*d3*c4
	dd2=c1*d4*b3+d1*b4*c3+c4*d3*b1-b1*d4*c3-d1*c4*b3-b4*d3*c1
	dd3=c1*d2*b4+d1*b2*c4+c2*d4*b1-b1*d2*c4-d1*c2*b4-b2*d4*c1

	dy(1)=1.d0
	dy(2)=dd1/dd
	dy(3)=dd2/dd
	dy(4)=dd3/dd
	if(dy(4).gt.0..or.dy(3).gt.0.) then
	endif

232	return
	end


	subroutine simp(a,b,n,s,hh,r,tau,rhorho,setae)
        implicit real*8(a-h,o-z)
        h=(b-a)/float(2*n)
        s=0.5*(f(a,hh,r,tau,rhorho,setae)-f(b,hh,r,tau,rhorho,setae))
        do 1 i=1,n
1       s=s+2.*f(a+float(2*i-1)*h,hh,r,tau,rhorho,setae)
     $          +f(a+2.*h*float(i),hh,r,tau,rhorho,setae)
        s=(b-a)*s/float(3*n)
        return
        end

        function f(x0,hh,r,tau,rhorho,setae)
c here small 'r' denotes the location in ADAF
        implicit real*8(a-h,o-z)
	dimension repr(100),reph(100),emis(100),tau2(200)
        double precision msun,m,mui,mue
	common alfa,beta,dotm0,rout,p0,gamai,ru,
     $          sl0,m,mui,mue,delta,rtran

        msun=1.989d33
        g=6.67d-8
        pi=3.1416d0
	rg=g*m*1.d6*msun/9.d20

	r0=r
        x=x0
c now r0 is r in units of cm.g.s

        dotm=dotm0/3.44d-15/m*1.39d18*m*1.d6

c       fvis=3./8./pi*g*m*1.d6*msun/x**3.d0*dotm*(1.-rtran*rg/x)
 	fvis=3./8./pi*g*m*1.d6*msun/x**3.d0*dotm

c now consider the reprocessing effect
c all the related quantites such as the location in ADAF, height
c and emissivity are started with 'rep'
c 'n' is the number of the integration zone of ADAF
c reph(i), repr(i): the height and radius of a point in ADAF

	repn=10

	if(kkrep.eq.0) then
	do 103 i=1,repn
	read(14,*) repr(i),reph(i),emis(i),tau2(i)
c emis is in units of cm.g.s 

	repr(i)=repr(i)*rg
	reph(i)=reph(i)*rg
103	continue
	kkrep=1
	endif

	frep=0.d0
	do 50 i=1,repn-1
c	ff1=reph(i)/2.*pi/(reph(i)*reph(i)/4.+(x-repr(i))**2.d0)**1.5

	ff1=reph(i)/2.*pi/(reph(i)*reph(i)/4.+repr(i)*repr(i)
     $          *sin(22.5/57.)*sin(22.5/57.)+(x-repr(i))**2.d0)**1.5

	ff2=reph(i)*pi/(reph(i)*reph(i)/4.+repr(i)*repr(i)+x*x)**1.5
	ff3=reph(i)/2.*pi/(reph(i)*reph(i)/4.+(repr(i)+x)**2.d0)**1.5

	if(tau2(i).gt.300.d0) then
	frep=frep
	else
	frep=frep+0.8d0*emis(i)*(ff1+ff2+ff3)
     $		*(repr(i)-repr(i+1))*repr(i)*reph(i)*exp(-tau2(i)*1.)
	endif

c11111111111111
c frep: corresponds to fvis, the flux produce per unit area due to
c the reprocessing of the central X-ray source
c here 0.8 is the assumed approximated value of the albedo rate, 
c and this approximation is not bad

50	continue

c       fff1=hh/4./(hh*hh/4.+(x-r0)**2.d0)**1.5d0
	
	fff1=hh/4./(hh*hh/4.+(r0*sin(22.5/57.))**2.d0+
     $          (x-r0)**2.d0)**1.5d0

        fff2=hh/2./(x*x+r0*r0+hh*hh/4.)**1.5d0
        fff3=hh/4./(hh*hh/4.+(x+r0)**2.d0)**1.5d0

c compared to teh code of calculating the spectrum, 'pi' is 
c absent in the above formula. This is due to the difference
c of spectrum and flux. Integration of bnu over nu is sigmaT**4/pi!!

c the following calculate the amplification factor
c the formular comes from Luo & Liang (1994; eqs. 11 & 12), or equivalently, 
c from Dermer et al. 1991 (eqs. 16 & 17) note there is a slight difference
        taoes=0.4*rhorho*hh
        delt=sqrt(9./4.+pi*pi/12./taoes/(taoes+2./3.)/setae)-1.5d0
c	delt=sqrt(9./4.+pi*pi/12./(taoes+2./3.d0)/(taoes+2./3.)/setae)-1.5d0
c note: can't use delta!!!!
c        sigmarep=5.67d-5
        ts=((frep+fvis)/5.67d-5)**0.25d0
c ts: the black boby temperature of the SSD
        te2=setae/1.686d-10
        ampl=(1.-(1.275*ts/te2)**(delt-1.d0))/(delt-1.)

        f=ampl*(fff1+fff2+fff3)/exp(tau)*x*(frep+fvis)

c	print*,'fvis=',fvis,frep,ampl,fff2*x
c	pause

        return
        end


	function fct(tt,setae,setai,rhorho,hh,qbr)
        implicit real*8(a-h,o-z)
        double precision m,nu,kapa,mui,mue,mbsl4,
     $		mgam2,jm
	common alfa,beta,dotm0,rout,p0,gamai,ru,
     $          sl0,m,mui,mue,delta,rtran
c nu denotes the frequency(transfered into tt),
c te for the electron temperature
c all the quan. are in cgs units


        nu=10.d0**tt
        te=setae/1.686d-10
        ti=setai/9.18d-14

	bnu=1.47d-47*nu*nu*nu/(exp(4.8d-11*nu/te)-1.d0)
c ----------------------------

        if(2.08d+10/nu*te.lt.1.d0) then
        gaunt=sqrt(3./3.1416d0*4.8d-11/nu/te)
        else
        gaunt=4.8d-11/te*sqrt(3.)/3.1416*log(4./1.781*2.08d10*te/nu)
        endif

        qiabr=qbr*gaunt*exp(-4.8d-11*nu/te)


c now the synchrotron emission
        xx=5.29d-12*nu/setae/setae/sqrt(rhorho)/sqrt((1.-beta)/beta)
     $          /sqrt(ti/mui+te/mue)


        fxx=4.0505/xx**(1./6.d0)*(1.+0.4d0/xx**0.25d0+
     $          0.5316d0/sqrt(xx))*exp(-1.8899d0*xx**0.33333333333d0)

        fk2=mbsl4(2,1.d0/setae)

c       qiasy=4.43d-6*4.*3.1416*rhorho/1.67*nu*fxx/fk2
	qiasy=33.33482155d-6*rhorho*nu*fxx/fk2/mue
c	qiasy=0.d0

c -----------------------------
        qianu=qiabr+qiasy

        kapa=qianu/4./3.1416/bnu

        taonux=sqrt(3.1416)/2.*kapa*hh
c here, taonux denotes taonu star

	if(taonux.gt.1.d-5) then
        fnu=2.*3.1416d0/sqrt(3.d0)*bnu*
     $          (1.-dexp(-2.d0*sqrt(3.)*taonux))
        else
        fnu=2.*3.1416d0/sqrt(3.d0)*bnu*(2.*sqrt(3.)*taonux-
     $          6.*taonux*taonux+4.*sqrt(3.)*taonux**3.d0)
        endif

        fnu=2.*3.1416d0/sqrt(3.d0)*bnu*(1.-exp(-2.*sqrt(3.)*taonux))
c now, the comptonization

	a=1.+4.*setae+16.d0*setae*setae
 	taoes1=0.8d0*rhorho*hh

	if(taoes1.gt.4.0) then
	taoes=taoes1
	else
	a1=0.552d0
	a2=0.713d0
	a3=0.182d0
	a4=0.315d0
	taoese=a1*taoes1**a2*setae**(-a3/taoes1**a4)
	if(taoese.gt.taoes1) taoes=taoes1
	if(taoese.le.taoes1) taoes=(taoese+taoes1)/2.
	endif

c	taoes=taoes1

        s=taoes+taoes*taoes
        etamax=3.d0*2.084d+10*te/nu
        if(etamax.gt.1.) then
        jm=log(etamax)/log(a)
        else
        print*,'# etamax < 1 !!!',etamax,te,nu
        stop
        jm=0.d0
        endif

        eta=exp(s*(a-1.d0))*(1.-mgam2(jm+1.,a*s))+etamax*mgam2(jm+1.,s)

        fct=eta*2.*fnu*nu*log(10.d0)

        end



	function mgam2(a,x)
	implicit real*8(a-h,o-z)
	DOUBLE PRECISION MGAM2,A,X
	DOUBLE PRECISION MGAM1,P,Q,D,S,S1,P0,Q0,P1,Q1,QQ
	IF((A.LE.0.0).OR.(X.LT.0.0)) THEN
	IF(A.LE.0.0) THEN
	WRITE(*,*) 'ERR** A <= 0!'
	stop
	ENDIF
	IF(X.LT.0.0) THEN
	WRITE(*,*) 'ERR * * X < 0!'
	stop
	ENDIF
	MGAM2=-1.0
	ENDIF
	IF(X+1.0.EQ.1.0) THEN
	MGAM2=0.0
	RETURN
	ENDIF
	IF(X.GT.1.0D+35) THEN
	MGAM2=1.0
	RETURN
	ENDIF
	Q=LOG(X)
	Q=A*Q
	QQ=EXP(Q)
	IF(X.LT.1.0+A) THEN
	P=A
	D=1.0/A
	S=D
	DO 10 N=1,100
	P=1.0+P
	D=D*X/P
	S=S+D
	IF(ABS(D).LT.ABS(S)*1.0D-7) THEN
	S=S*EXP(-X)*QQ/MGAM1(A)
	MGAM2=S
	RETURN
	ENDIF
10	CONTINUE
	ELSE
	S=1.0/X
	P0=0.0
	P1=1.0
	Q0=1.0
	Q1=X
	DO 20 N=1,100
	P0=P1+(N-A)*P0
	Q0=Q1+(N-A)*Q0
	P=X*P0+N*P1
	Q=X*Q0+N*Q1
	IF(ABS(Q)+1.0.NE.1.0) THEN
	S1=P/Q
	P1=P
	Q1=Q
	IF(ABS((S1-S)/S1).LT.1.0D-7) THEN
	S=S1*EXP(-X)*QQ/MGAM1(A)
	MGAM2=1.0-S
	RETURN
	ENDIF
	S=S1
	ENDIF
	P1=P
	Q1=Q
20 	CONTINUE
	ENDIF
	WRITE(*,*) 'A TOO LARGE'
	S=1.0-S*EXP(-X)*QQ/MGAM1(A)
	MGAM2=S
	RETURN
	END


	
	FUNCTION MGAM1(X)
	DOUBLE PRECISION MGAM1,X
	DOUBLE PRECISION Y,T,S,U,A(11)
	DATA A/0.0000677106,-0.0003442342,0.0015397681,-0.0024467480,
     $	0.0109736958,-0.0002109075,0.0742379071,0.0815782188,
     $	0.4118402518,
     $	0.4227843370,1.0/
	IF(X.LE.0.) THEN
	WRITE(*,*) 'ERR * * X< 0!'
	MGAM1=-1.0
	RETURN
	ENDIF
	Y=X
	IF(Y.LE.1.0) THEN
	T=1.0/(Y*(Y+1.0))
	Y=Y+2.0
	ELSE IF(Y.LE.2.0) THEN
	T=1.0/Y
	Y=Y+1.0
	ELSE IF(Y.LE.3.0) THEN
	T=1.0
	ELSE
	T=1.0
10 	IF(Y.GT.3.0) THEN
	Y=Y-1.0
	T=T*Y
	GOTO 10
	ENDIF
	ENDIF
	S=A(1)
	U=Y-2.0
	DO 20 I=1,10
20	S=S*U+A(I+1)
	S=S*T
	MGAM1=S
	RETURN
	END
	


	function fct3(x,a)
	implicit real*8(a-h,o-z)
	fct3=x**(a-1.)*exp(-x)
	end

	SUBROUTINE SIMP3(A,B,EPS,SUM,aa)
        implicit real*8(a-h,o-z)
        DIMENSION F(2,30),FM(2,30),E(2,30),KRTN(30)
        SUM=0.d0 
        T=1.d0
        ABSA=1.0 
        EST=1.d0 
        DA=B-A
        FA=FCT3(A,aa)
        FB=FCT3(B,aa)
        FP=4.0*FCT3(0.5*(A+B),aa)
        X=A
        L=0
 1      K=1
        L=L+1
        T=T*1.7
        DX=DA/3.0
        SX=DX/6.0
        FM1=4.0*FCT3(X+0.5*DX,aa)
        F1=FCT3(X+DX,aa)
	FM(1,L)=FP
        F(1,L)=FCT3(X+2.0*DX,aa)
        FM(2,L)=4.0*FCT3(X+2.5*DX,aa)
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




	SUBROUTINE SIMP2(A,B,EPS,SUM,setae,setai,rhorho,hh,qbr)
	implicit real*8(a-h,o-z)
	DIMENSION F(2,30),FM(2,30),E(2,30),KRTN(30)

	SUM=0.d0
	T=1.d0
	ABSA=1.0
	EST=1.d0
	DA=B-A
	FB=FCT(B,setae,setai,rhorho,hh,qbr)
	FP=4.0*FCT(0.5*(A+B),setae,setai,rhorho,hh,qbr)
	X=A
	L=0
 1	K=1
	L=L+1
	T=T*1.7
	DX=DA/3.0
	SX=DX/6.0
	FM1=4.0*FCT(X+0.5*DX,setae,setai,rhorho,hh,qbr)
	F1=FCT(X+DX,setae,setai,rhorho,hh,qbr)
	FM(1,L)=FP
	F(1,L)=FCT(X+2.0*DX,setae,setai,rhorho,hh,qbr)
	FM(2,L)=4.0*FCT(X+2.5*DX,setae,setai,rhorho,hh,qbr)
	F(2,L)=FB
	E1=SX*(FA+FM1+F1)
	E(1,L)=SX*(F1+FP+F(1,L))
	E(2,L)=SX*(F(1,L)+FM(2,L)+FB)
	S=E1+E(1,L)+E(2,L)
	ABSA=ABSA-ABS(EST)+ABS(E1)+ABS(E(1,L))+ABS(E(2,L))
	IF(EST.EQ.1.0) GOTO 5
	IF(T*ABS(EST-S).LE.EPS*ABSA) GOTO 2
	IF(L.LT.30) GOTO 5
2	SUM=SUM+S
3	L=L-1
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
4	DA=DX
	KRTN(L)=K
	GOTO 1
5	EST=E1
	FP=FM1
	FB=F1
	GOTO 4
6	RETURN
	END



	subroutine mers(n,inform,h,eps,y,dy,w,qn,qd,tau,sint)
        implicit real*8(a-h,o-z)
        dimension y(1),dy(1),w(4,1),a(4),b(4),c(4)
        data a(1),a(2),a(3),a(4)/0.3333333333,0.16666666667,0.375,2.0/,
     $    b(1),b(2),b(3),b(4)/1.d0,0.444444444444,0.9375d0,8.d0/,
     $    c(1),c(2),c(3),c(4)/0.d0,0.25d0,-0.2d0,-1.125d0/

        l=0
        k=w(3,1)
        if(inform.ge.10) goto 14
        ht=h/w(3,1)
        if(inform.ge.3) goto 3
        epsl=0.0781250*eps
 1      do 2 j=1,n
        w(3,j)=y(j)
 2      w(4,j)=dy(j)
  3     do 5 i=1,4
        do 4 j=1,n
        y(j)=y(j)+ht*a(i)*(dy(j)-w(2,j))
        w(2,j)=b(i)*(dy(j)+c(i)*w(1,j))
 4      if(i.eq.3) w(1,j)=dy(j)-0.22222222220*w(1,j)
 5      call dery(n,y,dy,qn,qd,tau,sint)
c ***
	print *,'Breakpoint 2'
	if(sig.eq.-1.) goto 17
        i=1
        do 9 j=1,n
        r=0.16666666670*ht*(dy(j)-w(2,j))
        y(j)=y(j)+r
        goto (7,6,9),inform
  6     if(y(j).ne.0.0) r=r/y(j)
 7      if(abs(r*0.2).gt.eps) if(y(1)-(0.125*ht+y(1))) 12,8,12
 8      if(abs(r).gt.epsl) i=3
 9      continue
        if(i.ne.l-l/2*2) goto 15
        k=k/2
        l=1+l/2
        ht=ht*2
 10     call dery(n,y,dy,qn,qd,tau,sint)
c ***
	print *,'Breakpoint 3'
        if(sig.eq.-1.) goto 17
        do 11 j=1,n
        w(1,j)=dy(j)
 11     w(2,j)=0.0
        if(l.ne.k) goto 1
        w(3,1)=k
17	inform=2
         return
 12     ht=0.5*ht


	cshu=cshu+1
        if(cshu.gt.150.) then
         inform=3
         print*,'###########################'
        endif

        l=l+l
        k=k+k
        do 13 j=1,n
        y(j)=w(3,j)
        dy(j)=w(4,j)
        w(1,j)=dy(j)
 13     w(2,j)=0.
        goto 3
 14     inform=inform-10
        k=1
 15     l=l+1
        goto 10
        end


	function cc(x)
	double precision c(7)
	data c/1.25331414,-0.07832358,0.02189568,-0.01062446,
     $          0.00587872,-0.0025154,0.00053208/
	p=c(7)
	do 100 i=6,1,-1
 100	p=p*x+c(i)
	cc=p
	return
	end
	
	function ddd(x)
	double precision d(7)
	data d/1.25331414,0.23498619,-0.0365562,0.01504268,
     $          -0.00780353,0.00325614,-0.00068245/
        p=d(7)
	do 101 i=6,1,-1
 101    p=p*x+d(i)
        ddd=p  
        return
        end
	
	
	function mbsl4(n,x)
	double precision mbsl4,x,mbsl3
	double precision y,p,b0,b1,a(7),b(7),c(7),d(7)
	data a/-0.57721566,0.4227842,0.23069756,0.0348859,
     $		0.00262698,0.0001075,0.0000074/
	data b/1.0,0.15443144,-0.67278579,-0.18156897,
     $		-0.01919402,-0.00110404,-0.00004686/
 	data c/1.25331414,-0.07832358,0.02189568,-0.01062446,
     $		0.00587872,-0.0025154,0.00053208/
	data d/1.25331414,0.23498619,-0.0365562,0.01504268,
     $		-0.00780353,0.00325614,-0.00068245/
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
 10	 p=p*y+a(i)
	p=p-mbsl3(0,x)*log(x/2.)
	else
	y=2.0/x
	p=c(7)
	do 20 i=6,1,-1
 20	p=p*y+c(i)
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
 30	p=p*y+b(i)
	p=p/x+mbsl3(1,x)*log(x/2.)
	else
	y=2./x
	p=d(7)
	do 40 i=6,1,-1
 40	p=p*y+d(i)
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
 50	continue
	mbsl4=p
	return
	end

C THIS PROGRAM IS TO CALCULATE THE MODIFIED BASEL FUNCTION OF THE FIRST KIND
C IT IS CORRECT!!!!!!!!!!!!!!!
	function mbsl3(n,x)
	double precision mbsl3,x
	double precision t,y,p,b0,b1,q,a(7),b(7),c(9),d(9)
	data a/1.0, 3.5156229,3.0899424,1.2067492,
     $		0.2659732,0.0360768,0.0045813/
	data b/0.5,0.87890594,0.51498869,0.15084934,0.02658773,
     $		0.00301532,0.00032411/
	data c/0.39894228,0.01328592,0.00225319,-0.00157565,
     $		0.00916281,-0.02057706,0.02635537,-0.01647663,0.00392377/
	data d/0.39894228,-0.03988024,-0.00362018,0.00163801,
     $		-0.01031555,0.02282967,-0.02895312,0.01787654,-0.00420059/

	if(n.lt.0) n=-n
	t=abs(x)
	if (n.ne.1) then
	 if(t.lt.3.75) then
	   y=(x/3.75)*(x/3.75)
	   p=a(7)
	   do 10 i=6,1,-1
  10	   p=p*y+a(i)
	 else
	   y=3.75/t
	   p=c(9)
	   do 20 i=8,1,-1
 20	   p=p*y+c(i)
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
 30	p=p*y+b(i)
	p=p*t
	else
	y=3.75/t
	p=d(9)
	do 40 i=8,1,-1
 40	p=p*y+d(i)
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
 50	continue
	p=t*q/b1
	if((x.lt.0.0).and.(mod(n,2).eq.1)) p=-p
	mbsl3=p
	return
	end


