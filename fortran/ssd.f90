! this code calculate the spectrum of a standard thin disk
! the code is based on the formula in the book "accretion
! power in astrophysics" by Frank et al. 1992
! in unit of cm.g.sec

! This is a modification of Feng's code ssd.f. Since the effect of the
! illumination of the thin disk by radiation from the ADAF seems
! to be negligible, I removed the reprocessing effect.


module globals
real(kind=8),save :: dotm,mass,rs,rin,dotme,rout,theta
real(kind=8),parameter :: h=6.626d-27, c=2.9979d+10, k=1.3807d-16, pi=3.1416d0, &
		pc=3.086d+18, solarmass=1.989d+33, g=6.672d-8, sigma=5.67d-5
end module globals




program ssd
use globals
implicit none
real :: y(4),dy(4),w(16),nu(200),fnu(200)
real(kind=8) :: nnu,s
integer :: i,n
	
! Output file
open(20,file='ssd_alone.dat')

read(*,*) mass
mass=mass*1.d6*solarmass

rs=2.*g*mass/c/c
! Note in the above defination, the factor of 0.1 is included,
! the unit is g/s
dotme=1.39d18*mass/solarmass

read(*,*) dotm
dotm=dotm*dotme
	
read(*,*) rin
rin=rin*rs
rout=1d5*rs

read(*,*) theta
theta=theta/57.2958 !--------------6-----------------!


do i=1,50
	nu(i)=10.d0**(11.+0.1d0*i)

	nnu=nu(i)
	s=0.d0
	call simpsn(rin,rout,1d-3,8000,n,s,nnu)

	fnu(i)=4.*pi*h*cos(theta)*nu(i)**3.d0/c/c*s 
	write(20,*) log10(nu(i)),log10(nu(i)*fnu(i)*4.*pi)
end do

end program
	
	
	
	
	
	
	
	
function f(x,nu)
use globals
implicit none
real(kind=8) :: nu
real :: x,tem,fvis,f

! Dissipation rate, see the first page of Chapter 8 in <<accretion>> book
fvis=3./8./pi/sigma*g*mass/x**3.d0*dotm*(1.-(3.*rs/x)**0.5d0)
	
! Temperature
tem=fvis**0.25d0

if (h*nu/(k*tem).gt.700.) then
	f=0.d0
else
    f=x/(dexp(h*nu/(k*tem))-1.d0)
endif

return
end function







subroutine simpsn(a,b,eps,k,n,s,nu)
implicit none
real(kind=8) :: nu
integer :: n,i,k
real :: a,b,c,h,t1,f,t2,x,s,s0,eps

n=1
c=abs(a)+abs(b)
h=0.5*(b-a)    
t1=h*(f(a,nu)+f(b,nu))

do while (c+h /= c)
	x=a-h
    t2=0.5*t1
    do i=1,n
		x=x+2.*h
		t2=t2+h*f(x,nu)
	end do

    s=0.33333333333d0*(4.*t2-t1)
    n=n+n
    
    if (abs(s-s0).le.(1.+abs(s))*eps) return
	s0=s
    t1=t2
    h=h*0.5
end do

n=-n
return
end subroutine
