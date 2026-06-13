program adv_yin_yang_CMC42_RK4
!********************************************************************************
! This program was written by Rasoul Mirzaei Shiri.
!
! It solves the Advection equations on a Rectangular Yin-Yang grid (RYY) using a
! fourth-order MacCormack scheme with Runge-Kutta time integration.
!
! The code is developed and employed to solve the 'Advection of the cosine bell over the pole' test case
! for validation and assessment of the numerical method.
!********************************************************************************
implicit none
integer,parameter::nx=256,ny=(nx/2)+1,n_epsilon_x=nx/32,n_epsilon_y=(ny-1)/16,mx=(3*nx/4)+1+2*n_epsilon_x&
,my=((ny-1)/2)+1+2*n_epsilon_y,days=12
integer::i,l,n,j,knex,kney
real*8,parameter::delt=900.0d0
real*8::a2,a3,a4,t,lx,ly,pi,twopi,dlanda,dphi&
,dtime,ta(2),bbbb,top_x_limit,top_y_limit
REAL*8::alpha,hene(nx,ny),pne(ny),lne(nx),norm_1,norm_2,norm_in,ut(nx,ny)
real*8::hee1(mx,my),hee2(mx,my),hen1(mx,my),hen2(mx,my),cnx(mx,my),cny(mx,my),cex(mx,my),cey(mx,my)&
,dccy(mx,2*(ny-1)),dccx(nx,my),hhn4(mx,my),hn1(mx,my),hn2(mx,my),hn3(mx,my),hhn1(mx,my),hhn2(mx,my)&
,hhn3(mx,my),hn4(mx,my),he1(mx,my),he2(mx,my)&
,he3(mx,my),hhe1(mx,my),hhe2(mx,my),hhe3(mx,my),hhe4(mx,my),he4(mx,my),dcx(mx,my),dcy(mx,my),ln(mx)&
,pn(my),le(mx),pe(my),he_zone(nx,my),he_meridion(mx,2*(ny-1))
character(5)::int_method
open(10,file='yin.plt')
open(20,file='yang.plt')
open(30,file='yin_yang.plt')
open(40,file='norm.dat')
open(50,file='comptime.dat')
!int_method='bilin'
int_method='bicub'
kney=(ny-my)/2
knex=(nx-mx+1)/2
print*,n_epsilon_x,n_epsilon_y
print*,knex,kney
print*,mx,my
print*,int_method
pi=2.0d0*dasin(1.0d0)
twopi=2.0d0*pi

lx=twopi;ly=pi
alpha=pi/2.0d0
dlanda=lx/dble(nx)
dphi=ly/dble(ny-1)
t=dble(days)*86400.0d0;l=int(t/delt)
top_x_limit=(pi*3.0d0/4.0d0)+dble(n_epsilon_x)*dlanda
top_y_limit=(pi/4.0d0)+dble(n_epsilon_y)*dphi
write(10,*)'zone            i=',mx,'   j=',my
write(20,*)'zone            i=',my,'   j=',mx
write(30,*)'zone            i=',nx,'   j=',ny
call initial_condition_yin(mx,my,cnx,cny,hen1,alpha,n_epsilon_x,n_epsilon_y,dlanda,dphi)
call initial_condition_yang(mx,my,cex,cey,hee1,alpha,n_epsilon_x,n_epsilon_y,dlanda,dphi)
call initial_condition_glob(nx,ny,ut,dlanda,dphi)
do i=1,mx
ln(i)=-top_x_limit+dble(i-1)*dlanda
enddo
do j=1,my
pn(j)=dble(j-1)*dphi-top_y_limit
enddo


do i=1,mx
le(i)=-top_x_limit+dble(i-1)*dlanda
enddo
do j=1,my
pe(j)=dble(j-1)*dphi-top_y_limit
enddo
do j=1,my
do i=1,mx
write(10,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,hen1(i,j)
enddo
enddo

do i=1,mx
do j=1,my
write(20,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,hee1(i,j)
enddo
enddo
do i=1,nx
lne(i)=-pi+dble(i-1)*dlanda
enddo
do j=1,ny
pne(j)=dble(j-1)*dphi-(pi/2.0d0)
enddo
do j=1,ny
do i=1,nx
write(30,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,ut(i,j)
enddo
enddo

a2=0.5d0;a3=0.5d0;a4=1.0d0

do n=1,l
!********************************************************************
!first step
!********************************************************************
call provide_Yin_scalar_to_derive(hen1,hee1,dlanda,dphi,nx,ny,mx,my&
,n_epsilon_x,n_epsilon_y,he_zone,he_meridion,int_method)

call dby42p(he_meridion,dccy,mx,2*(ny-1),dphi)
call dfx42p(he_zone,dccx,nx,my,dlanda)
do i=1,mx
do j=1,my
dcx(i,j)=dccx(i+knex,j)
dcy(i,j)=dccy(i,j+kney)
enddo
enddo
call RK_init(hen1,dcx,dcy,mx,my,cnx,cny,a2,hn1,hhn1,pn,delt)

call provide_Yang_scalar_to_derive(hen1,hee1,dlanda,dphi,nx,ny,mx,my&
,n_epsilon_x,n_epsilon_y,he_zone,he_meridion,int_method)
call dby42p(he_meridion,dccy,mx,2*(ny-1),dphi)
call dfx42p(he_zone,dccx,nx,my,dlanda)
do i=1,mx
do j=1,my
dcx(i,j)=dccx(i+knex,j)
dcy(i,j)=dccy(i,j+kney)
enddo
enddo
call RK_init(hee1,dcx,dcy,mx,my,cex,cey,a2,he1,hhe1,pe,delt)
!********************************************************************
!second step
!********************************************************************

call provide_Yin_scalar_to_derive(hhn1,hhe1,dlanda,dphi,nx,ny,mx,my&
,n_epsilon_x,n_epsilon_y,he_zone,he_meridion,int_method)

call dfy42p(he_meridion,dccy,mx,2*(ny-1),dphi)
call dbx42p(he_zone,dccx,nx,my,dlanda)
do i=1,mx
do j=1,my
dcx(i,j)=dccx(i+knex,j)
dcy(i,j)=dccy(i,j+kney)
enddo
enddo
call RK_init(hen1,dcx,dcy,mx,my,cnx,cny,a3,hn2,hhn2,pn,delt)

call provide_Yang_scalar_to_derive(hhn1,hhe1,dlanda,dphi,nx,ny,mx,my&
,n_epsilon_x,n_epsilon_y,he_zone,he_meridion,int_method)
call dfy42p(he_meridion,dccy,mx,2*(ny-1),dphi)
call dbx42p(he_zone,dccx,nx,my,dlanda)
do i=1,mx
do j=1,my
dcx(i,j)=dccx(i+knex,j)
dcy(i,j)=dccy(i,j+kney)
enddo
enddo
call RK_init(hee1,dcx,dcy,mx,my,cex,cey,a3,he2,hhe2,pe,delt)




!********************************************************************
!third step
!********************************************************************
call provide_Yin_scalar_to_derive(hhn2,hhe2,dlanda,dphi,nx,ny,mx,my&
,n_epsilon_x,n_epsilon_y,he_zone,he_meridion,int_method)

call dby42p(he_meridion,dccy,mx,2*(ny-1),dphi)
call dfx42p(he_zone,dccx,nx,my,dlanda)
do i=1,mx
do j=1,my
dcx(i,j)=dccx(i+knex,j)
dcy(i,j)=dccy(i,j+kney)
enddo
enddo
call RK_init(hen1,dcx,dcy,mx,my,cnx,cny,a4,hn3,hhn3,pn,delt)

call provide_Yang_scalar_to_derive(hhn2,hhe2,dlanda,dphi,nx,ny,mx,my&
,n_epsilon_x,n_epsilon_y,he_zone,he_meridion,int_method)
call dby42p(he_meridion,dccy,mx,2*(ny-1),dphi)
call dfx42p(he_zone,dccx,nx,my,dlanda)
do i=1,mx
do j=1,my
dcx(i,j)=dccx(i+knex,j)
dcy(i,j)=dccy(i,j+kney)
enddo
enddo
call RK_init(hee1,dcx,dcy,mx,my,cex,cey,a4,he3,hhe3,pe,delt)
!********************************************************************
!fourth step
!********************************************************************


call provide_Yin_scalar_to_derive(hhn3,hhe3,dlanda,dphi,nx,ny,mx,my&
,n_epsilon_x,n_epsilon_y,he_zone,he_meridion,int_method)

call dfy42p(he_meridion,dccy,mx,2*(ny-1),dphi)
call dbx42p(he_zone,dccx,nx,my,dlanda)
do i=1,mx
do j=1,my
dcx(i,j)=dccx(i+knex,j)
dcy(i,j)=dccy(i,j+kney)
enddo
enddo
call RK_init(hen1,dcx,dcy,mx,my,cnx,cny,a4,hn4,hhn4,pn,delt)
call RK4_final(hen1,hen2,hn1,hn2,hn3,hn4,mx,my)
call provide_Yang_scalar_to_derive(hhn3,hhe3,dlanda,dphi,nx,ny,mx,my&
,n_epsilon_x,n_epsilon_y,he_zone,he_meridion,int_method)
call dfy42p(he_meridion,dccy,mx,2*(ny-1),dphi)
call dbx42p(he_zone,dccx,nx,my,dlanda)
do i=1,mx
do j=1,my
dcx(i,j)=dccx(i+knex,j)
dcy(i,j)=dccy(i,j+kney)
enddo
enddo
call RK_init(hee1,dcx,dcy,mx,my,cex,cey,a4,he4,hhe4,pe,delt)
call RK4_final(hee1,hee2,he1,he2,he3,he4,mx,my)
















if(n*delt/21600.0d0==dble(floor(n*delt/21600.0d0)))then
write(10,*)'zone            i=',mx,' j=',my
write(20,*)'zone            i=',my,'   j=',mx
write(30,*)'zone            i=',nx,' j=',ny
do j=1,my
do i=1,mx
write(10,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,hen2(i,j)
enddo
enddo
do i=1,mx
do j=1,my
write(20,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,hee2(i,j)
enddo
enddo
call merging_Yin_Yang(hen2,hee2,hene,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
do j=1,ny
do i=1,nx
write(30,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,hene(i,j)
enddo
enddo
endif
print*,n*delt,dcx(1,1),dcx(2,1)
do i=1,mx
do j=1,my
hen1(i,j)=hen2(i,j)
hee1(i,j)=hee2(i,j)
enddo
enddo 
enddo
call norm_sphere(hene,ut,nx,ny,dlanda,dphi,pne,norm_1,norm_2,norm_in)
write(40,*)'norm1',norm_1
write(40,*)'norm2',norm_2
write(40,*)'norminfinite',norm_in
bbbb=dtime(ta)
write(50,*)'computational time=',bbbb
end program adv_yin_yang_CMC42_RK4
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  dbx42p(e,dbx,m,k,dx)
integer::i,j,m,k
real*8::e(m,k),dbx(m,k),a,dx,sum
a=0.5d0-0.5d0/dsqrt(3.0d0)
do j=1,k
sum=(e(1,j)-e(m,j))/dx
do i=m,2,-1
sum=sum+((e(i,j)-e(i-1,j))/dx)*((a/(a-1.0d0))**(m-i+1))
enddo
dbx(1,j)=sum/((1.0d0-a)+(((a/(a-1.0d0))**(m-1))*a))
do i=2,m
dbx(i,j)=(((e(i,j)-e(i-1,j))/dx)-(a*dbx(i-1,j)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
subroutine  dfx42p(e,dfx,m,k,dx)
integer::i,j,m,k
real*8::e(m,k),dfx(m,k),a,dx,sum
a=0.5d0-0.5d0/dsqrt(3.0d0)
do j=1,k
sum=(e(1,j)-e(m,j))/dx
do i=2,m
sum=sum+((e(i,j)-e(i-1,j))/dx)*((a/(a-1.0d0))**(i-1))
enddo
dfx(m,j)=sum/((1.0d0-a)+(((a/(a-1.0d0))**(m-1))*a))
do i=m-1,1,-1
dfx(i,j)=(((e(i+1,j)-e(i,j))/dx)-(a*dfx(i+1,j)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
subroutine  dfy42p(e,dfy,m,k,dy)
integer::i,j,m,k
real*8::e(m,k),dfy(m,k),a,dy,sum
a=0.5d0-0.5d0/dsqrt(3.0d0)
do i=1,m
sum=(e(i,1)-e(i,k))/dy
do j=2,k
sum=sum+((e(i,j)-e(i,j-1))/dy)*((a/(a-1.0d0))**(j-1))
enddo
dfy(i,k)=sum/((1.0d0-a)+(((a/(a-1.0d0))**(k-1))*a))
do j=k-1,1,-1
dfy(i,j)=(((e(i,j+1)-e(i,j))/dy)-(a*dfy(i,j+1)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
subroutine  dby42p(e,dby,m,k,dy)
integer::i,j,m,k
real*8::e(m,k),dby(m,k),a,dy,sum
a=0.5d0-0.5d0/dsqrt(3.0d0)
do i=1,m
sum=(e(i,1)-e(i,k))/dy
do j=k,2,-1
sum=sum+((e(i,j)-e(i,j-1))/dy)*((a/(a-1.0d0))**(k-j+1))
enddo
dby(i,1)=sum/((1.0d0-a)+(((a/(a-1.0d0))**(k-1))*a))
do j=2,k
dby(i,j)=(((e(i,j)-e(i,j-1))/dy)-(a*dby(i,j-1)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  yinvector_to_yangvector(phin,landan,phie,landae,vphin,vlandan,vphie,vlandae)
real*8::phin,landan,phie,landae,pi,phin1,phie1,vphin,vlandan,vphie,vlandae,cossaay,sinsaay
pi=dasin(1.0d0)*2.0d0
phin1=phin+(pi/2.0d0)
phie1=phie+(pi/2.0d0)
cossaay=-dsin(landae)*dsin(landan)
if(dabs(dsin(phin1))<1.0d-8)then
sinsaay=dcos(landan)/dsin(phie1)
else
sinsaay=-dcos(landae)/dsin(phin1)
endif
vphie=(vphin*cossaay)-(vlandan*sinsaay)
vlandae=(vphin*sinsaay)+(vlandan*cossaay)
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  yangvector_to_yinvector(phin,landan,phie,landae,vphin,vlandan,vphie,vlandae)
real*8::phin,landan,phie,landae,pi,phin1,phie1,vphin,vlandan,vphie,vlandae,cossaay,sinsaay
pi=dasin(1.0d0)*2.0d0
phin1=phin+(pi/2.0d0)
phie1=phie+(pi/2.0d0)
cossaay=-dsin(landae)*dsin(landan)
if(dabs(dsin(phin1))<1.0d-8)then
sinsaay=dcos(landan)/dsin(phie1)
else
sinsaay=-dcos(landae)/dsin(phin1)
endif
vphin=(vphie*cossaay)+(vlandae*sinsaay)
vlandan=-(vphie*sinsaay)+(vlandae*cossaay)
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine num(phi3,landa3,k1,k2,dlanda,dphi,n_epsilon_x,n_epsilon_y)
implicit none
integer::k1,k2,n_epsilon_x,n_epsilon_y
real*8::pi,landa3,phi3,dphi,dlanda,top_x_limit,top_y_limit
pi=2.d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)
k1=int(((landa3+top_x_limit+1.0d-10)/dlanda))+1
k2=int(((phi3+top_y_limit+1.0d-10)/dphi))+1
end subroutine  
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine norm_sphere(hene,ut,nx,ny,dlanda,dphi,phi,norm_1,norm_2,norm_in)
implicit none
integer::nx,ny,i,j
real*8::hene(nx,ny),ut(nx,ny),dlanda,dphi,phi(ny),norm_1,norm_2,norm_in,s1,s2
s1=0;s2=0
do i=1,nx
do j=1,ny
s1=s1+dabs(hene(i,j)-ut(i,j))*dabs(hene(i,j)-ut(i,j))*dcos(phi(j))*dlanda*dphi
s2=s2+dabs(ut(i,j))*dabs(ut(i,j))*dcos(phi(j))*dlanda*dphi
enddo
enddo
norm_2=dsqrt(s1)/dsqrt(s2)
s1=0;s2=0
do i=1,nx
do j=1,ny
s1=s1+dabs(hene(i,j)-ut(i,j))*dcos(phi(j))*dlanda*dphi
s2=s2+dabs(ut(i,j))*dcos(phi(j))*dlanda*dphi
enddo
enddo
norm_1=(s1)/(s2)
s1=0;s2=0
do i=1,nx
do j=1,ny
if(dabs(hene(i,j)-ut(i,j))>s1)then
s1=dabs(hene(i,j)-ut(i,j))
endif
if(dabs(ut(i,j))>s2)then
s2=dabs(ut(i,j))
endif
enddo
enddo
norm_in=(s1)/(s2)
return
end subroutine
!******************************************************************************************
!******************************************************************************************
!******************************************************************************************
!******************************************************************************************
subroutine  RK_init(he,dcx,dcy,nx,ny,cx,cy,a,h,hh,phi,delt)
implicit none
integer::i,j,nx,ny
real*8::he(nx,ny),dcx(nx,ny),dcy(nx,ny),cx(nx,ny),cy(nx,ny),h(nx,ny)&
,hh(nx,ny),phi(ny),a,delt
do i=1,nx
do j=1,ny
call Advection_sphere(h(i,j),dcx(i,j),dcy(i,j),cx(i,j),cy(i,j),phi(j),delt)
hh(i,j)=a*h(i,j)+he(i,j)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  RK4_final(he1,he2,h1,h2,h3,h4,nx,ny)
implicit none
integer::i,j,nx,ny
real*8::he1(nx,ny),he2(nx,ny),h1(nx,ny),h2(nx,ny),h3(nx,ny),h4(nx,ny)&
,b1,b2,b3,b4
b1=1.0d0/6.0d0;b2=1.0d0/3.0d0;b3=1.0d0/3.0d0;b4=1.0d0/6.0d0
do i=1,nx
do j=1,ny
he2(i,j)=he1(i,j)+b1*h1(i,j)+b2*h2(i,j)+b3*h3(i,j)+b4*h4(i,j)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  Advection_sphere(h,dcx,dcy,cx,cy,phi,delt)
implicit none
real*8::ux,uy,aa,dcx,dcy,cx,cy,phi,h,delt

aa=6.37122d6
ux=dcx*cx/(aa*dcos(phi))
uy=dcy*cy/aa
h=-delt*(ux+uy)
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  merging_Yin_Yang(Fn,Fe,Fne,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
implicit none
integer::nx,ny,my,mx,kney,knex,i,j,n_epsilon_x,n_epsilon_y
real*8::pi,top_x_limit,top_y_limit,dphi,dlanda,Fn(mx,my),Fe(mx,my),Fne(nx,ny)&
,a10,a11,a12,a13,hem,phi(ny),landa(nx)
character(5)::int_method
kney=(ny-my)/2
knex=(nx-mx+1)/2
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+dble(n_epsilon_x)*dlanda
top_y_limit=(pi/4.0d0)+dble(n_epsilon_y)*dphi

do j=1,ny
phi(j)=-(pi/2.0d0)+dble(j-1)*dphi
do i=1,nx
landa(i)=-pi+dble(i-1)*dlanda
if(phi(j)>=-top_y_limit.and.phi(j)<=top_y_limit.and.landa(i)>=-top_x_limit.and.landa(i)<=top_x_limit)then
Fne(i,j)=Fn(i-knex,j-kney)
else
a10=phi(j)
a11=landa(i)
call yin_to_yang(a10,a11,a12,a13)
call select_int_method(a12,a13,Fe,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
Fne(i,j)=hem
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  Filter_boundaries(Hen,Hee,Hen0,Hee0,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
implicit none
integer::my,mx,i,j,n_epsilon_x,n_epsilon_y
real*8::pi,top_x_limit,top_y_limit,dphi,dlanda,Hen(mx,my),Hee(mx,my),a10,a11,a12,a13,hem&
,Hen0(mx,my),Hee0(mx,my)
character(5)::int_method
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)

do j=1,my
do i=1,mx
if(i==1.or.i==mx.or.j==1.or.j==my)then
a10=-top_y_limit+dble(j-1)*dphi
a11=-top_x_limit+dble(i-1)*dlanda
call yin_to_yang(a10,a11,a12,a13)
call select_int_method(a12,a13,Hee0,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
Hen(i,j)=hem
a12=-top_y_limit+dble(j-1)*dphi
a13=-top_x_limit+dble(i-1)*dlanda
call yang_to_yin(a10,a11,a12,a13)
call select_int_method(a10,a11,Hen0,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
Hee(i,j)=hem

else
Hen(i,j)=Hen(i,j)
Hee(i,j)=Hee(i,j)
endif
enddo
enddo

return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  select_int_method(phi3,landa3,heb,nx,ny,hm,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
implicit none
integer::nx,ny,n_epsilon_x,n_epsilon_y
real*8::phi3,landa3,dphi,dlanda,hm,heb(nx,ny)
character(5)::int_method
if(int_method=='bilin')then
call bilinear(phi3,landa3,heb,nx,ny,hm,dlanda,dphi,n_epsilon_x,n_epsilon_y)
elseif(int_method=='bicub')then
call bicubic(phi3,landa3,heb,nx,ny,hm,dlanda,dphi,n_epsilon_x,n_epsilon_y)
endif
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine bicubic(phi3,landa3,heb,nx,ny,hm,dlanda,dphi,n_epsilon_x,n_epsilon_y)
implicit none
integer::nx,ny,n_epsilon_x,n_epsilon_y
integer::k1,k2
real*8::pi,landa3,phi3,dphi,dlanda,hm,hm1,hm2,phi,landa
real*8::heb(nx,ny),hm3,hm4,h1,h2,h3,h4

pi=2.d0*dasin(1.0d0)

k1=int(((landa3+(3.0d0*pi/4.0d0)+1.0d-10)/dlanda))+1+n_epsilon_x
k2=int(((phi3+(pi/4.0d0)+1.0d-10)/dphi))+1+n_epsilon_y
phi=(dble(k2-1-n_epsilon_y))*dphi-(pi/4.0d0)
landa=dble(k1-1-n_epsilon_x)*dlanda-3.0d0*pi/4.0d0

h1=heb(k1-1,k2-1)
h2=heb(k1,k2-1)
h3=heb(k1+1,k2-1)
h4=heb(k1+2,k2-1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm1)
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(k1+2,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm2)
h1=heb(k1-1,k2+1)
h2=heb(k1,k2+1)
h3=heb(k1+1,k2+1)
h4=heb(k1+2,k2+1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm3)
h1=heb(k1-1,k2+2)
h2=heb(k1,k2+2)
h3=heb(k1+1,k2+2)
h4=heb(k1+2,k2+2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm4)


call cubic(hm1,hm2,hm3,hm4,phi-dphi,phi,phi+dphi,phi+2.0d0*dphi,phi3,hm)
return
end subroutine 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine cubic(h1,h2,h3,h4,x1,x2,x3,x4,xm,hm)
implicit none
real*8::h1,h2,h3,h4,x1,x2,x3,x4,xm,hm
real*8::a,b,c,d,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11
c1=(h4-h1)/(x4-x1)
c2=(h2-h1)/(x2-x1)
c3=((h3-h1)*(x4-x2))/((x3-x1)*(x3-x2))
c4=((h2-h1)*(x4-x2))/((x2-x1)*(x3-x2))
c5=(x4*x4)-(x2*x4)-(x3*x4)+(x2*x3)
a=(c1/c5)-(c2/c5)-(c3/c5)+(c4/c5)
c6=(h3-h1)/((x3-x1)*(x3-x2))
c7=(h2-h1)/((x2-x1)*(x3-x2))
c8=x1+x2+x3
b=c6-c7-(c8*a)
c9=(h2-h1)/(x2-x1)
c10=(x2*x2)+(x1*x1)+(x2*x1)
c11=x1+x2
c=c9-(c10*a)-(b*c11)
d=h1-(a*x1*x1*x1)-(b*x1*x1)-(c*x1)
hm=(a*xm*xm*xm)+(b*xm*xm)+(c*xm)+d
return
end subroutine 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  yin_to_yang(phin,landan,phie,landae)
real*8::phin,landan,phie,landae,pi
pi=dasin(1.0d0)*2.0d0
if(phin==-(pi/2.0d0))then
phie=0.0d0
landae=pi/2.0d0
elseif(phin==(pi/2.0d0))then
phie=0.0d0
landae=-pi/2.0d0
elseif(phin==0.0d0.and.landan==(pi/2.0d0))then
phie=-pi/2.0d0
landae=0.0d0
elseif(phin==0.0d0.and.landan==-(pi/2.0d0))then
phie=pi/2.0d0
landae=0.0d0
elseif(phin==0.0d0.and.landan==0.0d0)then
phie=0.0d0
landae=-pi
elseif(phin==0.0d0.and.landan==-pi)then
phie=0.0d0
landae=0.0d0
elseif(phin==0.0d0.and.landan>-pi.and.landan<-(pi/2.0d0))then
phie=landan+pi
landae=0.0d0
elseif(phin==0.0d0.and.landan>-(pi/2.0d0).and.landan<0.0d0)then
phie=-landan
landae=-pi
elseif(phin==0.0d0.and.landan>0.0d0.and.landan<(pi/2.0d0))then
phie=-landan
landae=-pi
elseif(phin==0.0d0.and.landan>(pi/2.0d0))then
phie=landan-pi
landae=0.0d0

elseif(landan==0.0d0.and.phin>0.0d0)then
phie=0.0d0
landae=phin-pi
elseif(landan==0.0d0.and.phin<0.0d0)then
phie=0.0d0
landae=phin+pi
elseif(landan==-pi.and.phin==0.0d0)then
phie=0.0d0
landae=0.0d0
elseif(landan==-pi.and.phin>0.0d0)then
phie=0.0d0
landae=-phin
elseif(landan==-pi.and.phin<0.0d0)then
phie=0.0d0
landae=-phin
elseif(phin>0.0d0.and.landan==-pi/2.0d0)then
phie=dasin(dcos(phin))
landae=-pi/2.0d0
elseif(phin>0.0d0.and.landan==pi/2.0d0)then
phie=-dasin(dcos(phin))
landae=-pi/2.0d0
elseif(phin<0.0d0.and.landan==-pi/2.0d0)then
phie=dasin(dcos(phin))
landae=pi/2.0d0
elseif(phin<0.0d0.and.landan==pi/2.0d0)then
phie=-dasin(dcos(phin))
landae=pi/2.0d0
elseif(phin>0.0d0.and.landan>-pi.and.landan<-(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))
elseif(phin>0.0d0.and.landan>-(pi/2.0d0).and.landan<0.0d0)then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))-pi
elseif(phin>0.0d0.and.landan>0.0d0.and.landan<(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))-pi
elseif(phin>0.0d0.and.landan>(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))

elseif(phin<0.0d0.and.landan>-pi.and.landan<-(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))
elseif(phin<0.0d0.and.landan>-(pi/2.0d0).and.landan<0.0d0)then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))+pi
elseif(phin<0.0d0.and.landan>0.0d0.and.landan<(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))+pi
elseif(phin<0.0d0.and.landan>(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))
endif
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  yang_to_yin(phin,landan,phie,landae)
real*8::phin,landan,phie,landae,pi
pi=dasin(1.0d0)*2.0d0
if(phie==-(pi/2.0d0))then
phin=0.0d0
landan=pi/2.0d0
elseif(phie==(pi/2.0d0))then
phin=0.0d0
landan=-pi/2.0d0
elseif(phie==0.0d0.and.landae==(pi/2.0d0))then
phin=-pi/2.0d0
landan=0.0d0
elseif(phie==0.0d0.and.landae==-(pi/2.0d0))then
phin=pi/2.0d0
landan=0.0d0
elseif(phie==0.0d0.and.landae==0.0d0)then
phin=0.0d0
landan=-pi
elseif(phie==0.0d0.and.landae==-pi)then
phin=0.0d0
landan=0.0d0
elseif(phie==0.0d0.and.landae>-pi.and.landae<-(pi/2.0d0))then
phin=landae+pi
landan=0.0d0
elseif(phie==0.0d0.and.landae>-(pi/2.0d0).and.landae<0.0d0)then
phin=-landae
landan=-pi
elseif(phie==0.0d0.and.landae>0.0d0.and.landae<(pi/2.0d0))then
phin=-landae
landan=-pi
elseif(phie==0.0d0.and.landae>(pi/2.0d0))then
phin=landae-pi
landan=0.0d0
elseif(landae==0.0d0.and.phie>0.0d0)then
phin=0.0d0
landan=phie-pi
elseif(landae==0.0d0.and.phie<0.0d0)then
phin=0.0d0
landan=phie+pi
elseif(landae==-pi.and.phie==0.0d0)then
phin=0.0d0
landan=0.0d0
elseif(landae==-pi.and.phie>0.0d0)then
phin=0.0d0
landan=-phie
elseif(landae==-pi.and.phie<0.0d0)then
phin=0.0d0
landan=-phie
elseif(phie>0.0d0.and.landae==-pi/2.0d0)then
phin=dasin(dcos(phie))
landan=-pi/2.0d0
elseif(phie>0.0d0.and.landae==pi/2.0d0)then
phin=-dasin(dcos(phie))
landan=-pi/2.0d0
elseif(phie<0.0d0.and.landae==-pi/2.0d0)then
phin=dasin(dcos(phie))
landan=pi/2.0d0
elseif(phie<0.0d0.and.landae==pi/2.0d0)then
phin=-dasin(dcos(phie))
landan=pi/2.0d0




elseif(phie>0.0d0.and.landae>-pi.and.landae<-(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))

elseif(phie>0.0d0.and.landae>-(pi/2.0d0).and.landae<0.0d0)then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))-pi

elseif(phie>0.0d0.and.landae>0.0d0.and.landae<(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))-pi

elseif(phie>0.0d0.and.landae>(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))

elseif(phie<0.0d0.and.landae>-pi.and.landae<-(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))

elseif(phie<0.0d0.and.landae>-(pi/2.0d0).and.landae<0.0d0)then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))+pi

elseif(phie<0.0d0.and.landae>0.0d0.and.landae<(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))+pi

elseif(phie<0.0d0.and.landae>(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))
endif
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine bilinear(phi3,landa3,heb,nx,ny,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y)
implicit none
integer::nx,ny,n_epsilon_x,n_epsilon_y
integer::k1,k2
real*8::pi,landa3,phi3,dphi,dlanda,hem,dl,dp,phi,landa,hem1,hem2
real*8::heb(nx,ny),top_x_limit,top_y_limit


pi=2.d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)


k1=int(((landa3+top_x_limit+1.0d-10)/dlanda))+1
k2=int(((phi3+top_y_limit+1.0d-10)/dphi))+1

landa=dble(k1-1)*dlanda-top_x_limit
phi=(dble(k2-1))*dphi-top_y_limit
dl=landa3-landa
dp=phi3-phi
call linear(heb(k1,k2),heb(k1+1,k2),dlanda,dl,hem1)
call linear(heb(k1,k2+1),heb(k1+1,k2+1),dlanda,dl,hem2)
call linear(hem1,hem2,dphi,dp,hem)


return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine linear(h1,h2,d,dd,hem)
real*8::hem,h1,h2,d,dd
hem=((h1*(d-dd))+(h2*dd))/d
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  provide_Yin_scalar_to_derive(hen,hee,dlanda,dphi,nx,ny,mx,my&
,n_epsilon_x,n_epsilon_y,he_zone,he_meridion,int_method)
implicit none
integer::i,j,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,k1,k2,knex,kney
real*8::hen(mx,my),hee(mx,my),he_zone(nx,my),he_meridion(mx,2*(ny-1)),a10,a11,a12,a13,hem,pi,dlanda,dphi&
,top_x_limit,top_y_limit
character(5)::int_method
kney=(ny-my)/2
knex=(nx-mx+1)/2
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)
do j=1,my
a10=-top_y_limit+dble(j-1)*dphi
do i=1,nx
a11=-pi+dble(i-1)*dlanda
if(a11>=-top_x_limit.and.a11<=top_x_limit)then
he_zone(i,j)=hen(i-knex,j)
else

call yin_to_yang(a10,a11,a12,a13)
call select_int_method(a12,a13,hee,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
he_zone(i,j)=hem
endif
enddo
enddo
do i=1,mx
a11=-top_x_limit+dble(i-1)*dlanda
do j=1,ny
a10=-(pi/2.0d0)+dble(j-1)*dphi
if(a10>=-top_y_limit.and.a10<=top_y_limit)then
he_meridion(i,j)=hen(i,j-kney)
else

call yin_to_yang(a10,a11,a12,a13)
call select_int_method(a12,a13,hee,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
he_meridion(i,j)=hem
endif
enddo
do j=ny+1,2*(ny-1)
a10=(pi/2.0d0)-dble(j-ny)*dphi
if(a11<0.0d0)then
a11=a11+pi
elseif(a11>=0.0d0)then
a11=a11-pi
endif
if(a10>=-top_y_limit.and.a10<=top_y_limit.and.a11>=-top_x_limit.and.a11<=top_x_limit)then
call num(a10,a11,k1,k2,dlanda,dphi,n_epsilon_x,n_epsilon_y)
he_meridion(i,j)=hen(k1,k2)
else
call yin_to_yang(a10,a11,a12,a13)
call select_int_method(a12,a13,hee,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
he_meridion(i,j)=hem
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  provide_Yang_scalar_to_derive(hen,hee,dlanda,dphi,nx,ny,mx,my&
,n_epsilon_x,n_epsilon_y,he_zone,he_meridion,int_method)
implicit none
integer::i,j,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,k1,k2,knex,kney
real*8::hen(mx,my),hee(mx,my),he_zone(nx,my),he_meridion(mx,2*(ny-1)),a10,a11,a12,a13,hem,pi,dlanda,dphi&
,top_x_limit,top_y_limit
character(5)::int_method
kney=(ny-my)/2
knex=(nx-mx+1)/2
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)
do j=1,my
a12=-top_y_limit+dble(j-1)*dphi
do i=1,nx
a13=-pi+dble(i-1)*dlanda
if(a13>=-top_x_limit.and.a13<=top_x_limit)then
he_zone(i,j)=hee(i-knex,j)
else

call yang_to_yin(a10,a11,a12,a13)
call select_int_method(a10,a11,hen,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
he_zone(i,j)=hem
endif
enddo
enddo


do i=1,mx
a13=-top_x_limit+dble(i-1)*dlanda
do j=1,ny
a12=-(pi/2.0d0)+dble(j-1)*dphi
if(a12>=-top_y_limit.and.a12<=top_y_limit)then
he_meridion(i,j)=hee(i,j-kney)
else

call yang_to_yin(a10,a11,a12,a13)
call select_int_method(a10,a11,hen,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
he_meridion(i,j)=hem
endif
enddo
do j=ny+1,2*(ny-1)
a12=(pi/2.0d0)-dble(j-ny)*dphi
if(a13<0.0d0)then
a13=a13+pi
elseif(a13>=0.0d0)then
a13=a13-pi
endif
if(a12>=-top_y_limit.and.a12<=top_y_limit.and.a13>=-top_x_limit.and.a13<=top_x_limit)then
call num(a12,a13,k1,k2,dlanda,dphi,n_epsilon_x,n_epsilon_y)
he_meridion(i,j)=hee(k1,k2)
else
call yang_to_yin(a10,a11,a12,a13)
call select_int_method(a10,a11,hen,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
he_meridion(i,j)=hem
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine initial_condition_yin(mx,my,cnx,cny,hen1,alpha,n_epsilon_x,n_epsilon_y,dlanda,dphi)
implicit none
integer::i,j,mx,my,n_epsilon_x,n_epsilon_y
real*8::aa,rr,rrr,pi,twopi,u00,alpha,tetac,landac,dlanda,dphi,h00,top_x_limit,top_y_limit&
,ln(mx),pn(my),hen1(mx,my),cnx(mx,my),cny(mx,my)
aa=6.37122d6;rrr=aa/3.0d0
pi=2.0d0*dasin(1.0d0)
twopi=2.0d0*pi
u00=twopi*aa/(12.0d0*86400.0d0)
h00=1000.0d0
tetac=0.0d0
landac=3.0d0*pi/2.0d0
top_x_limit=(pi*3.0d0/4.0d0)+dble(n_epsilon_x)*dlanda
top_y_limit=(pi/4.0d0)+dble(n_epsilon_y)*dphi
do i=1,mx
ln(i)=-top_x_limit+dble(i-1)*dlanda
enddo
do j=1,my
pn(j)=dble(j-1)*dphi-top_y_limit
enddo
do j=1,my
do i=1,mx
rr=aa*dacos(dsin(pn(j))*dsin(tetac)+dcos(tetac)*dcos(pn(j))*dcos((ln(i)+pi)-landac))
if(rr<rrr)then
hen1(i,j)=(h00/2.0d0)*(1.0d0+dcos(pi*rr/rrr))
else
hen1(i,j)=0.0d0
endif
cnx(i,j)=u00*(dcos(pn(j))*dcos(alpha)+dcos(ln(i)+pi)*dsin(pn(j))*dsin(alpha))
cny(i,j)=-u00*dsin(ln(i)+pi)*dsin(alpha)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine initial_condition_yang(mx,my,cex,cey,hee1,alpha,n_epsilon_x,n_epsilon_y,dlanda,dphi)
implicit none
integer::i,j,mx,my,n_epsilon_x,n_epsilon_y
real*8::aa,rr,rrr,pi,twopi,u00,alpha,tetac,landac,dlanda,dphi,h00,top_x_limit,top_y_limit&
,le(mx),pe(my),hee1(mx,my),cex(mx,my),cey(mx,my),a10,a11,a12,a13,a14,a15,a16,a17
aa=6.37122d6;rrr=aa/3.0d0
pi=2.0d0*dasin(1.0d0)
twopi=2.0d0*pi
u00=twopi*aa/(12.0d0*86400.0d0)
h00=1000.0d0
tetac=0.0d0
landac=3.0d0*pi/2.0d0
top_x_limit=(pi*3.0d0/4.0d0)+dble(n_epsilon_x)*dlanda
top_y_limit=(pi/4.0d0)+dble(n_epsilon_y)*dphi
do i=1,mx
le(i)=-top_x_limit+dble(i-1)*dlanda
enddo
do j=1,my
pe(j)=dble(j-1)*dphi-top_y_limit
enddo
do i=1,mx
a13=le(i)
do j=1,my
a12=pe(j)
call yang_to_yin(a10,a11,a12,a13)
rr=aa*dacos(dsin(a10)*dsin(tetac)+dcos(tetac)*dcos(a10)*dcos((a11+pi)-landac))
if(rr<rrr)then
hee1(i,j)=(h00/2.0d0)*(1.0d0+dcos(pi*rr/rrr))
else
hee1(i,j)=0.0d0
endif
a15=u00*(dcos(a10)*dcos(alpha)+dcos(a11+pi)*dsin(a10)*dsin(alpha))
a14=-u00*dsin(a11+pi)*dsin(alpha)
call yinvector_to_yangvector(a10,a11,a12,a13,a14,a15,a16,a17)
cey(i,j)=a16;cex(i,j)=a17
write(20,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,hee1(i,j)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine initial_condition_glob(nx,ny,hene,dlanda,dphi)
implicit none
integer::i,j,nx,ny
real*8::aa,rr,rrr,pi,twopi,u00,tetac,landac,dlanda,dphi,h00&
,landa(nx),phi(ny),hene(nx,ny)
aa=6.37122d6;rrr=aa/3.0d0
pi=2.0d0*dasin(1.0d0)
twopi=2.0d0*pi
u00=twopi*aa/(12.0d0*86400.0d0)
h00=1000.0d0
tetac=0.0d0
landac=3.0d0*pi/2.0d0
do i=1,nx
landa(i)=-pi+dble(i-1)*dlanda
enddo
do j=1,ny
phi(j)=dble(j-1)*dphi-(pi/2.0d0)
enddo
do j=1,ny
do i=1,nx
rr=aa*dacos(dsin(phi(j))*dsin(tetac)+dcos(tetac)*dcos(phi(j))*dcos((landa(i)+pi)-landac))
if(rr<rrr)then
hene(i,j)=(h00/2.0d0)*(1.0d0+dcos(pi*rr/rrr))
else
hene(i,j)=0.0d0
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************