module func
use constant4pick
use typedef
implicit none
contains

function cal_dx(x1,x2)
	implicit none
	integer::x1,x2,cal_dx
	cal_dx=x2-x1
	if(cal_dx>x_num/2)then
		cal_dx=cal_dx-x_num
	else if(cal_dx<-x_num/2)then
		cal_dx=cal_dx+x_num
	end if
	return 
end function

function distance(x1,y1,x2,y2)
	implicit none
	integer::x1,x2,y1,y2
	real::dx,dy,distance
	dx=cal_dx(x1,x2)
	dy=y2-y1
	distance=sqrt(dx**2+dy**2)
	return
end function

subroutine coord_conv(x1,y1,x2,y2)
integer,intent(inout)::x1,y1,x2,y2
x2=nint((x0_1+(x1-1)*ddx1-x0_2)/ddx2+1)
y2=nint((y0_1+(y1-1)*ddx1-y0_2)/ddy2+1)
return
end subroutine

subroutine event_cal(p_st,p_end,p,flag)
implicit none
integer::i,j,k,x_max,x_min,y_max,y_min,cnt,dx,dy,tmpx,tmpy,x0,y0,x,y,dxx,dyy
integer::flag,step !flag==1 the open line is to the left, flag==2 open line to the right
real::x_sum,y_sum,cen_x,cen_y,left_tq,right_tq,tq_min,dll
real,allocatable::x_bk(:),y_bk(:)
type(isoline),pointer::p_st,p_end,p_iso,pr,pr_st,pr_end
type(event),pointer::p
type(point_record)::bndy(-x_num:x_num),bndx(-y_num:y_num)
integer::px,py
real::plc,plz,parea,plcl,plzl,pareal
!unt=pi*a/72 !unt is the unit distance of a degree of latitude
!===set the value pr pointing at, which uses relative cordinate, the initial point is (0,0)=== 
x0=p_st%x
y0=p_st%y
allocate(pr_st)
pr=>pr_st
pr%x=0
pr%y=0
pr%iso_type=p_st%iso_type
p_iso=>p_st
!=== 5th version ===================
!x=((x0_1+p_iso%x-1)*ddx1-x0_2)/ddx2+1
!y=((y0_1+p_iso%y-1)*ddx1-x0_2)/ddy2+1
call coord_conv(p_iso%x,p_iso%y,x,y)
!================================

if(p%ac_type==2.and.flag==2) p%coverage(x,y)=2
dyy=0
dll=0
do 
	dx=cal_dx(p_iso%x,p_iso%next%x)
	dy=p_iso%next%y-p_iso%y
	allocate(pr%next)
	pr%next%prev=>pr
	pr=>pr%next
	pr%x=pr%prev%x+dx
	pr%y=pr%prev%y+dy
!	print*,pr%x,pr%y
	pr%iso_type=p_iso%next%iso_type
	p_iso=>p_iso%next
!x=(p_iso%x-1)/ddx2+1
!y=(p_iso%y-1)/ddy2+1

	!=== 5th version ===================
	!x=((x0_1+p_iso%x-1)*ddx1-x0_2)/ddx2+1
	!y=((y0_1+p_iso%y-1)*ddx1-x0_2)/ddy2+1
	call coord_conv(p_iso%x,p_iso%y,x,y)
	!===========================
	if(p%ac_type==flag)then
		if(dx<=0)then
			dll=dll+sqrt((dx*1.0)**2+(dy*1.0)**2)
		end if
		if(p%ac_type==2.and.dx<=0)then
			p%coverage(x,y)=2
			dyy=dyy+dy
		end if
		if(p%ac_type==1.and.dx>=0.and.dy>=0) p%coverage(x,y)=2
	end if
	if(associated(p_iso,p_end))then
		exit
	end if
end do
if(dll>1.0) p%ll=dll
if(p%ac_type==2.and.flag==2.and.dyy<=0)then
	where(p%coverage==2) p%coverage=0
end if

!print*,pr%x,pr%y
pr%next=>null()
pr_st%prev=>null()
pr_end=>pr


!===calculate area===
!===1.calculate bndy===
bndy(:)%min=y_num+1
bndy(:)%max=-999
y_min=y_num+1
y_max=-999
pr=>pr_st
do 
	if(pr%y<bndy(pr%x)%min)then
		bndy(pr%x)%min=pr%y	
		if(pr%y<y_min)then
			y_min=pr%y
		end if
	end if
	if(pr%y>bndy(pr%x)%max)then
		bndy(pr%x)%max=pr%y
		if(pr%y>y_max)then
			y_max=pr%y
		end if
	end if
	if(associated(pr,pr_end))then
		exit
	end if
	pr=>pr%next
end do

!===2.calculate bndx
bndx(:)%max=-999
bndx(:)%min=x_num+1
x_max=-999
x_min=x_num+1

if(pr_end%y>pr_st%y)then
	step=1
else
	step=-1
end if
if(flag==1)then
	do i=pr_st%y,pr_end%y,step
		bndx(i)%min=pr_st%x
		x_min=pr_st%x
	end do
else if(flag==2)then
	do i=pr_st%y,pr_end%y,step
		bndx(i)%max=pr_st%x
		x_max=pr_st%x
	end do
end if

pr=>pr_st
do 
	if(pr%x>bndx(pr%y)%max)then
		bndx(pr%y)%max=pr%x
		if(pr%x>x_max)then
			x_max=pr%x
		end if
	end if
	if(pr%x<bndx(pr%y)%min)then
		bndx(pr%y)%min=pr%x
		if(pr%x<x_min)then
			x_min=pr%x
		end if
	end if
	if(associated(pr,pr_end))then
		exit
	end if
	pr=>pr%next
end do

!===3.provided bndx%bndy,calculate area===
parea=0
pareal=0
!if(any(p%coverage>50)==.true.) print*,'what a mess'
do i=x_min,x_max
	do j=y_min,y_max
		if(i>=bndx(j)%min.and.i<=bndx(j)%max.and.j>=bndy(i)%min.and.j<=bndy(i)%max)then
			!p%area1=p%area1+1
			parea=parea+1
			x=x0+i
			y=y0+j
			!p%areal1=p%areal1+cos((y-90)*1.0*pi/180.0)*unt**2
			pareal=pareal+cos((y-90)*1.0*pi/180.0)*unt**2
			!=== this line corrects the area ===
			if(x>x_num) x=x-x_num
			if(x<=0) x=x+x_num
!			x=(x-1)*ddx1/ddx2+1
!			y=(y-1)*ddy1/ddy2+1
	!		x=(x-1)/2.5+1
	!		y=(y-1)/2.5+1
			!=== 5th version ===================
			!x=((x0_1+x-1)*ddx1-x0_2)/ddx2+1
			!y=((y0_1+y-1)*ddx1-x0_2)/ddy2+1
			call coord_conv(x,y,x,y)
			!===================================
		!	print*,'...',x,y
!			if(y<37)then
!				print*,'occur in southern hemisphere'
!				print*,'p_st:',p_st%x,p_st%y
!			end if
			!print*,x,y
			!=== 5th version ====
			!if(p%ac_type==flag.and.p%coverage(x,y)/=2)then
			if(p%coverage(x,y)/=2)then
				p%coverage(x,y)=1
			end if
			!=======================
		end if
	end do
end do
!do i=1,nx
!	do j=1,ny
!		if(p%coverage(i,j)>50)then
!			print*,'value larger than normal'
!			print*,'x,y',i,j
!			print*,'value:',p%coverage(i,j)
!			print*,'p_st:',p_st%x,p_st%y
!			print*,'       '
!		end if
!	end do
!end do

!===use bndx&bndy to calculate the centroid===
allocate(x_bk(x_min:x_max))
allocate(y_bk(y_min:y_max))
x_bk=0
y_bk=0
do i=x_min,x_max
	do k=x_min,i-1
		x_bk(i)=x_bk(i)+abs(k-i)*abs(bndy(k)%max-bndy(k)%min)
	end do
	do k=i+1,x_max
		x_bk(i)=x_bk(i)-abs(k-i)*abs(bndy(k)%max-bndy(k)%min)
	end do
	x_bk(i)=abs(x_bk(i))
end do
tq_min=x_bk(x_min)
tmpx=x_min
do i=x_min+1,x_max
	if(x_bk(i)<tq_min)then
		tq_min=x_bk(i)
		tmpx=i
	end if
end do
cen_x=tmpx
do j=y_min,y_max
	do k=y_min,j-1
		y_bk(j)=y_bk(j)+abs(k-j)*abs(bndx(k)%max-bndx(k)%min)
	end do
	do k=j+1,y_max
		y_bk(j)=y_bk(j)-abs(k-j)*abs(bndx(k)%max-bndx(k)%min)
	end do
	y_bk(j)=abs(y_bk(j))
end do
tq_min=y_bk(y_min)
tmpy=y_min
do j=y_min+1,y_max
	if(y_bk(j)<tq_min)then
		tq_min=y_bk(j)
		tmpy=j
	end if
end do
cen_y=tmpy
px=cen_x+p_st%x
py=cen_y+p_st%y
!print*,p%x,p%y,x_min,x_max
if(px>x_num)then
	px=px-x_num
else if(px<=0)then
	px=px+x_num
end if
!print*,'centroid:',p%x,p%y,p%lc,p%lz
!=== calculate lc lz using area and centroid ===
plc=sqrt(parea/pi)
plz=(x_max-x_min+1)/2.0
plcl=sqrt(pareal/pi)
plzl=(x_max-x_min+1)*unt*cos((py-90)*1.0*pi/180.0)/2.0
!===============================================
if((x_max-x_min)/2.0>179)then
	print*,'wrong lz'
	print*,x_min,x_max,p%val
	print*,p_st%x,p_st%y,p_end%x,p_end%y
	print*,p%x1,p%y1,p%ac_type,p%area1
	stop
end if

!=== following is relatively new ===============
!=== to determine the pe_type ==================
if(flag==2)then
	p_iso=>p_end
	pr=>pr_end
	do 
		dx=cal_dx(p_iso%x,p_iso%next%x)
		dy=p_iso%next%y-p_iso%y
		allocate(pr%next)
		pr%next%prev=>pr
		pr=>pr%next
		pr%x=pr%prev%x+dx
		pr%y=pr%prev%y+dy
		pr%iso_type=p_iso%next%iso_type
		p_iso=>p_iso%next
		!x=(p_iso%x-1)/ddx2+1
		!y=(p_iso%y-1)/ddy2+1
		!=== 5th version ===================
		!x=((x0_1+p_iso%x-1)*ddx1-x0_2)/ddx2+1
		!y=((y0_1+p_iso%y-1)*ddx1-x0_2)/ddy2+1
		call coord_conv(p_iso%x,p_iso%y,x,y)
		!===================================
		!if(pr%x-pr_end%x<20)p%coverage(x,y)=2
		!cmt: to make capture less sloppy 
		if(dx<0)then
			p%right=0
			pr_end=>pr
			exit
		end if
		if(pr%y==pr_st%y)then
			p%right=pr%x-pr_st%x
			pr_end=>pr
			exit
		end if
	end do
end if

if(flag==1)then
	p_iso=>p_st
	pr=>pr_st
	!	print*,'mark'
	dyy=0
	do 
		dx=cal_dx(p_iso%x,p_iso%prev%x)
		dy=p_iso%prev%y-p_iso%y
		allocate(pr%prev)
		pr%prev%next=>pr
		pr=>pr%prev
		pr%x=pr%next%x+dx
		pr%y=pr%next%y+dy
	!	print*,'2',pr%x,pr%y
		pr%iso_type=p_iso%prev%iso_type
		p_iso=>p_iso%prev
		!x=(p_iso%x-1)/ddx2+1
		!y=(p_iso%y-1)/ddy2+1
		!=== 5th version ===================
		!x=((x0_1+p_iso%x-1)*ddx1-x0_2)/ddx2+1
		!y=((y0_1+p_iso%y-1)*ddx1-x0_2)/ddy2+1
		call coord_conv(p_iso%x,p_iso%y,x,y)
		!===================================
		!if(abs(pr%x-pr_st%x)<20)p%coverage(x,y)=2
		if(p%ac_type==1.and.abs(pr%x-pr_st%x)<2*plc)then
			p%coverage(x,y)=2
			dyy=dyy+dy
		end if
		!cmt: to make capture less sloppy
		!cmt: These two changes make the extension only be applied to awb left
		if(dx>0)then
			p%left=0
			pr_st=>pr
			exit
		end if
		if(pr%y==pr_end%y)then
			p%left=pr_end%x-pr%x
			pr_st=>pr
			exit
		end if
	end do
	if(p%ac_type==1.and.dyy>0)then
		where(p%coverage==2) p%coverage=0
	end if
end if

if(flag/=p%ac_type)then
p%pe_type=0
if(p%left/=0.and.p%right/=0)then
	if(p%left>p%right)then
		if(p%ac_type==1)then
			p%pe_type=2
		else
			p%pe_type=1
		end if
	else
		if(p%ac_type==1)then
			p%pe_type=1
		else
			p%pe_type=2
		end if
	end if
else
	if(p%left==0.and.p%right/=0)then
		if(p%ac_type==1)then
			p%pe_type=2
		else
			p%pe_type=1
		end if
	else if(p%right==0.and.p%left/=0)then
		if(p%ac_type==1)then
			p%pe_type=1
		else
			p%pe_type=2
		end if
	end if
end if
end if
!=== assign the value ==========================
if(flag==p%ac_type)then
	p%lz=plz
	p%lzl=plzl
	p%x1=px
	p%y1=py
	p%lc1=plc
	p%lcl1=plcl
	p%area1=parea
	p%areal1=pareal
else
	p%x2=px
	p%y2=py
	p%lc2=plc
	p%lcl2=plcl
	p%area2=parea
	p%areal2=pareal
	if(p%lz>0.and.p%lz/=plz)then
		!print*,plz,p%lz
	end if
end if
!=== deallocate the contour objects ============
pr=>pr_st
do
	pr=>pr%next
	deallocate(pr%prev)
	if(associated(pr,pr_end))exit
end do
deallocate(pr)
!print*,'mark2'
end subroutine

subroutine gen_grid(pv,grid,pv_val)
implicit none
real::pv(x_num,y_num),pv_val(1)
integer::grid(x_num,y_num),i,j,m,n,nx,ny,ex,ey,wx,wy,sx,sy,dx,dy,tmpx(0:7),tmpy(0:7),flag=0,r
integer::part_num,part(4),k
do i=1,x_num
	do j=3,y_num-2
		if(pv_val(1)<=pv(i,j))then  !If the value of the point is greater than pv_val(1),select it as a candidate.
			nx=i
			sx=i
			ey=j
			wy=j
			ny=j+1
			sy=j-1
			ex=i+1
			wx=i-1
			if(ex>x_num)then
				ex=ex-x_num
			end if
			if(wx<1)then
				wx=wx+x_num
			end if
			if(pv(nx,ny)<pv_val(1).or.pv(sx,sy)<pv_val(1).or.pv(ex,ey)<pv_val(1).or.pv(wx,wy)<pv_val(1))then  
			!And if there is at least one point imediately beside the candidate has value less than pv_val(1),recognize it as a grid point
				grid(i,j)=1
			end if
		end if
	end do
	if(pv_val(1)<=pv(i,2))then
		nx=i
		ey=2
		wy=2
		ny=2+1
		ex=i+1
		wx=i-1
		if(ex>x_num)then
			ex=ex-x_num
		end if
		if(wx<1)then
			wx=wx+x_num
		end if
		if(pv(nx,ny)<pv_val(1).or.pv(ex,ey)<pv_val(1).or.pv(wx,wy)<pv_val(1))then
			grid(i,j)=1
		end if
	end if
	if(pv_val(1)<=pv(i,y_num-1))then
		sx=i
		ey=y_num-1
		wy=y_num-1
		sy=y_num-1-1
		ex=i+1
		wx=i-1
		if(ex>x_num)then
			ex=ex-x_num
		end if
		if(wx<1)then
			wx=wx+x_num
		end if
		if(pv(sx,sy)<pv_val(1).or.pv(ex,ey)<pv_val(1).or.pv(wx,wy)<pv_val(1))then
			grid(i,j)=1
		end if
	end if
end do
do i=1,x_num
	do j=2,y_num-1
		if(grid(i,j)==0)then
			if(pv(i,j)>=pv_val(1))then
				grid(i,j)=2
			end if
		end if		
	end do
	grid(i,1)=-1
	grid(i,y_num)=-1 !这两行为后面的处理提供了方便
end do
do i=1,x_num
	do j=2,y_num-1
		if(grid(i,j)==1)then
			call get_cord(i,j,tmpx,tmpy)
			flag=0
			do k=0,7
!				print*,'clockwise:',tmpx(k),tmpy(k)			
				if(grid(tmpx(k),tmpy(k))==2)then
					flag=1
					exit
				end if
			end do
			if(flag==0)then
				grid(i,j)=0
			end if
		end if
	end do
end do
end subroutine

subroutine get_cord(x,y,en_x,en_y)
implicit none
integer::x,y,en_x(0:7),en_y(0:7),k,dx,dy
do k=0,7
	if(abs(cos(pi/4*k))>eps)then
		dx=cos(pi/4*k)/abs(cos(pi/4*k))
	else
		dx=cos(pi/4*k)
	end if
	if(abs(sin(pi/4*k))>eps)then
		dy=sin(pi/4*k)/abs(sin(pi/4*k))
	else
		dy=sin(pi/4*k)
	end if
	en_x(k)=x+dx
	en_y(k)=y+dy
	if(en_x(k)>x_num)then
		en_x(k)=en_x(k)-x_num
	else if(en_x(k)<1)then
		en_x(k)=en_x(k)+x_num
	end if
end do
end subroutine
end module func
