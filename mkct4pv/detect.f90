!子程序detect参数说明：
!1.grid是检测到的等值线格点场，是等值线的为1，大于等值线的为2，小于等值线的为0
!2.environ是一个标准空间场，每个点的整数值代表其周围8个格点中非等值线格点数
!3.trace_org是记录已经被接收的等值线格点的标准场，已接受的为1，未接受的为0
!4.p指向目前被接收的等值线变量的指针
!5.tmpx,tmpy传回去的检测结果（多个，通过要求的等值线格点），是数组
!6.cnt传回去的满足条件的格点数
!7.dr_mode 等值线走向类型，1为左高右低，2为左低右高。是调用函数指定的，即满足这种指定走向的格点才会被接收
!8.flag 在特殊情况下作特殊处理的参数
!  flag=-1 可以走prev格点 体现在if语句中满足.and.(i/=p%prev%x.or.j/=p%prev%y.or.flag==-1)的格点可以被接纳
!  flag=0 无视trace的限制，可以走之前走过的格点；也允许judge==0,即两边值相等的情况。
!  flag=1
!  flag=4 普通情况
!  flag=6 一种宽松方案，当judge/=dr_mode，judge=3即是因为等值线两边都是等值线格点无法判断才导致的时候，允许这样的点被接纳
!  flag=7 一种宽松方案，当judge==0，即两边相等的时候，允许。
!  flag=8 允许走之前画过的等值线格点
!  flag=9 无视trace的限制，和judge=3的限制

!  flag1=-1 允许走prev格点&&无视trace&&
!  flag1=0  无视trace
!  flag1=2  可以走已完成的等值线的点（grid=-1的点）
!  flag1=4  普通情况
!  flag2=6  允许judge==3的格点（两边都无法判断）
!  flag2=7  允许judge==0的格点（两边相等）
!  flag2=4  普通情况
!    flag2=9 allow previous grid point
!about direction: east is 0, and the number increase anticlockwisely.

module dtct
use typedef
use func

contains
!检测周围的等值线格点***************************************
subroutine detect(grid,environ,trace_org,p,tmpx,tmpy,cnt,dr_mode,flag1,flag2,stat1,stat2,p_node,mode,delta_x) !这里的stat1用来存放等值线方向，0--7
use constant
implicit none
type(contour),target::node
type(contour),pointer::p,p_node
integer::grid(x_num,y_num),environ(x_num,y_num),judgetmp(9),anything,flag,mode(x_num,y_num),left,right,delta_x
integer::tmpx(9),tmpy(9),cnt,flag1,flag2,i,j,k,temp,judge,judge_op,dr_mode,stat1(x_num,y_num),drct,stat2(x_num,y_num),m,n,en_x(0:7),en_y(0:7)
integer::gval(0:7),rcd(0:7),dlt(0:7),start,step,lo,hi,nb
integer::dx,dy,mini,idx(7)
integer::a,b,c,d
real::distance(0:7)
type(cord)::trace_org(x_num,y_num),trace(x_num,y_num)
if(dr_mode==1)then
    step=1
    lo=0
    hi=7
else if(dr_mode==2)then
    step=-1
    lo=7
    hi=0
else
    print*,'wrong dr_mode'
    stop
end if
!real::pv(x_num,y_num)
rcd=0
dlt=0
!print*,'flag1:',flag1
if(flag1==0.or.flag1==-1.or.flag1==2)then !flag=-1是允许走prev点的情况（到单等值线的尽头）
    do i=1,x_num
        do j=1,y_num
            trace(i,j)%x=0 !解除“不走回头路”的第二类限制
            trace(i,j)%y=0 !解除“不走回头路”的第二类限制
        end do
    end do
else
    do i=1,x_num
        do j=1,y_num
            trace(i,j)%x=trace_org(i,j)%x
            trace(i,j)%y=trace_org(i,j)%y
        end do
    end do
end if
call get_cord(p%x,p%y,en_x,en_y)
do i=0,7
    gval(i)=grid(en_x(i),en_y(i))
end do
if(stat1(p%x,p%y)<0)then
    start=0 !Need to be noticed
else
    start=stat1(p%x,p%y)
end if
cnt=0
idx=0
do k=lo,hi,step
    i=mod(start+k+7,8)
    if(en_x(i)==p%prev%x.and.en_y(i)==p%prev%y.and.flag2/=9)then
        cycle
    else
        nb=mod(i+step+8,8)
        if(gval(mod(i,8))==0.and.gval(nb)==2.and.stat1(en_x(i),en_y(i))/=i)then
            if(trace(en_x(i),en_y(i))%x/=en_x(nb).or.trace(en_x(i),en_y(i))%y/=en_y(nb))then
                cnt=cnt+1
                tmpx(cnt)=en_x(i)
                tmpy(cnt)=en_y(i)
                idx(cnt)=i
                rcd(i)=cnt
            end if
        end if
    end if
end do
if(cnt==0)then
    do k=0,7
        nb=mod(k+step+8,8)
        if(gval(k)==0.and.gval(nb)==2.and.((en_x(k)==p_node%x.and.en_y(k)==p_node%y)))then
            cnt=cnt+1
            tmpx(cnt)=en_x(k)
            tmpy(cnt)=en_y(k)
            rcd(k)=cnt
        end if 
    end do
end if
if(cnt>1)then
    do i=1,cnt
        nb=mod(idx(i)+step+8,8)
        dx=en_x(nb)-trace(p%x,p%y)%x
        if(dx>=(x_num-10))then
            dx=dx-x_num
        else if(dx<=-(x_num-10))then
            dx=dx+x_num
        end if
        dy=en_y(nb)-trace(p%x,p%y)%y
        distance(i)=sqrt(dx**2*1.0+dy**2*1.0)
    end do
    !print*,'before*************'
    !print*,tmpx(1),tmpy(1)
    !print*,tmpx(2),tmpy(2)
    !print*,'after*************'
    mini=1
    do i=1,cnt
        if(distance(i)<distance(mini))then
            mini=i
        end if
    end do
    tmpx(1)=tmpx(mini)
    tmpy(1)=tmpy(mini)
end if
if(cnt>1.and.p%x==p_node%x.and.p%y==p_node%y) flag1=8
end subroutine

subroutine gen_grid(pv,grid,wtf)
use func
real::pv(x_num,y_num),wtf
integer::grid(x_num,y_num),i,j,m,n,nx,ny,ex,ey,wx,wy,sx,sy,dx,dy,tmpx(0:7),tmpy(0:7),flag=0,r
integer::part_num,part(4)
do i=1,x_num
    do j=2,y_num-1
        if(pv(i,j)>wtf)then
            grid(i,j)=2
        else
            grid(i,j)=0
        end if
    end do
end do

where(abs(pv)>0.9*abs(undef)) grid=-2

do i=1,x_num
    do j=2,y_num-1
        if(grid(i,j)==1)then
            call get_cord(i,j,tmpx,tmpy)
            flag=0
            do k=0,7
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

subroutine gen_environ(grid,environ)
use func
implicit none
integer::grid(x_num,y_num),environ(x_num,y_num),i,j,m,n,sum
do i=2,x_num-1
    do j=2,y_num-1
        if(grid(i,j)==1)then
            sum=0
            do m=i-1,i+1
                do n=j-1,j+1
                    if(grid(m,n)/=1.and.(m/=i.or.n/=j))then
                        sum=sum+1
                    end if
                end do
            end do
            environ(i,j)=sum
        end if
    end do    
end do
do j=2,y_num-1
    if(grid(1,j)==1)then
        sum=0
        do m=1,2
            do n=j-1,j+1
                if(grid(m,n)/=1.and.(m/=1.or.n/=j))then
                    sum=sum+1
                end if
            end do
        end do
        do n=j-1,j+1
            if(grid(x_num,n)/=1)then
                sum=sum+1
            end if
        end do
        environ(1,j)=sum
    end if
    if(grid(x_num,j)==1)then
        sum=0
        do m=x_num-1,x_num
            do n=j-1,j+1
                if(grid(m,n)/=1.and.(m/=x_num.or.n/=j))then
                    sum=sum+1
                end if
            end do
        end do
        
        do n=j-1,j+1
            if(grid(1,n)/=1)then
                sum=sum+1
            end if
        end do
        environ(x_num,j)=sum
    end if
end do
end subroutine

function is_contour(grid,x,y)
integer::is_contour
integer::grid(x_num,y_num),x,y
integer::en_x(0:7),en_y(0:7)
integer::i
is_contour=0
call get_cord(x,y,en_x,en_y)
if(grid(x,y)==0)then
    do i=0,7,2
        if(grid(en_x(i),en_y(i))==2)then
            is_contour=1
            return
        end if
    end do
end if
end function

!given the point (x,y) and grid(x_num,y_num), get the patch including (x,y)
!using DFS
recursive subroutine gen_patch(grid,patch_map,x,y,val)
integer::x,y,val
integer::grid(x_num,y_num),patch_map(x_num,y_num)
integer::en_x(0:7),en_y(0:7)
patch_map(x,y)=val
call get_cord(x,y,en_x,en_y)
do i=0,7,2
    if(grid(en_x(i),en_y(i))==val.and.patch_map(en_x(i),en_y(i))/=val.and.y>y_num/2+2)then
        call gen_patch(grid,patch_map,en_x(i),en_y(i),val)
    end if
end do
!print*,'before return'
end subroutine

end module dtct


