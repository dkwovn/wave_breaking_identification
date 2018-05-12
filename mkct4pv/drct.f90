module func
use constant
integer::good

interface get_val
    module procedure get_val_real 
    module procedure get_val_int 
end interface
contains
subroutine llrr(x1,y1,x2,y2,leftx,lefty,rightx,righty,drct) !输入两个矢量格点（x1,y1,x2,y2），输出其左格点和右格点坐标(leftx,lefty,rightx,righty)
    integer::x1,y1,x2,y2,leftx,lefty,rightx,righty,drctx,drcty,drct
    drctx=x2-x1
    if(drctx<-100)then
        drctx=drctx+x_num
    else if(drctx>100)then
        drctx=drctx-x_num
    end if
    drcty=y2-y1 !设置方向
    if(drctx==0)then
        if(drcty>0)then
            drct=2
        else
            drct=6
        end if
    else if(drctx>0)then
        if(drcty>0)then
            drct=1
        else if(drcty<0)then
            drct=7
        else
            drct=0
        end if
    else 
        if(drcty>0)then
            drct=3
        else if(drcty<0)then
            drct=5
        else 
            drct=4
        end if
    end if
!对新格点进行方向检验
    if(abs(drctx)+abs(drcty)>1)then !leftx是左侧格点的x坐标
!        print*,'A'
        if(drctx>0)then
            if(drcty>0)then
                leftx=x1
                lefty=y1+1
                rightx=x1+1
                righty=y1
            else
                leftx=x1+1
                lefty=y1
                rightx=x1
                righty=y1-1
            end if
        else
            if(drcty>0)then
                leftx=x1-1
                lefty=y1
                rightx=x1
                righty=y1+1
            else
                leftx=x1
                lefty=y1-1
                rightx=x1-1
                righty=y1
            end if
        end if
    else
!        PRINT*,x1,y1,drctx,drcty
        leftx=x2-drcty
        lefty=y2+drctx
        rightx=x2+drcty
        righty=y2-drctx
    end if
    if(leftx<=0)then
        leftx=leftx+x_num
    else if(leftx>x_num)then
        leftx=leftx-x_num
    end if
    if(rightx<=0)then
        rightx=rightx+x_num
    else if(rightx>x_num)then
        rightx=rightx-x_num
    end if
end subroutine

subroutine get_val_real(array,x,y,val)
use constant
implicit none
integer::x,y
real::array(x_num,y_num),val(8)
return
end subroutine

subroutine get_val_int(array,x,y,val)
use constant
implicit none
integer::x,y,val(8),array(x_num,y_num)

return
end subroutine

subroutine get_cord(x,y,en_x,en_y)
!use constant
implicit none
integer::x,y,en_x(0:7),en_y(0:7),k
real::dx,dy
do k=0,7
    if(abs(cos(pi/4*k))>eps)then
        dx=cos(pi/4*k)/abs(cos(pi/4*k))
    else
        dx=0
    end if
    if(abs(sin(pi/4*k))>eps)then
        dy=sin(pi/4*k)/abs(sin(pi/4*k))
    else
        dy=0
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

function drct_judge(grid,x1,y1,x2,y2,drct) !判断两个点连成的等值线的方向类型(返回值：两边一样是0，左高右低是1，左低右高是2，两边都是等值线格点无法判断的为3)
    integer::grid(x_num,y_num),x1,x2,y1,y2,drctx,drcty,drct
    integer::drct_judge,leftx,lefty,rightx,righty
!            do i=1,cnt
    call llrr(x1,y1,x2,y2,leftx,lefty,rightx,righty,drct)
    if((grid(leftx,lefty)==1.or.grid(leftx,lefty)==-1.or.grid(leftx,lefty)==-2).and.(grid(rightx,righty)==1.or.grid(rightx,righty)==-1.or.grid(rightx,righty)==-2))then
        drct_judge=3
        good=0
    else if((grid(leftx,lefty)==2.or.grid(rightx,righty)==0).and.grid(leftx,lefty)/=grid(rightx,righty))then
        drct_judge=1
        if(grid(leftx,lefty)==2.and.grid(rightx,righty)==0)then
            good=1
        else
            good=0
        end if
    else if((grid(leftx,lefty)==0.or.grid(rightx,righty)==2).and.grid(leftx,lefty)/=grid(rightx,righty))then        
        drct_judge=2
        if(grid(leftx,lefty)==0.and.grid(rightx,righty)==2)then
            good=1
        else 
            good=0
        end if
    else if(grid(leftx,lefty)==grid(rightx,righty))then
        drct_judge=0
        good=0
    else
        print*,"!What's wrong?"
        good=0
        drct_judge=-5
        stop
    end if
    return
end function

function grad_judge(pv,x1,y1,x2,y2)
    real::pv(x_num,y_num),grad_judge
    integer::x1,y1,x2,y2,leftx,lefty,rightx,righty,drct
    call llrr(x1,y1,x2,y2,leftx,lefty,rightx,righty,drct)
    grad_judge=pv(leftx,lefty)-pv(rightx,righty)
end function
end module
