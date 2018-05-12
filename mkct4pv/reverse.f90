module reverse_mod
use constant
use typedef
contains
subroutine reverse(p,node,dr_mode,stat1,mode)
type(contour),pointer::p,p_node,p_tmp
type(contour),target::node
integer::dr_mode,stat1(x_num,y_num),mode(x_num,y_num),tt
print*,'in reverse'
nodex=p%x        !是界点
nodey=p%y
!***************将原格点的stat1转向（180度），头格点stat1清空
p=>node   
if(dr_mode==1)then
    dr_mode=2
else if(dr_mode==2)then
    dr_mode=1
else    
    print*,'dr_mode wrong!!'
    stop
end if                       
tt=0
do while(.TRUE.)
    tt=tt+1
    if(tt>2000)then
        print*,'dead circle in reverse:'
        stop
    end if
    stat1(p%x,p%y)=mod(stat1(p%next%x,p%next%y)+4,8)
    mode(p%x,p%y)=dr_mode
    if(stat1(p%x,p%y)==0)then
        stat1(p%x,p%y)=8
    end if
    if(p%next%x==nodex.and.p%next%y==nodey)then
        exit
    end if
    p=>p%next
end do
stat1(nodex,nodey)=0
mode(nodex,nodey)=dr_mode
!***************将原来的等值线反向
p=>p%next          !此后p指向界点
node%next=>p
tt=0
do while(.TRUE.)
    tt=tt+1
    if(tt>2000)then
        print*,'dead circle in reverse2:'
        stop
    end if
    p_tmp=>p%next
    p%next=>p%prev    
    p%prev=>p_tmp                
    if(p%next%x==node%x.and.p%next%y==node%y)then
        exit
    end if
    p=>p%next
end do
        
allocate(p%next)
p%next%x=node%x    !此时p%next指向新结点，并赋予原node点的值
p%next%y=node%y
p%next%prev=>p
p=>p%next          !p被赋值后继续新的搜寻（dr_mode==2的）
node%x=nodex
node%y=nodey
node%next=>node%next%next
node%next%prev=>node


end subroutine
end module
