!�ӳ���detect����˵����
!1.grid�Ǽ�⵽�ĵ�ֵ�߸�㳡���ǵ�ֵ�ߵ�Ϊ1�����ڵ�ֵ�ߵ�Ϊ2��С�ڵ�ֵ�ߵ�Ϊ0
!2.environ��һ����׼�ռ䳡��ÿ���������ֵ��������Χ8������зǵ�ֵ�߸����
!3.trace_org�Ǽ�¼�Ѿ������յĵ�ֵ�߸��ı�׼�����ѽ��ܵ�Ϊ1��δ���ܵ�Ϊ0
!4.pָ��Ŀǰ�����յĵ�ֵ�߱�����ָ��
!5.tmpx,tmpy����ȥ�ļ�����������ͨ��Ҫ��ĵ�ֵ�߸�㣩��������
!6.cnt����ȥ�����������ĸ����
!7.dr_mode ��ֵ���������ͣ�1Ϊ����ҵͣ�2Ϊ����Ҹߡ��ǵ��ú���ָ���ģ�����������ָ������ĸ��Żᱻ����
!8.flag ����������������⴦��Ĳ���
!  flag=-1 ������prev��� ������if���������.and.(i/=p%prev%x.or.j/=p%prev%y.or.flag==-1)�ĸ����Ա�����
!  flag=0 ����trace�����ƣ�������֮ǰ�߹��ĸ�㣻Ҳ����judge==0,������ֵ��ȵ������
!  flag=1
!  flag=4 ��ͨ���
!  flag=6 һ�ֿ��ɷ�������judge/=dr_mode��judge=3������Ϊ��ֵ�����߶��ǵ�ֵ�߸���޷��жϲŵ��µ�ʱ�����������ĵ㱻����
!  flag=7 һ�ֿ��ɷ�������judge==0����������ȵ�ʱ������
!  flag=8 ������֮ǰ�����ĵ�ֵ�߸��
!  flag=9 ����trace�����ƣ���judge=3������

!  flag1=-1 ������prev���&&����trace&&
!  flag1=0  ����trace
!  flag1=2  ����������ɵĵ�ֵ�ߵĵ㣨grid=-1�ĵ㣩
!  flag1=4  ��ͨ���
!  flag2=6  ����judge==3�ĸ�㣨���߶��޷��жϣ�
!  flag2=7  ����judge==0�ĸ�㣨������ȣ�
!  flag2=4  ��ͨ���
!    flag2=9 allow previous grid point
!about direction: east is 0, and the number increase anticlockwisely.

module dtct
use typedef
use func

contains
!�����Χ�ĵ�ֵ�߸��***************************************
subroutine detect(grid,environ,trace_org,p,tmpx,tmpy,cnt,dr_mode,flag1,flag2,stat1,stat2,p_node,mode,delta_x) !�����stat1������ŵ�ֵ�߷���0--7
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
if(flag1==0.or.flag1==-1.or.flag1==2)then !flag=-1��������prev��������������ֵ�ߵľ�ͷ��
    do i=1,x_num
        do j=1,y_num
            trace(i,j)%x=0 !��������߻�ͷ·���ĵڶ�������
            trace(i,j)%y=0 !��������߻�ͷ·���ĵڶ�������
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


