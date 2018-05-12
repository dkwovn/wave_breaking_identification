!Program:
!    get contour lon and lat for ipv.int1.0.????.bin
!History:
!2014/1/14    cliu    1st version

!主程序!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program mkct_for_pt
use typedef
use func
use dtct
use constant
use reverse_mod
implicit none
integer::dr_mode,stat1(x_num,y_num),stat2(x_num,y_num),drct,a,b,c,d,dy2,anything,fail,dx(9),dy(9),mode(x_num,y_num) !stat1数组用来记录方向
integer::q,i,j,k,l,m,n,year,t,cnt,tt,tmpx(9),tmpy(9),st_rec,nodex,nodey,max,good_cnt,ire_cnt&
&,drctx,drcty,environ(x_num,y_num),t_cnt,flag1,flag2,temp3,temp4,grd_cnt(10),line_cnt&
&,ire_x(100),ire_y(100),ttt,delta_x!记得每开始一条新的等值线，trace要清空,stat也要清空
integer::patch_map(x_num,y_num),grid(x_num,y_num)
real::var(x_num,y_num),val(11),temp,tempx,tempy,temp1,temp2
type(cord)::trace(x_num,y_num)
character(len=100)::dfile1,ofile,something,subfile,yr_str
integer::en_x(0:7),en_y(0:7)
!character(len=30),parameter::path='/newraid2/cjliu/ptpv/',path2='/newraid2/cjliu/contour/'
type(contour),target::node,head(10)
type(contour),pointer::p,p_node,p_head
data val /300,305,310,315,320,325,330,335,340,345,350/


loop_yr: do year=st_yr,end_yr  !最外（第一）层循环
!write(dfile1,"('/home/cliu/data/merra/6hourly/pv/','ipv.int1.0.',i2,'.6hr.bin')")year
write(yr_str,"(i2)")year
dfile1 = '/home/cliu/data/gfdl_dynamic_core/north_bl_ddt_A5n/ipv.int1.0.'//trim(adjustl(yr_str))//'.6hr.bin'
open(unit=10,file=dfile1,access='direct',recl=x_num*y_num*4,form='unformatted',status='old')
if(mod(year,4)==0.and.(mod(year,100)/=0.or.mod(year,400)==0))then
    !t_cnt=tlen+4
    t_cnt=tlen
else
    t_cnt=tlen
end if
loop_t: do t=1,tlen   !第二层循环
print*,'t=',t,'year=',year
loop_q: do q=1,11
read(unit=10,rec=(t-1)*11+q)var
write(ofile,"('/home/cliu/data/gfdl_dynamic_core/north_bl_ddt_A5n/rwbsm/isoline/','isoline',i3,'.',i4.4,'.',i2.2,'.txt')"),int(val(q)),t,year
write(subfile,"('/home/cliu/data/rwbsm/merra/','grid.bin')")
open(unit=12,file=ofile)
open(unit=14,file=subfile,access='stream',form='unformatted')
do i=1,x_num
    do j=1,y_num
        grid(i,j)=0
        environ(i,j)=0
    end do
end do

!标记等值线格点********************************************
call gen_grid(var,grid,2.0)
write(14) grid*1.0
close(14)

forall(i=1:x_num,j=1:y_num)
    trace(i,j)%x=0
    trace(i,j)%y=0
end forall!清空trace
stat2=0
line_cnt=0
loop_m: do m=1,100
loop_n: do n=y_num/2+1,y_num-3
if(m==1.or.n==2.or.n==y_num-1)then
if((grid(m,n)==0.and.(grid(m,n+1)==2)).and.trace(m,n)%x==0)then
    trace(m,n)%x=m
    trace(m,n)%y=n-1


!=== make the patch
patch_map=0
stat1=-1
call gen_patch(grid,patch_map,m,n+1,2)

!=== output patch ===
open(unit=31,file='/home/cliu/data/gfdl_dynamic_core/north_bl_ddt_A5n/rwbsm/patch.bin',access='stream',form='unformatted')
write(31)patch_map*1.0
close(31)
!==================

p_node=>node
p=>node

p%x=m
p%y=n
    p%val=val(q)
    temp=0
cnt=0
dr_mode=1
p%prev=>p              

ire_cnt=0
tt=0
st_rec=0
fail=0
delta_x=0
do while(.TRUE.)
    if(tt>2000)then
        print*,'dead circle'
        print*,node%x,node%y,trace(node%x,node%y)
        print*,'q=',q
        stop
    end if
    cnt=0
    temp=0
    do i=1,9
        tmpx(i)=0
        tmpy(i)=0
        dx(i)=0
        dy(i)=0
    end do    
    if(tt==2)then
        trace(node%x,node%y)%x=0
        trace(node%x,node%y)%y=0
        where(stat2==-9)
            trace%x=0
            trace%y=0
        end where
    end if
    flag1=4
    flag2=4
    call detect(patch_map,environ,trace,p,tmpx,tmpy,cnt,dr_mode,flag1,flag2,stat1,stat2,p_node,mode,delta_x)
    if(cnt==0)then
        
        if((p%y==2.or.p%y==y_num-1))then
            print*,'no candidate in boundary'
            stop
            if(p%x==node%x.and.p%y==node%y)then
                if(dr_mode==1)then
                    dr_mode=2
                    stop
                else if(dr_mode==2)then
                    dr_mode=1
                    stop
                end if
                cycle
            else
            st_rec=st_rec+1        
            if(st_rec==1)then   !若第一次碰到南北界点，清空之前的trace，把这个点设成头格点，颠倒等值线走向（dr_mode）
            call reverse(p,node,dr_mode,stat1,mode)
                cycle                
            else if(st_rec==2)then
                call reverse(p,node,dr_mode,stat1,mode)
                p%next=>node
                line_cnt=line_cnt+1
                exit
!                cycle
            else if(st_rec==3)then !若是第三次碰到南北界格点，完成等值线，跳出循环
                print*,'wrong:st_rec==3'
                stop
            else 
                print*,'st_rec wrong!!'
            end if    
            end if        
        else
        end if
    end if

    if(cnt==0)then
        if(tt==0)then
            fail=1
            exit
        else
            flag1=4
            flag2=9
            call detect(patch_map,environ,trace,p,tmpx,tmpy,cnt,dr_mode,flag1,flag2,stat1,stat2,p_node,mode,delta_x)
        end if
    end if
    if(cnt==0)then
        !close(12)
        print*,"cannot go ahead!",'q=',q,'t=',t,'year=',year
        print*,'dr_mode=',dr_mode
        print*,'node',node%x,node%y
        print*,'point',p%x,p%y
        !=== relax the retriction: skip dead lines ===
        fail=1
        exit
        !stop
    end if



        if(flag1==8)then
            exit loop_n
        end if

tt=tt+1
!Here flag1 carries message from sub detect and this number indicate that the
!head point has been updated and no need to allocate new grid
!设置结点值
        allocate(p%next)
        !stat2(tmpx(1),tmpy(1))=drct_judge(grid,p%x,p%y,tmpx(1),tmpy(1),drct)+1  !stat2值为1代表等值线两边值相等的点（对应尾点），在flag==0时应优先考虑        
        call llrr(p%x,p%y,tmpx(1),tmpy(1),a,b,c,d,drct)
        if(stat1(tmpx(1),tmpy(1))==-1)then
            stat1(tmpx(1),tmpy(1))=drct       !stat1 保存8个方向,如果方向相等，即剔除，所以选取尾点
        end if
        p%next%prev=>p
        !call get_cord(p%x,p%y,en_x,en_y)
        p=>p%next
        p%x=tmpx(1)
        p%y=tmpy(1)
        if(q>=8) print*,'contour point:',p%x,p%y
        !=== trace are a,b for dr_mode=1, are c,d for dr_mode=2
        trace(p%x,p%y)%x=a
        trace(p%x,p%y)%y=b
        p%drctx=p%x-p%prev%x
        if(p%drctx<-100)then
            p%drctx=p%drctx+x_num
        else if(p%drctx>100)then
            p%drctx=p%drctx-x_num
        end if
        p%drcty=p%y-p%prev%y !设置方向
        if(p%y<y_num-1)then
            delta_x=delta_x+p%drctx
        end if
        p%val=val(q)
        mode(p%x,p%y)=dr_mode
        if(p%x==node%x.and.p%y==node%y)then
            p%val=-541
            if(delta_x>0)then 
            !relax the restriction so that polar unclosed lines could be
            !included
                line_cnt=line_cnt+1
                fail=0
                exit
            else 
                fail=1
                print*,'fail',delta_x,q
                exit
            end if
        end if
end do !对应while的循环
if(fail==0)then
    grd_cnt(line_cnt)=tt
    p_head=>head(line_cnt)
    head(line_cnt)%x=node%x
    head(line_cnt)%y=node%y
    head(line_cnt)%next=>node%next
    head(line_cnt)%val=8849
    p=>head(line_cnt)%next
    ttt=0
    do while(.TRUE.)
        ttt=ttt+1
        if(ttt>2000)then
            print*,'dead cirle in tracing:',head(line_cnt)%x,head(line_cnt)%y
            stop
        end if
        if(p%next%x==head(line_cnt)%x.and.p%next%y==head(line_cnt)%y)then
            exit
        end if
        p=>p%next
    end do
    p%next=>head(line_cnt)
    head(line_cnt)%prev=>p
    exit loop_m
end if


    p=>node
    tt=0
    if(ire_cnt>0)then
        do i=1,ire_cnt
            grid(ire_x(ire_cnt),ire_y(ire_cnt))=1
        end do
    end if
end if
end if
end do loop_n
end do loop_m

max=1
if(line_cnt>1)then
    temp=grd_cnt(1)
    do i=2,line_cnt
        if(grd_cnt(i)>temp)then
            temp=grd_cnt(i)
            max=i
        end if
    end do
end if
if(line_cnt>=1.and.grd_cnt(max)>=x_num)then
    p=>head(max)
    tt=0
    do while(.TRUE.)
        tt=tt+1
        if(tt>2000)then
            print*,'dead circle in writing','q=',q,'t=',t,'year=',year
            print*,'head:',head(max)%x,head(max)%y
            stop
        end if
        write(12,"(2i6)")p%x,p%y
        print*,p%x,p%y
        
        p=>p%next
        if((p%x==head(max)%x).and.(p%y==head(max)%y))then
            exit
        end if
    end do
    do i=1,line_cnt
        nodex=head(i)%x
        nodey=head(i)%y
        p=>head(i)
        p=>p%next
        tt=0
        do while(.TRUE.)
            tt=tt+1
            if(tt>2000)then
                print*,'dead circle in deallocating:',nodex,nodey
                stop
            end if
            p=>p%next
            deallocate(p%prev)
            if(p%x==nodex.and.p%y==nodey)then
                exit
            end if
        end do
    end do
end if
close(12)
end do loop_q!对应q的循环

end do loop_t!对应t的循环
end do loop_yr!对应year的循环
close(10)

stop
end

