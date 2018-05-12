!program:
!    given the isoline record, identify the wave breaking, output the occurrence
!    frequency and relevant contour segments
!History:
!2014/02/10    cjliu    1st version
!2014/02/18    cjliu    2nd version
!    correctly output the RWB segments, extend the length within which two RWB are
!    considered one
!04/01/2014    cliu    3rd version
!    shorten the RWB segments, to make my capturing algorithm stricter
!04/17/2014 cliu     4th version
!    delete more small RWB, make the threshold larger
!08/14/2015    cliu    5th version
!    include the lower lobe into the coverage of RWB, make output of coverage
!    flexible in resolution

program pick_event
use func
implicit none
integer::year,leap,q,t,cnt,temp1,temp2,i,j,k,l1,l2,m,sum,sum2,flag,style,event_cnt,dx,something,space_cnt,rdc_cnt
integer::grid(x_num,y_num),dx1,dx2,lbin(l_num,l_num),k_mid,k_prev,k_next,circle,rcd(x_num),pq_cnt,valint,tt
real::cov1(nx,ny),cov2(nx,ny),seg1(nx,ny),seg2(nx,ny)
real::val(11),area_tmp,dll
type(isoline),pointer::p_iso,iso_hd,iso_end,p_st,p_end,p_tmp,p_restart
type(event),pointer::p,pq,p_hd,pq_hd,ps,pp_tmp
character(len=100)::file_src,file_obj,dfile1,fawbfile,fcwbfile,sawbfile,scwbfile
integer::error
data val /300,305,310,315,320,325,330,335,340,345,350/
lbin=0
space_cnt=0
rdc_cnt=0
do year=st_yr,end_yr
    print*,year
    if(mod(year,4)==0.and.(mod(year,100)/=0.or.mod(year,400)==0))then
        leap=1
    else
        leap=0
    end if
    write(file_obj,"('/home/cliu/data/rwbsm/era/rwbsm/','RWB.database.cmplt.',i2.2,'.bin')")year  !参数
    write(fawbfile,"('/home/cliu/data/rwbsm/era/','fawb.cmplt.',i2.2,'.bin')")year
    write(fcwbfile,"('/home/cliu/data/rwbsm/era/','fcwb.cmplt.',i2.2,'.bin')")year
    write(sawbfile,"('/home/cliu/data/rwbsm/era/','sawb.cmplt.',i2.2,'.bin')")year
    write(scwbfile,"('/home/cliu/data/rwbsm/era/','scwb.cmplt.',i2.2,'.bin')")year
    open(unit=16,file=fawbfile,access='stream',form='unformatted')
    open(unit=17,file=fcwbfile,access='stream',form='unformatted')
    open(unit=18,file=sawbfile,access='stream',form='unformatted')
    open(unit=19,file=scwbfile,access='stream',form='unformatted')

    open(unit=11,file=file_obj,access='stream',form='unformatted')
    do t=1,tlen
    print*,'t=',t

        cov1=0
        cov2=0
        seg1=0
        seg2=0
        allocate(p)      !modified
        p%coverage=0
        space_cnt=space_cnt+1
        event_cnt=0
        do q=1,11   !开始每一条等值线上breaking的检测
            !=== read contour ===
            valint=val(q)
            write(file_src,"('/home/cliu/data/rwbsm/era/isoline/','isoline',i3,'.',i4.4,'.',i2.2,'.txt')")valint,t,year
            open(unit=10,file=file_src,form='formatted',status='old')
                
            allocate(p_iso)
            iso_hd=>p_iso
            cnt=0
            do 
                read(unit=10,fmt="(2i6)",iostat=error)p_iso%x,p_iso%y
                if(error<0)then
                    if(cnt/=0)then
                        p_iso%prev%next=>null()
                    end if
                    deallocate(p_iso)
                    exit
                end if
                allocate(p_iso%next)
                p_iso%next%prev=>p_iso
                p_iso=>p_iso%next
                cnt=cnt+1
            end do
            close(10)

            if(cnt==0)then
                cycle
            end if
            !=== find the end of the contour
            p_iso=>iso_hd
            do 
                if(associated(p_iso%next))then
                    p_iso=>p_iso%next
                else
                    exit
                end if
            end do
            iso_end=>p_iso
!***********START breaking detect on one individual line*******************************
            !===if it is a closed contour, link the head and end===
            if(distance(p_iso%x,p_iso%y,iso_hd%x,iso_hd%y)<2)then
!                print*,'closed'
                p_iso%next=>iso_hd
                iso_hd%prev=>p_iso
                dx=0
                p_iso=>iso_hd
                do 
                    dx=dx+cal_dx(p_iso%x,p_iso%next%x)
                    if(associated(p_iso,iso_end))then
                        exit
                    end if
                    p_iso=>p_iso%next
                end do
                if(dx==0)then
!                    print*,'isolated contour'
                    style=2
                else if(dx<0)then
                    print*,'contour with wrong direction'
                    style=2
!                    stop                    
                else
                    style=1   !=== style=1 means circumpolar contour ===
                end if
            else
                style=2
            end if

            !===if it is not, release the space and start another detect===
            if(style==2)then        
                do 
                    if(associated(p_iso,iso_hd))then
                        exit
                    end if
                    p_iso=>p_iso%prev
                    deallocate(p_iso%next)
                end do
                deallocate(iso_hd)    
!                print*,'not qualified contour'
!                read(*,*)something
                cycle
            end if

            !===set the iso_type===
            p_iso=>iso_hd
            do i=1,cnt
                dx1=cal_dx(p_iso%prev%x,p_iso%x)
                dx2=cal_dx(p_iso%x,p_iso%next%x)
                if(dx2<=0.and.dx1<=0)then
                    p_iso%iso_type=2
                else if((dx2<0.and.dx1>0).or.(dx2>0.and.dx1<0))then
                    p_iso%iso_type=4
                else
                    p_iso%iso_type=0
                end if
                p_iso=>p_iso%next
            end do

            !===start from where there is no breaking===
            do 
                p_iso=>p_iso%next
                if(p_iso%iso_type/=2)then
                    iso_hd=>p_iso
                    exit
                end if
            end do

            !===detect and set the breaking events===
            p_iso=>p_iso%next
            circle=0
            do 
                !print*,'p_iso:',p_iso%x,p_iso%y,p_iso%iso_type
                circle=circle+1
                if(circle>3200)then
                    print*,'dead circle'
                    stop
                end if
                if(associated(p_iso,iso_hd))then
                    exit
                end if
                if(p_iso%iso_type==2.and.p_iso%prev%iso_type/=2)then
                    p_st=>p_iso
                    dll=0
                end if
                if(p_iso%iso_type==2.and.p_iso%next%iso_type/=2)then
                    p_end=>p_iso
            !        p_tmp=>p_iso
                    !=== set the restart point as the end of type=2 segment ====
                    p_restart=>p_end
                    !===recognize this ac_type==2 segment as an breaking event，calculate relative variables===
                    if(cal_dx(p_st%x,p_end%x)/=0)then
                        allocate(p%next)
                        p%next%coverage=0
                        space_cnt=space_cnt+1
                        if(event_cnt/=0)then
                            p%next%prev=>p
                            p=>p%next
                        else
                            p_hd=>p%next
                            pp_tmp=>p
                            p=>p%next
                            deallocate(pp_tmp)
                            rdc_cnt=rdc_cnt-1
                        end if
                        event_cnt=event_cnt+1
                    !=== move back and decide the breaking type ===
                        p_iso=>p_st
                        do 
                            if(p_iso%x==p_end%x)then
                                exit
                            end if
                            p_iso=>p_iso%prev
                        end do
                        if(p_iso%y>p_end%y)then
                            p%ac_type=1
                            p_tmp=>p_st
                            p_st=>p_iso
                        else if(p_iso%y<=p_end%y)then
                            p%ac_type=2
                        else
                            print*,'error in deciding ac_type!'
                    !        read(*,*)something
                    !        exit
                            cycle
                        end if
                        if(p%ac_type==2)then
                            p_iso=>p_end
                            p_tmp=>p_end
                            do 
                                if(p_iso%x==p_st%x)then
                                    p_end=>p_iso
                                    exit
                                end if
                                p_iso=>p_iso%next
                            end do
                        !=== not sure if the following part is redundant ===
                            p_iso=>p_end
                            do 
                                p_iso=>p_iso%prev
                                if(p_iso%x==p_end%x)then
                                    p_st=>p_iso
                                    exit
                                end if
                            end do
                        else if(p%ac_type==1)then
                            do 
                                p_iso=>p_iso%next
                                if(p_iso%x==p_end%x)then
                                    p_end=>p_iso
                                    exit
                                end if
                            end do
                        end if
                        !=== finish locate a segment of contour having area ===
                        p%val=val(q)                            
                        p%lz=0
                        call event_cal(p_st,p_end,p,p%ac_type)
                        if(p%ac_type==1)then
                            p_st=>p_tmp
                            p_iso=>p_end
                            do 
                                if(p_iso%x==p_st%x)then
                                    p_end=>p_iso
                                    exit
                                end if
                                p_iso=>p_iso%next
                            end do
                            call event_cal(p_st,p_end,p,2)
                        else
                            p_end=>p_tmp
                            p_iso=>p_st
                            do 
                                if(p_iso%x==p_end%x)then
                                    p_st=>p_iso
                                    exit
                                end if
                                p_iso=>p_iso%prev
                            end do
                            call event_cal(p_st,p_end,p,1)
                        end if
                    end if
                    !===end detecting and move on from where it stop===
                    p_iso=>p_restart
                end if
                p_iso=>p_iso%next
            end do
!***********FINISH breaking detecting in one individual line************************************            
            p_iso=>iso_hd
            do
                if(t==210.and.q==11) print*,p_iso%x,p_iso%y
                p_iso=>p_iso%next
                deallocate(p_iso%prev)
                if(associated(p_iso,iso_end))exit
            end do
            deallocate(p_iso)
        end do
!*******START FILTERING*******************************************************************************
        if(event_cnt==0)then
            print*,'no event'
            write(16)cov1
            write(17)cov2
            cycle
        end if
        p%next=>null()
        !====delete events not big enough====
        p=>p_hd
        pq_cnt=0
        do 
            if(associated(p).eqv..FALSE.)then
                exit
            end if
            if(p%lz<=lz_lmt.or.p%ll<=ll_lmt.or.p%lc1<=lc_lmt)then
                if(associated(p,p_hd))then
                    if(associated(p%next))then
                        p_hd=>p%next
                        pp_tmp=>p
                        p=>p%next
                        deallocate(pp_tmp)
                        rdc_cnt=rdc_cnt-1
                        cycle
                    else
                        pq_hd=>null()
                        exit
                    end if
                else if(associated(p%next).eqv..FALSE.)then
                    p%prev%next=>null()
                    deallocate(p)
                    rdc_cnt=rdc_cnt-1
                    exit
                else
                    p%prev%next=>p%next
                    p%next%prev=>p%prev
                    pp_tmp=>p
                    p=>p%next
                    deallocate(pp_tmp)
                    rdc_cnt=rdc_cnt-1
                    cycle
                end if
            end if
            p=>p%next
        end do
        !=====================================
        pq_cnt=0

        tt=0        
        do 
            tt=tt+1
            if(.not.associated(p_hd))then
                exit
            end if
            !*****pick the biggest breaking, also record the frequency of various scale*****
            ps=>p_hd
            cnt=0
            p=>p_hd     
            do 
                if(associated(p).eqv..FALSE.)then
                    exit
                end if
                l1=(log(p%lz)/log(10.0)+0.5)/lbin_len+1
                if(l1>l_num.or.l1<=0)then
                    print*,'lz=',p%lz
                end if
                l2=(log(p%lc1)/log(10.0)+0.5)/lbin_len+1
                if(l2>l_num.or.l2<=0)then
                    print*,'lc=',p%lc1
                end if
                lbin(l1,l2)=lbin(l1,l2)+1
                if(ps%area1<p%area1)then
                    ps=>p
                end if
                p=>p%next
                cnt=cnt+1
            end do
            !====move the biggest one from p to pq====
            if(pq_cnt==0)then
                pq=>ps
                pq_hd=>pq
            else
                pq%next=>ps
                pq=>pq%next
            end if
            pq_cnt=pq_cnt+1
            if(associated(ps,p_hd))then  !delete the biggest one in this if block
                p_hd=>ps%next
                if(associated(p_hd))then
                    p_hd%prev=>null()
                end if
            else if(associated(ps%next).eqv..FALSE.)then
                ps%prev%next=>null()
            else
                ps%prev%next=>ps%next
                ps%next%prev=>ps%prev
            end if        
            ps%next=>null()
            ps%prev=>null()

            !=== count event number ==============
            p=>p_hd
            cnt=0
            do 
                if(associated(p).eqv..FALSE.)then
                    exit
                end if
                if(ps%area1<p%area1)then
                    ps=>p       !=== ??? ===
                end if
                p=>p%next
                cnt=cnt+1
            end do
            !====delete the events overlaping with(within a certain distance from) the biggest one====
            p=>p_hd
            do 
                if(associated(p).eqv..FALSE.)then
                    exit
                end if
                !=== A: added 03/02/2014 ===
                flag=0
                do i=1,nx
                    do j=ny/2,ny
                        if(p%coverage(i,j)/=0.and.pq%coverage(i,j)/=0.and.p%ac_type==pq%ac_type) flag=1
                    end do
                end do
                !=== end A ===
                    !=== reserve the sawb/scwb of p ===
                    if(p%ac_type==1)then
                        do i=1,nx
                            do j=1,ny
                                if(p%coverage(i,j)==1) cov1(i,j)=1
                                if(p%coverage(i,j)==2) seg1(i,j)=1
                            end do
                        end do
                    else if(p%ac_type==2)then
                        do i=1,nx
                            do j=1,ny
                                if(p%coverage(i,j)==1) cov2(i,j)=1
                                if(p%coverage(i,j)==2) seg2(i,j)=1
                            end do
                        end do
                    else
                        print*,'ac_type wrong?'
                    end if
                    !==================================
                    if(associated(p,p_hd))then
                        if(associated(p%next))then
                            p_hd=>p%next
                        else
                            deallocate(p_hd)
                            exit
                        end if
                        if(associated(p_hd%prev))then
                            pp_tmp=>p
                            p=>p%next
                            deallocate(pp_tmp)
                            rdc_cnt=rdc_cnt-1
                            cycle
                        end if
                    else if(associated(p%next).eqv..FALSE.)then
                        deallocate(p%prev%next) !modified: add this line
                        rdc_cnt=rdc_cnt-1
                        exit                !modified: add this line
                    else
                        p%prev%next=>p%next
                        p%next%prev=>p%prev
                        pp_tmp=>p                !modified
                        p=>p%next                !modified
                        deallocate(pp_tmp)        !modified
                        rdc_cnt=rdc_cnt-1
                        cycle                    !modified
                    end if
                !end if
                p=>p%next
            end do
            !=== count event number ==============
            p=>p_hd
            cnt=0
            do 
                if(associated(p).eqv..FALSE.)then
                    exit
                end if
                if(ps%area1<p%area1)then
                    ps=>p       !=== ??? ===
                end if
                p=>p%next
                cnt=cnt+1
            end do
            event_cnt=cnt
        end do

!        print*,'===================='
        pq=>pq_hd
        do 
            if(.not.associated(pq))then
                exit
            end if
            write(11)t,pq%ac_type,pq%left,pq%right
            write(11)pq%x1,pq%y1,pq%lc1,pq%areal1
            write(11)pq%x2,pq%y2,pq%lc2,pq%areal2
            write(11)pq%pe_type,pq%lz,pq%lzl,int(pq%val)
            if(pq%ac_type==1)then
                do i=1,nx
                    do j=1,ny
                        if(pq%coverage(i,j)==1) cov1(i,j)=1
                        if(pq%coverage(i,j)==2) seg1(i,j)=1
                    end do
                end do
            else if(pq%ac_type==2)then
                do i=1,nx
                    do j=1,ny
                        if(pq%coverage(i,j)==1) cov2(i,j)=1
                        if(pq%coverage(i,j)==2) seg2(i,j)=1
                    end do
                end do
            else
                print*,'ac_type wrong?'
            end if
            pp_tmp=>pq    !modified
            pq=>pq%next
            deallocate(pp_tmp) !modified
            rdc_cnt=rdc_cnt-1
        end do
        write(16)cov1
        write(17)cov2
        write(18)seg1
        write(19)seg2
    end do !===for t===
    close(11)
    close(16)
    close(17)
end do !===for year===
open(unit=15,file='/home/cliu/data/rwbsm/era/lbin.bin',form='unformatted')
write(15)lbin
close(15)
end 
