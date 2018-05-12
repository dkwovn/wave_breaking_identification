module constant4pick
    integer,parameter::x_num=360,y_num=181,tlen=1440,st_yr=4,end_yr=25,default=-999,l_num=25
    integer,parameter::nx=128,ny=64
    real,parameter::pi=3.1415927,eps=1e-6,lc_lmt=3.0,lz_lmt=3.5,ll_lmt=15,e=2.718282,a=6370949,unt=pi*a/180

    real,parameter::logl_min=-0.5,logl_max=2.0,lbin_len=(logl_max-logl_min)/l_num
    real,parameter::x0_1=0,y0_1=-90,x0_2=0,y0_2=-87.8638
    real,parameter::ddx1=360.0/x_num,ddy1=ddx1,ddx2=360.0/nx,ddy2=2.7893
    character(len=50),parameter::path='/newraid2/cjliu/contour/',path2='/newraid2/cjliu/rwbSM/'
end module

module typedef
use constant4pick
implicit none


type::isoline
    integer::x,y,year,t,iso_type  !iso_type=1为从左向右走（上高下低），iso_type=2为从右向左走（上低下高）
    type(isoline),pointer::next,prev
end type isoline

type::event
    integer::ac_type    !ac_type=1为AWB,ac_type=2为CWB
    integer::pe_type         !pe_type=1 poleward, pe_type=2 equatorward
    integer::x1,y1,x2,y2 !x1,y1 are upstream bulb, x2,y2 are downstream bulb
    integer::left,right
    real::area1,areal1,lc1,lcl1,lz,lzl
    real::area2,areal2,lc2,lcl2,ll
    real::val
    integer::coverage(nx,ny)
    type(event),pointer::next,prev
end type event

type::point_record
    integer::max,min
end type point_record
end module typedef
