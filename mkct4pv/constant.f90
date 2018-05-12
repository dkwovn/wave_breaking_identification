module constant
    integer,parameter::x_num=360,y_num=181,tlen=1440,st_yr=4,end_yr=25
    real,parameter::pi=3.1415927,eps=1e-6,undef=-9.99e8
!    integer,save::grid(x_num,y_num),patch_map(x_num,y_num)
end module

module typedef
implicit none
type::contour
    integer::x,y,ct_type            !1代表绕纬圈的封闭等值线，2代表局地的封闭等值线，3代表不封闭的等值线
    integer::drctx,drcty,dr_type  !drctx,drcty代表以此点为端点的单位长度等值线的走向，以矢量表示 ; dr_type为1代表，左高右低，为2代表右高左低  
    real::val
    type(contour),pointer::prev
    type(contour),pointer::next
end type contour
type::cord
    integer::x,y
end type cord
end module typedef
