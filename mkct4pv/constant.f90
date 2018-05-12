module constant
    integer,parameter::x_num=360,y_num=181,tlen=1440,st_yr=4,end_yr=25
    real,parameter::pi=3.1415927,eps=1e-6,undef=-9.99e8
!    integer,save::grid(x_num,y_num),patch_map(x_num,y_num)
end module

module typedef
implicit none
type::contour
    integer::x,y,ct_type            !1������γȦ�ķ�յ�ֵ�ߣ�2����ֵصķ�յ�ֵ�ߣ�3������յĵ�ֵ��
    integer::drctx,drcty,dr_type  !drctx,drcty�����Դ˵�Ϊ�˵�ĵ�λ���ȵ�ֵ�ߵ�������ʸ����ʾ ; dr_typeΪ1��������ҵͣ�Ϊ2�����Ҹ����  
    real::val
    type(contour),pointer::prev
    type(contour),pointer::next
end type contour
type::cord
    integer::x,y
end type cord
end module typedef
