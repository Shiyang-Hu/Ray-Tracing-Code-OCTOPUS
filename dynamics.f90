subroutine hamiltonian(vector,hami)
    implicit none
    real*8,dimension(8)::vector
    real*8 t,r,theta,phi,p_t,p_r,p_theta,p_phi
    real*8 g00,g11,g22,g33
    real*8 hami
    call anti_dimensional(vector,t,r,theta,phi,p_t,p_r,p_theta,p_phi)
    call contravariant_metric(r,theta,g00,g11,g22,g33)
    hami=(1d0/2d0)*(p_t*p_t*g00+p_r*p_r*g11+p_theta*p_theta*g22+p_phi*p_phi*g33)
    return
    end subroutine
    
subroutine event_horizon(r_e)
    implicit none
    real*8 r_e
    real*8 r_0,r
    real*8 f,df
    real*8::error=1d-15
    real*8,external::df_1st
    integer i,max_n
    i=0
    max_n=1000
    r=2.5d0
    do while(.true.)
        r_0=r
        call metric_potential(r,f)
        r=r-f/df_1st(r)
        i=i+1
        if(i.gt.max_n)exit
        if(dabs(r-r_0).le.error)exit
    end do
    r_e=dabs(r)
    return
    end subroutine    
    
subroutine photon_sphere(r_ps)
    implicit none
    real*8 r_ps
    real*8 r_0,r
    real*8 f,df
    real*8::error=1d-15
    integer i,max_n
    i=0
    max_n=1000
    r=r_ps
    do while(.true.)
        r_0=r
        call ps_function(r,f)
        call diff_ps_function(r,df)
        r=r-f/df
        i=i+1
        if(i.gt.max_n)exit
        if(dabs(r-r_0).le.error)exit
    end do
    r_ps=dabs(r)
    return
    end subroutine    
    
subroutine impact_parameter(r_ps,b_ps)
    implicit none
    real*8 r_ps,b_ps
    real*8 f
    call metric_potential(r_ps,f)
    b_ps=r_ps/dsqrt(f)
    return
    end subroutine    
    
subroutine isco_orbit(r_isco)
    implicit none
    real*8 r_isco
    real*8 r_0,r
    real*8 f,df
    real*8::error=1d-15
    integer i,max_n
    i=0
    max_n=1000
    r=r_isco
    do while(.true.)
        r_0=r
        call isco_function(r,f)
        call diff_isco_function(r,df)
        r=r-f/df
        i=i+1
        if(i.gt.max_n)exit
        if(dabs(r-r_0).le.error)exit
    end do
    r_isco=dabs(r)
    return
    end subroutine        

subroutine angular_obs(b_p,angular_diameter)   
    use black_hole_parameters
    use metric_parameters
    implicit none
    real*8 b_p,angular_diameter
    angular_diameter=6.191165d0*1d-8*b_p*mass_ratio/distance_obs/math_pi
    return
    end subroutine

subroutine infalling_velocity(r,theta,e_isco,l_isco,r_dot)
    implicit none
    real*8 e_isco,l_isco,r_dot,r,theta
    real*8 g00,g11,g22,g33
    real*8 p_r
    call contravariant_metric(r,theta,g00,g11,g22,g33)
    p_r=dsqrt((-1d0-g00*e_isco*e_isco-g33*l_isco*l_isco)/g11+1d-16)
    r_dot=-(p_r*g11)
    return
    end subroutine    
    
subroutine isco_function(r,f)
    implicit none
    real*8 r,f,f_r
    real*8,external::df_1st,df_2nd
    call metric_potential(r,f_r)
    f=(3d0*f_r*df_1st(r))/(2d0*df_1st(r)*df_1st(r)-f_r*df_2nd(r))-r    
    return
    end subroutine   
    
subroutine diff_isco_function(r,df)
    implicit none
    real*8 r,df,f_r
    real*8,external::df_1st,df_2nd,df_3rd
    call metric_potential(r,f_r)
    df=((3d0*df_1st(r)*df_1st(r)+3d0*f_r*df_2nd(r))*(2d0*df_1st(r)*df_1st(r)-f_r*df_2nd(r))-3d0*f_r*df_1st(r)*(4d0*df_2nd(r)*df_1st(r)-(df_1st(r)*df_2nd(r)+f_r*df_3rd(r))))/((2d0*df_1st(r)*df_1st(r)-f_r*df_2nd(r))**2d0)-1d0
    return
    end subroutine

subroutine ps_function(r,f)
    implicit none
    real*8 r,f,f_r
    real*8,external::df_1st
    call metric_potential(r,f_r)
    f=((1d0/2d0)*f_r**(-1d0/2d0)*df_1st(r)*r-f_r**(1d0/2d0))/r/r
    return
    end subroutine

subroutine diff_ps_function(r,df)
    implicit none
    real*8 r,df,f_r
    real*8,external::df_1st,df_2nd
    call metric_potential(r,f_r)
    df=((((-1d0/4d0)*f_r**(-3d0/2d0)*df_1st(r)*df_1st(r)+(1d0/2d0)*f_r**(-1d0/2d0)*df_2nd(r))*r+(1d0/2d0)*f_r**(-1d0/2d0)*df_1st(r)-(1d0/2d0)*f_r**(-1d0/2d0)*df_1st(r))*r*r-((1d0/2d0)*f_r**(-1d0/2d0)*df_1st(r)*r-f_r**(1d0/2d0))*2d0*r)/r/r/r/r
    return
    end subroutine
    
subroutine difffunction(u,diff_u)
implicit none
real*8,dimension(8)::u,diff_u
real*8 t,r,theta,phi,p_t,p_r,p_theta,p_phi
real*8 g00_r,g11_r,g22_r,g33_r,g33_theta
real*8 g00,g11,g22,g33
call anti_dimensional(u,t,r,theta,phi,p_t,p_r,p_theta,p_phi)
!if(theta.lt.1d-8)theta=1d-8
call diff_contravariant_metric(r,theta,g00_r,g11_r,g22_r,g33_r,g33_theta)
call contravariant_metric(r,theta,g00,g11,g22,g33)
diff_u(1)=-(g00*p_t)
diff_u(2)=-(g11*p_r)
diff_u(3)=-(g22*p_theta)
diff_u(4)=-(g33*p_phi)
diff_u(5)=0d0
diff_u(6)=-(-(1d0/2d0)*(g00_r*p_t*p_t+g11_r*p_r*p_r+g22_r*p_theta*p_theta+g33_r*p_phi*p_phi))
diff_u(7)=-(-(1d0/2d0)*(g33_theta*p_phi*p_phi))
diff_u(8)=0d0
return
end subroutine  
    
subroutine covariant_metric(r,theta,g_00,g_11,g_22,g_33)
    implicit none
    real*8 g_00,g_11,g_22,g_33
    real*8 r,theta
    real*8 f
    call metric_potential(r,f)
    g_00=-f
    g_11=1d0/f
    g_22=r*r
    g_33=r*r*dsin(theta)*dsin(theta)
    return
    end subroutine 
   
subroutine contravariant_metric(r,theta,g00,g11,g22,g33)
    implicit none
    real*8 g00,g11,g22,g33
    real*8 g_00,g_11,g_22,g_33
    real*8 r,theta
    call covariant_metric(r,theta,g_00,g_11,g_22,g_33)
    g00=1d0/g_00
    g11=1d0/g_11
    g22=1d0/g_22
    g33=1d0/g_33
    return
    end subroutine    
        
subroutine diff_contravariant_metric(r,theta,g00_r,g11_r,g22_r,g33_r,g33_theta)
    use metric_parameters
    implicit none
    real*8 r,theta
    real*8 f
    real*8 g00_r,g11_r,g22_r,g33_r
    real*8 g33_theta
    real*8,external::df_1st
    call metric_potential(r,f)
    g00_r=df_1st(r)/f/f
    g11_r=df_1st(r)
    g22_r=-2d0/r**3d0
    g33_r=-2d0/(r**3d0*dsin(theta)**2d0)
    g33_theta=-(2d0*dcos(theta))/(r**2d0*dsin(theta)**3d0)
    return
    end subroutine     
    
subroutine motion_constant(r_s,energy,angular_momentum)
    implicit none
    real*8 r_s,energy,angular_momentum
    real*8 f_r
    real*8,external::df_1st
    call metric_potential(r_s,f_r)
    energy=f_r/dsqrt(f_r-(1d0/2d0)*r_s*df_1st(r_s))
    angular_momentum=r_s*dsqrt(r_s*df_1st(r_s))/dsqrt(2d0*f_r-r_s*df_1st(r_s))
    return
    end subroutine     