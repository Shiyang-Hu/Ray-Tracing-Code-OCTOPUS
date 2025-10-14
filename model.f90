!Please enter the metric potential f(r).
subroutine metric_potential(r,f)
    use metric_parameters
    implicit none
    real*8 r,f
    f=1d0-2d0/r-4d0*math_pi*(r_s+2d0*r)*r_s*r_s*r_s*rho_s/3d0/(r_s+r)/(r_s+r)
    return
    end subroutine
!Please enter the first-order derivative of the metric potential f(r) with respect to r.    
function df_1st(r)
    use metric_parameters
    implicit none
    real*8 df_1st
    real*8 r
    df_1st=2d0/r**2d0-(8d0*math_pi*r_s**3d0*rho_s)/(3d0*(r+r_s)**2d0)+(8d0*math_pi*r_s**3d0*rho_s*(2d0*r+r_s))/(3d0*(r+r_s)**3d0)
    return
    end 
!Please enter the second-order derivative of the metric potential f(r) with respect to r.
function df_2nd(r)
    use metric_parameters
    implicit none
    real*8 df_2nd
    real*8 r
    df_2nd=(32d0*math_pi*r_s**3d0*rho_s)/(3d0*(r+r_s)**3d0)-4d0/r**3d0-(8d0*math_pi*r_s**3d0*rho_s*(2d0*r+r_s))/(r+r_s)**4d0
    return
    end
!Please enter the third-order derivative of the metric potential f(r) with respect to r.    
function df_3rd(r)
    use metric_parameters
    implicit none
    real*8 df_3rd
    real*8 r
    df_3rd=12d0/r**4d0-(48d0*math_pi*r_s**3d0*rho_s)/(r+r_s)**4d0+(32d0*math_pi*r_s**3d0*rho_s*(2d0*r+r_s))/(r+r_s)**5d0
    return
    end    