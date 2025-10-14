subroutine energy_flux(r,r_eh,flux230,flux86)
    implicit none
    real*8 r,r_eh,flux230,flux86
    real*8 z
    z=dlog(r/r_eh)
    flux230=dexp(-2d0*z+(-1d0/2d0)*z*z)
    flux86=dexp((-3d0/4d0)*z*z)
    !flux86=1d0/r/r
    return
    end subroutine
    
subroutine source_profile(x_p,y_p,z_p,x_m,y_m,z_m,j_v,alpha_v) 
    use lightsource_parameters
    implicit none
    real*8 x_p,y_p,z_p,x_m,y_m,z_m,j_v,alpha_v
    real*8 distance,sigma
    sigma=source_radius/2d0
    distance=dsqrt((x_m-x_p)**2d0+(y_m-y_p)**2d0+(z_m-z_p)**2d0)
    j_v=j_0*dexp(-distance*distance/2d0/sigma/sigma)
    alpha_v=alpha_0*(1d0-distance/source_radius)**2d0
    return
    end subroutine
    
subroutine gravitational_waves(r,theta,phi,h_plus,h_cross)
    use black_hole_parameters
    use gravitational_parameters
    implicit none
    real*8 r,theta,phi,h_plus,h_cross
    h_plus=-2d0*(1d0+dcos(iota)*dcos(iota))*dcos(2d0*phi+2d0*zeta)/r
    h_cross=-4d0*dcos(iota)*dsin(2d0*phi+2d0*zeta)/r
    return
    end subroutine