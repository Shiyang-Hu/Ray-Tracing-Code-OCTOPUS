subroutine initial_values(x_screen,y_screen,u)
use observation_parameters
implicit none
real*8,dimension(8)::u
real*8 x_screen,y_screen
real*8 x_bh,y_bh,z_bh
real*8 radius,t,r,theta,phi
real*8 r_dot,theta_dot,phi_dot
real*8 g_00,g_11,g_22,g_33
real*8 p_t,p_r,p_theta,p_phi,energy
u=0d0
!---------------------------------------------------------------------------------------------!
x_bh=(r_obs*dsin(theta_obs)-y_screen*dcos(theta_obs))*dcos(phi_obs)-x_screen*dsin(phi_obs)
y_bh=(r_obs*dsin(theta_obs)-y_screen*dcos(theta_obs))*dsin(phi_obs)+x_screen*dcos(phi_obs)
z_bh=r_obs*dcos(theta_obs)+y_screen*dsin(theta_obs)
!---------------------------------------------------------------------------------------------!
radius=x_bh*x_bh+y_bh*y_bh+z_bh*z_bh
r=dsqrt(radius)                                                
theta=dacos(z_bh/r)
phi=datan2(y_bh,x_bh)
t=t_obs
!---------------------------------------------------------------------------------------------!
r_dot=-(-dcos(theta_obs)*dcos(theta)-dsin(theta_obs)*dcos(phi-phi_obs)*dsin(theta))
theta_dot=(dsin(theta_obs)*dcos(phi-phi_obs)*dcos(theta)-dcos(theta_obs)*dsin(theta))/r
phi_dot=-dsin(theta_obs)*dsin(phi-phi_obs)/r/dsin(theta)
!---------------------------------------------------------------------------------------------!
call covariant_metric(r,theta,g_00,g_11,g_22,g_33)
p_r=g_11*r_dot
p_theta=g_22*theta_dot
p_phi=g_33*phi_dot
call initial_photon_energy(r,theta,p_r,p_theta,p_phi,energy)
!---------------------------------------------------------------------------------------------!
p_r=p_r/energy
p_theta=p_theta/energy
p_phi=p_phi/energy
p_t=-1d0
call dimensional(u,t,r,theta,phi,p_t,p_r,p_theta,p_phi)
return
end subroutine

subroutine initial_photon_energy(r,theta,p_r,p_theta,p_phi,energy)
implicit none
real*8 r,theta
real*8 p_r,p_theta,p_phi
real*8 energy
real*8 g00,g11,g22,g33
call contravariant_metric(r,theta,g00,g11,g22,g33)
energy=dsqrt((-g11*p_r*p_r-g22*p_theta*p_theta-g33*p_phi*p_phi)/g00)
return
end subroutine
    
subroutine particle_initialization(particle)
    use particle_parameters
    implicit none
    real*8,dimension(8)::particle
    real*8 hami
    hami=0d0
    select case(path_class)
    case(1)
        call initial_pr(particle)
    case(2)
        call initial_ptheta(particle)
    case(3)
        call initial_pt(particle)
    case(4)
        call circular_orbit(particle)
    case(5)
        call dimensional(particle,t_ini,r_ini,theta_ini,phi_ini,pt_ini,0d0,0d0,pphi_ini)
    end select
    write(*,*)"*******************************************************************************"
    write(*,*)"Initial condition of massive particle:"
    write(*,*)"t:",particle(1)
    write(*,*)"r:",particle(2)
    write(*,*)"theta:",particle(3)
    write(*,*)"phi:",particle(4)
    write(*,*)"Energy:",-particle(5)
    write(*,*)"Pr:",particle(6)
    write(*,*)"Ptheta:",particle(7)
    write(*,*)"Angular Momentum:",particle(8)
    call hamiltonian(particle,hami)
    write(*,*)"Particle's Hamiltonian:",hami
    return
    end subroutine
  
subroutine initial_pr(particle)
    use particle_parameters
    implicit none
    real*8,dimension(8)::particle
    real*8 p_theta,p_r
    real*8 g00,g11,g22,g33
    p_theta=0d0
    call contravariant_metric(r_ini,theta_ini,g00,g11,g22,g33)
    p_r=dsqrt((-1d0-g00*pt_ini*pt_ini-g33*pphi_ini*pphi_ini)/g11)
    call dimensional(particle,t_ini,r_ini,theta_ini,phi_ini,pt_ini,p_r,p_theta,pphi_ini)
    return
    end subroutine
    
subroutine initial_ptheta(particle)
    use particle_parameters
    implicit none
    real*8,dimension(8)::particle
    real*8 p_theta,p_r
    real*8 g00,g11,g22,g33
    p_r=0d0
    call contravariant_metric(r_ini,theta_ini,g00,g11,g22,g33)
    p_theta=dsqrt((-1d0-g00*pt_ini*pt_ini-g33*pphi_ini*pphi_ini)/g22)
    call dimensional(particle,t_ini,r_ini,theta_ini,phi_ini,pt_ini,p_r,p_theta,pphi_ini)
    return
    end subroutine    
    
subroutine initial_pt(particle)
    use particle_parameters
    implicit none
    real*8,dimension(8)::particle
    real*8 p_theta,p_r,p_t
    real*8 g00,g11,g22,g33
    p_r=0d0
    p_theta=0d0
    call contravariant_metric(r_ini,theta_ini,g00,g11,g22,g33)
    p_t=-dsqrt((-1d0-g33*pphi_ini*pphi_ini)/g00)
    call dimensional(particle,t_ini,r_ini,theta_ini,phi_ini,p_t,p_r,p_theta,pphi_ini)
    return
    end subroutine    
    
subroutine circular_orbit(particle)
    use particle_parameters
    implicit none
    real*8,dimension(8)::particle
    real*8 energy,angular_momentum
    call motion_constant(r_ini,energy,angular_momentum)
    call dimensional(particle,t_ini,r_ini,theta_ini,phi_ini,-energy,0d0,0d0,angular_momentum)
    return
    end subroutine
    
    
    
    
    
    
    
    
    
    
    
    
    