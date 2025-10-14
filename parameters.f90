!-------------------------------------------------------------------------------------------!
!In module metric_parameters, please enter the required parameters for the metric potential. 
!Please note that the variable math_pi should not be deleted.
!-------------------------------------------------------------------------------------------!
module metric_parameters
    implicit none
    real*8::r_s=0.5d0
    real*8::rho_s=0.5d0
    real*8::math_pi=dacos(-1d0)
    end module metric_parameters   
      
module black_hole_parameters
    implicit none
    real*8::mass_ratio=4.14d0*1d6
    real*8::distance_obs=8.127d0*1d-3
    end module black_hole_parameters    
    
module gravitational_parameters
    implicit none
    real*8::zeta=dacos(-1d0)/4d0
    real*8::iota=dacos(-1d0)/4d0
    real*8::eta=1d0
    real*8::light_speed=299792458d0
    real*8::gravitational_constant=6.69d0*1d-11
    real*8::sun_mass=1.989d0*1d30
    end module gravitational_parameters

module observation_parameters
    implicit none
    real*8::t_obs=0d0
    real*8::r_obs=1000d0
    real*8::theta_obs=50d0*dacos(-1d0)/180d0
    real*8::phi_obs=0d0*dacos(-1d0)/180d0
    integer::resolution_x=1000
    integer::resolution_y=1000
    real*8::x_ini=-15,x_end=15
    real*8::y_ini=-15,y_end=15
    end module observation_parameters

module raytracing_parameters
    implicit none
    real*8::hit_error=5d-4
    real*8::observer=1500d0
    real*8::ini_step=5d-4
    real*8::step_ratio=1.8d0
    end module raytracing_parameters

module accretion_disk_parameters
    implicit none
    real*8::r_inner=6d0
    real*8::r_outer=200d0
    integer::disk_model=1
    integer::max_hit=3
    end module accretion_disk_parameters
    
module lightsource_parameters
    implicit none
    real*8::source_x=-8d0,source_y=0d0,source_z=0d0
    real*8::source_radius=0.5d0
    real*8::j_0=1d0,alpha_0=1d0
    end module lightsource_parameters

module particle_parameters
    implicit none
    real*8::t_ini=0d0
    real*8::r_ini=15d0
    real*8::theta_ini=(90d0/180d0)*dacos(-1d0)
    real*8::phi_ini=(0d0/180d0)*dacos(-1d0)
    real*8::pt_ini=-0.972679d0
    real*8::pphi_ini=3.887002d0
    real*8::interval=-0.1d0
    integer::sample=50
    integer::counts=500
    integer::path_class=4
    end module particle_parameters
    
module task_parameters
    implicit none
    integer::task_model=8
    integer::cpu=128
    end module task_parameters