subroutine black_hole_features(r_eh,r_p,b_p,r_isco,angular_diameter)
    implicit none
    real*8 r_eh,r_p,b_p,r_isco,angular_diameter
    call event_horizon(r_eh)
    r_p=r_eh+0.5d0
    call photon_sphere(r_p)
    r_isco=r_p+1d0
    call isco_orbit(r_isco)
    call impact_parameter(r_p,b_p)
    call angular_obs(b_p,angular_diameter)
    write(*,*)"*******************************************************************************"
    write(*,*)"Radius of the Event Horizon:",r_eh,"M"
    write(*,*)"Radius of the Photon Sphere:",r_p,"M"
    write(*,*)"Radius of the Black Hole Shadow (critical curve):",b_p,"M"
    write(*,*)"Diameter of the Black Hole Shadow in uas:",angular_diameter,"uas"
    write(*,*)"Radius of the Innermost Stable Circular Orbit:",r_isco,"M"
    return
    end subroutine
    
subroutine shadow(x_screen,y_screen,r_eh,r_p)
    use observation_parameters
    use raytracing_parameters
    implicit none
    real*8 r_eh,r_p
    real*8 x_screen,y_screen
    real*8,dimension(8)::vector
    real*8 step
    call initial_values(x_screen,y_screen,vector)
    do while(.true.)
        step=ini_step*(vector(2)/r_eh)**step_ratio
        !step=ini_step*(vector(2)/r_p)**step_ratio
        call rkf(vector,step)
        !call adjust_angle(vector(3),vector(4))
        if(vector(2).le.(r_eh+hit_error))then
        !if(vector(2).le.(r_p+hit_error))then
            write(1,*)x_screen,y_screen
            exit
        else
        end if
        if(vector(2).ge.observer)exit
    end do
    return
    end subroutine
      
subroutine hamiltonian_error(x_screen,y_screen,r_eh,r_p)
    use observation_parameters
    use raytracing_parameters
    implicit none
    real*8 r_eh,r_p
    real*8 x_screen,y_screen
    real*8,dimension(8)::vector,vector0
    real*8 step
    real*8 hami
    call initial_values(x_screen,y_screen,vector)
    do while(.true.)
        step=ini_step*(vector(2)/r_eh)**step_ratio
        call record_values(vector,vector0)
        call rkf(vector,step)
        !call adjust_angle(vector(3),vector(4))
        if(vector(2).le.(r_eh+hit_error))then
            write(1,*)x_screen,y_screen
            vector=vector0
            exit
        else
        end if
        if(vector(2).ge.observer)exit
    end do
    call hamiltonian(vector,hami)
    write(15,*)x_screen,y_screen,dlog10(dabs(hami)+1d-16)
    return
    end subroutine
    
subroutine disk_reconstruction(x_screen,y_screen,r_eh,r_p,r_isco)
    use observation_parameters
    use raytracing_parameters
    use accretion_disk_parameters
    use metric_parameters
    implicit none
    real*8 r_eh,r_p,r_isco
    real*8 x_screen,y_screen
    real*8,dimension(8)::vector,vector0
    real*8 step
    call initial_values(x_screen,y_screen,vector)
    select case(disk_model)
    case(1)
        r_inner=r_eh
    case(2)
        r_inner=r_p
    case(3)
        r_inner=r_isco
    case(4)
        r_inner=r_inner
    end select
    do while(.true.)
        step=ini_step*(vector(2)/r_eh)**step_ratio
        call record_values(vector,vector0)
        call rkf(vector,step)
        call adjust_angle(vector(3),vector(4))
        if(vector(2).le.(r_eh+hit_error))then
            write(1,*)x_screen,y_screen
            exit
        else if(vector(2).ge.observer)then
            exit
        else
        end if
        if(dcos(vector(3))*dcos(vector0(3)).le.0d0.and.vector(2).le.r_outer.and.vector(2).ge.r_inner)then
            write(2,*)x_screen,y_screen
            exit
        else
        end if
    end do
    return
    end subroutine
        
subroutine lensing(x_screen,y_screen,r_eh)
    use observation_parameters
    use raytracing_parameters
    use lightsource_parameters
    implicit none
    real*8 x_screen,y_screen,r_eh
    real*8,dimension(8)::vector,vector0
    real*8 step
    real*8 radius,radius0
    real*8 x,y,z,x0,y0,z0
    call initial_values(x_screen,y_screen,vector)
    do while(.true.)
        step=ini_step*(vector(2)/r_eh)**step_ratio
        call record_values(vector,vector0)
        call rkf(vector,step)
        call adjust_angle(vector(3),vector(4))
        if(vector(2).le.(r_eh+hit_error))then
            write(1,*)x_screen,y_screen
            write(3,*)x_screen,y_screen,0
            exit
        else 
        end if
        call spherical_transfer(vector,x,y,z)
        call spherical_transfer(vector0,x0,y0,z0)
        radius=dsqrt((x-source_x)**2d0+(y-source_y)**2d0+(z-source_z)**2d0)   
        radius0=dsqrt((x0-source_x)**2d0+(y0-source_y)**2d0+(z0-source_z)**2d0)
        if(radius.ge.source_radius.and.radius0.le.source_radius)then
            write(3,*)x_screen,y_screen,5
            exit
        else
        end if
        if(vector(2).ge.observer)then
            if(z.gt.0d0.and.y.gt.0d0)write(3,*)x_screen,y_screen,1
            if(z.gt.0d0.and.y.lt.0d0)write(3,*)x_screen,y_screen,2
            if(z.lt.0d0.and.y.gt.0d0)write(3,*)x_screen,y_screen,3
            if(z.lt.0d0.and.y.lt.0d0)write(3,*)x_screen,y_screen,4
            exit
        else
        end if
    end do
    return
    end subroutine
        
subroutine hit_times(x_screen,y_screen,r_eh,r_p,r_isco)
    use raytracing_parameters
    use accretion_disk_parameters
    implicit none
    real*8 x_screen,y_screen,r_eh,r_p,r_isco
    real*8,dimension(8)::vector,vector0
    real*8 step
    integer hit_time
    hit_time=0
    call initial_values(x_screen,y_screen,vector)
    select case(disk_model)
    case(1)
        r_inner=r_eh
    case(2)
        r_inner=r_p
    case(3)
        r_inner=r_isco
    case(4)
        r_inner=r_inner
    end select
    do while(.true.)
        step=ini_step*(vector(2)/r_eh)**step_ratio
        call record_values(vector,vector0)
        call rkf(vector,step)
        !call adjust_angle(vector(3),vector(4))
        if(vector(2).le.(r_eh+hit_error))exit
        if(vector(2).ge.observer)exit
        if(dcos(vector(3))*dcos(vector0(3)).le.0d0.and.vector(2).le.r_outer.and.vector(2).ge.r_inner)then
            hit_time=hit_time+1
        else
        end if
        if(hit_time.ge.max_hit)exit
    end do
    write(4,*)x_screen,y_screen,hit_time
    return
    end subroutine
    
subroutine redshift_factor(x_screen,y_screen,r_eh,r_p,r_isco)
    use metric_parameters
    use raytracing_parameters
    use accretion_disk_parameters
    implicit none
    real*8 x_screen,y_screen,r_eh,r_p,r_isco
    real*8,dimension(8)::vector,vector0
    real*8,dimension(max_hit)::rs,redshift
    real*8 g_00,g_11,g_22,g_33
    real*8 e,l,e_isco,l_isco
    real*8 t_dot,r_dot,phi_dot
    real*8 step
    integer i,hit
    integer file_name
    hit=0
    call initial_values(x_screen,y_screen,vector)
    select case(disk_model)
    case(1)
        r_inner=r_eh
    case(2)
        r_inner=r_p
    case(3)
        r_inner=r_isco
    case(4)
        r_inner=r_inner
    end select
    do while(.true.)
        step=ini_step*(vector(2)/r_eh)**step_ratio
        call record_values(vector,vector0)
        call rkf(vector,step)
        call adjust_angle(vector(3),vector(4))
        if(vector(2).le.(r_eh+hit_error))exit
        if(vector(2).ge.observer)exit
        if(dcos(vector(3))*dcos(vector0(3)).le.0d0.and.vector(2).le.r_outer.and.vector(2).gt.r_inner)then
            hit=hit+1
            rs(hit)=(vector(2)+vector0(2))/2d0
            if(rs(hit).ge.r_isco)then
                call motion_constant(rs(hit),e,l)
                call covariant_metric(rs(hit),math_pi/2d0,g_00,g_11,g_22,g_33)
                t_dot=-e/g_00
                phi_dot=l/g_33
                redshift(hit)=vector(5)/(vector(5)*t_dot+vector(8)*phi_dot)
            else
                call motion_constant(r_isco,e_isco,l_isco)
                if(dabs(r_isco-rs(hit)).ge.0.001d0)then
                    call infalling_velocity(rs(hit),math_pi/2d0,e_isco,l_isco,r_dot)
                else
                    r_dot=0d0
                end if
                call covariant_metric(rs(hit),math_pi/2d0,g_00,g_11,g_22,g_33)
                t_dot=-e_isco/g_00
                phi_dot=l_isco/g_33
                redshift(hit)=vector(5)/(vector(5)*t_dot+(vector(6)+vector0(6))*r_dot/2d0+vector(8)*phi_dot)
            end if
        else
        end if
        if(hit.ge.max_hit)exit
    end do
    if(hit==0)then
        write(1,*)x_screen,y_screen,-1
    else
        do i=1,hit,1
            file_name=i+4
            write(file_name,*)x_screen,y_screen,redshift(i)
        end do
    end if
    return
    end subroutine
        
subroutine image(x_screen,y_screen,r_eh,r_p,r_isco)
    use metric_parameters
    use raytracing_parameters
    use accretion_disk_parameters
    implicit none
    real*8 x_screen,y_screen,r_eh,r_p,r_isco
    real*8,dimension(8)::vector,vector0
    real*8,dimension(max_hit)::rs,redshift,flux1,flux2
    real*8 intensity1,intensity2
    real*8 g_00,g_11,g_22,g_33
    real*8 e,l,e_isco,l_isco
    real*8 t_dot,r_dot,phi_dot
    real*8 step
    integer i,hit
    integer file_name
    hit=0
    flux1=0d0
    flux2=0d0
    call initial_values(x_screen,y_screen,vector)
    select case(disk_model)
    case(1)
        r_inner=r_eh
    case(2)
        r_inner=r_p
    case(3)
        r_inner=r_isco
    case(4)
        r_inner=r_inner
    end select
    do while(.true.)
        step=ini_step*(vector(2)/r_eh)**step_ratio
        call record_values(vector,vector0)
        call rkf(vector,step)
        call adjust_angle(vector(3),vector(4))
        if(vector(2).le.(r_eh+hit_error))exit
        if(vector(2).ge.observer)exit
        if(dcos(vector(3))*dcos(vector0(3)).le.0d0.and.vector(2).le.r_outer.and.vector(2).gt.r_inner)then
            hit=hit+1
            rs(hit)=(vector(2)+vector0(2))/2d0
            call energy_flux(rs(hit),r_eh,flux1(hit),flux2(hit))
            if(rs(hit).ge.r_isco)then
                call motion_constant(rs(hit),e,l)
                call covariant_metric(rs(hit),math_pi/2d0,g_00,g_11,g_22,g_33)
                t_dot=-e/g_00
                phi_dot=l/g_33
                redshift(hit)=vector(5)/(vector(5)*t_dot+vector(8)*phi_dot)
            else
                call motion_constant(r_isco,e_isco,l_isco)
                if(dabs(r_isco-rs(hit)).ge.0.001d0)then
                    call infalling_velocity(rs(hit),math_pi/2d0,e_isco,l_isco,r_dot)
                else
                    r_dot=0d0
                end if
                call covariant_metric(rs(hit),math_pi/2d0,g_00,g_11,g_22,g_33)
                t_dot=-e_isco/g_00
                phi_dot=l_isco/g_33
                redshift(hit)=vector(5)/(vector(5)*t_dot+(vector(6)+vector0(6))*r_dot/2d0+vector(8)*phi_dot)
            end if
        else
        end if
        if(hit.ge.max_hit)exit
    end do
    if(hit==0)then
        write(1,*)x_screen,y_screen,-1
    else
        do i=1,hit,1
            file_name=i+4
            write(file_name,*)x_screen,y_screen,redshift(i)
        end do
    end if
    select case(hit)
    case(0)
        intensity1=0d0
        intensity2=0d0
    case(1)
        intensity1=flux1(1)*redshift(1)**3d0
        intensity2=flux2(1)*redshift(1)**3d0
    case(2)
        intensity1=flux1(1)*redshift(1)**3d0+(2d0/3d0)*flux1(2)*redshift(2)**3d0
        intensity2=flux2(1)*redshift(1)**3d0+(2d0/3d0)*flux2(2)*redshift(2)**3d0
    case(3)
        intensity1=flux1(1)*redshift(1)**3d0+(2d0/3d0)*flux1(2)*redshift(2)**3d0+(2d0/3d0)*flux1(3)*redshift(3)**3d0
        intensity2=flux2(1)*redshift(1)**3d0+(2d0/3d0)*flux2(2)*redshift(2)**3d0+(2d0/3d0)*flux2(3)*redshift(3)**3d0
    case(4)
        intensity1=flux1(1)*redshift(1)**3d0+(2d0/3d0)*flux1(2)*redshift(2)**3d0+(2d0/3d0)*flux1(3)*redshift(3)**3d0+(2d0/3d0)*flux1(4)*redshift(4)**3d0
        intensity2=flux2(1)*redshift(1)**3d0+(2d0/3d0)*flux2(2)*redshift(2)**3d0+(2d0/3d0)*flux2(3)*redshift(3)**3d0+(2d0/3d0)*flux2(4)*redshift(4)**3d0
    end select
    write(9,*)x_screen,y_screen,intensity1
    write(10,*)x_screen,y_screen,intensity2
    write(4,*)x_screen,y_screen,hit
    return
    end subroutine                
            
subroutine light_curves(x,y,delta_x,delta_y,r_eh)
    use particle_parameters
    use observation_parameters
    use raytracing_parameters
    implicit none
    real*8,dimension(counts,4)::path_data
    real*8,dimension(counts)::radius
    real*8,dimension(8)::particle,particle_0
    real*8 x,y,delta_x,delta_y,r_eh
    real*8 t_int
    real*8 x_m,y_m,z_m
    real*8 h_plus,h_cross
    real*8 r_max,r_min
    real*8 period
    integer orbit
    integer i,j,k,l
    h_plus=0d0
    h_cross=0d0
    period=0d0
    orbit=0
    path_data=0d0
    radius=0d0
    particle_0=0d0
    t_int=0d0
    call particle_initialization(particle)
    !if(path_class==4)then
    !    period=2d0*dacos(-1d0)*r_ini*r_ini/particle(8)
    !    interval=-2d0*period/(float(counts)-1d0)
    !    sample=1
    !else
    !end if    
    do i=1,counts,1
        do j=1,sample,1
            call record_values(particle,particle_0)
            call timestep(t_int,interval)
            call rksixth(particle,interval)
            if(particle(2).le.(r_eh+hit_error))then
                write(*,*)"*********************************************************************"
                write(*,*)"WARNING"
                write(*,*)"The hot-spot is approaching event horizon!"
                write(*,*)"Calculation stop!"
                write(*,*)"Please change the initial conditions of the hot-spot and try again."
                stop
            else
            end if
            call spherical_transfer(particle,x_m,y_m,z_m)
        end do
        call gravitational_waves(particle(2),particle(3),particle(4),h_plus,h_cross)
        write(13,*)x_m,y_m,z_m
        write(16,*)dabs(t_int),h_plus,h_cross
        path_data(i,1)=particle(1)
        path_data(i,2)=x_m
        path_data(i,3)=y_m
        path_data(i,4)=z_m
        radius(i)=particle(2)
        orbit=orbit+1
    end do
    r_max=maxval(radius)+1d0
    r_min=minval(radius)-1d0
!$omp parallel do schedule(dynamic,1) default(none) firstprivate(x_ini,y_ini,delta_x,delta_y,resolution_x,resolution_y,x,y,r_eh) shared(path_data,r_max,r_min)
    do k=1,resolution_y,1
        y=y_ini+delta_y*(float(k)-1d0)
        do l=1,resolution_x,1
            x=x_ini+delta_x*(float(l)-1d0)
            call geodesic_tracing(x,y,r_eh,path_data,r_max,r_min)
        end do
    end do
!$omp end parallel do
    return
    end subroutine
        
subroutine geodesic_tracing(x,y,r_eh,path_data,r_max,r_min) 
    use raytracing_parameters
    use lightsource_parameters
    use particle_parameters
    implicit none
    real*8,dimension(counts,4)::path_data
    real*8,dimension(8)::photon,photon_0
    real*8 x,y,r_eh,r_max,r_min
    real*8 step
    real*8 x_p,y_p,z_p
    real*8 x_p0,y_p0,z_p0
    real*8 r_1,r_0
    integer hits
    integer m
    hits=0
    photon_0=0d0
    call initial_values(x,y,photon)
    do while(.true.)
        step=ini_step*(photon(2)/r_eh)**step_ratio
        call record_values(photon,photon_0)
        call rkf(photon,step)
        if(photon(2).le.(r_eh+hit_error))exit
        if(photon(2).ge.observer)exit
        call spherical_transfer(photon,x_p,y_p,z_p)
        call spherical_transfer(photon_0,x_p0,y_p0,z_p0)
        if(photon(2).le.r_max.and.photon(2).ge.r_min)then
            do m=1,counts,1
                r_0=dsqrt((x_p0-path_data(m,2))**2d0+(y_p0-path_data(m,3))**2d0+(z_p0-path_data(m,4))**2d0)
                r_1=dsqrt((x_p-path_data(m,2))**2d0+(y_p-path_data(m,3))**2d0+(z_p-path_data(m,4))**2d0)
                if(r_0.le.source_radius.and.r_1.ge.source_radius)then
                    hits=hits+1
                    !write(*,*)"photon + 1"
                    if(hits==1)write(11,*)x,y,dabs(path_data(m,1))+dabs(photon(1))
                    write(12,*)dabs(path_data(m,1))+dabs(photon(1))
                else
                end if
            end do
        else
        end if
    end do
    write(4,*)x,y,hits
    return
    end subroutine
    
subroutine lensing_image(x_screen,y_screen,r_eh)
    use observation_parameters
    use raytracing_parameters
    use lightsource_parameters
    implicit none
    real*8 x_screen,y_screen,r_eh
    real*8,dimension(8)::vector,vector0
    real*8,dimension(5000,4)::path_data
    real*8 step,step_0
    real*8 radius
    real*8 x_p,y_p,z_p
    real*8 luminosity
    real*8 j_v,alpha_v
    integer orbit
    integer i
    orbit=0
    luminosity=0d0
    path_data=0d0
    call initial_values(x_screen,y_screen,vector)
    do while(.true.)
        step_0=step
        step=ini_step*(vector(2)/r_eh)**step_ratio
        call record_values(vector,vector0)
        call rkf(vector,step)
        if(vector(2).le.(r_eh+hit_error))exit
        call spherical_transfer(vector,x_p,y_p,z_p)
        radius=dsqrt((x_p-source_x)**2d0+(y_p-source_y)**2d0+(z_p-source_z)**2d0)
        if(radius.le.source_radius)then
            orbit=orbit+1
            path_data(orbit,1)=x_p
            path_data(orbit,2)=y_p
            path_data(orbit,3)=z_p
            path_data(orbit,4)=step_0
        else
        end if
        if(vector(2).ge.observer)exit   
    end do
    !write(*,*)x_screen,y_screen,orbit
    if(orbit.ge.1)then
        do i=1,orbit,1
            call source_profile(path_data(i,1),path_data(i,2),path_data(i,3),source_x,source_y,source_z,j_v,alpha_v)
            luminosity=luminosity+path_data(i,4)*(j_v-alpha_v*luminosity)
        end do
        write(14,*)x_screen,y_screen,luminosity
    else
        write(14,*)x_screen,y_screen,0d0
    end if
    return
    end subroutine
    
subroutine lensing_animation(x,y,delta_x,delta_y,r_eh)
    use observation_parameters
    use particle_parameters
    use raytracing_parameters
    implicit none
    real*8,dimension(8)::particle,particle_0
    real*8 x,y,delta_x,delta_y,r_eh
    real*8 t_int
    real*8 x_m,y_m,z_m
    integer i,j,k,l
    t_int=0d0
    call particle_initialization(particle)
    do i=1,counts,1
        do j=1,sample,1
            call record_values(particle,particle_0)
            call timestep(t_int,interval)
            call rksixth(particle,interval)
            if(particle(2).le.(r_eh+hit_error))then
                write(*,*)"*********************************************************************"
                write(*,*)"WARNING"
                write(*,*)"The hot-spot is approaching event horizon!"
                write(*,*)"Calculation stop!"
                write(*,*)"Please change the initial conditions of the hot-spot and try again."
                stop
            else
            end if
            call spherical_transfer(particle,x_m,y_m,z_m)
        end do
        write(13,*)x_m,y_m,z_m
!$omp parallel do schedule(dynamic,1) default(none) firstprivate(x_ini,y_ini,delta_x,delta_y,resolution_x,resolution_y,x,y,r_eh) shared(particle,x_m,y_m,z_m)
    do k=1,resolution_y,1
        y=y_ini+delta_y*(float(k)-1d0)
        do l=1,resolution_x,1
            x=x_ini+delta_x*(float(l)-1d0)
            call tracing(x,y,r_eh,particle,x_m,y_m,z_m)
        end do
    end do
!$omp end parallel do
write(*,*)i
!pause
    end do
    return
    end subroutine
    
subroutine tracing(x,y,r_eh,particle,x_m,y_m,z_m)
    use raytracing_parameters
    use lightsource_parameters
    implicit none
    real*8,dimension(1000,9)::path_data   
    real*8,dimension(8)::vector,vector0,particle
    real*8 x,y,r_eh
    real*8 x_p,y_p,z_p,x_m,y_m,z_m
    real*8 radius
    real*8 g_00,g_11,g_22,g_33
    real*8 r_s,theta_s
    real*8 t_velocity,r_velocity,theta_velocity,phi_velocity
    real*8 j_v,alpha_v
    real*8 redshift
    real*8 step
    real*8 luminosity
    integer i,orbit
    orbit=0
    luminosity=0d0
    path_data=0d0
    call initial_values(x,y,vector)
    do while(.true.)
        step=ini_step*(vector(2)/r_eh)**step_ratio
        call record_values(vector,vector0)
        call rkf(vector,step)
        if(vector(2).le.(r_eh+hit_error))exit
        call spherical_transfer(vector,x_p,y_p,z_p)
        radius=dsqrt((x_p-x_m)**2d0+(y_p-y_m)**2d0+(z_p-z_m)**2d0)
        if(radius.le.source_radius)then
            orbit=orbit+1
            path_data(orbit,1)=x_p
            path_data(orbit,2)=y_p
            path_data(orbit,3)=z_p
            path_data(orbit,4)=step
            path_data(orbit,5)=vector(5)
            path_data(orbit,6)=vector(6)
            path_data(orbit,7)=vector(7)
            path_data(orbit,8)=vector(8)
            path_data(orbit,9)=dabs(vector(1))
        else
        end if
        if(vector(2).ge.observer)exit
    end do
    if(orbit.ge.1)then
        do i=1,orbit,1
            call source_profile(path_data(i,1),path_data(i,2),path_data(i,3),x_m,y_m,z_m,j_v,alpha_v)
            r_s=particle(2)
            theta_s=particle(3)
            call covariant_metric(r_s,theta_s,g_00,g_11,g_22,g_33)
            t_velocity=particle(5)/g_00
            r_velocity=particle(6)/g_11
            theta_velocity=particle(7)/g_22
            phi_velocity=particle(8)/g_33
            redshift=path_data(i,5)/(path_data(i,5)*t_velocity+path_data(i,6)*r_velocity+path_data(i,7)*theta_velocity+path_data(i,8)*phi_velocity)
            luminosity=luminosity+path_data(i,4)*(j_v-alpha_v*luminosity)*redshift*redshift*redshift
        end do
        write(14,*)x,y,luminosity
    else
        write(14,*)x,y,0d0
    end if   
    return
    end subroutine
    
subroutine massive_particle_path(r_eh)    
    use particle_parameters
    use raytracing_parameters
    implicit none
    real*8,dimension(8)::particle
    real*8 r_eh
    real*8 period
    real*8 t_int
    real*8 x_m,y_m,z_m
    real*8 h_plus,h_cross
    integer i,j
    period=0d0
    t_int=0d0
    call particle_initialization(particle)    
    do i=1,counts,1
        do j=1,sample,1
            call timestep(t_int,interval)
            call rksixth(particle,interval)
            if(particle(2).le.(r_eh+hit_error))then
                write(*,*)"*********************************************************************"
                write(*,*)"WARNING"
                write(*,*)"The hot-spot is approaching event horizon!"
                write(*,*)"Calculation stop!"
                write(*,*)"Please change the initial conditions of the hot-spot and try again."
                stop
            else
            end if
            call spherical_transfer(particle,x_m,y_m,z_m)
        end do
        call gravitational_waves(particle(2),particle(3),particle(4),h_plus,h_cross)
        write(13,*)x_m,y_m,z_m
        write(16,*)dabs(t_int),h_plus,h_cross
    end do
    return
    end subroutine