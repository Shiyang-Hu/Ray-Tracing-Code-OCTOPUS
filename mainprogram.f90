program RayTracing_Octopus
use omp_lib
use task_parameters
use observation_parameters
implicit none
character(len=30)::file1="shadow.dat",file2="disk profile.dat",file3="lensing.dat",file4="hit.dat",file5="redshift1.dat"
character(len=30)::file6="redshift2.dat",file7="redshift3.dat",file8="redshift4.dat",file9="image230.dat",file10="image86.dat"
character(len=30)::file11="hot-spot time delay.dat",file12="time sequence.dat",file13="particle path.dat",file14="hot-spot image.dat"
character(len=30)::file15="hamiltonian error.dat",file16="GW emission.dat"
real*8::r_eh=0d0,r_p=0d0,b_p=0d0,r_isco=0d0,angular_diameter=0d0
real*8::ini_t=0d0,over_t=0d0
real*8::delta_x=0d0,delta_y=0d0,x=0d0,y=0d0
integer i,j
open(1,file=file1)
open(2,file=file2)
open(3,file=file3)
open(4,file=file4)
open(5,file=file5)
open(6,file=file6)
open(7,file=file7)
open(8,file=file8)
open(9,file=file9)
open(10,file=file10)
open(11,file=file11)
open(12,file=file12)
open(13,file=file13)
open(14,file=file14)
open(15,file=file15)
open(16,file=file16)
ini_t=omp_get_wtime()
call black_hole_features(r_eh,r_p,b_p,r_isco,angular_diameter)
delta_x=dabs(x_end-x_ini)/(float(resolution_x)-1d0)
delta_y=dabs(y_end-y_ini)/(float(resolution_y)-1d0)
call omp_set_num_threads(cpu)
select case(task_model)
case(1)
    call light_curves(x,y,delta_x,delta_y,r_eh)
    goto 1
case(2)
    call lensing_animation(x,y,delta_x,delta_y,r_eh)
    goto 1
case(3)
    call massive_particle_path(r_eh)
    goto 1
end select
!$omp parallel do schedule(dynamic,1) default(none) firstprivate(x_ini,y_ini,delta_x,delta_y,resolution_x,resolution_y,x,y,r_eh,r_isco,r_p,task_model)
do i=1,resolution_y,1
    y=y_ini+(float(i)-1d0)*delta_y
    do j=1,resolution_x,1
        x=x_ini+(float(j)-1d0)*delta_x
        select case(task_model)
        case(4)
            call shadow(x,y,r_eh,r_p)
        case(5)
            call disk_reconstruction(x,y,r_eh,r_p,r_isco)
        case(6)
            call hit_times(x,y,r_eh,r_p,r_isco)
        case(7)
            call redshift_factor(x,y,r_eh,r_p,r_isco)
        case(8)
            call image(x,y,r_eh,r_p,r_isco)
        case(9)
            call lensing(x,y,r_eh)
        case(10)
            call lensing_image(x,y,r_eh)
        case(11)
            call hamiltonian_error(x,y,r_eh,r_p)
        end select
    end do
end do
!$omp end parallel do
1 over_t=omp_get_wtime()
write(*,*)"*******************************************************************************"
write(*,*)"Cpu Cost:",dabs(over_t-ini_t),"seconds"
write(*,*)"Program Completed"
close(1)
close(2)
close(3)
close(4)
close(5)
close(6)
close(7)
close(8)
close(9)
close(10)
close(11)
close(12)
close(13)
close(14)
close(15)
close(16)
stop
end 