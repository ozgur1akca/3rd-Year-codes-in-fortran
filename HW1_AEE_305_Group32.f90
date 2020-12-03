!  hw1allinone.f90 
!
!  FUNCTIONS:
!  hw1allinone - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: hw1allinone
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
Module data
  parameter ( g = 9.81 )
  parameter ( u = 0.02 ) !friction
  parameter ( w_span = 16.25 )
  parameter ( w_parea = 29.24 )
  parameter ( Wto = 88250. )
  parameter ( Ta_Max = 16256. )
  parameter ( Ne = 2. )
  parameter ( cl_max = 1.792 )
  parameter (p0 = 1.225 )
  
    End module
    
    program hw1allinone
    use data
    character*40 fname
    real:: n,cd,vf,p

   print*, "1-estimating minimum time of taking off with respect to altitude and stepsize with using euler method (for question 1 and 2)" 
   print*, ""
   print*, "!if you want to calculate first question, then you have to type 1.225 for air density and change stepsize(0.125->0.250->0.50). "
   print*, ""
   print*, "!if you want to calculate second question, then you have to choose one stepsize(0.125) and type air density at 1000m(1.112) and 2000m (1.007)"
   print*, ""
   print*, "2-estimating minimum time of taking off with respect to altitude and stepsize with using runge kutta order of second method(for question 3)" 
   print*, ""
   print*, "!if you want to calculate third question, then you have to type 0.125 for stepsize and 1.225 for air density"
   n=0
   do while (n .lt. 3)
       
       print*,""
       print*,""
       print*, "type the number of the method"
       print*, ""
       print*, "if you want to stop the program type 3 or more"
       n=0
       read*, n
     
    if(n==1) then
           
        !..open the output file for euler method

        write(*,'(/,(a))',advance='no')'  Enter stepsize and air density  :> '
        read(*,*) stepsize,p
        write(*,'(a)',advance='no')'  Enter the output file name [velocity.dat]: '
        read(*,"(a)") fname
        if( fname .eq. " ") fname = "velocity.dat"
        open(1,file=fname,form="formatted")
   
        !..Set the Initial Conditions and output them
   
        time      = 0.
        velocity  = 0.
        bonus1 = 0.
        
        !find  final velocity 
   
        vf = (Wto*2/(cl_max*p*w_parea))**(0.5)   
        write(1,"(5f12.3)") time, velocity,bonus1,vf,p

        !.. Euler Solution loop for given altitude and stepsize
        do while ( velocity .lt. vf )
        time  = time + stepsize
        bonus = (velocity+(velocity + ODE(time,velocity,p)*stepsize))*stepsize*0.5 
        bonus1 = bonus1 + bonus
        velocity = velocity + ODE(time,velocity,p)*stepsize
        
        write(1,"(3f12.3)") time, velocity , bonus1
        enddo
        !..Close the output file for euler method
   
        close(1)
         
        elseif(n==2) then
         
        !..open the output file for RK2 method

        write(*,'(/,(a))',advance='no')'  Enter stepsize and air density  :> '
        read(*,*) stepsize,p
        write(*,'(a)',advance='no')'  Enter the output file name [velocity.dat]: '
        read(*,"(a)") fname
        if( fname .eq. " ") fname = "velocity.dat"
        open(1,file=fname,form="formatted")
         
        !..Set the Initial Conditions and output them
   
        time      = 0.
        velocity  = 0.
   
        !find  final velocity 
   
        vf = (Wto*2/(cl_max*p*w_parea))**(0.5)
        write(1,"(4f12.3)") time, velocity,vf,p

        !.. loop for velocity
   
        do while ( velocity .lt. vf )   
        velocity = rk2(time,velocity,stepsize,p)
        time  = time + stepsize
        write(1,"(2f12.3)") time, velocity 
        enddo
   
        !..Close the output file for RK2
        close(1)      
        endif
    enddo
   
    stop   
end program hw1allinone
    
    
    
    
    
!ode function
Function ODE(time,velocity,p)
  use data
  real :: T,cd,ode,L,D
  cd = 0.0207 + 0.0605* cl_max**2 !drag coefficient
  L = cl_max * 0.5 * p * (velocity**2) * w_parea !lift equation
  D = cd * 0.5 * p * (velocity**2) * w_parea !drag equation
  T = 2 * (Ta_max*p/p0)  !thrust with respect to air density
  ODE  = ((T-D-u*(Wto-L))*g)/Wto
  return
    End
    
!runge kutta second order method
    Function rk2(x1,y1,dx,p)
    real :: k1,k2,p1,a1,a2
    p1=0.5
    a1=0
    a2=1

    k1 = ode(x1,y1,p)
    k2 = ode(x1+p1*dx,y1+p1*dx*k1,p)

    rk2 = y1 + (k1*a1+k2*a2)*dx
    return
end