!------------------------------------------------------------------------
!  RK4 SOLVER for a system of ODEs                                      |
!  Course:  AE305                                                       |
!------------------------------------------------------------------------
module data
integer, parameter :: neq=2



real*8,parameter :: rho_p=1140.
real*8,parameter ::  a=0.0000555
real*8,parameter ::  n=0.305
real*8,parameter :: r0=0.05
real*8,parameter :: rf=0.15
real*8,parameter :: L=1.25
real*8,parameter :: Tece=2810.
real*8,parameter :: g_c=365.
real*8,parameter :: sp_hr=1.25
real*8,parameter :: r_star=0.03
real*8,parameter :: p_a=101325.
real*8,parameter ::  g = 9.81d0
real*8,parameter :: pi = 3.141592


end module

    
program sysRK4
	use data
    integer tip
	real*8 :: y(neq),f(neq),rt,A_star,ximpals,E0,err,e_allowed,y_temp(neq),nstep,xama(neq),dtmin,dtmax,ytemp2(neq),ytemp(neq)
	character*40 :: fname
    write(*,'(/(a))',advance='no') "if you want to compile this code for first question type 1"
    write(*,'(/(a))',advance='no') "if you want to compile this code for second question type 2"
    write(*,'(/(a))',advance='no') "if you want to compile this code for third question type 3"
    write(*,'(/(a))',advance='no') ""
    
read*, tip
    
!for the first question
!----------------------------------------------------------------------------------- 
    if(tip.eq.1) then
!..read the input data
	write(*,'(/(a))',advance='no')' Enter the step size :>'
	read(*,*) dt
!..open the output file 
	write(*,'(a)',advance='no')' Enter the output file name [solution.dat]:>'
	read(*,'(a)') fname
	if( fname .eq. ' ') fname = 'solution.dat'
	open(1,file=fname,form='formatted')

!..set the initial conditions
	time = 0.
	y(1) = r0
 	y(2) = p_a
    ximpals = I_sp(y(2))
    A_star = 0.0028274
    nstep = 0
    
    call ODEs(A_star,time,y,f)
    
	write(1,'(6E14.5)') time,(y(i),i=1,neq),I,f(1)

!..solution loop
	DO WHILE (y(2).ge.p_a)   
        
		call SRK4(A_star,time,dt,y)
        
        call ODEs(A_star,time,y,f)
        
        ximpals = I_sp(y(2))
        
		time = time + dt 
        nstep = nstep + 1
        
        if(y(2).lt.p_a) exit
        
		write(1,'(7E14.5)') time,(y(i),i=1,neq),ximpals,f(1),nstep
	ENDDO
    
	close(1)
    
    !for the second question
    !-----------------------------------------------------------------------------------
    elseif(tip.eq.2) then
    !..read the input data
	write(*,'(/(a))',advance='no')' Enter the throat radius(0.02 ->0.03 ->0.04 ->0.05 ->0.06) :>'
	read(*,*) rt
!..open the output file 
	write(*,'(a)',advance='no')' Enter the output file name [solution.dat]:>'
	read(*,'(a)') fname
	if( fname .eq. ' ') fname = 'solution.dat'
	open(1,file=fname,form='formatted')

!..set the initial conditions
    dt = 0.0003
	time = 0.
	y(1) = r0
 	y(2) = p_a
    ximpals = I_sp(y(2))
    A_star = pi*rt**2
    call ODEs(A_star,time,y,f)
    
	write(1,'(6E14.5)') time,(y(i),i=1,neq),xI,f(1)

!..solution loop
	DO WHILE (y(2).ge.p_a) 
        
          
		call SRK4(A_star,time,dt,y)
        
        call ODEs(A_star,time,y,f)
        
        
        ximpals = I_sp(y(2))
                  
		time = time + dt 
        
        if(y(2).lt.p_a) exit
        
		write(1,'(6E14.5)') time,(y(i),i=1,neq),ximpals,f(1)
    ENDDO

    
	close(1)

!for the third question
!-------------------------------------------------------------------------------
    elseif(tip.eq.3) then
!..read the input data
!..open the output file 
	write(*,'(a)',advance='no')' Enter the output file name [solution.dat]:>'
	read(*,'(a)') fname
	if( fname .eq. ' ') fname = 'solution.dat'
	open(1,file=fname,form='formatted')


!..set the initial conditions
    e_allowed =1e-5
    dt = 0.0003
	time = 0.
	y(1) = r0
 	y(2) = p_a
    ximpals = I_sp(y(2))
    A_star = pi*0.03**2
    nstep = 0
    err =0
    E0 = 0
    dtmax = 0.025
    dtmin = 0.000003

    call ODEs(A_star,time,y,f)

	write(1,'(6E14.5)') time,(y(i),i=1,neq),I,f(1)


!..solution loop
    DO WHILE (y(2).ge.p_a)
        
        do i = 1,neq
        ytemp(i) = y(i)
        enddo

		call SRK4(A_star,time,dt,y)
        call ODEs(A_star,time,y,f)
        
        do j=1,neq
        ytemp2(j) = y(j)
        enddo
        
        
        
        ximpals = I_sp(y(2))
        
		time = time + dt 
        
        nstep = nstep + 1
        
        if(y(2).lt.p_a) exit
        
		write(1,'(7E14.5)') time,(y(i),i=1,neq),ximpals,f(1),nstep,dt
        
!----------adaptive stepping with using rk2 and rk4 method
               
        do k = 1,neq
        y(k) = ytemp(k)
        enddo
        
   
        call SRK2(a_star,time,dt,y)
        
        do m =1,neq
        xama(m) = y(m)
        enddo
        
        do i = 1,neq 
        y(i) = ytemp(i)
        enddo
        
        call SRK4(a_star,time,dt,y)
        
        err = abs(e_allowed/((y(2)-xama(2))/y(2))) 
       
        
        dt = dt*err**(0.2)
        
        if(dt.lt.dtmin) dt = dtmin
        if(dt.gt.dtmax) dt = dtmax
        
        do i = 1,neq
        y(i) = ytemp2(i)
        enddo


	ENDDO

	close(1)
    
        
        
!---------------------------------------------------------------------------
    else
        
        
        endif
	stop
    end
	
!------------------------------------------------------------------------
!subroutine for rk4 method
subroutine SRK4(A_star,time,dt,y)
      use data
	real*8 :: y(neq),ytmp(neq),k1(neq),k2(neq),k3(neq),k4(neq),f(neq),A_star
    
	dt2 = 0.5*dt
!by using ODEs we calculate K1(2),K2(2),K3(2),K4(2)
	call odes(A_star,time,y,k1)
	do i = 1,neq
         ytmp(i)  = y(i) + k1(i)*dt2
    enddo
    
 	call ODES(A_star,time+dt2,ytmp,k2)
    	do i = 1,neq
         ytmp(i)  = y(i) + k2(i)*dt2
        enddo
        
    call ODES(A_star,time+dt2,ytmp,k3)
    	do i = 1,neq
         ytmp(i)  = y(i) + k3(i)*dt
        enddo
        
        call ODES(A_star,time+dt2,ytmp,k4)
!..obtain the solution at t+dt and update y for the next step
	do i = 1,neq
		phi  = (k1(i) + 2*(k2(i)+k3(i)) + k4(i))/6.
		y(i) = y(i) + phi*dt
	enddo

	return
    end subroutine
    
    
!--------------------------------------------------------------------
! subroutine for rk2 method
    subroutine SRK2(A_star,time,dt,y)
    use data
    real*8 :: y(neq),ytmp(neq),k1(neq),k2(neq),A_star,p1,a1,a2
    p1 = 0.5
    a1 = 0.
    a2 = 1.0
    !by using ODEs we calculate K1(2),K2(2)
    call odes(A_star,time,y,k1)
    
	do i = 1,neq
         ytmp(i)  = y(i) + k1(i)*p1*dt
    enddo
 	call ODES(A_star,time+dt2,ytmp,k2)
    !..obtain the solution at t+dt and update y for the next step
    do i= 1,neq
    y(i) = y(i) +k2(i)*dt
    enddo
    
    return
    end subroutine
    
!------------------------------------------------------------------------
!subroutine for ODEs
subroutine ODEs(A_star,time,y,f)
      use data
	real*8 :: y(neq),f(neq),A_star

!..define the ODE's & return the slopes in the "f" array
    if(y(1).gt.rf)then  !I wrote an "if statment" here because we have limited radius and after reaching the final radius there is no burn rate 
        y(1) = rf
        f(1) = 0
    else
    f(1) = a*(y(2)**n)
    endif
    
 	f(2) = g_c*Tece*((f_cor(y(1))*2*a*(y(2)**n)/y(1)*(rho_p-(y(2)/(g_c*Tece))))-((y(2)*A_star)/(pi*(y(1)**2)*L)*sqrt(sp_hr/(g_c*Tece))*(((sp_hr+1)/2)**(-(sp_hr+1)/(2*(sp_hr-1))))))

	return
    end
!---------------------------------------------------------------------------------------------------
!function for correction factor
    function f_cor(r)
   ! perimeter factor
    Use  data
    real*8 :: eta,r
        eta =  (rf-r)/(rf-r0)
         if( eta .le. 0. )then
         f_cor = 0.
         elseif( eta .gt. 0.15)then
         f_cor=1.
        else
         f_cor=1- exp(-7*eta)
         endif
    end function
!--------------------------------------------------------------------------------------------------------    
!function for specific impulse
 function I_sp(p_c)
  use data
  real*8 :: p_c
   ! specific impulse
       I_sp = (1/g)*sqrt(((2*sp_hr*g_c*Tece)/(sp_hr-1))*(1-(p_a/p_c)**((sp_hr-1)/sp_hr)))
    end function
    