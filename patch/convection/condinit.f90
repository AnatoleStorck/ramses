!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,pert,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp)::pert                          ! perturbations on/off
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::i,id,iu,iv,iw,ip,is                     ! Indices Euler
  real(dp)::p1,p2,p3,rho1,rho2,rho3,x1,x2,g,T0,drho ! Variables
  real(dp)::gammainit1,gammainit2,gammainit3
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  id=1; iu=2; iv=3; iw=4; ip=ndim+2; is=ndim+3
  x1=x_center(1); x2=x_center(2)
  rho1=d_region(1); rho2=d_region(2); rho3=d_region(3)
  p1=p_region(1); p2=p_region(2); p3=p_region(3)
  gammainit1=gamma_region(1)
  gammainit2=gamma_region(2)
  gammainit3=gamma_region(3)
  g = abs(gravity_params(1))

  do i=1,nn
    if(x(i,1) .le. x1)then
      T0 = 1 - ((gammainit1-1.0d0)/gammainit1)*g*(rho1/p1)*(x(i,1))
      q(i,id)=rho1*((T0)**(1.0/(gammainit1-1.0d0)))
      q(i,ip)=p1*((T0)**(gammainit1/(gammainit1-1.0d0)))
    else if ((x(i,1) .gt. x1) .and. (x(i,1) .lt. x2)) then
      T0 = 1 - ((gammainit2-1.0d0)/gammainit2)*g*(rho2/p2)*((x(i,1)-x1))
      drho = 0.0d0
      if ((pert_r .gt. 0.) .and. (pert_dx .gt. 0.)) then
        if (x(i,1) .lt. x1+pert_dx) then
          !! produce density perturbation in small layer of convection zone!
          call random_number(drho)
          drho = pert_r*pert*2.0*(drho-0.5)*10.0**(-2.0)
        end if 
      end if
      q(i,id)=rho2*((T0)**(1.0/(gammainit2-1.0d0)))*(1.0-drho)
      q(i,ip)=p2*((T0)**(gammainit2/(gammainit2-1.0d0)))
    else 
      T0 = 1 - ((gammainit3-1.0d0)/gammainit3)*g*(rho3/p3)*((x(i,1)-x2))
      q(i,id)=rho3*((T0)**(1.0/(gammainit3-1.0d0)))
      q(i,ip)=p3*((T0)**(gammainit3/(gammainit3-1.0d0)))
    endif
    q(i,iu)=0.0d0
    if(ndim>1)q(i,iv)=0.0d0
    if(ndim>2)q(i,iw)=0.0d0
#if NVAR>NDIM+2
    ! Set entropy as a passive scalar
    q(i,is)=q(i,ip)/q(i,id)**gamma
#endif
 end do

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif

  ! pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)

  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================

subroutine eneana(x,e,dx,t,ncell)

  use amr_parameters
  use hydro_parameters

  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector)::e        ! Energy 
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! Heat bottom of convection zone with 
  ! e = 2*10**10 erg/g/s * rho0       for HeFlash
  ! e = 2.489*10.0**12 erg/g/s * rho0 for Model S
  ! where rho0 is the density at the bottom of convection zone
  !
  ! We heat the first layer of 
  !   2000 km for Model S
  !   500 km  for HeFlash
  ! of the convection zone
  !================================================================
  integer :: i, idim
  real(dp):: rho0, e0, x1, x2, dxq, pi

  pi = 4.d0*ATAN(1.d0)

  x1 = x_center(1)
  x2 = x_center(2)
  rho0 = d_region(2)

  e0 = heating_r
  dxq = heating_dx
  ! HeFlash
  !e0 = 2.0e-6 ! erg/g/s (Normalized 10**16)
  !dxq = 0.5d0
  ! Model S
  !e0 = 2.489e-4 ! erg/g/s (Normalized 10**16)
  ! dxq = 2.0d0
  
  ! Initialize
  do i=1,ncell
    e(i) = 0.0d0
  end do 

  ! Heating loop
  do i=1,ncell
    if ((x(i,1) .gt. x1) .and. (x(i,1) .lt. x1+dxq)) then
      ! heating
      !e(i) = e0*(1.0 + cos(2.0*pi*(x(i,1)-x1-dxq/2.0)/dxq))/dxq
      e(i) = e0*rho0 ! erg/s/cm^3
      !e(i) = 0.d0
    else if ((x(i,1) .gt. x2-dxq) .and. (x(i,1) .lt. x2)) then
      ! cooling
      !e(i) = e0*(-1.0 - cos(2.0*pi*(x(i,1)-x2+dxq/2.0)/dxq))/dxq
      e(i) = (-e0)*rho0 ! erg/s/cm^3
      !e(i) = 0.d0
    end if
  end do 


end subroutine eneana

!================================================================
!================================================================
!================================================================
!================================================================

subroutine spongelayers(x,u,req,peq,t,ncell)

  use amr_parameters
  use hydro_parameters

  implicit none
  integer ::ncell                             ! Size of input arrays
  real(dp)::t                                 ! Current time
  real(dp),dimension(1:nvector,1:ndim)::x     ! Cell center position.
  real(dp),dimension(1:nvector,1:nvar)::u     ! Conservative variables 
  real(dp),dimension(1:nvector)::req,peq      ! Equilibrium profiles
  !================================================================
  integer :: i, irad
  real(dp):: x1, x2, crho, cp, cv
  
  ! c: factor of "damping"
  crho=0.1d0
  cv=0.1d0
  cp=0.1d0

  x1 = x_center(1)
  x2 = x_center(2)
  
!   ! Sponge loop 
!   do i=1,ncell
!     if ((x(i,1) .lt. x1-0.5) .or. (x(i,1) .gt. x2+0.5)) then
       
!       ! Damp the internal energy
!       u(i,ndim+2) = u(i,ndim+2)*(1.0d0-cp) + cp*peq(i)/(gamma-1.0d0)
      
!       ! damp density
!       u(i,1) = u(i,1)*(1.0d0 - crho) + crho*req(i)
      
!       ! damp momentum components
!       u(i,2) = u(i,2)*(1.0d0-cv)
! #if NDIM>1
!       u(i,3) = u(i,3)*(1.0d0-cv)
! #endif
! #if NDIM>2
!       u(i,4) = u(i,4)*(1.0d0-cv)
! #endif

!     ! else
! !       if (t .lt. 500.0) then
! !         u(i,2) = u(i,2)*(t/500.0)
! ! #if NDIM>1
! !         u(i,3) = u(i,3)*(t/500.0)
! ! #endif
! ! #if NDIM>2
! !         u(i,4) = u(i,4)*(t/500.0)
! ! #endif
! !       end if 
!     end if

!   end do 

end subroutine spongelayers
