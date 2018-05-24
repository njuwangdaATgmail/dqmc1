! I am designing a new/simpler version:
! + multi-flavor fermions are supportted
! + both discrete and continuous fields are treated in the same way 
!   so that the core module is simpler
! + each field has a unique form
! + BB...B are stored as Bstring to accelerate. A cheap realizaion:
!     Bstring_0(:,:,i)=B(i*nscratch)*B(i*nscratch-1)*...*B(1) 
!     Bstring_L(:,:,i)=B(ntime)*...*B(ntime-i*nscratch+1)
! + checkerboard algorithm for expk? maybe not so emergent
MODULE dqmc_complex

#ifdef MPI
  USE mpi
#endif
  USE randomlib
  USE matrixlib
  USE measpool

  LOGICAL :: restart                = .false.   ! restart mode or not
  INTEGER :: ith_start              = 0         ! start from a given bin
  INTEGER :: mcs_start              = 1         ! start from a given MC step
  LOGICAL :: proj                   = .false.   ! true for T=0, false for T>0
  logical :: scratch_global_useful  = .false.   ! whether update_scratch() has been performed in global update
  
  INTEGER nsite    ! number of sites
  INTEGER nelec    ! number of filled electrons, only used for T=0
  INTEGER ntime    ! number of times
  INTEGER nsp      ! do measurements from ntime/2+1-nsp to ntime/2+1+nsp for T=0
                   ! measure G(tau) from time=1 to time=nsp for T>0
  
  COMPLEX(8), ALLOCATABLE :: expk(:,:)            ! exp(-dtau*K), (nsite,nsite)
  COMPLEX(8), ALLOCATABLE :: expk_half(:,:)       ! exp(-dtau*K/2)
  COMPLEX(8), ALLOCATABLE :: inv_expk(:,:)        ! exp(dtau*K)
  COMPLEX(8), ALLOCATABLE :: inv_expk_half(:,:)   ! exp(dtau*K/2)
  COMPLEX(8), ALLOCATABLE :: slater(:,:)          ! (nsite,nelec), trial wave function
  COMPLEX(8), ALLOCATABLE :: g(:,:)               ! Green's function, (nsite,nsite)
  COMPLEX(8), ALLOCATABLE :: Bstring_Q(:,:,:,:)       ! (nsite,nsite/nelec,nblock-1,nflv)
  COMPLEX(8), ALLOCATABLE :: Bstring_D(:,:,:)       ! (nsite/nelec,nblock-1,nflv)
  COMPLEX(8), ALLOCATABLE :: Bstring_R(:,:,:,:)       ! (nsite/nelec,nsite/nelec,nblock-1,nflv)

  INTEGER nising                                  ! descrete auxiliary fields
  LOGICAL, ALLOCATABLE :: mask_ising(:)           ! whether the ising field is simulated. Then we can write a general input file:
                                                  ! e.g. containing different boson fields and only choose a few of them in each run.
  INTEGER, ALLOCATABLE :: ising(:,:,:)            ! (nsite,ntime,nising), I define 1<=ising<=isingmax.
  INTEGER, ALLOCATABLE :: isingmax(:)             ! (nising)
  INTEGER, ALLOCATABLE :: isingflip(:,:)          ! (maxval(isingmax)-1,maxval(isingmax))
  INTEGER, ALLOCATABLE :: form_ising(:)           ! (nising)
  COMPLEX(8), ALLOCATABLE :: lam_ising(:,:,:)       ! (maxval(isingmax),nsite,nising)
  COMPLEX(8), ALLOCATABLE :: expflam(:,:,:,:,:)     ! (maxval(dim_form),maxval(dim_form),maxval(isingmax),nsite,nising)
  COMPLEX(8), ALLOCATABLE :: inv_expflam(:,:,:,:,:) ! (maxval(dim_form),maxval(dim_form),maxval(isingmax),nsite,nising)
  COMPLEX(8), ALLOCATABLE :: diff_ef(:,:,:,:,:,:)   ! (maxval(dim_form),maxval(dim_form),maxval(isingmax),maxval(isingmax),nsite,nising)
  
  INTEGER nphi                                    ! continuous auxiliary fields
  LOGICAL, ALLOCATABLE :: mask_phi(:)             ! whether the phi field is simulated
  COMPLEX(8), ALLOCATABLE :: phi(:,:,:)           ! (nsite,ntime,nphi), I define -phimax<=phi<=phimax
  COMPLEX(8), ALLOCATABLE :: dphi(:)              ! (nphi)
  INTEGER, ALLOCATABLE :: form_phi(:)             ! (nphi)

  INTEGER nform                                   ! number of form factors in c'*F*c
  INTEGER, ALLOCATABLE :: ndim_form(:)            ! (nform), dimension of the form
  LOGICAL, ALLOCATABLE :: mask_form(:,:)          ! (nsite,nform), whether the site belongs to the form
  INTEGER, ALLOCATABLE :: nb_form(:,:,:)          ! (nsite,maxval(ndim_form),nform), neighbour sites of the form
  COMPLEX(8), ALLOCATABLE :: fmat(:,:,:)          ! (maxval(ndim_form),maxval(ndim_form),nform)
  REAL(8), ALLOCATABLE :: expf_E(:,:)             ! (maxval(ndim_form),nform)
  COMPLEX(8), ALLOCATABLE :: expf_U(:,:,:)        ! (maxval(ndim_form),maxval(ndim_form),nform)
  COMPLEX(8), ALLOCATABLE :: expf_Udag(:,:,:)     ! (maxval(ndim_form),maxval(ndim_form),nform)

  COMPLEX(8) :: currentphase  = (1d0,0d0)            ! the phase of the current configuration
  REAL(8)    :: newMetro      = (0d0,0d0)            ! correction of Metropolis ratio
  INTEGER nbin              ! number of data bins. The covariance is calculated between different bins.
  INTEGER nwarmup           ! warmup steps, in unit of Monte Carlo steps (MCS)
  INTEGER nmeasure          ! MCS in each bin
  INTEGER ninterval         ! do measurement every ninterval MCSs
  INTEGER ntmpout           ! output temporarily every ntmpout MCSs
  INTEGER nscratch          ! update G from scratch every nscratch steps
  INTEGER ngroup            ! every ngroup matrices are producted directly before performing QDR
  INTEGER randomseed        ! seed of random number generator
  INTEGER poolsize_r            ! how many real(8) type observables to measure
  INTEGER poolsize_z            ! how many complex(8) type observables to measure
                                ! Should we move poolsize_* outside?
  INTEGER, ALLOCATABLE :: nglobal_ising(:)      ! do global update for ising
  INTEGER, ALLOCATABLE :: nglobal_phi(:)        ! do global update for phi
  
  INTEGER :: id=0           ! the index of the current node
  INTEGER :: nd=1           ! number of total nodes

  REAL(8), PRIVATE :: err_fast      ! difference between scratch and fast-update
  REAL(8), PRIVATE :: t0            ! used to save starting time, in order to obtain running time
  
  INTEGER(8), ALLOCATABLE, PRIVATE :: Ntotal_field(:)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Naccept_field(:)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Ntotal_field_global(:)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Naccept_field_global(:)

  PRIVATE :: dqmc_update_local,dqmc_update_global,runtime_,tmpout_

CONTAINS

  !> This is the entrance of this module. It should be called from outside.
  SUBROUTINE dqmc_driver()
    IMPLICIT NONE
    INTEGER ith,mcs,ierr,ifield

#ifdef MPI
    ! set MPI mode
    CALL mpi_init(ierr)
    CALL mpi_comm_size(mpi_comm_world,nd,ierr)
    CALL mpi_comm_rank(mpi_comm_world,id,ierr)
    IF(id==0)t0=mpi_wtime()  ! It seems mpi_wtime() may fail on very few machines.
#else
    ! set sequential mode
    nd=1
    id=0
    CALL cpu_time(t0)
#endif

    ! external initialization subroutine, set required parameters listed above
    CALL init()
    
    ! internal initialization, set remaining parameters
    CALL init_()  

    DO ith=ith_start,nbin

      IF(ith==0)THEN

        ! do warmup
        DO mcs=mcs_start,nwarmup
          
          ! do global update
          DO ifield=1,nfield; IF(.not.mask_field(ifield))CYCLE
            IF(mod(mcs,ninterval_global(ifield))==0) CALL dqmc_update_global(ifield)
          END DO
          
          ! do local update without doing measurement
          CALL dqmc_update_local(.false.)  
          
          ! output information to screen and save work sapce into hardware
          IF(mod(mcs,ntmpout)==0)THEN
            CALL tmpout_(ith,mcs)
            CALL tmpout()
          END IF

        END DO ! DO mcs=mcs_start,nwarmup

        ! reset mcs_start=1 in case the program starts from a break point with mcs_start>1
        mcs_start=1
    
      ELSE  ! IF(ith==0)
      
        ! do measurement
        DO mcs=mcs_start,nmeasure*ninterval

          ! do global update
          DO ifield=1,nfield; IF(.not.mask_field(ifield))CYCLE
            IF(mod(mcs,ninterval_global(ifield))==0) CALL dqmc_update_global(ifield)
          END DO
          
          ! only do measurement every ninterval space-time sweeps
          IF(mod(mcs,ninterval)==0)THEN
            CALL dqmc_update_local(.true.)
          ELSE
            CALL dqmc_update_local(.false.)
          END IF
          
          IF(mod(mcs,ntmpout)==0)THEN
            CALL tmpout_(ith,mcs)
            CALL tmpout()
          END IF
        
        END DO  ! DO mcs=mcs_start,nmeasure*ninterval
        
        ! reset mcs_start=1 in case the program starts from a break point with mcs_start>1
        mcs_start=1

        ! prepare the pool for the next data bin: move bin to bin+1 and reset head
        CALL move_pool()

      END IF  ! IF(ith==0)

    END DO  ! DO ith=ith_start,nbin

    ! do average of the measured observables in the pool
    CALL average_pool()
    
    ! external subroutine, data analysis and output
    CALL postprocess()

#ifdef MPI
    ! close MPI
    CALL mpi_barrier(mpi_comm_world,ierr)
    CALL mpi_finalize(ierr)
#endif

    IF(id==0)THEN
      PRINT*,'Job is done with',nd,' cores.'
    END IF

  END SUBROUTINE dqmc_driver

  ! internal initial subroutine, to set remaining parameters
  SUBROUTINE init_()
    IMPLICIT NONE
    
    ! set up the random number seed on different cores
    randomseed=randomseed+701703*id
    CALL init_rng(randomseed)
    
    IF(nising>0) ALLOCATE(ising(nsite,ntime,nising))
    IF(nphi>0) ALLOCATE(phi(nsite,ntime,nphi))
    ALLOCATE(g(nsite,nsite))
    
    CALL set_expf()

    IF(nising>0)THEN
      CALL set_isingflip()
      CALL set_expflam()
      CALL set_diff_ef()
    END IF
    
    IF(nising>0)THEN
      ALLOCATE(Ntotal_ising(nising),Naccept_ising(nising))
      ALLOCATE(Ntotal_ising_global(nising),Naccept_ising_global(nising))
      Ntotal_ising(:)=0
      Naccept_ising(:)=0
      Ntotal_ising_global(:)=0
      Naccept_ising_global(:)=0
    END IF

    IF(nphi>0)THEN
      ALLOCATE(Ntotal_phi(nphi),Naccept_phi(nphi))
      ALLOCATE(Ntotal_phi_global(nphi),Naccept_phi_global(nphi))
      Ntotal_phi(:)=0
      Naccept_phi(:)=0
      Ntotal_phi_global(:)=0
      Naccept_phi_global(:)=0
    END IF

    IF(restart)THEN
      !CALL set_pool(nbin,poolsize_r,poolsize_z)
      call loadtmp_()
    ELSE
      CALL set_pool(nbin,poolsize_r,poolsize_z)
      currentphase=(1d0,0d0)
      IF(nising>0) CALL init_ising_random()
      IF(nphi>0) CALL init_phi_random()
    END IF

  END SUBROUTINE init_
  
  ! internal subroutine to output running status to screen 
  ! and save variables temporarily in case of unexpected breakdown (e.g. power off)
  SUBROUTINE tmpout_(ith,mcs)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ith,mcs
    CHARACTER*4 cid
            
    ! temporirally output the running status
    IF(id==0)THEN
      PRINT'(2i8,1a,1e13.6,1a,1e13.6,1a,1f15.2,1A,1e13.6)',ith,mcs,'     (', &
      &  real(currentphase),',',aimag(currentphase),')',runtime_(),'s     ',err_fast
      err_fast=0d0
      PRINT'(1a50,100f7.4)','acceptance rate of local update:', &
      & Naccept_field(:)*1d0/(Ntotal_field(:)+1d-30)
      PRINT'(1a50,100f7.4)','acceptance rate of global update:', &
      & Naccept_field_global(:)*1d0/(Ntotal_field_global(:)+1d-30)
      PRINT*,'field range:'
      DO ifield=1,nfield
        IF(mask_field(ifield)) PRINT'(1i4,1a46,4e15.4)',ifield,'  max(Re), min(Re), max(Im), min(Im):', &
        & maxval(real(field(:,:,ifield))),minval(real(field(:,:,ifield))), &
        & maxval(aimag(field(:,:,ifield))),minval(aimag(field(:,:,ifield)))
      END DO
      PRINT*,'-----------------------------------'
    END IF

    ! save QMC variables to files
    WRITE(cid,'(1i4)') id
    OPEN(81,FILE='tmpout_.dat'//trim(adjustl(cid)),FORM='UNFORMATTED')
    WRITE(81) ith,mcs
    WRITE(81) currentphase
    WRITE(81) g
    IF(nising>0) WRITE(81) ising,Ntotal_ising,Naccept_ising,Ntotal_ising_global,Naccept_ising_global
    IF(nphi>0) WRITE(81) phi,Ntotal_phi,Naccept_phi,Ntotal_phi_global,Naccept_phi_global
    CLOSE(81)
    
    IF(ith>0)CALL save_pool()

  END SUBROUTINE tmpout_

  ! internal subroutine to load variables from previous break point.
  SUBROUTINE loadtmp_()
    IMPLICIT NONE
    CHARACTER*4 cid
    WRITE(cid,'(1i4)') id
    OPEN(81,FILE='tmpout_.dat'//trim(adjustl(cid)),FORM='UNFORMATTED')
    READ(81) ith_start,mcs_start
    READ(81) currentphase
    READ(81) g
    IF(nising>0) READ(81) ising,Ntotal_ising,Naccept_ising,Ntotal_ising_global,Naccept_ising_global
    IF(nphi>0) READ(81) phi,Ntotal_phi,Naccept_phi,Ntotal_phi_global,Naccept_phi_global
    CLOSE(81)

    IF(ith_start>0)THEN
      CALL load_pool()
    ELSE
      CALL set_pool(nbin,poolsize_r,poolsize_z)
    END IF

    mcs_start=mcs_start+1

  END SUBROUTINE loadtmp_
  
  ! internal function to get the running time
  FUNCTION runtime_()
    IMPLICIT NONE
    REAL(8) runtime_,t1
#ifdef MPI
    IF(id==0)t1=mpi_wtime()
#else
    CALL cpu_time(t1)
#endif
    runtime_=t1-t0
  END FUNCTION runtime_

  !
  !------------------------------------------------------------------------------
  !  The following 8 subroutines are provided to perform evolutions of a matrix.
  !       evolve_left_K   : exp(K)*matrix
  !       evolve_right_K  : matrix*exp(K)
  !       evolve_left_V   : exp(V)*matrix
  !       evolve_right_V  : matrix*exp(V)
  !       evolve_left     : B*matrix
  !       evolve_right    : matrix*B
  !       evolve_left_2nd : B_2nd*matrix
  !       evolve_right_2nd: matrix*B_2nd
  !------------------------------------------------------------------------------
  !
  !> exp( +- (dtau, dtau/2)*K ) * matrix
  SUBROUTINE evolve_left_K(matrix,d,flv,inv,half)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: d,flv
    COMPLEX(8) matrix(nsite,d)
    LOGICAL, INTENT(IN) :: inv,half
    IF(inv.and.half)THEN
      matrix=matmul(inv_expk_half(:,:,flv),matrix)
    ELSE IF(.not.inv.and.half)THEN
      matrix=matmul(expk_half(:,:,flv),matrix)
    ELSE IF(inv.and..not.half)THEN
      matrix=matmul(inv_expk(:,:,flv),matrix)
    ELSE
      matrix=matmul(expk(:,:,flv),matrix)
    END IF
  END SUBROUTINE evolve_left_K
  
  !> matrix * exp( +- (dtau, dtau/2)*K )
  SUBROUTINE evolve_right_K(matrix,d,flv,inv,half)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: d,flv
    COMPLEX(8) matrix(d,nsite)
    LOGICAL, INTENT(IN) :: inv,half
    IF(inv.and.half)THEN
      matrix=matmul(matrix,inv_expk_half(:,:,flv))
    ELSE IF(.not.inv.and.half)THEN
      matrix=matmul(matrix,expk_half(:,:,flv))
    ELSE IF(inv.and..not.half)THEN
      matrix=matmul(matrix,inv_expk(:,:,flv))
    ELSE
      matrix=matmul(matrix,expk(:,:,flv))
    END IF
  END SUBROUTINE evolve_right_K

  !> exp( +- V) * matrix
  SUBROUTINE evolve_left_V(time,ifield,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,ifield,d,flv
    LOGICAL, INTENT(IN) :: inv
    INTEGER a,b,sitea,siteb
    COMPLEX(8) matrix(nsite,d),expV(ndim_field(ifield),ndim_field(ifield))
    COMPLEX(8) gtmp(ndim_field(ifield),d)
    DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE
      CALL get_expV(site,time,ifield,flv,inv,expV)
      IF(ndim_field(ifield)==1)THEN
        matrix(site,:)=expV(1,1)*matrix(site,:)
      ELSE
        DO a=1,ndim_field(ifield); sitea=nb_field(site,a,ifield)
          gtmp(a,:)=matrix(sitea,:)
        END DO
        gtmp=matmul(expV,gtmp)
        DO a=1,ndim_field(ifield); sitea=nb_field(site,a,ifield)
          matrix(sitea,:)=gtmp(a,:)
        END DO
      END IF
    END DO
  END SUBROUTINE evolve_left_V
  
  !> matrix * exp( +- V)
  SUBROUTINE evolve_right_V(time,ifield,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,ifield,d,flv
    LOGICAL, INTENT(IN) :: inv
    INTEGER a,b,sitea,siteb
    COMPLEX(8) matrix(d,nsite),expV(ndim_field(ifield),ndim_field(ifield))
    COMPLEX(8) gtmp(d,ndim_field(ifield))
    DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE
      CALL get_expV(site,time,ifield,flv,inv,expV)
      IF(ndim_field(ifield)==1)THEN
        matrix(:,site)=matrix(:,site)*expV(1,1)
      ELSE
        DO a=1,ndim_field(ifield); sitea=nb_field(site,a,ifield)
          gtmp(:,a)=matrix(:,sitea)
        END DO
        gtmp=matmul(gtmp,expV)
        DO a=1,ndim_field(ifield); sitea=nb_field(site,a,ifield)
          matrix(:,sitea)=gtmp(:,a)
        END DO
      END IF
    END DO
  END SUBROUTINE evolve_right_V
  
  ! B(time)*matrix
  SUBROUTINE evolve_left(time,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d,flv
    LOGICAL, INTENT(IN) :: inv
    COMPLEX(8) matrix(nsite,d)
    INTEGER ifield
    DO ifield=1,nfield
      CALL evolve_left_V(time,ifield,matrix,d,flv,inv)
    END DO
    CALL evolve_left_K(matrix,d,flv,inv,.false.)
  END SUBROUTINE evolve_left

  ! matrix*B(time)
  SUBROUTINE evolve_right(time,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d,flv
    LOGICAL, INTENT(IN) :: inv
    COMPLEX(8) matrix(d,nsite)
    INTEGER ifield
    CALL evolve_right_K(matrix,d,flv,inv,.false.)
    DO ifield=nfield,1,-1
      CALL evolve_right_V(time,ifield,matrix,d,flv,inv)
    END DO
  END SUBROUTINE evolve_right

  ! B_2nd(time)*matrix
  SUBROUTINE evolve_left_2nd(time,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d,flv
    LOGICAL, INTENT(IN) :: inv
    COMPLEX(8) matrix(nsite,d)
    INTEGER ifield
    CALL evolve_left_K(matrix,d,flv,inv,.true.)
    DO ifield=1,nfield
      CALL evolve_left_V(time,ifield,matrix,d,flv,inv)
    END DO
    CALL evolve_left_K(matrix,d,flv,inv,.true.)
  END SUBROUTINE evolve_left_2nd

  ! matrix*B_2nd(time)
  SUBROUTINE evolve_right_2nd(time,matrix,d,flv,inv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d,flv
    LOGICAL, INTENT(IN) :: inv
    COMPLEX(8) matrix(d,nsite)
    INTEGER ifield
    CALL evolve_right_K(matrix,d,flv,inv,.true.)
    DO ifield=nfield,1,-1
      CALL evolve_right_V(time,ifield,matrix,d,flv,inv)
    END DO
    CALL evolve_right_K(matrix,d,flv,inv,.true.)
  END SUBROUTINE evolve_right_2nd

  !
  !--------------------------------------------------------------------
  ! several usually used discrete Hubbard-Stratonovich transformation
  ! HS1:
  !   exp(-dtau*g*op^2) = ga*exp(lam*op) + ga*exp(-lam*op)
  ! HS2, HSgeneral:
  !   exp(-dtau*g*op^2) = ga1*exp(lam1*op) + ga1*exp(-lam1*op)
  !                     + ga2*exp(lam2*op) + ga2*exp(-lam2*op)
  !--------------------------------------------------------------------
  !
  ! |eig(op)|=0,1, exact
  SUBROUTINE HS1(ga,lam,a)  ! a=exp(-dtau*g)
    IMPLICIT NONE
    COMPLEX(8) ga,lam
    REAL(8) a
    ga=0.5d0
    IF(a>1d0)THEN
      lam=acosh(a)
    ELSE
      lam=dcmplx(0d0,acos(a))
    END IF
  END SUBROUTINE HS1

  ! |eig(op)|=0,1,2,3, exact
  SUBROUTINE HS2(ga1,lam1,ga2,lam2,a)  ! a=exp(-dtau*g)
    IMPLICIT NONE
    COMPLEX(8) ga1,ga2,lam1,lam2
    REAL(8) a,d
    d=sqrt(8d0+a**2*(3+a**2)**2)
    ga1=(-a*(3+a**2)+d)/(4*d)
    ga2=(a*(3+a**2)+d)/(4*d)
    IF(a>1d0)THEN
      lam1=acosh((a+2*a**3+a**5+(a**2-1)*d)/4)
      lam2=acosh((a+2*a**3+a**5-(a**2-1)*d)/4)
    ELSE
      lam1=DCMPLX(0d0,acos((a+2*a**3+a**5+(a**2-1)*d)/4))
      lam2=DCMPLX(0d0,acos((a+2*a**3+a**5-(a**2-1)*d)/4))
    END IF
  END SUBROUTINE HS2

  ! |eig(op)|=any value, but approximate.
  SUBROUTINE HSgeneral(ga1,lam1,ga2,lam2,a) ! a=-dtau*g
    IMPLICIT NONE
    COMPLEX(8) ga1,lam1,ga2,lam2
    REAL(8) a
    ga1=(3d0+sqrt(6d0))/12d0
    ga2=(3d0-sqrt(6d0))/12d0
    IF(a>0d0)THEN
      lam1=sqrt(a*2*(3d0-sqrt(6d0)))
      lam2=sqrt(a*2*(3d0+sqrt(6d0)))
    ELSE
      lam1=DCMPLX(0d0,sqrt(-a*2*(3d0-sqrt(6d0))))
      lam2=DCMPLX(0d0,sqrt(-a*2*(3d0+sqrt(6d0))))
    END IF
  END SUBROUTINE HSgeneral

END MODULE mod_dqmc_complex


