! I am designing a new/simpler version:
! + multi-flavor fermions are supportted
! + both discrete and continuous fields are treated in the same way 
!   so that the core module is simpler
! + each field has a unique form
! + BB...B are stored as Bstring to accelerate. A cheap realization::q

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

  COMPLEX(8), PARAMETER, PRIVATE :: zone=(1d0,0d0)
  COMPLEX(8), PARAMETER, PRIVATE :: zmone=(-1d0,0d0)
  COMPLEX(8), PARAMETER, PRIVATE :: zzero=(0d0,0d0)

  LOGICAL :: restart=.false.  ! restart mode or not
  INTEGER :: ith_start=0      ! start from a given bin
  INTEGER :: mcs_start=1      ! start from a given MC step
  LOGICAL :: proj             ! true for T=0, false for T>0
  
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

  COMPLEX(8) :: currentphase=(1d0,0d0)            ! the phase of the current configuration
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
  !INTEGER(8), PRIVATE :: Ntotal=0   ! total number of tries to update locally
  !INTEGER(8), PRIVATE :: Naccept=0  ! total number of accepted tries
  INTEGER(8), ALLOCATABLE, PRIVATE :: Ntotal_ising(:)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Ntotal_phi(:)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Naccept_ising(:)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Naccept_phi(:)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Ntotal_ising_global(:)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Ntotal_phi_global(:)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Naccept_ising_global(:)
  INTEGER(8), ALLOCATABLE, PRIVATE :: Naccept_phi_global(:)

  PRIVATE :: dqmc_update, runtime

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
    CALL init_internal()  

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
          IF(mod(mcs,ntmpout)==0) CALL tmpout(ith,mcs)

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
          
          IF(mod(mcs,ntmpout)==0) CALL tmpout(ith,mcs)
        
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

  END SUBROUTINE

  !> initialize ising configuration randomly
  SUBROUTINE init_ising_random()
    IMPLICIT NONE
    INTEGER ifield,time,site
    ising=0
    DO ifield=1,nising; IF(.not.mask_ising(ifield))CYCLE
      DO time=1,ntime
        DO site=1,nsite
          IF(mask_form(site,form_ising(ifield))) ising(site,time,ifield)=irand(isingmax(ifield))+1
        END DO
      END DO
    END DO
  END SUBROUTINE

  !> initialize phi configuration randomly
  SUBROUTINE init_phi_random()
    IMPLICIT NONE
    INTEGER ifield,time,site
    phi=0d0
    DO ifield=1,nphi; IF(.not.mask_phi(ifield))CYCLE
      DO time=1,ntime
        DO site=1,nsite
          IF(mask_form(site,form_phi(ifield))) phi(site,time,ifield)=drand_sym()*dphi(ifield)
        END DO
      END DO
    END DO
  END SUBROUTINE

  !> set isingflip
  SUBROUTINE set_isingflip()
    IMPLICIT NONE
    INTEGER n,oldising,i
    n=maxval(isingmax)
    IF(.not.allocated(isingflip))ALLOCATE(isingflip(n-1,n))
    DO oldising=1,n
      DO i=1,n
        IF(i<oldising)THEN
          isingflip(i,oldising)=i
        ELSEIF(i>oldising)THEN
          isingflip(i-1,oldising)=i
        END IF
      END DO
    END DO
  END SUBROUTINE

  !> set expf_E, expf_U, expf_Udag, ie U*exp(E)*Udag=exp(fmat) from fmat
  SUBROUTINE set_expf()
    IMPLICIT NONE
    INTEGER iform,ndim
    IF(.not.allocated(expf_E))ALLOCATE(expf_E(maxval(ndim_form),nform))
    IF(.not.allocated(expf_U))ALLOCATE(expf_U(maxval(ndim_form),maxval(ndim_form),nform))
    IF(.not.allocated(expf_Udag))ALLOCATE(expf_Udag(maxval(ndim_form),maxval(ndim_form),nform))
    DO iform=1,nform
      ndim=ndim_form(iform)
      IF(ndim==1)THEN
        expf_E(1,iform)=exp(1d0)
        expf_U(1,1,iform)=1d0
        expf_Udag(1,1,iform)=1d0
      ELSE
        expf_U(1:ndim,1:ndim,iform)=fmat(1:ndim,1:ndim,iform)
        CALL eigen(ndim,expf_U(1:ndim,1:ndim,iform),expf_E(1:ndim,iform))
        expf_E(1:ndim,iform)=exp(expf_E(1:ndim,iform))
        expf_Udag(1:ndim,1:ndim,iform)=conjg(transpose(expf_U(1:ndim,1:ndim,iform)))
      END IF
    END DO
  END SUBROUTINE

  !> set exp(lam*fmat) and exp(-lam*fmat)
  SUBROUTINE set_expflam()
    IMPLICIT NONE
    INTEGER ifield,iform,ndim,site,s,a,b
    IF(.not.allocated(expflam))THEN
      ALLOCATE(expflam(maxval(ndim_form),maxval(ndim_form),maxval(isingmax),nsite,nising))
    END IF
    IF(.not.allocated(inv_expflam))THEN
      ALLOCATE(inv_expflam(maxval(ndim_form),maxval(ndim_form),maxval(isingmax),nsite,nising))
    END IF
    DO ifield=1,nising
      iform=form_ising(ifield)
      ndim=ndim_form(iform)
      DO site=1,nsite
        DO s=1,isingmax(ifield)
          DO a=1,ndim
            DO b=1,ndim
              expflam(a,b,s,site,ifield)=sum(expf_U(a,1:ndim,iform) &
              & *expf_E(1:ndim,iform)**lam_ising(s,site,ifield)*expf_Udag(1:ndim,b,iform))
              inv_expflam(a,b,s,site,ifield)=sum(expf_U(a,1:ndim,iform) &
              & /expf_E(1:ndim,iform)**lam_ising(s,site,ifield)*expf_Udag(1:ndim,b,iform))
            END DO
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE
  
  !> set exp(lam'*fmat)/exp(lam*fmat)-1
  SUBROUTINE set_diff_ef()
    IMPLICIT NONE
    INTEGER ifield,ndim,site,a,b,d
    IF(.not.allocated(diff_ef))THEN
      ALLOCATE(diff_ef(maxval(ndim_form),maxval(ndim_form),maxval(isingmax),maxval(isingmax),nsite,nising))
    END IF
    DO ifield=1,nising
      ndim=ndim_form(form_ising(ifield))
      DO site=1,nsite
        DO a=1,isingmax(ifield)
          DO b=1,isingmax(ifield)
            diff_ef(1:ndim,1:ndim,a,b,site,ifield)=matmul(expflam(1:ndim,1:ndim,a,site,ifield),&
            & inv_expflam(1:ndim,1:ndim,b,site,ifield))
            DO d=1,ndim
              diff_ef(d,d,a,b,site,ifield)=diff_ef(d,d,a,b,site,ifield)-1d0
            END DO
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE
  
  !> set slater for T=0 dqmc
  SUBROUTINE set_slater(kmat)
    IMPLICIT NONE
    COMPLEX(8) kmat(nsite,nsite)
    REAL(8) eval(nsite)
    IF(.not.proj)THEN
      IF(id==0)PRINT*,'Slater is not needed in finite-T DQMC.'
      RETURN
    END IF
    IF(.not.allocated(slater))ALLOCATE(slater(nsite,nelec))
    CALL eigen(nsite,kmat,eval)
    slater(:,:)=kmat(:,1:nelec)
    slater=matmul(expk_half,slater)
  END SUBROUTINE

  !> set exp(K)=exp(-dtau*H0)
  SUBROUTINE set_expk(kmat,dtau)
    IMPLICIT NONE
    COMPLEX(8) kmat(nsite,nsite)
    REAL(8) eval(nsite),dtau
    INTEGER i,j
    IF(.not.allocated(expk))ALLOCATE(expk(nsite,nsite))
    IF(.not.allocated(inv_expk))ALLOCATE(inv_expk(nsite,nsite))
    IF(.not.allocated(expk_half))ALLOCATE(expk_half(nsite,nsite))
    IF(.not.allocated(inv_expk_half))ALLOCATE(inv_expk_half(nsite,nsite))
    CALL eigen(nsite,kmat,eval)
    DO i=1,nsite
      DO j=1,nsite
        expk(i,j)=sum(kmat(i,:)*exp(-dtau*eval(:))*conjg(kmat(j,:)))
        inv_expk(i,j)=sum(kmat(i,:)*exp(dtau*eval(:))*conjg(kmat(j,:)))
        expk_half(i,j)=sum(kmat(i,:)*exp(-dtau*eval(:)/2)*conjg(kmat(j,:)))
        inv_expk_half(i,j)=sum(kmat(i,:)*exp(dtau*eval(:)/2)*conjg(kmat(j,:)))
      END DO
    END DO
  END SUBROUTINE

  !> internal init subroutine, to set remaining parameters
  SUBROUTINE init_internal()
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
      call loadtmp()
    ELSE
      CALL set_pool(nbin,poolsize_r,poolsize_z)
      currentphase=(1d0,0d0)
      IF(nising>0) CALL init_ising_random()
      IF(nphi>0) CALL init_phi_random()
    END IF

  END SUBROUTINE

  
  !> save variables temporarily in case of unexpected break off.
  SUBROUTINE tmpout(ith,mcs)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ith,mcs
    CHARACTER*4 cid
            
            IF(id==0)THEN
              ! temporirally output the running status
              PRINT'(2i8,1a,1e13.6,1a,1e13.6,1a,1f15.2,1A,1e13.6)',ith,mcs,'     (', &
              &  real(currentphase),',',aimag(currentphase),')',runtime(),'s     ',err_fast
              err_fast=0d0
              IF(nising>0)THEN
                PRINT'(1a50,12f7.4)','ising acceptance rate of local update:', &
                        & Naccept_ising(:)*1d0/(Ntotal_ising(:)+1d-30)
                PRINT'(1a50,12f7.4)','ising acceptance rate of global update:', &
                        & Naccept_ising_global(:)*1d0/(Ntotal_ising_global(:)+1d-30)
              END IF
              IF(nphi>0)THEN
                PRINT'(1a50,12f7.4)','phi acceptance rate of local update:', &
                        & Naccept_phi(:)*1d0/(Ntotal_phi(:)+1d-30)
                PRINT'(1a50,12f7.4)','phi acceptance rate of global update:', &
                        & Naccept_phi_global(:)*1d0/(Ntotal_phi_global(:)+1d-30)
                PRINT*,'phi range:'
                DO ifield=1,nphi
                  IF(mask_phi(ifield)) PRINT'(1i4,1a46,4e15.4)',ifield,'  max(Re), min(Re), max(Im), min(Im):', &
                  & maxval(real(phi(:,:,ifield))),minval(real(phi(:,:,ifield))), &
                  & maxval(aimag(phi(:,:,ifield))),minval(aimag(phi(:,:,ifield)))
                END DO
              END IF
              PRINT*,'-----------------------------------'
            END IF

    WRITE(cid,'(1i4)') id
    
    OPEN(81,FILE='tmpout.dat'//trim(adjustl(cid)),FORM='UNFORMATTED')
    WRITE(81) ith,mcs
    WRITE(81) currentphase
    WRITE(81) g
    IF(nising>0) WRITE(81) ising,Ntotal_ising,Naccept_ising,Ntotal_ising_global,Naccept_ising_global
    IF(nphi>0) WRITE(81) phi,Ntotal_phi,Naccept_phi,Ntotal_phi_global,Naccept_phi_global
    CLOSE(81)

    IF(ith>0)CALL save_pool()

  END SUBROUTINE

  !> load variables from previous break point.
  SUBROUTINE loadtmp()
    IMPLICIT NONE
    CHARACTER*4 cid
    WRITE(cid,'(1i4)') id

    OPEN(81,FILE='tmpout.dat'//trim(adjustl(cid)),FORM='UNFORMATTED')
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

  END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! These following evolve subroutines are most time consuming.        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  SUBROUTINE evolve_left_ising(time,i,matrix,d)
  ! exp(lam*F)*B
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,i,d
    INTEGER iform,ndim,site,a,b,sitea,siteb,newising
    COMPLEX(8) matrix(nsite,d),gtmp(nsite,d)
    iform=form_ising(i)
    ndim=ndim_form(iform)
    IF(ndim==1)THEN
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
        newising=ising(site,time,i)
        matrix(site,:)=expflam(1,1,newising,site,i)*matrix(site,:)
      END DO
    ELSE
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
      gtmp=matrix
        newising=ising(site,time,i)
        DO a=1,ndim; sitea=nb_form(site,a,iform)
          DO b=1,ndim; siteb=nb_form(site,b,iform)
            gtmp(sitea,:)=gtmp(sitea,:)+(expflam(a,b,newising,site,i)-transfer(a==b,a))*matrix(siteb,:)
          END DO
        END DO
      matrix=gtmp
      END DO
    END IF
  END SUBROUTINE

  SUBROUTINE evolve_right_ising(time,i,matrix,d)
  ! matrix*exp(lam*F)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,i,d
    INTEGER iform,ndim,site,a,b,sitea,siteb,newising
    COMPLEX(8) matrix(d,nsite),gtmp(d,nsite)
    iform=form_ising(i)
    ndim=ndim_form(iform)
    IF(ndim==1)THEN
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
        newising=ising(site,time,i)
        matrix(:,site)=matrix(:,site)*expflam(1,1,newising,site,i)
      END DO
    ELSE
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
      gtmp=matrix
        newising=ising(site,time,i)
        DO a=1,ndim; sitea=nb_form(site,a,iform)
          DO b=1,ndim; siteb=nb_form(site,b,iform)
            gtmp(:,siteb)=gtmp(:,siteb)+matrix(:,sitea)*(expflam(a,b,newising,site,i)-transfer(a==b,a))
          END DO
        END DO
      matrix=gtmp
      END DO
    END IF
  END SUBROUTINE

  SUBROUTINE inv_evolve_left_ising(time,i,matrix,d)
  ! exp(-lam*F)*B
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,i,d
    INTEGER iform,ndim,site,a,b,sitea,siteb,newising
    COMPLEX(8) matrix(nsite,d),gtmp(nsite,d)
    iform=form_ising(i)
    ndim=ndim_form(iform)
    IF(ndim==1)THEN
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
        newising=ising(site,time,i)
        matrix(site,:)=matrix(site,:)/expflam(1,1,newising,site,i)
      END DO
    ELSE
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
      gtmp=matrix
        newising=ising(site,time,i)
        DO a=1,ndim; sitea=nb_form(site,a,iform)
          DO b=1,ndim; siteb=nb_form(site,b,iform)
            gtmp(sitea,:)=gtmp(sitea,:)+(inv_expflam(a,b,newising,site,i)-transfer(a==b,a))*matrix(siteb,:)
          END DO
        END DO
      matrix=gtmp
      END DO
    END IF
  END SUBROUTINE

  SUBROUTINE inv_evolve_right_ising(time,i,matrix,d)
  ! matrix*exp(-lam*F)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,i,d
    INTEGER iform,ndim,site,a,b,sitea,siteb,newising
    COMPLEX(8) matrix(d,nsite),gtmp(d,nsite)
    iform=form_ising(i)
    ndim=ndim_form(iform)
    IF(ndim==1)THEN
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
        newising=ising(site,time,i)
        matrix(:,site)=matrix(:,site)/expflam(1,1,newising,site,i)
      END DO
    ELSE
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
      gtmp=matrix
        newising=ising(site,time,i)
        DO a=1,ndim; sitea=nb_form(site,a,iform)
          DO b=1,ndim; siteb=nb_form(site,b,iform)
            gtmp(:,siteb)=gtmp(:,siteb)+matrix(:,sitea)*(inv_expflam(a,b,newising,site,i)-transfer(a==b,a))
          END DO
        END DO
      matrix=gtmp
      END DO
    END IF
  END SUBROUTINE

  SUBROUTINE evolve_left_field(time,ifield,matrix,d,flv)
    IMPLICIT NONE
    
    ndim=ndim_field(iform)
  
  END SUBROUTINE


  
  SUBROUTINE evolve_left_phi(time,i,matrix,d)
  ! exp(phi*F)*B
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,i,d
    INTEGER iform,ndim,site,a,b,sitea,siteb
    COMPLEX(8) matrix(nsite,d),gtmp(nsite,d),newphi,delta
    iform=form_phi(i)
    ndim=ndim_form(iform)
    IF(ndim==1)THEN
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
        newphi=phi(site,time,i)
        matrix(site,:)=expf_E(1,iform)**newphi*matrix(site,:)
      END DO
    ELSE
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
      gtmp=matrix
        newphi=phi(site,time,i)
        DO a=1,ndim; sitea=nb_form(site,a,iform)
          DO b=1,ndim; siteb=nb_form(site,b,iform)
            delta=sum(expf_U(a,1:ndim,iform)*expf_E(1:ndim,iform)**newphi*expf_Udag(1:ndim,b,iform))-transfer(a==b,a)
            gtmp(sitea,:)=gtmp(sitea,:)+delta*matrix(siteb,:)
          END DO
        END DO
        matrix=gtmp
      END DO
      !matrix=gtmp
    END IF
  END SUBROUTINE

  SUBROUTINE evolve_right_phi(time,i,matrix,d)
  ! matrix*exp(phi*F)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,i,d
    INTEGER iform,ndim,site,a,b,sitea,siteb
    COMPLEX(8) matrix(d,nsite),gtmp(d,nsite),newphi,delta
    iform=form_phi(i)
    ndim=ndim_form(iform)
    IF(ndim==1)THEN
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
        newphi=phi(site,time,i)
        matrix(:,site)=matrix(:,site)*expf_E(1,iform)**newphi
      END DO
    ELSE
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
      gtmp=matrix
        newphi=phi(site,time,i)
        DO a=1,ndim; sitea=nb_form(site,a,iform)
          DO b=1,ndim; siteb=nb_form(site,b,iform)
            delta=sum(expf_U(a,1:ndim,iform)*expf_E(1:ndim,iform)**newphi*expf_Udag(1:ndim,b,iform))-transfer(a==b,a)
            gtmp(:,siteb)=gtmp(:,siteb)+matrix(:,sitea)*delta
          END DO
        END DO
      matrix=gtmp
      END DO
    END IF
  END SUBROUTINE

  SUBROUTINE inv_evolve_left_phi(time,i,matrix,d)
  ! exp(-phi*F)*B
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,i,d
    INTEGER iform,ndim,site,a,b,sitea,siteb
    COMPLEX(8) matrix(nsite,d),gtmp(nsite,d),newphi,delta
    iform=form_phi(i)
    ndim=ndim_form(iform)
    IF(ndim==1)THEN
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
        newphi=phi(site,time,i)
        matrix(site,:)=matrix(site,:)/expf_E(1,iform)**newphi
      END DO
    ELSE
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
      gtmp=matrix
        newphi=phi(site,time,i)
        DO a=1,ndim; sitea=nb_form(site,a,iform)
          DO b=1,ndim; siteb=nb_form(site,b,iform)
            delta=sum(expf_U(a,1:ndim,iform)/expf_E(1:ndim,iform)**newphi*expf_Udag(1:ndim,b,iform))-transfer(a==b,a)
            gtmp(sitea,:)=gtmp(sitea,:)+delta*matrix(siteb,:)
          END DO
        END DO
      matrix=gtmp
      END DO
    END IF
  END SUBROUTINE

  SUBROUTINE inv_evolve_right_phi(time,i,matrix,d)
  ! matrix*exp(-phi*F)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,i,d
    INTEGER iform,ndim,site,a,b,sitea,siteb
    COMPLEX(8) matrix(d,nsite),gtmp(d,nsite),newphi,delta
    iform=form_phi(i)
    ndim=ndim_form(iform)
    IF(ndim==1)THEN
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
        newphi=phi(site,time,i)
        matrix(:,site)=matrix(:,site)/expf_E(1,iform)**newphi
      END DO
    ELSE
      DO site=1,nsite; IF(.not.mask_form(site,iform))CYCLE
      gtmp=matrix
        newphi=phi(site,time,i)
        DO a=1,ndim; sitea=nb_form(site,a,iform)
          DO b=1,ndim; siteb=nb_form(site,b,iform)
            delta=sum(expf_U(a,1:ndim,iform)/expf_E(1:ndim,iform)**newphi*expf_Udag(1:ndim,b,iform))-transfer(a==b,a)
            gtmp(:,siteb)=gtmp(:,siteb)+matrix(:,sitea)*delta
          END DO
        END DO
      matrix=gtmp
      END DO
    END IF
  END SUBROUTINE

  SUBROUTINE evolve_left(time,matrix,d)
  ! B(time)*matrix
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d
    INTEGER i
    COMPLEX(8) matrix(nsite,d),gtmp(nsite,d)
    gtmp=matrix
    DO i=1,nising;
      IF(mask_ising(i))CALL evolve_left_ising(time,i,gtmp,d)
    END DO
    DO i=1,nphi
      IF(mask_phi(i))CALL evolve_left_phi(time,i,gtmp,d)
    END DO
    matrix=matmul(expk,gtmp)
    !CALL zgemm('n','n',nsite,d,nsite,zone,expk,nsite,gtmp,nsite,zzero,matrix,nsite)
  END SUBROUTINE
  
  SUBROUTINE evolve_right(time,matrix,d)
  ! matrix*B(time)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d
    INTEGER i
    COMPLEX(8) matrix(d,nsite),gtmp(d,nsite)
    gtmp=matmul(matrix,expk)
    !CALL zgemm('n','n',d,nsite,nsite,zone,matrix,d,expk,nsite,zzero,gtmp,d)
    DO i=nphi,1,-1
      IF(mask_phi(i))CALL evolve_right_phi(time,i,gtmp,d)
    END DO
    DO i=nising,1,-1
      IF(mask_ising(i))CALL evolve_right_ising(time,i,gtmp,d)
    END DO
    matrix=gtmp
  END SUBROUTINE
  
  SUBROUTINE inv_evolve_left(time,matrix,d)
  ! B(time)^{-1}*matrix
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d
    INTEGER i
    COMPLEX(8) matrix(nsite,d),gtmp(nsite,d)
    gtmp=matmul(inv_expk,matrix)
    !CALL zgemm('n','n',nsite,d,nsite,zone,inv_expk,nsite,matrix,nsite,zzero,gtmp,nsite)
    DO i=nphi,1,-1
      IF(mask_phi(i))CALL inv_evolve_left_phi(time,i,gtmp,d)
    END DO
    DO i=nising,1,-1
      IF(mask_ising(i))CALL inv_evolve_left_ising(time,i,gtmp,d)
    END DO
    matrix=gtmp
  END SUBROUTINE
  
  SUBROUTINE inv_evolve_right(time,matrix,d)
  ! matrix*B(time)^{-1}
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d
    INTEGER i
    COMPLEX(8) matrix(d,nsite),gtmp(d,nsite)
    gtmp=matrix
    DO i=1,nising
      IF(mask_ising(i))CALL inv_evolve_right_ising(time,i,gtmp,d)
    END DO
    DO i=1,nphi
      IF(mask_phi(i))CALL inv_evolve_right_phi(time,i,gtmp,d)
    END DO
    matrix=matmul(gtmp,inv_expk)
    !CALL zgemm('n','n',d,nsite,nsite,zone,gtmp,d,inv_expk,nsite,zzero,matrix,d)
  END SUBROUTINE

  ! xxx_evolve_xxx_2nd() are used ONLY in calculations of time evolved Green's functions  
  SUBROUTINE evolve_left_2nd(time,matrix,d)
  ! B(time)*matrix
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d
    INTEGER i
    COMPLEX(8) matrix(nsite,d),gtmp(nsite,d)
    gtmp=matmul(expk_half,matrix)
    !CALL zgemm('n','n',nsite,d,nsite,zone,expk_half,nsite,matrix,nsite,zzero,gtmp,nsite)
    DO i=1,nising
      IF(mask_ising(i))CALL evolve_left_ising(time,i,gtmp,d)
    END DO
    DO i=1,nphi
      IF(mask_phi(i))CALL evolve_left_phi(time,i,gtmp,d)
    END DO
    matrix=matmul(expk_half,gtmp)
    !CALL zgemm('n','n',nsite,d,nsite,zone,expk_half,nsite,gtmp,nsite,zzero,matrix,nsite)
  END SUBROUTINE
  
  SUBROUTINE evolve_right_2nd(time,matrix,d)
  ! matrix*B(time)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d
    INTEGER i
    COMPLEX(8) matrix(d,nsite),gtmp(d,nsite)
    gtmp=matmul(matrix,expk_half)
    !CALL zgemm('n','n',d,nsite,nsite,zone,matrix,d,expk_half,nsite,zzero,gtmp,d)
    DO i=nphi,1,-1
      IF(mask_phi(i))CALL evolve_right_phi(time,i,gtmp,d)
    END DO
    DO i=nising,1,-1
      IF(mask_ising(i))CALL evolve_right_ising(time,i,gtmp,d)
    END DO
    matrix=matmul(gtmp,expk_half)
    !CALL zgemm('n','n',d,nsite,nsite,zone,gtmp,d,expk_half,nsite,zzero,matrix,d)
  END SUBROUTINE
  
  SUBROUTINE inv_evolve_left_2nd(time,matrix,d)
  ! B(time)^{-1}*matrix
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d
    INTEGER i
    COMPLEX(8) matrix(nsite,d),gtmp(nsite,d)
    gtmp=matmul(inv_expk_half,matrix)
    !CALL zgemm('n','n',nsite,d,nsite,zone,inv_expk_half,nsite,matrix,nsite,zzero,gtmp,nsite)
    DO i=nphi,1,-1
      IF(mask_phi(i))CALL inv_evolve_left_phi(time,i,gtmp,d)
    END DO
    DO i=nising,1,-1
      IF(mask_ising(i))CALL inv_evolve_left_ising(time,i,gtmp,d)
    END DO
    matrix=matmul(inv_expk_half,gtmp)
    !CALL zgemm('n','n',nsite,d,nsite,zone,inv_expk_half,nsite,gtmp,nsite,zzero,matrix,nsite)
  END SUBROUTINE
  
  SUBROUTINE inv_evolve_right_2nd(time,matrix,d)
  ! matrix*B(time)^{-1}
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time,d
    INTEGER i
    COMPLEX(8) matrix(d,nsite),gtmp(d,nsite)
    gtmp=matmul(matrix,inv_expk_half)
    !CALL zgemm('n','n',d,nsite,nsite,zone,matrix,d,inv_expk_half,nsite,zzero,gtmp,d)
    DO i=1,nising
      IF(mask_ising(i))CALL inv_evolve_right_ising(time,i,gtmp,d)
    END DO
    DO i=1,nphi
      IF(mask_phi(i))CALL inv_evolve_right_phi(time,i,gtmp,d)
    END DO
    matrix=matmul(gtmp,inv_expk_half)
    !CALL zgemm('n','n',d,nsite,nsite,zone,gtmp,d,inv_expk_half,nsite,zzero,matrix,d)
  END SUBROUTINE

  FUNCTION runtime()
    IMPLICIT NONE
    REAL(8) runtime
    REAL(8) :: t1=0d0
#ifdef MPI
    IF(id==0)t1=mpi_wtime()
#else
    CALL cpu_time(t1)
#endif
    runtime=t1-t0
  END FUNCTION

  !--------------------------------------------------------------------!
  ! several usually used discrete Hubbard-Stratonovich transformation. !
  !   exp(-dtau*g*op^2)=...                                            !
  !--------------------------------------------------------------------!
  SUBROUTINE HS1(ga,lam,a)  ! a=exp(-dtau*g)
  ! |eig(op)|=0,1, exact
    IMPLICIT NONE
    COMPLEX(8) ga,lam
    REAL(8) a
    ga=0.5d0
    IF(a>1d0)THEN
      lam=acosh(a)
    ELSE
      lam=dcmplx(0d0,acos(a))
    END IF
  END SUBROUTINE

  SUBROUTINE HS2(ga1,lam1,ga2,lam2,a)  ! a=exp(-dtau*g)
  ! |eig(op)|=0,1,2,3, exact
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
  END SUBROUTINE 

  SUBROUTINE HSgeneral(ga1,lam1,ga2,lam2,a) ! a=-dtau*g
  ! |eig(op)|=any value, but approximate.
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
  END SUBROUTINE

  !> sort the array a(n) from large to small absulute values
  SUBROUTINE ordering(n,a)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER i,j
    COMPLEX(8) a(n),tmp
    DO j=n-1,1,-1
      DO i=1,j
        IF(abs(a(i))<abs(a(i+1)))THEN
          tmp=a(i);a(i)=a(i+1);a(i+1)=tmp
        END IF
      END DO
    END DO
  END SUBROUTINE

END MODULE


