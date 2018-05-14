SUBROUTINE init()
  USE dqmc_complex, lam_ising_=>lam_ising
  USE model
  IMPLICIT NONE
  
  INTEGER orb,site,i,a,b,da,db,orb2,nhop
  REAL(8) re,im
  COMPLEX(8) ga,lam,ga2,lam2
  
  OPEN(10,file='job.in')
  READ(10,*) restart
  READ(10,*) proj
  READ(10,*) sun
  READ(10,*) beta
  READ(10,*) dtau
 
  ntime=ceiling(beta/dtau)
  IF(mod(ntime,2)==1)ntime=ntime+1
  IF(ntime<2)ntime=2
  dtau=beta/ntime
  IF(id==0) PRINT*, 'In this program: dtau=',dtau,' and ntime=',ntime
  
  READ(10,*) nsp
  IF(proj.and.(nsp>ntime/2-1))THEN
    IF(id==0)PRINT*, 'ERROR: nsp should be less than ntime/2 in T=0 version'
    CALL exit(0)
  END IF
  IF((.not.proj).and.(nsp>ntime))THEN
    nsp=ntime
    IF(id==0)PRINT*,'PLEASE NOTICE: nsp is reset to ntime.'
  END IF
  READ(10,*) nbin
  READ(10,*) nwarmup
  READ(10,*) nmeasure
  READ(10,*) ninterval
  READ(10,*) ntmpout
  READ(10,*) nscratch
  READ(10,*) ngroup
  READ(10,*) nising,nphi
  IF(nising>0)ALLOCATE(nglobal_ising(nising))
  IF(nphi>0)ALLOCATE(nglobal_phi(nphi))
  READ(10,*) nglobal_ising,nglobal_phi
  READ(10,*) randomseed
  READ(10,*) norb
  READ(10,*) La,Lb;    nsite=La*Lb*norb
  READ(10,*) nelec;    
  IF(nelec>nsite)THEN
    IF(id==0)PRINT*, 'ERROR: nelec should not be larger than nsite'
    CALL exit(0)
  END IF
  READ(10,*) pbca,pbcb
  READ(10,*) twista,twistb
  READ(10,*) a0r(1:2)
  READ(10,*) b0r(1:2)
  
  ALLOCATE(rorb(2,norb))
  DO orb=1,norb
    READ(10,*) rorb(1:2,orb)
  END DO

  ALLOCATE(U(norb))
  DO orb=1,norb
    READ(10,*) U(orb)
  END DO
  
  READ(10,*) cuta,cutb
  ALLOCATE(hop(-cuta:cuta,-cutb:cutb,norb,norb))
  hop=0d0

  READ(10,*) nhop
  DO i=1,nhop
    READ(10,*) da,db,orb,orb2,re,im
    hop(da,db,orb,orb2)=dcmplx(re,im)
    hop(-da,-db,orb2,orb)=dcmplx(re,-im)
  END DO 

  CLOSE(10)
  !=========================================================
  
  ALLOCATE(kmat(nsite,nsite))
  
  ! set up expk
  CALL init_kmat(); CALL set_expk(kmat,dtau)

  ! set up slater
  IF(proj)THEN
    re=twista; im=twistb
    twista=3d0/17; twistb=3d0/13;   ! add this twisted boundary condition to break degeneraty
    CALL init_kmat(); CALL set_slater(kmat)
    twista=re; twistb=im
  END IF

  ! set up kmat for measurement possibly
  CALL init_kmat()
  OPEN(10,FILE='kmat.dat'); WRITE(10,'(2f18.10)') kmat; CLOSE(10)

  ! set form factor
  nform=1
  ALLOCATE(ndim_form(1)); ndim_form(:)=1
  ALLOCATE(nb_form(nsite,1,1))
  ALLOCATE(fmat(1,1,1))
  DO site=1,nsite
    nb_form(site,1,1)=site
  END DO

  ! set up ising
  nising=1  
  ALLOCATE(mask_form(nsite,1)); mask_form(:,:)=.true.
  ALLOCATE(isingmax(1))
  IF(sun==2)THEN
    isingmax(1)=2
  ELSE
    isingmax(1)=4
  END IF
  ALLOCATE(form_ising(1)); form_ising(:)=1

  ! set lam and gamma for ising fields
  ALLOCATE(gamma_ising(isingmax(1),nsite))
  ALLOCATE(lam_ising(isingmax(1),nsite))
  ALLOCATE(lam_ising_(isingmax(1),nsite,1))
  
  ! get lam and gamma
  IF(sun==2)THEN
    DO orb=1,norb
      CALL HS1(ga,lam,exp(-dtau*U(orb)/2))
      DO a=1,La; DO b=1,Lb; i=label(a,b,orb)
        gamma_ising(1,i)=ga*exp(-lam*sun/2)
        gamma_ising(2,i)=ga*exp(lam*sun/2)
        lam_ising(1,i)=lam
        lam_ising(2,i)=-lam
      END DO; END DO
    END DO
  ELSE
    DO orb=1,norb
      IF(sun==4.or.sun==6)THEN
        CALL HS2(ga,lam,ga2,lam2,exp(-dtau*U(orb)/2))
      ELSE
        CALL HSgeneral(ga,lam,ga2,lam2,-dtau*U(orb)/2)
      END IF
      DO a=1,La; DO b=1,Lb; I=label(a,b,orb)
        gamma_ising(1,i)=ga*exp(-lam*sun/2)
        gamma_ising(2,i)=ga*exp(lam*sun/2)
        gamma_ising(3,i)=ga2*exp(-lam2*sun/2)
        gamma_ising(4,i)=ga2*exp(lam2*sun/2)
        lam_ising(1,i)=lam
        lam_ising(2,i)=-lam
        lam_ising(3,i)=lam2
        lam_ising(4,i)=-lam2
      END DO; END DO
    END DO
  END IF

  lam_ising_(:,:,1)=lam_ising(:,:)

  ! set phi
  nphi=0

  ALLOCATE(mask_ising(1))
  mask_ising(1)=.true.

  ! poolsize_r and poolsize_z, depends on measurement
  poolsize_r=0
  poolsize_z=6+nsite**2

END SUBROUTINE
