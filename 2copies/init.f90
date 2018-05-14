!> This is an realization of the SU(2)-symmetric Hubbard-Holstein model on SQUARE lattice.
!> The Hubbard interaction is decoupled in the magnetic channel.
!> Spin up and Spin down are viewed as two layers.
!!
!! el-el interactions considered: V_U
!! 1.     -(U/2)(nup-ndn)^2, form 1
!! 
!! el-ph interactions considered: Vph_hol
!! 1.     onsite, Holstein: phi*(nup+ndn-1), form 2
!!
!! Forms:
!! 1.     rung-bond, sigma_3
!! 2.     rung-bond, sigma_0
!!
SUBROUTINE init()
  USE dqmc_complex, lam_ising_=>lam_ising
  USE model
  IMPLICIT NONE
  
  INTEGER orb,site,i,a,b,da,db,orb2,nhop,j,iform,ifield,ndim
  REAL(8) re,im,re2,re3,re4
  COMPLEX(8) ga,lam,ga2,lam2
  
  nising=1
  nphi=1
  nform=2

  OPEN(10,file='job.in')
  READ(10,*) restart
  READ(10,*) proj
  READ(10,*) sun;  sun=1
  
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

  READ(10,*) a,b;
  ALLOCATE(nglobal_ising(nising))
  ALLOCATE(nglobal_phi(nphi))
  IF(nising>0) nglobal_ising(:)=a
  IF(nphi>0) nglobal_phi(:)=b
  
  READ(10,*) randomseed
  
  READ(10,*) norb; norb=2

  READ(10,*) La,Lb;    nsite=La*Lb*norb
  READ(10,*) pbca,pbcb
  IF(La<=2)pbca=.false.
  IF(Lb<=2)pbcb=.false.
  READ(10,*) twista,twistb
  READ(10,*) a0r(1:2)
  READ(10,*) b0r(1:2)
  
  ALLOCATE(rorb(2,norb))
  DO orb=1,norb
    READ(10,*) rorb(1:2,orb)
  END DO

  READ(10,*) V_U, V_V, V_J
  READ(10,*) Vph_hol, Vph_brea, Vph_den_buck, Vph_hop_buck
  READ(10,*) debye_hol,debye_brea, debye_den_buck, debye_hop_buck
  ALLOCATE(dphi(nphi))
  READ(10,*) re,re2,re3,re4
  dphi(1)=re
  
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
  CALL init_kmat()
  CALL set_expk(kmat,dtau)

  ! set up kmat for measurement possibly
  CALL init_kmat()
  OPEN(10,FILE='kmat.dat'); WRITE(10,'(2f18.10)') kmat; CLOSE(10)
  
  !=========================================================

  ! there are 2 kinds of form

  ALLOCATE(form_ising(nising),form_phi(nphi),ndim_form(nform))
  form_ising=(/1/)
  form_phi=(/2/)

  !------------------------------------------------------------
  ALLOCATE(fmat(2,2,nform),mask_form(nsite,nform),nb_form(nsite,2,nform))
  fmat=0d0
  mask_form=.false.
  nb_form=0

  ! form 1
  ndim_form(1)=2
  fmat(1:2,1,1)=(/1d0,0d0/)
  fmat(1:2,2,1)=(/0d0,-1d0/)
  mask_form(1:nsite,1)=.false.
  nb_form(:,:,1)=0
  DO a=1,La
    DO b=1,Lb
      i=label(a,b,1)
      j=label(a,b,2)
      mask_form(i,1)=.true.
      nb_form(i,1:2,1)=(/i,j/)
    END DO
  END DO

  ! form 2
  ndim_form(2)=2
  fmat(1:2,1,2)=(/1d0,0d0/)
  fmat(1:2,2,2)=(/0d0,1d0/)
  mask_form(:,2)=mask_form(:,1)
  nb_form(:,:,2)=nb_form(:,:,1)
  
  !----------------------------------------------------------------------------
  ALLOCATE(isingmax(nising),lam_ising(2,nising),gamma_ising(2,nising))
  
  ! V_U
    isingmax(1)=2
      CALL HS1(ga,lam,exp(dtau*V_U/2))
      DO a=1,La; DO b=1,Lb; i=label(a,b,1)
        gamma_ising(1,1)=ga
        gamma_ising(2,1)=ga
        lam_ising(1,1)=lam
        lam_ising(2,1)=-lam
      END DO; END DO
 
  ! save lam_ising to lam_ising_ from the module dqmc_complex
  IF(nising>0)THEN
    ALLOCATE(lam_ising_(maxval(isingmax),nsite,nising))
    DO site=1,nsite
      DO ifield=1,nising
        lam_ising_(1:isingmax(ifield),site,ifield)=lam_ising(1:isingmax(ifield),ifield)
      END DO
    END DO
  END IF

  !--------------------------------------------------------------
  ALLOCATE(gamma_phi(nphi))
  gamma_phi(1)=exp(-1d0)

  IF(nising>0)ALLOCATE(mask_ising(nising))
  IF(nphi>0)ALLOCATE(mask_phi(nphi))

  mask_ising(:)=.false.
  mask_phi(:)=.false.
  IF(abs(V_U)>1d-6)mask_ising(1)=.true.  ! U
  !IF(abs(V_V)>1d-6)mask_ising(2:3)=.true.  ! V
  !IF(abs(V_J)>1d-6)mask_ising(4:11)=.true.  ! J
  IF(abs(Vph_hol)>1d-6)mask_phi(1)=.true.
  !IF(abs(Vph_brea)>1d-6)mask_phi(2:3)=.true.
  !IF(abs(Vph_den_buck)>1d-6)mask_phi(4:5)=.true.
  !IF(abs(Vph_hop_buck)>1d-6)mask_phi(6:9)=.true.


  IF(nphi>0)THEN
    ALLOCATE(gph_x2(nphi),gph_p2(nphi))
    IF(abs(Vph_hol)>1d-6)gph_x2(1)=1d0/2/dtau/Vph_hol
    IF(abs(Vph_hol)>1d-6)gph_p2(1)=1d0/dtau**3/Vph_hol/debye_hol**2
  !  IF(abs(Vph_brea)>1d-6)gph_x2(2:3)=1d0/2/dtau/Vph_brea
  !  IF(abs(Vph_brea)>1d-6)gph_p2(2:3)=1d0/dtau**3/Vph_brea/debye_brea**2
  !  IF(abs(Vph_den_buck)>1d-6)gph_x2(4:5)=1d0/2/dtau/Vph_den_buck
  !  IF(abs(Vph_den_buck)>1d-6)gph_p2(4:5)=1d0/dtau**3/Vph_den_buck/debye_den_buck**2
  !  IF(abs(Vph_hop_buck)>1d-6)gph_x2(6:9)=1d0/2/dtau/Vph_hop_buck
  !  IF(abs(Vph_hop_buck)>1d-6)gph_p2(6:9)=1d0/dtau**3/Vph_hop_buck/debye_hop_buck**2
  END IF
  
  !----------------------------------------------------------
  poolsize_r=0
  poolsize_z=10

  !------------------------------------------------------------
  IF(proj)THEN
    DO i=1,nphi
      IF(mask_phi(i))THEN
        IF(id==0)PRINT*,'phonon can only be added when T>0?'
        CALL exit
      END IF
    END DO
  END IF

  IF(id==0) CALL print_input()

END SUBROUTINE

SUBROUTINE print_input()
  USE dqmc_complex, lam_ising_=>lam_ising
  USE model
  IMPLICIT NONE
  INTEGER i,j,site,a,b

  OPEN(10,FILE='input.dat')
  WRITE(10,*),'nform,nising,nphi=',nform,nising,nphi

  WRITE(10,'(1A12,100I3)') 'form_ising=',form_ising
  WRITE(10,'(1A12,100I3)') 'form_phi=',form_phi
  WRITE(10,'(1A12,100I3)') 'isingmax=',isingmax
  WRITE(10,'(1A12,100L2)') 'mask_ising=',mask_ising
  WRITE(10,'(1A12,100L2)') 'mask_phi=',mask_phi

  ! print form factors
  WRITE(10,*) '------ mask_form -----------------'
  DO i=1,nform
    WRITE(10,*),'form',i
    DO a=1,La
      DO b=1,Lb
        site=label(a,b,1)
        WRITE(10,'(1L2)',advance='no') mask_form(site,i)
      END DO
      WRITE(10,*)
    END DO
  END DO
  
  DO i=1,nform
    WRITE(10,*),'form',i
    DO a=1,La
      DO b=1,Lb
        site=label(a,b,2)
        WRITE(10,'(1L2)',advance='no') mask_form(site,i)
      END DO
      WRITE(10,*)
    END DO
  END DO


  WRITE(10,*) '----------- nb_form ---------------'
  DO i=1,nform
    DO j=1,ndim_form(i)
      WRITE(10,*) 'form',i,'dim',j
      DO a=1,La
        DO b=1,Lb
          site=label(a,b,1)
          WRITE(10,'(1I4)',advance='no') nb_form(site,j,i)
        END DO
        WRITE(10,*)
      END DO
    END DO
  END DO

  WRITE(10,*) '---------- fmat -------------------'
  DO i=1,nform
    WRITE(10,*) 'form',i
    DO j=1,ndim_form(i)
      WRITE(10,'(100f6.1)') fmat(j,1:ndim_form(i),i)
    END DO
  END DO

  WRITE(10,*) '---------- lam_ising ---------------'
  DO i=1,nising
    WRITE(10,*) 'ising field',i
    WRITE(10,'(100F12.6)') lam_ising(1:isingmax(i),i)
  END DO

  WRITE(10,*) '----------- gamma_ising ---------------'
  DO i=1,nising
    WRITE(10,*) 'ising field',i
    WRITE(10,'(100F12.6)') gamma_ising(1:isingmax(i),i)
  END DO

  WRITE(10,*) '----------- gamma_phi --------------'
  DO i=1,nphi
    WRITE(10,'(1A12,1I4,2F12.6)') 'phi field',i,gamma_phi(i)
  END DO

  CLOSE(10)
END SUBROUTINE
