MODULE model
! SU(N) Hubbard model

  INTEGER sun

  INTEGER norb

  REAL(8) a0r(2),b0r(2)

  REAL(8) a0k(2),b0k(2)
  
  REAL(8), ALLOCATABLE :: rorb(:,:)  ! (2,norb)
 
  REAL(8) V_U,V_V,V_J
  REAL(8) Vph_hol, Vph_brea, Vph_den_buck, Vph_hop_buck
  REAL(8) debye_hol, debye_brea, debye_den_buck, debye_hop_buck
 
  ! about phonon, Vph=g^2/(M*debye^2), see Johnston etal, 2013
  
  REAL(8), ALLOCATABLE :: gph_x2(:),gph_p2(:)  !
  ! gph_x2=1/2/dtau/V
  ! gph_p2=1/dtau^3/Omega^2/V

  INTEGER cuta,cutb
  
  COMPLEX(8), ALLOCATABLE :: hop(:,:,:,:)  ! (cuta,cutb,norb,norb)

  INTEGER La,Lb

  REAL(8) beta

  REAL(8) dtau

  LOGICAL pbca,pbcb

  REAL(8) twista,twistb

  COMPLEX(8), ALLOCATABLE :: lam_ising(:,:), gamma_ising(:,:)

  COMPLEX(8), ALLOCATABLE :: gamma_phi(:)  ! size of nphi

  COMPLEX(8), ALLOCATABLE :: kmat(:,:)

CONTAINS

  INTEGER FUNCTION label(a,b,orb)
    IMPLICIT NONE
    INTEGER a,b,orb
    label=(a-1)*Lb*norb+(b-1)*norb+orb
  END FUNCTION

  SUBROUTINE init_kmat()
    IMPLICIT NONE
    REAL(8), PARAMETER :: twopi=acos(-1d0)*2
    INTEGER a,b,orb,i
    INTEGER a2,b2,orb2,i2
    INTEGER da,db
    COMPLEX(8) boundary

    kmat=0d0

    DO a=1,La; DO b=1,Lb; DO orb=1,norb; i=label(a,b,orb)
      DO da=-cuta,cuta; DO db=-cutb,cutb; DO orb2=1,norb
        a2=a+da; b2=b+db

        IF(a2>=1.and.a2<=La.and.b2>=1.and.b2<=Lb)THEN ! without crossing the boundary

          i2=label(a2,b2,orb2)
          kmat(i,i2)=hop(da,db,orb,orb2)
          !kmat(i2,i)=conjg(kmat(i,i2))

        ELSE ! crossing the boundary
          
          IF(a2<1.and.(.not.pbca))CYCLE
          IF(b2<1.and.(.not.pbcb))CYCLE

          IF(a2>La.and.(.not.pbca))CYCLE
          IF(b2>Lb.and.(.not.pbcb))CYCLE

          boundary=1d0
          
          IF(a2<1)THEN
            boundary=exp(dcmplx(0d0,-twista*twopi))
            a2=a2+La
          END IF

          IF(a2>La)THEN
            boundary=exp(dcmplx(0d0,twista*twopi))
            a2=a2-La
          END IF

          IF(b2<1)THEN
            boundary=boundary*exp(dcmplx(0d0,-twistb*twopi))
            b2=b2+Lb
          END IF

          IF(b2>Lb)THEN
            boundary=boundary*exp(dcmplx(0d0,twistb*twopi))
            b2=b2-Lb
          END IF

          IF(a2==a.and.b2==b)CYCLE

          i2=label(a2,b2,orb2)
          kmat(i,i2)=hop(da,db,orb,orb2)*boundary
          !kmat(i2,i)=conjg(kmat(i,i2))

        END IF
 
      END DO; END DO; END DO
    END DO; END DO; END DO

  END SUBROUTINE

END MODULE

SUBROUTINE acceptprob_ising(ratio,newising,site,time,ifield)
  USE model
  USE dqmc_complex
  IMPLICIT NONE
  COMPLEX(8) ratio
  INTEGER newising,site,time,ifield,oldising
  INTEGER newising_array(nsite,ntime)
  ratio=ratio**sun*gamma_ising(newising,ifield)/gamma_ising(ising(site,time,ifield),ifield)
!  print*,newising,site,time,ifield,ratio
!  read*
  oldising=ising(site,time,ifield)
  print*,'site,time,ifield=',site,time,ifield
  print*,'ising:    ratio=',ratio
  newising_array(:,:)=ising(:,:,ifield)
  newising_array(site,time)=newising
!read*;return
  open(21,file='newising.dat',form='unformatted')
  write(21)newising_array
  close(21)
  call dqmc_update_global(ifield,'ising')

  ising(site,time,ifield)=oldising
  

END SUBROUTINE

SUBROUTINE acceptprob_phi(ratio,newphi,site,time,ifield)
  USE model
  USE dqmc_complex
  IMPLICIT NONE
  COMPLEX(8) ratio,newphi,diffphi,sumphi,xx
  INTEGER site,time,ifield
  COMPLEX(8) newphi_array(nsite,ntime),oldphi
  diffphi=newphi-phi(site,time,ifield)
  sumphi=newphi+phi(site,time,ifield)
  xx=phi(site,mod(time,ntime)+1,ifield)+phi(site,mod(time-2+ntime,ntime)+1,ifield)
  ratio=ratio**sun*exp( -diffphi*( gph_x2(ifield)*sumphi + gph_p2(ifield)*(sumphi-xx) ) ) &
  & *gamma_phi(ifield)**diffphi
!  print*,'phimax,phimin=',maxval(real(phi)),minval(real(phi))
!  print*,'newphi,site,time,ifield=',newphi,site,time,ifield
!  print*,'ratio=',ratio
!  read*
  
  oldphi=phi(site,time,ifield)
  print*,'phi:      ratio=',ratio
  newphi_array(:,:)=phi(:,:,ifield)
  newphi_array(site,time)=newphi

  open(21,file='newphi.dat',form='unformatted')
  write(21)newphi_array
  close(21)
  call dqmc_update_global(ifield,'phi  ')

  phi(site,time,ifield)=oldphi
 
END SUBROUTINE

SUBROUTINE generate_newising_global(ifield)
  USE model
  USE dqmc_complex
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield
  INTEGER time,site
  open(21,file='newising.dat',form='unformatted')
  read(21)ising(:,:,ifield)
  close(21)
  return

  
  DO site=1,nsite; IF(.not.mask_form(site,form_ising(ifield)))CYCLE
    DO time=1,ntime
      ising(site,time,ifield)=irand(isingmax(ifield))+1
    END DO
  END DO
END SUBROUTINE

SUBROUTINE generate_newphi_global(ifield)
! global update phi(i,t)->phi(i,t)+dphi_ according to Johnston et al. 2013
  USE model
  USE dqmc_complex
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield
  COMPLEX(8) dphi_
  INTEGER site
  open(21,file='newphi.dat',form='unformatted')
  read(21)phi(:,:,ifield)
  close(21)
  return

  
  dphi_=dphi(ifield)*drand_sym()
  DO site=1,nsite; IF(.not.mask_form(site,form_phi(ifield)))CYCLE
    phi(site,:,ifield)=phi(site,:,ifield)+dphi_
  END DO
  !print*,'dphi_=',dphi_,maxval(real(phi(:,:,ifield))),minval(real(phi(:,:,ifield)))
END SUBROUTINE

SUBROUTINE acceptprob_ising_global(ratio,newising,ifield)
  USE model
  USE dqmc_complex
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield,newising(nsite,ntime)
  COMPLEX(8) ratio
  INTEGER site,time
  ratio=ratio**sun
  DO site=1,nsite; IF(.not.mask_form(site,form_ising(ifield)))CYCLE
    DO time=1,ntime
      ratio=ratio*gamma_ising(newising(site,time),ifield)/gamma_ising(ising(site,time,ifield),ifield)
    END DO
  END DO
  print*,'ratio by global:',ratio
  read*
END SUBROUTINE

SUBROUTINE acceptprob_phi_global(ratio,newphi,ifield)
! global update phi(i,t)->phi(i,t)+dphi_ according to Johnston et al. 2013
! Such kind of global update is quite simple since the phonon kinetic energy does not change.
  USE model
  USE dqmc_complex
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield
  COMPLEX(8), INTENT(IN) :: newphi(nsite,ntime)
  COMPLEX(8) ratio
  INTEGER site,time
  ratio=ratio**sun
  DO site=1,nsite; IF(.not.mask_form(site,form_phi(ifield)))CYCLE
    DO time=1,ntime
      ratio=ratio*exp( -gph_x2(ifield)*(newphi(site,time)+phi(site,time,ifield))*(newphi(site,time)-phi(site,time,ifield)) ) &
      &          *gamma_phi(ifield)**(newphi(site,time)-phi(site,time,ifield))
      !===============
      
      ratio=ratio*exp( -0.5*gph_p2(ifield)*( (newphi(site,mod(time,ntime)+1)-newphi(site,time))**2 &
      &                                     -(phi(site,mod(time,ntime)+1,ifield)-phi(site,time,ifield))**2 ))

      !================
    END DO
  END DO
  print*,'ratio by global:',ratio
  read*
END SUBROUTINE
