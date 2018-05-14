MODULE model
! SU(N) Hubbard model

  INTEGER sun

  INTEGER norb

  REAL(8) a0r(2),b0r(2)

  REAL(8) a0k(2),b0k(2)
  
  REAL(8), ALLOCATABLE :: rorb(:,:)  ! (2,norb)
  
  REAL(8), ALLOCATABLE :: U(:)    ! (norb)
  
  INTEGER cuta,cutb
  
  COMPLEX(8), ALLOCATABLE :: hop(:,:,:,:)  ! (cuta,cutb,norb,norb)

  INTEGER La,Lb

  REAL(8) beta

  REAL(8) dtau

  LOGICAL pbca,pbcb

  REAL(8) twista,twistb

  COMPLEX(8), ALLOCATABLE :: lam_ising(:,:), gamma_ising(:,:)

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
            boundary=exp(cmplx(0d0,-twista*twopi))
            a2=a2+La
          END IF

          IF(a2>La)THEN
            boundary=exp(cmplx(0d0,twista*twopi))
            a2=a2-La
          END IF

          IF(b2<1)THEN
            boundary=boundary*exp(cmplx(0d0,-twistb*twopi))
            b2=b2+Lb
          END IF

         IF(b2>Lb)THEN
            boundary=boundary*exp(cmplx(0d0,twistb*twopi))
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
  INTEGER newising,site,time,ifield
  ratio=ratio**sun*gamma_ising(newising,site)/gamma_ising(ising(site,time,ifield),site)
END SUBROUTINE

SUBROUTINE acceptprob_phi(ratio,newphi,site,time,ifield)
  USE model
  USE dqmc_complex
  IMPLICIT NONE
  COMPLEX(8) ratio,newphi
  INTEGER site,time,ifield
END SUBROUTINE

SUBROUTINE generate_newising_global(ifield)
  USE model
  USE dqmc_complex
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield
  INTEGER time,site,tmp(nsite)
  DO time=1,ntime-1
    DO site=1,nsite
      ising(site,time,ifield)=irand(isingmax(ifield))+1
    END DO
  END DO
END SUBROUTINE

SUBROUTINE generate_newphi_global(ifield)
  USE model
  USE dqmc_complex
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield
  IF(id==0)PRINT*,'undefined subroutine'
  CALL exit(0)
END SUBROUTINE

SUBROUTINE acceptprob_ising_global(ratio,newising,ifield)
  USE model
  USE dqmc_complex
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield,newising(nsite,ntime)
  COMPLEX(8) ratio
  INTEGER site,time
  ratio=ratio**sun
  DO site=1,nsite
    DO time=1,ntime
      ratio=ratio*gamma_ising(newising(site,time),site)/gamma_ising(ising(site,time,ifield),site)
    END DO
  END DO
END SUBROUTINE

SUBROUTINE acceptprob_phi_global(ratio,newphi,ifield)
  USE model
  USE dqmc_complex
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield
  COMPLEX(8), INTENT(IN) :: newphi(nsite,ntime)
  COMPLEX(8) ratio
  IF(id==0)PRINT*,'undefined subroutine'
  CALL exit(0)
END SUBROUTINE
