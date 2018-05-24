!> do global update
SUBROUTINE dqmc_update_global(ifield)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield
  INTEGER flv,ifield
  INTEGER, ALLOCATABLE :: ising_old(:,:),ising_new(:,:)
  REAL(8) paccept
  COMPLEX(8) oldfield(nsite,ntime),newfield(nsite,ntime)
  COMPLEX(8) ratio(nflv),rtot
  COMPLEX(8), ALLOCATABLE :: dvec(:,:),dvec_old(:,:)
  COMPLEX(8), ALLOCATABLE :: dvec(:,:),dvec_old(:,:)
  
 
  ! calculate old determinant from scratch
  IF(proj)THEN
    CALL update_scratch_T0(1)
  ELSE
    CALL update_scratch(1)
  END IF

  ! save eigenvalues of the old determinant
  DO flv=1,nflv

    IF(proj)THEN

      ratio(flv)=1d0/det(nelec,Bstring_Q(1:nelec,1:nelec,1,flv))
      dvec_old(1:nelec,flv)=Bstring_D(1:nelec,1,flv)
      CALL sort(nelec,dvec_old(1:nelec,flv))

    ELSE

      ratio(flv)=1d0/det(nsite,Bstring_Q(1:nsite,1:nsite,1,flv))
      dvec_old(1:nsite,flv)=Bstring_D(1:site,1,flv)
      CALL sort(nsite,dvec_old(1:nsite,flv))

    END IF

  END DO

  ! save old field
  oldfield=field(:,:,ifield)

  ! generate new field
  CALL generate_newfield_global(ifield)

  ! calculate new determinant from scratch
  IF(proj)THEN
    CALL update_scratch_T0(1)
  ELSE
    CALL update_scratch(1)
  END IF
  
  ! get eigenvalues of the new determinant
  DO flv=1,nflv

    IF(proj)THEN

      ratio(flv)=ratio(flv)*det(nelec,Bstring_Q(1:nelec,1:nelec,1,flv))
      dvec(1:nelec,flv)=Bstring_D(1:nelec,1,flv)
      CALL sort(nelec,dvec(1:nelec,flv))

      DO i=1,nelec
        ratio(flv)=ratio(flv)*dvec(i,flv)/dvec_old(i,flv)
      END DO

    ELSE

      ratio(flv)=ratio(flv)*det(nsite,Bstring_Q(1:nsite,1:nsite,1,flv))
      dvec(1:nsite,flv)=Bstring_D(1:site,1,flv)
      CALL sort(nsite,dvec(1:nsite,flv))

      DO i=1,nsite
        ratio(flv)=ratio(flv)*dvec(i,flv)/dvec_old(i,flv)
      END DO

    END IF

  END DO

  ! get acceptance ratio
  newfield=field(:,:,ifield)
  field(:,:,ifield)=oldfield
  CALL acceptprob_global(ratio,newfield,ifield,rtot)

  ! correction of Metropolis ratio
  IF(abs(rtot)<1d0)THEN
    paccept=abs(rtot)/(1d0+newMetro*abs(rtot))
  ELSE
    paccept=abs(rtot)/(newMetro+abs(rtot))
  END IF

  ! record the total number of tryings
  Ntotal_field_global(ifield)=Ntotal_field_global(ifield)+1

  ! whether update_scratch is useful
  scratch_global_useful=.false.
  
  ! accept the new configuration with probability paccept
  IF(drand()>paccept)RETURN
  
  ! record the number of accepted tryings
  Naccept_field_global(ifield)=Naccept_ising_global(ifield)+1

  ! obtain the phase of the current configuration
  currentphase=currentphase*rtot/abs(rtot)

  ! update field
  field(:,:,ifield)=newfield

  ! whether update_scratch useful
  scratch_global_useful=.true.

CONTAINS
  
  !> sort the array a(n) from large to small absulute values
  SUBROUTINE sort(n,a)
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
  
END SUBROUTINE
