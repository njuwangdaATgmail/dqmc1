!> do global update
SUBROUTINE dqmc_update_global(ifield)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifield
  INTEGER i,p,flag,j
  INTEGER, ALLOCATABLE :: ising_old(:,:),ising_new(:,:)
  COMPLEX(8), ALLOCATABLE :: phi_old(:,:),phi_new(:,:)
  COMPLEX(8) ratio(nflv),dvec(nsite),dvec_old(nsite)
  COMPLEX(8) qmat(nsite,nsite),rmat(nsite,nsite),Bmat(nsite,nsite)
  
  DO flv=1,nflv

    ratio=1d0

  ! get det(P'BP) or det(1+B). In practice, we save their eigenvalues
    dvec=1d0
    qmat=0d0
    Bmat=0d0
    DO i=1,nsite
      qmat(i,i)=1d0
      Bmat(i,i)=1d0
    END DO
    
    flag=0
    DO p=1,ntime
      CALL evolve_left(p,qmat,nsite,flv)
      flag=flag+1
      IF(flag/=ngroup)CYCLE
      flag=0
      DO i=1,nsite
        qmat(:,i)=qmat(:,i)*dvec(i)
      END DO
      CALL zqdr(nsite,nsite,qmat,rmat,dvec)
      Bmat=matmul(rmat,Bmat)
    END DO
    ! B*...*B=qmat*dvec*Bmat

    IF(proj)THEN  !det(P'*B*P)
    
      qmat(1:nelec,1:nsite)=matmul(conjg(transpose(slater)),qmat)
      rmat(1:nsite,1:nelec)=matmul(Bmat,slater)
      DO i=1,nelec
        DO j=1,nelec
          Bmat(i,j)=sum(qmat(i,:)*dvec(:)*rmat(:,j))
        END DO
      END DO
      CALL qdr(nelec,nelec,Bmat(1:nelec,1:nelec),rmat(1:nelec,1:nelec),dvec_old)
      ratio(flv)=ratio(flv)/det(nelec,Bmat(1:nelec,1:nelec))
      CALL ordering(nelec,dvec_old(1:nelec))
  
    ELSE   !det(1+B)
    
      ratio(flv)=ratio(flv)/det(nsite,qmat)
      CALL inverse(nsite,qmat)
      DO i=1,nsite
        qmat(i,:)=qmat(i,:)+dvec(i)*Bmat(i,:)
      END DO
      CALL qdr(nsite,nsite,qmat,rmat,dvec_old)
      ratio(flv)=ratio(flv)/det(nsite,qmat)
      CALL ordering(nsite,dvec_old)

    END IF
  
  ! generate newising or newphi externally
  IF(job=='ising')THEN
    ALLOCATE(ising_old(nsite,ntime),ising_new(nsite,ntime))
    ising_old=ising(:,:,ifield)
    CALL generate_newising_global(ifield)
  ELSE
    ALLOCATE(phi_old(nsite,ntime),phi_new(nsite,ntime))
    phi_old=phi(:,:,ifield)
    CALL generate_newphi_global(ifield)
  END IF

  ! get new det(P'BP) or det(1+B)
    dvec=1d0
    qmat=0d0
    Bmat=0d0
    DO i=1,nsite
      qmat(i,i)=1d0
      Bmat(i,i)=1d0
    END DO
    
    flag=0
    DO p=1,ntime
      CALL evolve_left(p,qmat,nsite)
      flag=flag+1
      IF(flag/=ngroup)CYCLE
      flag=0
      DO i=1,nsite
        qmat(:,i)=qmat(:,i)*dvec(i)
      END DO
      CALL zqdr(nsite,nsite,qmat,rmat,dvec)
      Bmat=matmul(rmat,Bmat)
    END DO
    ! B*...*B=qmat*dvec*Bmat

  IF(proj)THEN  !det(P'*B*P)
    
    qmat(1:nelec,1:nsite)=matmul(conjg(transpose(slater)),qmat)
    rmat(1:nsite,1:nelec)=matmul(Bmat,slater)
    DO i=1,nelec
      DO j=1,nelec
        Bmat(i,j)=sum(qmat(i,:)*dvec(:)*rmat(:,j))
      END DO
    END DO
    CALL qdr(nelec,nelec,Bmat(1:nelec,1:nelec),rmat(1:nelec,1:nelec),dvec)
    ratio=ratio*det(nelec,Bmat(1:nelec,1:nelec))
    CALL ordering(nelec,dvec(1:nelec))
  
  ELSE   !det(1+B)
    
    ratio=ratio*det(nsite,qmat)
    CALL inverse(nsite,qmat)
    DO i=1,nsite
      qmat(i,:)=qmat(i,:)+dvec(i)*Bmat(i,:)
    END DO
    CALL qdr(nsite,nsite,qmat,rmat,dvec)
    ratio=ratio*det(nsite,qmat)
    CALL ordering(nsite,dvec)

  END IF

  ! get accept ratio from the determinant ratio externally
  IF(proj)THEN
    DO i=1,nelec
      ratio=ratio*dvec(i)/dvec_old(i)
    END DO
  ELSE
    DO i=1,nsite
      ratio=ratio*dvec(i)/dvec_old(i)
    END DO
  END IF
  IF(job=='ising')THEN
    ising_new(:,:)=ising(:,:,ifield)
    ising(:,:,ifield)=ising_old(:,:)
    CALL acceptprob_ising_global(ratio,ising_new,ifield)
  ELSE  ! job='phi'
    phi_new(:,:)=phi(:,:,ifield)
    phi(:,:,ifield)=phi_old(:,:)
    CALL acceptprob_phi_global(ratio,phi_new,ifield)
  END IF
  
  ! if accepted, update ising or phi
  ! Here, we don't need to update Green's function since we will 
  ! always do that in the following step of local update
  IF(job=='ising')THEN
    Ntotal_ising_global(ifield)=Ntotal_ising_global(ifield)+1
  ELSE
    Ntotal_phi_global(ifield)=Ntotal_phi_global(ifield)+1
  END IF

  IF(drand()>abs(ratio))RETURN
  
  IF(job=='ising')THEN
    Naccept_ising_global(ifield)=Naccept_ising_global(ifield)+1
  ELSE
    Naccept_phi_global(ifield)=Naccept_phi_global(ifield)+1
  END IF

  currentphase=currentphase*ratio/abs(ratio)

  IF(job=='ising')THEN
    ising(:,:,ifield)=ising_new
  ELSE
    phi(:,:,ifield)=phi_new
  END IF

END SUBROUTINE
