!> calculate T=0 Green's function from definition, using QR decomposition stabilization algorithm.
SUBROUTINE update_scratch_T0(time)
  USE mod_dqmc_complex
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: time
  INTEGER flag,p,i,flv
  COMPLEX(8) t(nelec,nelec),lmat(nelec,nsite),rmat(nsite,nelec),gfast(nsite,nsite)
   
  DO flv=1,nflv
    
    gfast=g(:,:,flv)
    
    lmat=conjg(transpose(slater(:,:,flv)))
    CALL zlq(nelec,nsite,lmat,t)
    flag=0
    DO p=ntime,time,-1
      CALL evolve_right(p,lmat,nelec,flv)
      flag=flag+1
      IF(flag/=ngroup)CYCLE
      flag=0
      CALL zlq(nelec,nsite,lmat,t)
    END DO

    rmat=slater(:,:,flv)
    CALL zqr(nsite,nelec,rmat,t)
    flag=0
    DO p=1,time-1
      CALL evolve_left(p,rmat,nelec,flv)
      flag=flag+1
      IF(flag/=ngroup)CYCLE
      flag=0
      CALL zqr(nsite,nelec,rmat,t)
    END DO

    t=matmul(lmat,rmat)
    CALL inverse(nelec,t)

    g(:,:,flv)=-matmul(rmat,matmul(t,lmat))

    DO i=1,nsite
      g(i,i,flv)=g(i,i,flv)+1d0
    END DO
    
    IF(time>1)THEN
      gfast=matmul(matmul(expk(:,:,flv),gfast),inv_expk(:,:,flv))
      err_fast=max(maxval(abs(gfast-g(:,:,flv))),err_fast)
    END IF

  END DO

END SUBROUTINE

!> calculate T>0 Green's function from definition, using QDR decomposition stabilization algorithm.
SUBROUTINE update_scratch(time)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: time
  INTEGER d,flag,p0,p,i,flv
  COMPLEX(8) dvec(nsite),qmat(nsite,nsite),rmat(nsite,nsite),gtmp(nsite,nsite),gfast(nsite,nsite)
  
  DO flv=1,nflv

    gfast=g(:,:,flv)

    d=nsite
    
    dvec=1d0
    qmat=0d0
    g(:,:,flv)=0d0
    DO i=1,d
      qmat(i,i)=1d0
      g(i,i,flv)=1d0
    END DO
    
    flag=0
    DO p0=time,ntime+time-1
      p=p0;IF(p>ntime)p=p-ntime
      CALL evolve_left(p,qmat,d,flv)
      flag=flag+1
      IF(flag/=ngroup)CYCLE
      flag=0
      DO i=1,d
        qmat(:,i)=qmat(:,i)*dvec(i)
      END DO
      CALL zqdr(d,d,qmat,rmat,dvec)
      g(:,:,flv)=matmul(rmat,g(:,:,flv))
    END DO
    
    CALL inverse(d,qmat)
    DO i=1,d
      g(:,i,flv)=dvec(:)*g(:,i,flv)
    END DO
    g(:,:,flv)=g(:,:,flv)+qmat
    CALL inverse(d,g(:,:,flv))
    g(:,:,flv)=matmul(g(:,:,flv),qmat)

  IF(time>1)THEN  ! when time=1, we may come through a global update, then the following fast update has no meaning
    gfast=matmul(matmul(expk(:,:,flv),gfast),inv_expk(:,:,flv))
    err_fast=max(maxval(abs(gfast-g(:,:,flv))),err_fast)
  END IF

END SUBROUTINE
