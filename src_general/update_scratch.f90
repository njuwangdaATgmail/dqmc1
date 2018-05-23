!> evaluate Green's function from definition directly using the Bstring technique
!  to accelerate (a cheap realization, may not be the most efficient)
!
! For T>0 algorithm:
! define Bstring(nsite,nsite,nblock)
! when time=1 or block=1,
!   used  :   none
!   update:   Bstring(:,:,nblock)=BBB(nblock)
!             Bstring(:,:,nblock-1)=BBB(nblock)*BBB(nblock-1)
!             ...
!             Bstring(:,:,2)=BBB(nblock)*...*BBB(2)
!             Bstring(:,:,1)=BBB(nblock)*...*BBB(1)+1 which further gives the determinant
! when time=nscratch+1 or block=2,
!   used  :   Bstring(:,:,block)
!   update:   Bstring(:,:,1)=BBB(1)
!
! when nblock>block>2,
!   used  :   Bstring(:,:,block)
!             Bstring(:,:,block-2)
!   update:   Bstring(:,:,block-1)=BBB(block-1)*Bstring(block-2)
! ...
! when block=nblock
!   used  :   Bstring(:,:,nblock)
!             Bstring(:,:,nblock-2)
!   update:   none
!
! For T=0 algorithm:
! define Bstring(nsite,nelec,nblock)
! when time=1 or block=1,
!   used  :   none
!   update:   Bstring(:,:,nblock)=transpose( P'*BBB(nblock) )
!             Bstring(:,:,nblock-1)=transpose( P'*BBB(nblock)*BBB(nblock-1) )
!             ...
!             Bstring(:,:,2)=transpose( P'*BBB(nblock)*...*BBB(2) )
!             Bstring(1:nelec,:,1)=P'*BBB(nblock)*...*BBB(1)*P which further gives the determinant
! when time=nscratch+1 or block=2,
!   used  :   Bstring(:,:,block)
!   update:   Bstring(:,:,1)=BBB(1)*P
!
! when nblock>block>2,
!   used  :   Bstring(:,:,block)
!             Bstring(:,:,block-2)
!   update:   Bstring(:,:,block-1)=BBB(block-1)*Bstring(block-2)
! ...
! when block=nblock
!   used  :   Bstring(:,:,nblock)
!             Bstring(:,:,nblock-2)
!   update:   none
!

!> calculate T=0 Green's function from definition, using QR decomposition stabilization algorithm.
SUBROUTINE update_scratch_T0(time)
  USE mod_dqmc_complex
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: time
  INTEGER flag,p,i,flv
  COMPLEX(8) t(nelec,nelec),lmat(nelec,nsite),rmat(nsite,nelec),gfast(nsite,nsite)

  ! get the block index: [1,nscratch]->1, [nscratch+1,2*nscratch]->2,...
  ! we come here only when mod(time-1,nscratch)==0
  block=(time-1)/nscratch+1 

  DO flv=1,nflv
    
    gfast=g(:,:,flv)

    
    IF(block==1)THEN
      
      ! correspondingly, time=1
      ! L=B(ntime)...B(1)
      ! save Bstring(Nblock-1), BString(Nblock-2),...,Bstring(1)


      lmat=conjg(transpose(slater(:,:,flv)))
      CALL zldq(nelec,nsite,lmat,tri,dvec)     ! we can save ldq of slater
      flag=0
      DO p=ntime,1,-1
        CALL evolve_right(p,lmat,nelec,flv,.false.)
        flag=flag+1
        IF(flag==ngroup)THEN
          flag=0
          DO i=1,nelec
            lmat(i,:)=dvec(i)*lmat(i,:)
          END DO
          CALL zldq(nelec,nsite,lmat,t,dvec)
          tri=matmul(tri,t)
        END IF
        IF(mod(p-1,nscratch)=0.and.p>1)THEN
          pblock=(p-1)/nscratch
          Bstring_Q(:,:,pblock)=transpose(lmat)
          Bstring_D(:,pblock)=dvec
          Bstring_R(:,:,pblock)=transpose(tri)
        END IF
      END DO

      rmat=slater(:,:,flv)
      CALL zqdr(nsite,nelec,rmat,tri,devc)    


    ELSE

      ! L=Bstring(block-1)
      ! R=B(time-1)...B(time-nscratch)*[Bstring(1) or 1(if block==2)]
      ! Bstring(1)=R

    END IF

    
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
