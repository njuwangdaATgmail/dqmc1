!> evaluate Green's function from definition directly 
! Bstring technique is used to accelerate (a cheap realization, may not be the most efficient)
!
!------------------------------------------------------------------------------------------------
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
!   update:   Bstring(:,:,1)=BBB+1
!-----------------------------------------------------------------------------------------------
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
!   update:   Bstring(:,:,1)=P'BBB*P
!-------------------------------------------------------------------------------------------------

!> calculate T=0 Green's function from definition, using QDR decomposition stabilization algorithm.
SUBROUTINE update_scratch_T0(time)
  USE mod_dqmc_complex
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: time
  INTEGER block,flv,flag,p,i
  COMPLEX(8) dvec(nelec),tri(nelec,nelec),tmat(nelec,nelec),gfast(nsite,nsite)
  COMPLEX(8), ALLOCATABLE :: qmat(:,:)
  
  IF(scratch_global_useful)THEN
    scratch_global_useful=.false.
    RETURN
  END IF
 
  ! get the block index: [1,nscratch]->1, [nscratch+1,2*nscratch]->2,...
  ! we come here only when mod(time-1,nscratch)==0
  block=(time-1)/nscratch+1 

  DO flv=1,nflv
    
    ! save g to evaluate the error of fast-update
    gfast=g(:,:,flv)
    
    IF(block==1)THEN

      ALLOCATE(qmat(nelec,nsite))
      
      ! tmat*dvec*qmat=P'
      qmat=conjg(transpose(slater_Q(:,:,flv)))
      dvec=conjg(slater_D(:,flv))
      tmat=conjg(transpose(slater_R(:,:,flv)))

      flag=0
      DO p=ntime,1,-1
        CALL evolve_right(p,qmat,nelec,flv,.false.)
        flag=flag+1
        IF(flag==ngroup.or.mod(p-1,nscratch)==0)THEN
          DO i=1,nsite
            qmat(:,i)=dvec(:)*qmat(:,i)
          END DO
          CALL ldq(nelec,nsite,qmat,tri,dvec)
          tmat=matmul(tmat,tri)
          flag=0
        END IF
        IF(mod(p-1,nscratch)==0.and.p>1)THEN
          pblock=(p-1)/nscratch+1
          Bstring_Q(:,:,pblock,flv)=transpose(qmat)
          Bstring_D(:,pblock,flv)=dvec
          Bstring_T(:,:,pblock,flv)=tmat
        END IF
      END DO
     
      !   g = 1 - Q_R * (Q_L*Q_R)^{-1} * Q_L
      tri=matmul(qmat,slater_Q(:,:,flv))
      CALL inverse(nelec,tri)
      g(:,:,flv)=-matmul(matmul(slater_Q(:,:,flv),tri),qmat)
      DO i=1,nsite
        g(i,i,flv)=g(i,i,flv)+1d0
      END DO

      ! Bstring_T(1)*Bstring_D(1)*Bstring_Q(1) = P'*BBB*P
      DO i=1,nsite
        qmat(:,i)=dvec(:)*qmat(:,i)
      END DO
      qmat(1:nelec,1:nelec)=matmul(qmat,slater(:,:,flv))
      CALL ldq(nelec,nelec,qmat(1:nelec,1:nelec),tri,dvec)
      Bstring_Q(1:nelec,1:nelec,1,flv)=qmat(1:nelec,1:nelec)
      Bstring_D(:,1,flv)=dvec
      Bstring_T(:,:,1,flv)=matmul(tmat,tri)

      DEALLOCATE(qmat)

    ELSE

      ALLOCATE(qmat(nsite,nelec))
      
      IF(block==2)THEN
        ! qmat*dvec*tmat = slater
        qmat=slater_Q(:,:,flv)
        dvec=slater_D(:,flv)
        tmat=slater_R(:,:,flv)
      ELSE
        ! qmat*dvec*tmat = Bstring(block-2)
        qmat=Bstring_Q(:,:,block-2,flv)
        dvec=Bstring_D(:,block-2,flv)
        tmat=Bstring_T(:,:,block-2,flv)
      END IF

      flag=0
      DO p=(block-2)*nscratch+1,(block-1)*nscratch
        CALL evolve_left(p,qmat,nelec,flv,.false.)
        flag=flag+1
        IF(flag==ngroup.or.mod(p,nscratch)==0)THEN
          DO i=1,nelec
            qmat(:,i)=qmat(:,i)*dvec(i)
          END DO
          CALL qdr(nsite,nelec,qmat,tri,dvec)
          tmat=matmul(tri,tmat)
          flag=0
        END IF
      END DO
      Bstring_Q(:,:,block-1,flv)=qmat
      Bstring_D(:,block-1,flv)=dvec
      Bstring_T(:,:,block-1,flv)=tmat

      ! g = 1 - qmat*dvec*tmat * (BL_T*BL_D*BL_Q*qmat*dvec*tmat)^{-1} * BL_T*BL_D*BL_Q
      !   = 1 - qmat * (BL_Q*qmat)^{-1} * BL_Q

      tri=matmul(transpose(Bstring_Q(:,:,block,flv)),qmat)
      CALL inverse(nelec,tri)
      g(:,:,flv)=-matmul(matmul(qmat,tri),transpose(Bstring_Q(:,:,block,flv)))
      DO i=1,nsite
        g(i,i,flv)=g(i,i,flv)+1d0
      END DO

      DEALLOCATE(qmat)

    END IF
   
    IF(time>1)THEN
      CALL evolve_left_K(gfast,nsite,flv,.false.,.false.)
      CALL evolve_right_K(gfast,nsite,flv,.true.,.false.)
      err_fast=max(maxval(abs(gfast-g(:,:,flv))),err_fast)
    END IF

  END DO

END SUBROUTINE

!> calculate T>0 Green's function from definition, using QDR decomposition stabilization algorithm.
SUBROUTINE update_scratch(time)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: time
  INTEGER block,flv,flag,p,i
  COMPLEX(8) dvec(nsite),qmat(nsite,nsite),tmat(nsite,nsite),tri(nsite,nsite),gfast(nsite,nsite)
  
  IF(scratch_global_useful)THEN
    scratch_global_useful=.false.
    RETURN
  END IF
 
  ! get the block index: [1,nscratch]->1, [nscratch+1,2*nscratch]->2,...
  ! we come here only when mod(time-1,nscratch)==0
  block=(time-1)/nscratch+1

  DO flv=1,nflv

    gfast=g(:,:,flv)

    IF(block==1)THEN

      ! tmat*dvec*qmat=1
      qmat=0d0
      dvec=1d0
      tmat=0d0
      DO i=1,nsite
        qmat(i,i)=1d0
        tmat(i,i)=1d0
      END DO

      flag=0
      DO p=ntime,1,-1
        CALL evolve_right(p,qmat,nsite,flv,.false.)
        flag=flag+1
        IF(flag==ngroup.or.mod(p-1,nscratch)==0)THEN
          DO i=1,nsite
            qmat(:,i)=dvec(:)*qmat(:,i)
          END DO
          CALL ldq(nsite,nsite,qmat,tri,dvec)
          tmat=matmul(tmat,tri)
          flag=0
        END IF
        IF(mod(p-1,nscratch)==0.and.p>1)THEN
          pblock=(p-1)/nscratch+1
          Bstring_Q(:,:,pblock,flv)=qmat
          Bstring_D(:,pblock,flv)=dvec
          Bstring_T(:,:,pblock,flv)=tmat
        END IF
      END DO

      ! tmat*dvec*qmat = BBB
      ! g^{-1}= BBB + 1 = ( tmat*dvec + qmat^{-1} ) * qmat

      Bstring_Q(:,:,1,flv)=qmat
      CALL inverse(nsite,qmat)
      g(:,:,flv)=qmat
      DO i=1,nsite
        qmat(:,i)=qmat(:,i)+tmat(:,i)*dvec(i)
      END DO
      CALL ldq(nsite,nsite,qmat,tri,dvec)
      Bstring_Q(:,:,1,flv)=matmul(qmat,Bstring_Q(:,:,1,flv))
      Bstring_D(:,1,flv)=dvec
      Bstring_T(:,:,1,flv)=tri
      CALL inverse(nsite,tri)
      CALL inverse(nsite,qmat)
      DO i=1,nsite
        tri(:,i)=tri(:,i)/dvec(:)
      END DO
      g(:,:,flv)=matmul(matmul(g(:,:,flv),qmat),tri)

    ELSE

      IF(block==2)THEN
        ! tmat*dvec*qmat=1
        qmat=0d0
        dvec=1d0
        tmat=0d0
        DO i=1,nsite
          qmat(i,i)=1d0
          tmat(i,i)=1d0
        END DO
      ELSE
        ! tmat*dvec*qmat=Bstring(block-2)
        tmat=Bstring_T(:,:,block-2,flv)
        qmat=Bstring_Q(:,:,block-2,flv)
        dvec=Bstring_D(:,:,block-2,flv)
      END IF

      flag=0
      DO p=(block-2)*nscratch+1,(block-1)*nscratch
        CALL evolve_left(p,qmat,nsite,flv,.false.)
        flag=flag+1
        IF(flag==ngroup.or.mod(p,nscratch)==0)THEN
          DO i=1,nsite
            qmat(:,i)=qmat(:,i)*dvec(i)
          END DO
          CALL qdr(nsite,nsite,qmat,tri,dvec)
          tmat=matmul(tri,tmat)
          flag=0
        END IF
      END DO
      Bstring_Q(:,:,block-1,flv)=qmat
      Bstring_D(:,block-1,flv)=dvec
      Bstring_T(:,:,block-1,flv)=tmat

      ! g^{-1} = qmat*dvec*tmat*BL_T*BL_D_BL_Q + 1

      tmat=matmul(tmat,Bstring_T(:,:,block,flv))
      DO i=1,nsite
        tmat(:,i)=dvec(:)*tmat(:,i)*Bstring_D(i,block,flv)
      END DO
      CALL qdr(nsite,nsite,tmat,tri,dvec)
      qmat=matmul(qmat,tmat)
      tri=matmul(tri,Bstring_Q(:,:,block,flv))

      ! g^{-1} = qmat*dvec*tri + 1

      CALL inverse(nsite,qmat)
      g(:,:,flv)=qmat
      DO i=1,nsite
        qmat(:,i)=qmat(:,i)+dvec(:)*tri(:,i)
      END DO
      CALL qdr(nsite,nsite,qmat,tri,dvec)
      CALL inverse(nsite,qmat)
      CALL inverse(nsite,tri)
      DO i=1,nsite
        tri(:,i)=tri(:,i)/dvec(i)
      END DO
      g(:,:,flv)=matmul(tri,matmul(qmat,g(:,:,flv)))

    END IF

    IF(time>1)THEN  ! when time=1, we may come through a global update, then the following fast update has no meaning
      CALL evolve_left_K(gfast,nsite,flv,.false.,.false.)
      CALL evolve_right_K(gfast,nsite,flv,.true.,.false.)
      err_fast=max(maxval(abs(gfast-g(:,:,flv))),err_fast)
    END IF

  END DO

END SUBROUTINE
