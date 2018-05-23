!> do local update trying to update all auxiliary fields at each site and time
SUBROUTINE dqmc_update_local(do_meas)
  USE mod_dqmc_complex
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: do_meas
  INTEGER time,site,ifield,ndim,j,a,b,sitea,siteb,flv
  COMPLEX(8) ratio,newfield,oldfield,delta_ja,gtmp(nsite,nsite)
  COMPLEX(8) delta(maxval(ndim_field),maxval(ndim_field),nflv)
  COMPLEX(8) ratiomat(maxval(ndim_field),maxval(ndim_field),nflv)
  
  DO time=1,ntime
    
    IF(mod(time-1,nscratch)==0)THEN
      
      ! update Green's function from scratch every nscratch steps
      IF(proj)THEN
        CALL update_scratch_T0(time)
      ELSE
        CALL update_scratch(time)
      END IF

    ELSE

      ! update Green's function by fast-update algorithm
      DO flv=1,nflv
        CALL evolve_left_K(g(:,:,flv),nsite,flv,.false.,.false.)
        CALL evolve_right_K(g(:,:,flv),nsite,flv,.true.,.false.)
        !g(:,:,flv)=matmul(matmul(expk(:,:,flv),g(:,:,flv)),inv_expk(:,:,flv))
      END DO

    END IF
    
    IF(proj)THEN

      ! at T=0, do measurement only during [ntime/2+1-nsp,ntime/2+1+nsp]
      IF(do_meas.and.abs(time-ntime/2-1)<=nsp)THEN
        ! prepare the pool for the following measurement
        CALL add_pool()
        ! do measurement by calling external subroutine
        CALL measurement(time)
      END IF

    ELSE
      
      ! at T>0, do measurement only when time<=nsp
      IF(do_meas.and.time<=nsp)THEN
        ! prepare the pool for the following measurement
        CALL add_pool()
        ! do measurement by calling external subroutine
        CALL measurement(time)
      END IF

    ENDIF

    ! update the i-th field locally
    DO ifield=1,nfield; IF(.not.mask_field(ifield))CYCLE  
      
      ndim=ndim_field(ifield)
      
      ! try to update i-th field on each (block-)site
      DO site=1,nsite; IF(.not.mask_field_site(site,ifield))CYCLE
        
        ! record old configuration
        oldfield=field(site,time,i)
        
        ! generate new configuration by calling an external subroutine which also returns delta=exp(V'-V)-1
        CALL generate_newfield(newfield,site,time,delta(1:ndim,1:ndim,1:nflv))

        ! calculate determinant ratio for each flavor
        IF(ndim==1)THEN
          ratio=1d0
          DO flv=1,nflv
            ratiomat(1,1,flv)=1d0+(1d0-g(site,site,flv))*delta(1,1,flv)
            ratio=ratio*ratiomat(1,1,flv)
          END DO
        ELSE
          ratio=1d0
          DO flv=1,nflv
            DO a=1,ndim; sitea=nb_field(site,a,ifield)
              DO b=1,ndim; siteb=nb_field(site,b,ifield)
                gtmp(a,b)=g(sitea,siteb,flv)
              END DO
            END DO
            ratiomat(1:ndim,1:ndim,flv)=delta(1:ndim,1:ndim,flv)-matmul(gtmp(1:ndim,1:ndim),delta(1:ndim,1:ndim,flv))
            DO a=1,ndim
              ratiomat(a,a,flv)=ratiomat(a,a,flv)+1d0
            END DO
            ratio=ratio*det(ndim,ratiomat(1:ndim,1:ndim,flv))
          END DO
        END IF

        ! calculate total accept probability externally
        CALL acceptprob_local(ratio,newfield,site,time,ifield)

        ! correction of Methopolis ratio
        IF(abs(ratio)<1d0)THEN
          paccept=abs(ratio)/(1+newMetro*abs(ratio))
        ELSE
          paccept=abs(ratio)/(newMetro+abs(ratio))
        END IF
        
        ! record the total number of tryings
        Ntotal_field(ifield)=Ntotal_field(ifield)+1

        ! accept the new configuration with probability paccept
        IF(drand()>paccept)CYCLE

        ! record the number of accepted tryings
        Naccept_field(ifield)=Naccept_field(ifield)+1
        
        ! obtain the phase of the current configuration
        currentphase=currentphase*ratio/abs(ratio)     
        
        ! update Green's function based on Dyson equation
        IF(ndim==1)THEN
          DO flv=1,nflv
            gtmp=g(:,:,flv)
            DO j=1,nsite
              delta_ja=0d0
              IF(j==site)delta_ja=1d0
              gtmp(j,:)=gtmp(j,:)+(g(j,site,flv)-delta_ja)*delta(1,1,flv)/ratiomat(1,1,flv)*g(site,:,flv)
            END DO
            g(:,:,flv)=gtmp
          END DO
        ELSE
          DO flv=1,nflv
            CALL inverse(ndim,ratiomat(1:ndim,1:ndim,flv))
            delta(1:ndim,1:ndim,flv)=matmul(delta(1:ndim,1:ndim,flv),ratiomat(1:ndim,1:ndim,flv))
            gtmp=g(:,:,flv)
            DO j=1,nsite
              DO a=1,ndim; sitea=nb_field(site,a,ifield)
                delta_ja=0d0; IF(j==sitea) delta_ja=1d0
                DO b=1,ndim; siteb=nb_field(site,b,ifield)
                  gtmp(j,:)=gtmp(j,:)+(g(j,sitea,flv)-delta_ja)*delta(a,b,flv)*g(siteb,:,flv)
                END DO
              END DO
            END DO
            g(:,:,flv)=gtmp
          END DO
        END IF

        ! update field
        field(site,time,ifield)=newfield

      END DO ! site-loop

      ! udpate Green's function by exp(phi*F)*G*exp(-phi*F)
      DO flv=1,nflv
        CALL evolve_left_V(time,ifield,g(:,:,flv),nsite,flv,.false.)
        CALL evolve_right_V(time,ifield,g(:,:,flv),nsite,flv,.true.)
      END DO

    END DO ! ifield-loop

  END DO ! time-loop

END SUBROUTINE