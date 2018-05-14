SUBROUTINE postprocess
  USE dqmc_complex
  USE model
  IMPLICIT NONE
  COMPLEX(8) mean_phase,mean_phase_err,g3(nsite,nsite),g3_err(nsite,nsite)
  COMPLEX(8) kinetic,kinetic_err
  IF(id==0)THEN
  CALL get_pool(mean_phase,mean_phase_err)
  PRINT*,'sign average:',mean_phase,mean_phase_err
  CALL get_pool(kinetic,kinetic_err)
  PRINT*,'kinetic energy:',kinetic,kinetic_err
  CALL get_pool(kinetic,kinetic_err)
  PRINT*,'double occupation:',kinetic,kinetic_err
  CALL get_pool(kinetic,kinetic_err)
  PRINT*,'structure factor:',kinetic,kinetic_err
  CALL get_pool(2,g3(1:2,1),g3_err(1:2,1))
  PRINT*,'g3(1:2,1):',g3(1:2,1),g3_err(1:2,1)
  CALL zget_array_err(nsite*nsite,g3,g3_err)
  PRINT*,'g3:',g3(1:2,1),g3_err(1:2,1)
  PRINT*,'Job is done with',nd,' cores'
  END IF
END SUBROUTINE
