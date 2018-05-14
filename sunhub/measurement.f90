SUBROUTINE measurement(time)
  USE dqmc_complex
  USE model
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: time
  COMPLEX(8) g3(nsite,nsite)
  INTEGER i,j,a,b,a2,b2
  COMPLEX(8) kinetic,double,Sq,deltaij
  CALL put_pool(currentphase)
  g3=-matmul(inv_expk_half,matmul(g,expk_half))
  DO i=1,nsite
    g3(i,i)=g3(i,i)+1d0
  END DO
  kinetic=0d0
  double=0d0
  DO i=1,nsite
    DO j=1,nsite
      kinetic=kinetic+kmat(i,j)*g3(j,i)*currentphase/(La*Lb)*sun
    END DO
    double=double+g3(i,i)**2*currentphase/(La*Lb)
  END DO
  CALL put_pool(kinetic)
  CALL put_pool(double)
  Sq=0d0
  DO a=1,La; DO b=1,Lb; i=label(a,b,1)
    DO a2=1,La; DO b2=1,Lb; j=label(a2,b2,1)
      deltaij=0d0
      IF(i==j)deltaij=1d0
      Sq=Sq+2*g3(i,j)*(deltaij-g3(j,i))*(-1)**(a-a2+b-b2)*currentphase/(La*Lb)
    END DO; END DO
  END DO; END DO
  CALL put_pool(Sq)
  
  g3=g3*currentphase
  CALL put_pool(2,g3(1:2,1))
  CALL zput_array(nsite*nsite,g3(:,:))
END SUBROUTINE
