SUBROUTINE measurement(time)
  USE dqmc_complex
  USE model
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: time
  COMPLEX(8) g3(nsite,nsite)
  INTEGER i,j,a,a2,b,b2
  COMPLEX(8) kinetic,double,Sqz,Sqx,deltaij,Sqc,Pph
  
  CALL put_pool(currentphase)
  
  g3=-matmul(inv_expk_half,matmul(g,expk_half))
  DO i=1,nsite
    g3(i,i)=g3(i,i)+1d0
  END DO
  
  kinetic=0d0
  double=0d0
  DO i=1,nsite
    DO j=1,nsite
      kinetic=kinetic+kmat(i,j)*g3(j,i)*currentphase/nsite*2
    END DO
    IF(mod(i,2)==1) double=double+g3(i,i)*g3(i+1,i+1)*currentphase/nsite*2
  END DO
  CALL put_pool(kinetic)
  CALL put_pool(double)
  
  Sqz=0d0
  Sqx=0d0
  Sqc=0d0
  DO a=1,La; DO b=1,Lb; i=label(a,b,2)
    DO a2=1,La; DO b2=1,Lb; j=label(a2,b2,2)
      deltaij=0d0
      IF(i==j)deltaij=1d0
      Sqx=Sqx+2*g3(i,j)*(deltaij-g3(j-1,i-1))*(-1)**(a-a2+b-b2)*currentphase/(La*Lb)
      Sqz=Sqz+((g3(i,i)-g3(i-1,i-1))*(g3(j,j)-g3(j-1,j-1))+g3(i,j)*(deltaij-g3(j,i))+ &
              & g3(i-1,j-1)*(deltaij-g3(j-1,i-1)) )*(-1)**(a-a2+b-b2)*currentphase/(La*Lb)
      Sqc=Sqc+((g3(i,i)+g3(i-1,i-1))*(g3(j,j)+g3(j-1,j-1))+g3(i,j)*(deltaij-g3(j,i))+ &
              & g3(i-1,j-1)*(deltaij-g3(j-1,i-1)) )*(-1)**(a-a2+b-b2)*currentphase/(La*Lb)
    END DO; END DO
  END DO; END DO
  CALL put_pool(Sqz)
  CALL put_pool(Sqx)
  CALL put_pool(Sqc)

  Pph=gph_x2(1)/beta*sum(phi(:,:,1)**2)*currentphase/(La*Lb)
  CALL put_pool(Pph)

  g3=g3*currentphase
  CALL put_pool(g3(1,1))
  CALL put_pool(g3(3,1))
END SUBROUTINE
