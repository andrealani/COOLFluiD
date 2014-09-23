C-----------------------------------------------------------------------
C     This subroutone is used when the NASA Ames vibrational collisional 
C     model for N2-N2 system is being used. It interpolates by means of 
C     bi-cubic spline approximation the following quantities:
C     - Rotational partition functions of vibrational levels Qrot_v ; v = 0,..,v_max 
C     - Constant volume and constant pressure rotational specific heats
C     - Rotational energies of vibrational levels 
      SUBROUTINE N2LOOKUP ()

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'

      INTEGER V, J, N  
      DOUBLE PRECISION YP1, YPN, ERV, QRV, CPV, T  

C     Maximum and minimum limit of temperature
      TMIN = 100.D0
      TMAX = 60000.D0
      STEP = 10.D0

C     Number of components in the temperature vector
      N = INT((TMAX - TMIN)/STEP) + 1

      TEMP(1) = TMIN
      DO J = 2,N
         TEMP(J) = TEMP(J-1) + STEP 
      ENDDO

C     Data inizialization
      YP1 = 0.D0   ;   YPN = 0.D0
      YQRV(:,:) = 0.D0   ;  Y2QRV(:,:) = 0.D0
      YERV(:,:) = 0.D0   ;  Y2ERV(:,:) = 0.D0
      YCPV(:,:) = 0.D0   ;  Y2CPV(:,:) = 0.D0

C     Spline interpolation
      WRITE(*,*)
      WRITE(*,*)'Starting interpolation loop'
      WRITE(*,*)

      DO V = 0,VMAX-1

         J = 1
         T = TEMP(J)
         DO WHILE (J.LE.N) 

            ! Thermodynamic properties
            CALL N2AMESPROP (V, T, QRV, ERV, CPV)

            YQRV(V+1,J) = QRV
            YERV(V+1,J) = ERV
            YCPV(V+1,J) = CPV  
             
            J = J + 1
            T = TEMP(J-1) + STEP 
 
         ENDDO
 
         ! Interpolation algorithm
         CALL SPLINEN2(N, TEMP, YQRV(V+1,1:N), YP1, YPN, Y2QRV(V+1,1:N))
         CALL SPLINEN2(N, TEMP, YERV(V+1,1:N), YP1, YPN, Y2ERV(V+1,1:N))
         CALL SPLINEN2(N, TEMP, YCPV(V+1,1:N), YP1, YPN, Y2CPV(V+1,1:N))

      ENDDO
      WRITE(*,*)
      WRITE(*,*)'Finished interpolation loop'
      WRITE(*,*)

      END SUBROUTINE N2LOOKUP
C-----------------------------------------------------------------------
C     Subroutine for interpolation: here coeffiecients for spline
C     interpolation are computed and stored
      SUBROUTINE SPLINEN2 (N, X, Y, YP1, YPN, Y2)

      IMPLICIT NONE

      INTEGER I, J, N
      DOUBLE PRECISION X(1:N), Y(1:N), YP1, YPN, Y2(1:N), U(1:N), SIG, 
     &                 P, QN, UN

C     Data inizialization
      U(:) = 0.D0

      IF (YP1.GT.99D30) THEN 
         Y2(1) = 0.D0
         U(1) = 0.D0
      ELSE
         Y2(1) = -0.5D0
         U(1) = (3.D0/(X(2)-X(1)))*((Y(2) - Y(1))/(X(2)-X(1))-YP1)
      ENDIF

C     Computation of second derivatives: a small band linear system is
C     solved     
      DO I =2,N-1         
         SIG = (X(I) - X(I-1))/(X(I+1)-X(I-1))
         P = SIG*Y2(I-1)+2.D0
         Y2(I) = (SIG-1.D0)/P
         U(I) = (6.D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/
     &          (X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      ENDDO 

      IF (YPN.GT.99D30) THEN
         QN = 0.D0
         UN = 0.D0
      ELSE
         QN = 0.5D0
         UN = (3.D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
              ENDIF

      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.D0)

      DO J = N-1,1,-1
         Y2(J) = Y2(J)*Y2(J+1)+U(J)
      ENDDO

      END SUBROUTINE SPLINEN2
C-----------------------------------------------------------------------
