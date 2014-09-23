C-----------------------------------------------------------------------
C     This subroutine computes the rotational partition function, 
C     rotational energy and rotational specific heat for a given
C     vibrational level when eploying the NASA Ames vibrational
C     collisional model
      SUBROUTINE N2AMESPROP (V, T, QR, ER, CVR)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'

      INTEGER I, J, V
      DOUBLE PRECISION SUM1, SUM2, SUM3, T, DELTAEVJ, DEXP1, QR, ER, CVR

C     Data inizialization
      QR = 0.D0    ;  ER = 0.D0    ;  CVR = 0.D0
      SUM1 = 0.D0  ;  SUM2 = 0.D0  ;  SUM3 = 0.D0

C     Thermodynamic property computation
      DO J = 0,JMAX(V+1)
         I = IND(V+1,J+1)
         DELTAEVJ = E(I,3)
         DEXP1 = DEXP(-DELTAEVJ/UKB/T) 
         SUM1 = SUM1 + (2.D0*J + 1.D0)*G(I)*DEXP1
         SUM2 = SUM2 + (2.D0*J + 1.D0)*G(I)*DELTAEVJ*DEXP1
         SUM3 = SUM3 + (2.D0*J + 1.D0)*G(I)*DELTAEVJ**2.D0*DEXP1
      ENDDO

C     Rotational partition function, energy and specific heat for the 
C     vibrational level v
      QR = SUM1
      ER = UNA/MM(2)*SUM2/SUM1
      CVR = UNA/MM(2)/UKB/T**2.D0*(-(SUM2/SUM1)**2.D0 +SUM3/SUM1)       

      END SUBROUTINE N2AMESPROP 
C----------------------------------------------------------------------
C     This subroutine computes the translation partition function
C     for a given species. The volume of the system is fixed to
C     1.D0 for sake of simplicity  
      SUBROUTINE PARTTRANS (T, MASS, VOL, QTR)

      IMPLICIT NONE
 
      INCLUDE 'commondata.cmn'
      DOUBLE PRECISION MASS, T, VOL, QTR       

C     Translational partition function 
      QTR = ((2.D0*UPI*MASS*UKB*T/UH**2.D0/UNA)**1.5D0)*VOL

      END SUBROUTINE PARTTRANS
C-----------------------------------------------------------------------
C     This subroutine computes the rotational partition function for a 
C     given vibrational level v. The computatation is performed under
C     the assumption of rotation in equilibrium at T
      SUBROUTINE QROT (V, T, QR)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER I, J, JPLUS1, V
      DOUBLE PRECISION A, B, H, HH, T, QR

C     Kinetic model selection
      SELECT CASE(FLAGVIBCR)

        CASE('NASA_Ames')

           ! Useful data
           H = STEP
           HH = STEP*STEP

           ! Rotational partition function
           J = INT((T-TMIN)/H)+1
           JPLUS1 = J+1           
           A = (TEMP(JPLUS1)-T)/H
           B = (T-TEMP(J))/H  

           QR = A*YQRV(V+1,J) + B*YQRV(V+1,JPLUS1) + 
     &          ((A*A*A-A)*Y2QRV(V+1,J) + (B*B*B-B)*Y2QRV(V+1,JPLUS1))
     &          *HH/6.D0

         CASE('Capitelli')

           ! Rotational partition function
           QR = (T/THETAR/2.D0)

      END SELECT 

      END SUBROUTINE QROT
C-----------------------------------------------------------------------
C     This subroutine computes the rotaional energy when 
C     rotational levels of a given vibrational quamtum state v are in 
C     equilibrium at T (Boltzmann distribution)
      SUBROUTINE EROT (V, T, ER)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER I, J, JPLUS1, V
      DOUBLE PRECISION A, B, H, HH, T, ER

C     Kinetic model selection
      SELECT CASE(FLAGVIBCR)

        CASE('NASA_Ames')

           ! Useful data
           H = STEP
           HH = STEP*STEP

           ! Rotational energy per unit mass of corresponfing to vibrational level v
           J = INT((T-TMIN)/H)+1
           JPLUS1 = J+1           
           A = (TEMP(JPLUS1)-T)/H
           B = (T-TEMP(J))/H  

           ER = A*YERV(V+1,J) + B*YERV(V+1,JPLUS1) + 
     &          ((A*A*A-A)*Y2ERV(V+1,J) + (B*B*B-B)*Y2ERV(V+1,JPLUS1))*
     &          HH/6.D0

        CASE('Capitelli')

           ! Rotational energy per unit mass of corresponfing to vibrational level v
           ER = URG/MM(2)*T

      END SELECT 


      END SUBROUTINE EROT
C-----------------------------------------------------------------------
C     This subroutine computes the internal partition function and enery
C     per unit mass in equilibrium conditions (Boltzmann distribution)
      SUBROUTINE INTENE (T, TV, QINT, EINT)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER I, J, V 
      DOUBLE PRECISION DELTAEVJ, EVJ0, EINT, DEXP1, DEXP2, QINT, SUM1, 
     &                 SUM2, SUM3, SUM4, T, TV

C     Internal partition function and energy (per unit mass)
      SUM1 = 0.D0  ;  SUM2 = 0.D0  ; SUM3 = 0.D0  ;  SUM4 = 0.D0

      SELECT CASE(NT)

         ! Vibrational collisional model 
         CASE(1)

           ! Kinetic model selection
           SELECT CASE(FLAGVIBCR)

             CASE('NASA_Ames')

                DO V = 0,VMAX-1 

                   I = IND(V+1,1)
                   EVJ0 = E(I,2)
                   DEXP1 = DEXP(-EVJ0/UKB/T)
                   SUM2 = 0.D0 ; SUM3 = 0.D0 

                   DO J = 0,JMAX(V+1)
                      I = IND(V+1,J+1)
                      DELTAEVJ = E(I,3)
                      DEXP2 = DEXP(-DELTAEVJ/UKB/T)
                      SUM2 = SUM2 + (2.D0*J + 1.D0)*G(I)*DEXP2 
                      SUM3 = SUM3 + (2.D0*J + 1.D0)*G(I)*DELTAEVJ*DEXP2
                   ENDDO

                   SUM1 = SUM1 + SUM2*DEXP1
                   SUM4 = SUM4 + DEXP1*(EVJ0*SUM2 + SUM3) 

                ENDDO

                QINT = SUM1
                EINT = UNA/MM(2)/QINT*SUM4

             CASE('Capitelli')

                DO V = 0,VMAX-1
                   EVJ0 = EVIBR(V+2)
                   DEXP1 = DEXP(-EVJ0/UKB/TV) 
                   SUM1 = SUM1 + DEXP1
                   SUM2 = SUM2 + EVJ0*DEXP1
                ENDDO

                QINT = SUM1*(T/THETAR/2.D0)
                EINT = URG/MM(2)*T + UNA/MM(2)*SUM2/SUM1
                
           END SELECT 

         ! Park model
         CASE(2)

           QINT = (T/THETAR/2.D0)*(1.D0/(1.D0 - DEXP(-THETAV/TV)))
           EINT = URG/MM(2)*T + URG/MM(2)*THETAV/(DEXP(THETAV/TV) 
     &            - 1.D0)

      END SELECT 

      END SUBROUTINE INTENE
C-----------------------------------------------------------------------
C     This subroutine computes the specific heat at constant volume 
C     per unit mass of all components of N2-N mixtures according to the 
C     NASA Ames database energy data 
      SUBROUTINE N2CV (T, CV) 

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER I, J, JPLUS1, V 
      DOUBLE PRECISION A, B, H, HH, T, CV(1:NS)

C     Data inizialization
      CV(:) = 0.D0

C     Atomic nitrogen N (Translation only)
      CV(1) = 1.5D0*URG/MM(1)
 
C     N2 molecule 
      SELECT CASE(NT)

        ! Vibrational collisional model 
        CASE(1)

          ! Kinetic model selection
          SELECT CASE(FLAGVIBCR)

            CASE('NASA_Ames')

               ! Vibrational levels of N2 
               H = STEP  ;  HH = STEP*STEP
               J = INT((T-TMIN)/H)+1
               JPLUS1 = J+1           
               A = (TEMP(JPLUS1)-T)/H
               B = (T-TEMP(J))/H  
 
               ! Translation and rotation 
               DO V = 0,VMAX-1

                  CV(V+2) = 1.5D0*URG/MM(2) 

                  CV(V+2) = CV(V+2) + A*YCPV(V+1,J) + B*YCPV(V+1,JPLUS1)
     &                      + ((A*A*A-A)*Y2CPV(V+1,J) + (B*B*B-B)*
     &                      Y2CPV(V+1,JPLUS1))*HH/6.D0

               ENDDO

            CASE('Capitelli')

               ! Translation and rotation 
               DO V = 0,VMAX-1
                  CV(V+2) = 2.5D0*URG/MM(2)
               ENDDO

          END SELECT

        ! Park multi-temperature model
        CASE(2)

          CV(2) = 2.5D0*URG/MM(2)

      END SELECT 

      END SUBROUTINE N2CV
C-----------------------------------------------------------------------
C     This subroutine computes the specific heat at constant pressure 
C     per unit mass of all components of N2-N mixtures according to the 
C     NASA Ames database energy data. 
      SUBROUTINE N2CP (T, CP) 

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER I  
      DOUBLE PRECISION T, CP(1:NS)

C     Data inizialization
      CP(:) = 0.D0

C     Constant volume species specific heats
      CALL N2CV (T, CP)  

C     Additional contribution in order to obtain constant pressure 
C     species specific heats
      DO I = 1,NS
         CP(I) = CP(I) + URG/MM(I)
      ENDDO    

      END SUBROUTINE N2CP
C-----------------------------------------------------------------------
C     This subroutine computes the internal energy per unit mass of 
C     N2-N mixtures according to the NASA Ames database energy data. 
      SUBROUTINE N2ENERGY (T, TV, ET, ER, EV, HF, ES) 

      IMPLICIT NONE
       
      INCLUDE 'commondata.cmn' 
      INTEGER I, V
      DOUBLE PRECISION T, TV, ERR, ES(1:NS), ET(1:NS), ER(1:NS), 
     &                 EV(1:NS), HF(1:NS)

C     Data inizialization
      ES(:) = 0.D0 ; ET(:) = 0.D0 
      ER(:) = 0.D0 ; EV(:) = 0.D0 
      HF(:) = 0.D0

C     N atom (formation enthalpy must be added)
      ET(1) = 1.5D0*URG/MM(1)*T 
      HF(1) = HFORN

C     N2 molecule
      SELECT CASE(NT) 

        ! NASA Ames database
        CASE(1)

          DO V = 0,VMAX-1
             CALL EROT (V, T, ERR)
             ET(V+2) = 1.5D0*URG/MM(2)*T
             ER(V+2) = ERR
             EV(V+2) = EVIBRUM(V+2)
          ENDDO
        
        ! Park multi-temperature model
        CASE(2)
         
         ET(2) = 1.5D0*URG/MM(2)*T
         ER(2) = URG/MM(2)*T
         EV(2) = URG/MM(2)*THETAV/(DEXP(THETAV/TV) - 1.D0) 

      END SELECT 

C     Internal energy per unit mass of mixture species 
      DO I = 1,NS
         ES(I) = ET(I) + ER(I) + EV(I) + HF(I) 
      ENDDO

      END SUBROUTINE N2ENERGY 
C------------------------------------------------------------------------- 
C     This subroutine computes the internal energy per unit mass of 
C     N2-N mixtures according to the NASA Ames database energy data 
      SUBROUTINE N2ENTHALPY (T, TV, HT, HR, HV, HF, HS) 

      IMPLICIT NONE
       
      INCLUDE 'commondata.cmn' 
      INTEGER I, V
      DOUBLE PRECISION T, TV, ER, HS(1:NS), HT(1:NS), HR(1:NS), 
     &                 HV(1:NS), HF(1:NS)

C     Data inizialization
      HS(:) = 0.D0 ; HT(:) = 0.D0 
      HR(:) = 0.D0 ; HV(:) = 0.D0 
      HF(:) = 0.D0

C     N atom (formation enthalpy must be added)
      HT(1) = 2.5D0*URG/MM(1)*T 
      HF(1) = HFORN

C     N2 molecule
      SELECT CASE(NT) 

        ! NASA Ames database
        CASE(1)

          DO V = 0,VMAX-1
             CALL EROT (V, T, ER)
             HT(V+2) = 2.5D0*URG/MM(2)*T
             HR(V+2) = ER
             HV(V+2) = EVIBRUM(V+2)
          ENDDO

        ! Park multi-temperature model
        CASE(2) 

           HT(2) = 2.5D0*URG/MM(2)*T
           HR(2) = URG/MM(2)*T
           HV(2) = URG/MM(2)*THETAV/(DEXP(THETAV/TV) - 1.D0)

      END SELECT 

C     Static enthalpy per unit mass of mixture species 
      DO I = 1,NS
         HS(I) = HT(I) + HR(I) + HV(I) + HF(I) 
      ENDDO

      END SUBROUTINE N2ENTHALPY 
C-----------------------------------------------------------------------
C     This subroutine computes the frozen specific heat ratio for the
C     NASA Ames vibrational collisional model and the Park
C     multi-temperature model for N2-N mistures
      SUBROUTINE FROZENGAMMAN2 (YI, T, GAMMA)

      IMPLICIT NONE
 
      INCLUDE 'commondata.cmn'     
      INTEGER IS
      DOUBLE PRECISION SUM1, SUM2, YI(1:NS), CV(1:NS), T, KI, 
     &                 GAMMA

C     Constant volume specific heat
      CALL N2CV (T, CV) 

C     Data inizialization (only degrees of freedom in equilibrium
C     with translation are taken into account)
      SUM1 = 0.D0 ; SUM2 = 0.D0 ; KI = 0.D0

      DO IS = 1,NS
         SUM1 = SUM1 + YI(IS)*URG/MM(IS)
         SUM2 = SUM2 + YI(IS)*CV(IS)
      ENDDO

      KI = SUM1/SUM2

C     Frozen specific heat
      GAMMA = 1.D0 + KI

      END SUBROUTINE FROZENGAMMAN2
C-----------------------------------------------------------------------
C     This subroutine computes the frozen speed of sound for the
C     NASA Ames vibrational collisional model and the Park 
C     multi-temperature model for N2-N mistures
      SUBROUTINE SOUNDSPEED (RHO, YI, T, C)

      IMPLICIT NONE
  
      INCLUDE 'commondata.cmn'     
      INTEGER IS
      DOUBLE PRECISION SUM1, SUM2, YI(1:NS), RHO, P, CV(1:NS), T,  
     &                 KI, GAMMA, C

C     Constant volume specific heat (only degrees of freedom in equilibrium
C     with translation are taken into account)
      CALL N2CV (T, CV) 

C     Data inizialization
      SUM1 = 0.D0 ; SUM2 = 0.D0 
      KI = 0.D0   ; RHO  = 0.D0 ; P = 0.D0

      DO IS = 1,NS
         P = P + YI(IS)*URG/MM(IS)
         SUM1 = SUM1 + YI(IS)*URG/MM(IS)
         SUM2 = SUM2 + YI(IS)*CV(IS)
      ENDDO
      P = P*RHO*T

      KI = SUM1/SUM2

C     Frozen speed of sound
      C = DSQRT((1.D0 + KI)*P/RHO)

      END SUBROUTINE SOUNDSPEED
C-----------------------------------------------------------------------
C     This subroutine computes the equilibrium composition of N-N2
C     system in equilibrium conditions in terms of molar fractions
      SUBROUTINE EQCOMPN2 (T, P, XTOL, X)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER I, J, V
      DOUBLE PRECISION T, P, PN2, QTN, QTN2, QINT, RHS, EINT, VOL,
     &                 DELTA2, MASS, NBN2, RHO, EVJ, XTOL, EV0, QRV, 
     &                 DEXP0, X(1:NS), NI(1:NS), YI(1:NS), XX(1:2)

C     Partition function evaluation
      VOL = 1.D0
      CALL PARTTRANS (T, MM(1), VOL, QTN)
      CALL PARTTRANS (T, MM(2), VOL, QTN2)
      CALL INTENE (T, T, QINT, EINT)      
 
      RHS = (1.D0/P)*(UKB*T/VOL)*DEXP(-EDN2/UKB/T)*((GN*QTN)**2.D0)/
     &      (QTN2*QINT)

      DELTA2 = RHS**2.D0 + 4.D0*RHS
      
C     Molar fraction of N and N2  
      X(:) = 0.D0
      XX(1) = (- RHS + DSQRT(DELTA2))/2.D0
      X(1) = XX(1)
      XX(2) = 1.D0 - XX(1)

      MASS = XX(1)*MM(1) + XX(2)*MM(2)
      RHO = P/URG/T*MASS
      PN2 = XX(2)*P
      NBN2 = PN2/UKB/T

      NI(:) = 0.D0
      SELECT CASE(NT)

        ! Vibrational collisional model
        CASE(1)

          ! Ro-vibrational levels of N2 molecule 
          ! (Boltzmann distribution at temperature T)  
          DO V = 0,VMAX-1

             CALL QROT (V, T, QRV)
             EV0 = EVIBR(V+2)
             DEXP0 = DEXP(-EV0/UKB/T)
             NI(V+2) = NBN2*QRV/QINT*DEXP0

         ENDDO

         ! Conversion from number densities to molar fractions
         DO V = 0,VMAX-1
            X(V+2) = NI(V+2)*MASS/RHO/UNA
         ENDDO

        ! Park multi-temperature model 
        CASE(2)

          X(2) = XX(2)

      END SELECT 

C     Tolerance on molar fractions in order to avoi numerical problems
      DO I = 1,NS
         IF (X(I).LT.XTOL) THEN
             X(I) = XTOL
         ENDIF
      ENDDO

      CALL MOLFRACTOSPECFRAC (X, YI)
c      WRITE(*,*)
c      DO I = 1,NS
c         WRITE(*,*)YI(I)
c      ENDDO
c      WRITE(*,*)
c      DO I = 1,NS
c         WRITE(*,*)RHO*YI(I)
c      ENDDO
c      PAUSE
      
      END SUBROUTINE EQCOMPN2
C-----------------------------------------------------------------------
C     This subroutine stores (for sake of simplicity) the species molar 
C     masses for the N-N2 system in the vector masses
      SUBROUTINE MOLARMASSES (MASSES)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER IS
      DOUBLE PRECISION MASSES(1:NS)
     
      DO IS = 1,NS
         MASSES(IS) = MM(IS) 
      ENDDO

      END SUBROUTINE MOLARMASSES
C-----------------------------------------------------------------------
C     This subroutine assigns to the variable RGAS the value of the 
C     universal gas constant
      SUBROUTINE UNGASCONST (RGAS)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      DOUBLE PRECISION RGAS

      RGAS = URG

      END SUBROUTINE UNGASCONST
C-----------------------------------------------------------------------
C     This subroutine converts mass fractions ys to molar fraction xs
      SUBROUTINE SPECMASSTOMOLFRAC (YI, XI)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER IS
      DOUBLE PRECISION MASS, YI(1:NS), XI(1:NS)
  
C     Molar mass
      MASS = 0.D0
      DO IS = 1,NS
         MASS = MASS + YI(IS)/MM(IS)
      ENDDO    
      MASS = 1.D0/MASS

C     Conversion from mass fractions to molar fractions
      DO IS = 1,NS 
         XI(IS) = YI(IS)*MASS/MM(IS)
      ENDDO

      END SUBROUTINE SPECMASSTOMOLFRAC
C-----------------------------------------------------------------------
C     This subroutine converts molar fractions xs to mass fractions ys
      SUBROUTINE MOLFRACTOSPECFRAC (XI, YI)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER IS
      DOUBLE PRECISION MASS, YI(1:NS), XI(1:NS)

C     Molar mass
      MASS = 0.D0
      DO IS = 1,NS
         MASS = MASS + XI(IS)*MM(IS)
      ENDDO 

C     Conversion from molar fractions to mass fractions
      DO IS = 1,NS
         YI(IS) = XI(IS)*MM(IS)/MASS
      ENDDO

      END SUBROUTINE MOLFRACTOSPECFRAC
C-----------------------------------------------------------------------
C     This subroutine computes the global number density given the compo
C     sition expressed in terms of molar fractions  
      SUBROUTINE NUMBERD (P, T, ND)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      DOUBLE PRECISION P, T, ND

C     Mixture number density
      ND = 0.D0
      ND = P/UKB/T

      END SUBROUTINE NUMBERD 
C-----------------------------------------------------------------------
C     This subroutine computes the mixture density given temperature, 
C     pressure and number density 
      SUBROUTINE DENSITY (XI, ND, RHO)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER IS
      DOUBLE PRECISION XI(1:NS), ND, RHO

      RHO = 0.D0
      DO IS = 1, NS
        RHO = RHO + XI(IS)*MM(IS)/UNA
      ENDDO
      
      RHO = RHO*ND

      END SUBROUTINE DENSITY
C-----------------------------------------------------------------------
C     This subroutine computes the pressure according to Dalton's law 
C     valid for a mixture of perfect gases
      SUBROUTINE PRESSURE (YI, RHO, T, P)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER IS
      DOUBLE PRECISION YI(1:NS), RHO, T, P

C     Mixture pressure
      P = 0.D0
      DO IS = 1,NS 
         P = P + YI(IS)/MM(IS)      
      ENDDO
      P = P*RHO*URG*T

      END SUBROUTINE PRESSURE
C-----------------------------------------------------------------------
C     This subroutine computes the mixture molar mass given in input 
C     density and number density
      SUBROUTINE MOLARMASS (RHO, ND, MASS)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      DOUBLE PRECISION ND, RHO, MASS

      MASS = RHO*UNA/ND

      END SUBROUTINE MOLARMASS
C-----------------------------------------------------------------------
