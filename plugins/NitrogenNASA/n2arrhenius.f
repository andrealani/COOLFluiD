C-----------------------------------------------------------------------
C     This subroutine computes the mass production terms for
C     vibrationally specific collisional model derived from
C     NASA Ames database.
C     The output is stored in the vector OMEGA in unit of
C     mol/m^3/s 
      SUBROUTINE OMEGAN2 (RHO, YI, T, TV, OMEGA, OMEGAVIB)

      IMPLICIT NONE
 
      INCLUDE 'commondata.cmn'

      INTEGER I, II, IS, V, VP, DV
      DOUBLE PRECISION RHO, T, TV, TTV, LNTTV, KFE, EV, EVP, QTN, QTN2,
     &                 FIX, KFD, KFDEQ, KEQ, LNT, EXP1, EXP2, OMEGAVT, 
     &                 OMEGACV, SIGMA, MU, MILLIKAN, PARK, TAU, SUM1, 
     &                 SUM2, P, ET, ETV, T2, T3, T4
      DOUBLE PRECISION YI(1:NS), RHOIT(1:NS), MOLEI(1:NS), OMEGA(1:NS),
     &                 QRV(1:VMAX), TAUVT(1:NS), OMEGAVIB(1:1)

C     Logarithm of the temperature
      LNT  = DLOG(T)
      TTV = DSQRT(T*TV)
      LNTTV = DLOG(TTV) 

C     Other useful data
      T2 = T**2.D0
      T3 = T**3.D0
      T4 = T**4.D0

C     Conversion in concentrations [mol/m^3]
      MOLEI(:) = 0.D0

      DO IS = 1,NS
         IF (YI(IS).LT.YTOL) THEN 
             RHOIT(IS) = YTOL*RHO
         ELSE
             RHOIT(IS) = YI(IS)*RHO
         ENDIF
         MOLEI(IS) = RHOIT(IS)/MM(IS)
      ENDDO

C     Data inizialization
      OMEGA(:) = 0.D0 ; QRV(:) = 0.D0 

C     Translational partition function of N atoms and vibrational 
C     levels of N2
      CALL PARTTRANS (T, MM(1), 1.D0, QTN)
      CALL PARTTRANS (T, MM(2), 1.D0, QTN2)

      SELECT CASE(NT)

        ! Vibrational collisional model
        CASE(1)

          ! Kinetic model selection
          SELECT CASE(FLAGVIBCR)

            
            CASE('NASA_Ames')

               ! Dissociation reactions
               ! N2(v) + N = N + N + N
               DO V = 0,VMAX-1

                  KFD = 1.D-6*UNA*AD(V+1)*DEXP(BD(V+1)*LNT - CD(V+1)/T)

                  CALL KEQDISS (V, T, QTN, QTN2, UKB, UNA, QRV(V+1), 
     &                          KEQ)

                  EXP1 = MOLEI(1)*(KFD*MOLEI(V+2) - KFD*MOLEI(1)**2/KEQ)
                  OMEGA(1) = OMEGA(1) + 2.D0*EXP1
                  OMEGA(V+2) = - EXP1

               ENDDO

               ! Excitation reactions
               ! N2(v) + N = N2(v') + N ; v' > v 
               DO V = 0,VMAX-2

                  I = IND(V+1,1)
                  EV = E(I,2) 

                  DO VP = V+1,VMAX-1
 
                     II = IND(VP+1,1)
                     EVP = E(II,2)

                     KFE = 1.D-6*UNA*AE(V+1,VP+1)*DEXP(BE(V+1,VP+1)
     &                     *LNT - CE(V+1,VP+1)/T)        
                     CALL KEQEXC (T, EV, EVP, QRV(V+1), QRV(VP+1), UKB, 
     &                            KEQ) 

                     EXP1 = KFE*MOLEI(1)*(MOLEI(V+2) - MOLEI(VP+2)/KEQ) 
                     OMEGA(V+2) = OMEGA(V+2) - EXP1 
                     OMEGA(VP+2) = OMEGA(VP+2) + EXP1   

                  ENDDO

               ENDDO

            CASE('Capitelli')

               ! Dissociation reactions
               ! N2(v) + N = N + N + N
               DO V = 0,VMAX-1

                  KFD = 1.D-6*UNA*(DEXP(D1(V+1) + D2(V+1)/T + D3(V+1)/T2
     &                  + D4(V+1)/T3 + D5(V+1)*LNT))

                  CALL KEQDISS (V, T, QTN, QTN2, UKB, UNA, QRV(V+1), 
     &                          KEQ)

                  EXP1 = MOLEI(1)*(KFD*MOLEI(V+2) - KFD*MOLEI(1)**2/KEQ)
                  OMEGA(1) = OMEGA(1) + 2.D0*EXP1
                  OMEGA(V+2) = - EXP1

               ENDDO

               ! De-excitation reactions
               ! N2(v') + N = N2(v) + N ; v' > v
               DO V = 0,VMAX-2
                 
                  DO VP = V+1,VMAX-1

                     DV = VP - V

                     ! Maximum quantum jump is 50 (see the original
                     ! paper of Capitelli et al.)
                     IF (DV.LE.50) THEN

                         IF (VP.LT.(VMAX-1)) THEN

                             IF (DV.LE.30) THEN

                                 KFE = 1.D-6*UNA*DEXP(A1(VP+1,DV) + 
     &                                 A2(VP+1,DV)/T + A3(VP+1,DV)/T2 +
     &                                 A4(VP+1,DV)/T3 + A5(VP+1,DV)*LNT)

                             ELSE

                                 KFE = 1.D-6*UNA*DEXP(A1(VP+1,DV) + 
     &                                 A2(VP+1,DV)/T + A3(VP+1,DV)/T4 +
     &                                 A4(VP+1,DV)/LNT)

                             ENDIF

                         ! De-excitatation from the last vibrational level
                         ELSE

                             KFE = 1.D-6*UNA*DEXP(A1(VP+1,DV) + 
     &                             A2(VP+1,DV)/T + A3(VP+1,DV)/T2 + 
     &                             A4(VP+1,DV)/T3 + A5(VP+1,DV)*LNT) 

                         ENDIF

                         EVP = EVIBR(VP+2)
                         EV = EVIBR(V+2) 

                         CALL KEQEXC (T, EVP, EV, 1.D0, 1.D0, UKB, KEQ) 

                         EXP1 = KFE*MOLEI(1)*(MOLEI(VP+2) - MOLEI(V+2)/
     &                          KEQ) 
                         OMEGA(VP+2) = OMEGA(VP+2) - EXP1 
                         OMEGA(V+2) = OMEGA(V+2) + EXP1

                     ENDIF

                  ENDDO
                  
               ENDDO

          END SELECT 

        ! Park multi-temperature model
        CASE(2)

          KFD = 1.D-6*UNA*AP(1)*DEXP(AP(2)*LNTTV - AP(3)/TTV)
          KFDEQ = 1.D-6*UNA*AP(1)*DEXP(AP(2)*LNT - AP(3)/T)

          CALL KEQDISSMULTI (NT, GN, T, QTN, QTN2, UKB, UNA, EDN2, KEQ)

          EXP1 = MOLEI(1)*(KFD*MOLEI(2) - MOLEI(1)**2.D0*KFDEQ/KEQ)
          OMEGA(1) = 2.D0*EXP1
          OMEGA(2) = -EXP1

      END SELECT 

C     Mass production terms [mol/m^3/s]
      DO IS = 1,NS
         OMEGA(IS) = OMEGA(IS)*MM(IS)
      ENDDO

C     Vibrational energy conservation equation source terms
C     (This must be used just for the Park multi-temperature model)
      IF (NT.EQ.2) THEN 

         OMEGAVT = 0.D0 ; OMEGACV = 0.D0

         ET = THETAV*URG/MM(2)/(DEXP(THETAV/T) - 1.D0)
         ETV = THETAV*URG/MM(2)/(DEXP(THETAV/TV) - 1.D0)

         ! Chemistry-vibration coupling term
         OMEGACV = OMEGA(2)*ETV

         ! Limiting cross section for Park's correction [m^2]
         IF (T.GT.20000.D0) THEN
            SIGMA = 1.D-21*(50000.D0/20000.D0)**2.D0
         ELSE
            SIGMA = 1.D-21 *(50000.D0/T)**2.D0
         ENDIF

         CALL PRESSURE (YI, RHO, T, P)

         ! Relaxation times for VT interactions
         DO IS = 1, NS
            MU = MM(IS)*MM(2)/(MM(IS) + MM(2))
            MILLIKAN  = DEXP(MWA(IS) *(T**(-0.3333333333) -MWB(IS))
     &                  -18.42D0)*101325.D0/P
            PARK      = DSQRT(UPI*MU*UKB*T/(8.D0*UNA))/(SIGMA*P)
            TAUVT(IS) = MILLIKAN + PARK
         ENDDO

         ! Frequency average of relaxation times
         ! Beware, for consistency with the Ames data base, N2-N2 excitation is neglected
         SUM1 = 0.D0 ; SUM2 = 0.D0 
c         DO IS = 1,NS
         DO IS = 1,1  
            SUM1 = SUM1 + RHOIT(IS)/MM(IS)
            SUM2 = SUM2 + RHOIT(IS)/(MM(IS)*TAUVT(IS))
         ENDDO
         TAU = SUM1 /SUM2

         OMEGAVT   =  RHOIT(2)*(ET - ETV)/TAU
         OMEGAVIB(1) = OMEGAVT + OMEGACV

      ENDIF

      END SUBROUTINE OMEGAN2
C-----------------------------------------------------------------------
C     This subroutine computes the equilibrium constant for the reaction
C     N2 + N = N + N + N when using the multitemperature approach 
      SUBROUTINE KEQDISSMULTI (NT, GN, T, QTN, QTN2, KB, NA, EDN2, KEQ)

      IMPLICIT NONE

      INTEGER GN, NT 
      DOUBLE PRECISION EDN2, EINT, NA, QTN, QTN2, QINT, KB, T, KEQ

C     Internal partion funnction of N2 molecule
      CALL INTENE (T, T, QINT, EINT) 
     
C     Equilibrium constant 
      KEQ = 1.D0/NA*(GN*QTN)**2.D0/(QTN2*QINT)*DEXP(-EDN2/KB/T)

      END SUBROUTINE KEQDISSMULTI
C-----------------------------------------------------------------------
C     This subroutine computes the equilibrium constant for the reacion
C     N2(v) + N = N + N + N at the temperature T for the vibrationally 
C     specific collisional model derived from NASA Ames database.
      SUBROUTINE KEQDISS (V, T, QTN, QTN2, KB, NA, QR, KEQ)

      IMPLICIT NONE

      INCLUDE 'commondata.cmn'
      INTEGER V
      DOUBLE PRECISION EV, T, KEQ, KB, NA, QTN2, QTN, QR, EPS

C     Fix to avoid numerical problems
      EPS = 1.D-25

C     Rotational partition function and vibrational energy
      EV = EVIBR(V+2)
      CALL QROT (V, T, QR)

C     Equilibrium constant based on micro-reversibility
      KEQ = 1.D0/NA*(GN*QTN)**2.D0/(QTN2*QR)
     &      *DEXP(-(2.D0*HFN - EV)/(KB*T)) + EPS  
    
      END SUBROUTINE KEQDISS
C----------------------------------------------------------------------- 
C     This subroutine computes the equilibrium constants for the reaction
C     N2(v) + N = N2(v') + N at the temperature T for the vibrationally 
C     specific collisional model derived from NASA Ames database.
      SUBROUTINE KEQEXC (T, EV, EVP, QV, QVP, KB, KEQ)

      IMPLICIT NONE
      DOUBLE PRECISION EV, EVP, KB, QV, QVP, T, KEQ, EPS 

C     Fix to avoid numerical problems
      EPS = 1.D-25      

C     Equilibrium constant based on micro-reversibility
      KEQ = QVP/QV*EXP(-(EVP - EV)/(KB*T)) + EPS

      END SUBROUTINE KEQEXC
C-----------------------------------------------------------------------
