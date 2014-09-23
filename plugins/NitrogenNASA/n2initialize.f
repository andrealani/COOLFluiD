C     This subroutine   initialize common arrays for computations 
C     based on NASA Ames database data for N2-N system
      SUBROUTINE N2INITIALIZE (NBS, NBT) 

      IMPLICIT NONE
      
      INCLUDE 'commondata.cmn'
  
      INTEGER NBS, NBT, DUM1, DUM2, DUM3, DUMOLD, I, IN1, IN2, IN3, J, 
     &        JOLD, LCHAR, KK, V, VP, DV, LPATH1, LPATH2
      PARAMETER (IN1 = 10, IN2 = 11, IN3 = 12)
      DOUBLE PRECISION MMN, MMN2, CONV, EDISGR, ER(1:LEV), EVV(1:LEV), 
     &                 EV(1:VAMES)      
      CHARACTER(100) PATH1, PATH2 

C     Number of species 
      NS = NBS

C     Number of temperatures
      NT = NBT + 1

C     Flag for vibrationally specific collisional model
      IF (NT.EQ.1) THEN
         IF(NS.EQ.62) THEN 
            FLAGVIBCR = 'NASA_Ames'
         ELSEIF (NS.EQ.69) THEN
            FLAGVIBCR = 'Capitelli'
         ENDIF 
      ENDIF

C     Physical contants
      UKB  = 1.380658D-23   
      UNA  = 6.0221367D23   
      URG  = UKB*UNA
      UE   = 1.602191D-19
      UH   = 6.626075D-34
      UPI  = 3.14159265D0      

C     Tolerance on species mass fractions
      YTOL = 1.D-20

C     Molar masses of N and N2 [kg/mol]
      MMN = 14.0067D-3
      MMN2 = 28.0134D-3  

C     Molar masses stored in vector MM
      MM(:) = 0.D0
      MM(1) = MMN
      DO I = 2,NS
         MM(I) = MMN2
      ENDDO

C     Conversion factor from eV to J
      CONV = UE

C     Data inizialization
      EVIBR(:) = 0.D0  ;  EVIBRUM(:) = 0.D0

C     Common data allocation. This depends on the physical model
C     being used:
C     - NASA Ames database vibrational collisional model (1 temp.)
C     - Capitelli vibrational collisional model          (1 temp.)
C     - Park multi-temperature model                     (2 temp.)

C     Vibrational collisional models
      IF (NT.EQ.1) THEN 

         ! Path to the library data         
         PATH1 = 'data'
         LPATH1 = LCHAR(PATH1)

         SELECT CASE(FLAGVIBCR)

            ! NASA Ames database
            CASE('NASA_Ames')

               ! Nuclear spin of N
               GN = 12

               ! Number of vibrational levels       
               VMAX = VAMES
     
               ! Maximum rotational quantum for each vibrational state 
               WRITE(*,*)'NitrogenNASALibrary:: Reading levels from 
     & file'
               JMAX(:) = 0.D0
               OPEN(UNIT=IN1,FILE=''//PATH1(1:LPATH1)//'/Lev.dat',
     &              STATUS='OLD')
               DUM1 = 0 ; DUM2 = 0 ; DUM3 = 0 ; V = 0
               DUMOLD = DUM2 
        
               ! Loop over energy levels. The nuclear spin is stored
               DO I = 1,LEV
                  READ(IN1,*)DUM1,DUM2,DUM3
                  IF (MOD(DUM3,2).EQ.0) THEN  
                     G(I) = 6
                  ELSE
                     G(I) = 3
                  ENDIF 
                  IF(DUM2.EQ.DUMOLD) THEN
                     IF (V.EQ.0) THEN  
                         JMAX(1) = DUM3
                     ELSE
                         JMAX(DUM2+1) = DUM3
                     ENDIF
                  ELSE 
                     DUMOLD = DUM2
                     V = V + 1
                  ENDIF         
               ENDDO
   
               ! Dissociation energy of the ground state (v=0, J=0) 
               EDISGR = -9.75377D0
             
               ! Energy values converted in Joule and including the
               ! ground state dissociation energy
               REWIND(IN1)
        
               DO I = 1,LEV
                  READ(IN1,*)DUM1,DUM2,DUM3,ENE(I)
                  ENE(I) = (ENE(I) - EDISGR)*CONV
               ENDDO  
              
               CLOSE(IN1)
               
               ! Energy splitting in vibrational and rotational contributions
               ! E(i) = E(v,J) = E(v,J=0) + deltaE(v,J)
         
               ! Computation of E(v,J=0) for each vibrational state
               J = 0 ; EV(:) = 0.D0 
               DO I = 1,VMAX
                  EV(I) = ENE(J+1)
                  EVIBR(I+1) = EV(I)
                  EVIBRUM(I+1) = EVIBR(I+1)*UNA/MM(2)
                  J = J + (JMAX(I) + 1)
               ENDDO

               ! Computation of rotational contributions deltaE(v,J)
               JOLD = 0
               DO I = 1,VMAX
                  J = JMAX(I) + 1
                  DO KK = 1,J
                     ER(KK+JOLD) = ENE(KK+JOLD) - EV(I)
                     EVV(KK+JOLD) = EV(I)        
                  ENDDO 
                  JOLD = JOLD + J
               ENDDO
        
               ! Common vector to store
               ! - Total ro-vibrational energy E(v,J)
               ! - Vibrational energy          E(v,J=0) 
               ! - Rotational correction       deltaE(v,J)
               DO I = 1,LEV
                  E(I,1) = ENE(I)
                  E(I,2) = EVV(I)
                  E(I,3) = ER(I)         
               ENDDO
      
               ! Matrix for energy level computation given quantum numbers v and J
               IND(:,:) = 0
               DO V = 0,VMAX-1
                  DO J = 0,JMAX(V+1)
                     CALL ENLEVEL (V, J, I)
                     IND(V+1,J+1) = I
                  ENDDO
               ENDDO
         
               ! Formation enthalpy of nitrogen atom N in [J] and [J/kg]
               EDN2 = -CONV*EDISGR
               HFN = -CONV*EDISGR/2.D0
               HFORN = -CONV*EDISGR/2.D0*UNA/MMN 
        
               ! Fit coefficients for vibrational dissociation and excitation
         
               ! Dissociation N2(v) + N = N + N + N
               OPEN(UNIT=IN2,FILE=''//PATH1(1:LPATH1)//'/VibDis.dat',
     &              STATUS='OLD')
     
               DO V = 0,VMAX-1 
                  READ(IN2,*)DUM1,AD(V+1),BD(V+1),CD(V+1)
               ENDDO
               CLOSE(IN2)

               ! Excitation N2(v) + N = N2(v') + N
               OPEN(UNIT=IN3,FILE=''//PATH1(1:LPATH1)//'/VibExc.dat',
     &              STATUS='OLD')
     
               DO I = 1,1830
                  READ(IN3,*)V,VP,AE(V+1,VP+1),BE(V+1,VP+1),CE(V+1,VP+1)
               ENDDO
     
               CLOSE(IN3)

            ! Capitelli model
            CASE('Capitelli')

               ! Data inizialization
               D1(:) = 0.D0 ; D2(:) = 0.D0 ; D3(:) = 0.D0 ; D4(:) = 0.D0

               A1(:,:) = 0.D0 ; A2(:,:) = 0.D0 ; A3(:,:) = 0.D0 
               A4(:,:) = 0.D0 ; A5(:,:) = 0.D0

               ! Nuclear spin of N
               GN = 4

               ! Number of vibrational levels
               VMAX = VCAP

               ! Characteristic rotational temperature of N2 molecule 
               THETAR = 2.88D0
             
               ! Dissociation energy of the ground state (v=0, J=0) 
               EDISGR = 9.7538D0  

               OPEN(UNIT=IN1,FILE=''//PATH1(1:LPATH1)//
     &              '/Lev_Capitelli.dat',STATUS='OLD')
               
               DO V = 0,VMAX-1
                  READ(IN1,*)DUM1,EVIBR(V+2)
                  EVIBR(V+2) = EVIBR(V+2)*CONV
                  EVIBRUM(V+2) = EVIBR(V+2)*UNA/MM(2)  
               ENDDO              

               CLOSE(IN1)

               ! Formation enthalpy of nitrogen atom N in [J] and [J/kg]
               EDN2 = CONV*EDISGR
               HFN = CONV*EDISGR/2.D0
               HFORN = CONV*EDISGR/2.D0*UNA/MMN 

               ! Dissociation N2(v) + N = N + N + N
               OPEN(UNIT=IN2,FILE=''//PATH1(1:LPATH1)//
     &              '/VibDis_Capitelli.dat',STATUS='OLD')
     
               DO V = 0,VMAX-1 
                  READ(IN2,*)DUM1,D1(V+1),D2(V+1),D3(V+1),D4(V+1), 
     &                       D5(V+1)
               ENDDO
               CLOSE(IN2)


               ! De-excitation N2(v) + N = N2(v-dv) + N
               OPEN(UNIT=IN3,FILE=''//PATH1(1:LPATH1)//
     &              '/VibExc_Capitelli.dat',STATUS='OLD')
     
               DO V = 0,VMAX-2
                  DO VP = V+1,VMAX-1
                     DV = VP-V
                     READ(IN3,*)DUM1,DUM2,A1(VP+1,DV),A2(VP+1,DV),
     &                          A3(VP+1,DV),A4(VP+1,DV),A5(VP+1,DV)
                  ENDDO
               ENDDO
    
               CLOSE(IN3)

         END SELECT 

C     Park multi-temperature model   
      ELSE

         ! Nuclear spin of N
         GN = 4

         ! Characteristic vibrational and rotational temperatures of N2 molecule 
         THETAR = 2.88D0
         THETAV = 3392.7D0

         ! Formation enthalpy and dissociation energy of N2 molecule     
         EDN2 = 113272.D0*UKB
         HFORN = (113272.D0*UKB/2.D0)*UNA/MM(1) 

         !  Dissociation reaction N2 + N = N + N + N
         AP(1) = 0.049816D0
         AP(2) = -1.60D0
         AP(3) = 113200.D0 

         ! VT interactions (parameters for Millikan and White formula)

         ! N2-N  interaction
         MWA(1) = 180.D0
         MWB(1) = 0.0262D0

         ! N2-N2 interaction
         MWA(2) = 221.D0 
         MWB(2) = 0.0290D0

      ENDIF

      END SUBROUTINE N2INITIALIZE
C-----------------------------------------------------------------------
C     This subroutine computes the global enrgy index i given the
C     vibrational and rotational quantum numbers v and J 
      SUBROUTINE ENLEVEL (V, J, I)

      IMPLICIT NONE
      INCLUDE 'commondata.cmn'
      INTEGER K, V, J, POS, I

      POS = 0
      DO K = 1,V
         POS = POS + (JMAX(K) + 1)
      ENDDO

C     Global energy index 
      I = POS + J + 1

      END SUBROUTINE ENLEVEL
C-----------------------------------------------------------------------
C     This function determines the effective length (non blank) 
C     of a character 
      FUNCTION LCHAR (SYMBOL)

      IMPLICIT NONE

      CHARACTER(*) SYMBOL
      INTEGER LCHAR, C, I

      C = 0 
      DO I = 1, LEN(SYMBOL)
         IF ( SYMBOL(I:I) /= " " ) THEN
            C = C +1
         ENDIF
      ENDDO
      
      LCHAR = C 

      END FUNCTION LCHAR
C-----------------------------------------------------------------------
