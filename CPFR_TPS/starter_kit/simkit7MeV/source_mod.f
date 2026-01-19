*$ CREATE SOURCE.FOR
*COPY SOURCE
*
*=== source ===========================================================*
*
      SUBROUTINE SOURCE ( NOMORE )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2020      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     New source for FLUKA9x-FLUKA20xy:                                *
*                                                                      *
*     Created on 07 January 1990   by    Alfredo Ferrari & Paola Sala  *
*                                                Infn - Milan          *
*                                                                      *
*     Last change on  03-Jan-20    by             Paola Sala           *
*                                                Infn - Milan          *
*                                                                      *
*     Exemple of source sampling from a spectrum                       *
*     the spectrum is read in by the rdmysp routine                    *
*     the file name is read from the SDUM in the source card           *
*     the format of the spectrum file is chosen by the first what      *
*     in the source card                                               *
*     1.0 : emin emax content (default)                                *
*     2.0 : emid content (uniform bin spacing)                         *
*     3.0 : first line=number of bins, then all lower bin limits,      *
*           then all contents                                          *
*     WHASOU (1) .LE. 0 : Spectrum is assumed in kinetic energgy (def.)*
*     WHASOU (1) .GT. 0 : Spectrum is assumed in momentum              *
*     WHASOU (2) .GT. 0 : The absolute normalization is preserved      *
*                         (def.: the cum. spectrum is normalized to 1.)*
*     All other beam parameters from BEAM and BEAMPOS cards            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*       Output variables:                                              *
*                                                                      *
*              Nomore = if > 0 the run will be terminated              *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(BEAMCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(IOIOCM)'
      INCLUDE '(LTCLCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SOURCM)'
      INCLUDE '(SUMCOU)'
*
      LOGICAL LFIRST, LISNUT
*
*     max number of energy bins
      PARAMETER ( MXEFIL = 1000 )
*
      DIMENSION EKNMIN (MXEFIL), EKNMAX (MXEFIL), CUMSPE(MXEFIL),
     &          EMED(MXEFIL), DELTAE (MXEFIL), SPECTR  (MXEFIL)
      LOGICAL LEKINE, LSNORM
      COMMON /MYSPEC/ EKNMIN, CUMSPE, SPECTR, EMED, DELTAE, TOTSPE, 
     &                ESUPSP, NOBINS, LEKINE, LSNORM
      SAVE LFIRST
      DATA LFIRST / .TRUE. /
*  Statement function:
      LISNUT (IJ) = INDEX ( PRNAME (IJ), 'NEUTRI' ) .GT. 0
*======================================================================*
*                                                                      *
*                 BASIC VERSION                                        *
*                                                                      *
*======================================================================*
      NOMORE = 0
*  +-------------------------------------------------------------------*
*  |  First call initializations:
*      WRITE (LUNOUT, *) 'ciaooooooooooooooooooooo'
*      WRITE (*, *) 'ciaooooooooooooooooooooo'
*      CALL FLUSH(LUNOUT)
      IF ( LFIRST ) THEN
*  |  *** The following 3 cards are mandatory ***
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.
*  |  *** User initialization ***
         IF ( WHASOU(1) .GT. ZERZER ) THEN
            LEKINE = .FALSE. 
         ELSE
            LEKINE = .TRUE. 
         END IF
         IF ( WHASOU(2) .GT. ZERZER ) THEN
            LSNORM = .TRUE. 
         ELSE
            LSNORM = .FALSE. 
         END IF
*         WRITE (LUNOUT, *) 'start source init'
*         CALL FLUSH(LUNOUT)
         CALL RDMYSP 
*         WRITE ( LUNOUT, *) ' SOURCE succesfully initialized'
*         WRITE ( LUNOUT, *) ' Cumulative spectrum is: ',TOTSPE
*        CALL FLUSH(LUNOUT)
         IF ( LEKINE ) THEN
            WRITE ( LUNOUT, *)' Max kinetic energy: ', ESUPSP
         ELSE
            WRITE ( LUNOUT, *)' Max momentum: ', ESUPSP
         END IF
         IF ( LSNORM ) THEN
            WRITE ( LUNOUT, *) ' Absolute normalization is preserved: ',
     $           TOTSPE 
         ELSE
            WRITE ( LUNOUT, *) ' Cumulative spectrum is norm. to 1. '
         ENDIF
         DO IG = 1,NOBINS
            EKNMAX(IG) = EKNMIN(IG) + DELTAE(IG)
            EMED(IG) = SQRT(EKNMIN(IG)*EKNMAX(IG))
            WRITE(LUNOUT,*) EKNMIN(IG),EKNMAX(IG),DELTAE(IG),SPECTR(IG)
     $           ,CUMSPE(IG) 
         END DO
      END IF
*  |
*  +-------------------------------------------------------------------*
*  Push one source particle to the stack. Note that you could as well
*  push many but this way we reserve a maximum amount of space in the
*  stack for the secondaries to be generated
*  Npflka is the stack counter: of course any time source is called it
*  must be =0
      NPFLKA = NPFLKA + 1
*  Wt is the weight of the particle
      WTFLK  (NPFLKA) = ONEONE
      WEIPRI = WEIPRI + WTFLK (NPFLKA)
*  Particle type (1=proton.....). Ijbeam is the type set by the BEAM
*  card
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope:
      IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
         IARES  = IPROA
         IZRES  = IPROZ
         IISRES = IPROM
         CALL STISBM ( IARES, IZRES, IISRES )
         IJHION = IPROM  * 100000 + MOD ( IPROZ, 100 ) * 1000 + IPROA
         IJHION = IJHION * 100    + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
         EEXSTK (NPFLKA) = EXENRG (IONID)
         TMNSTK (NPFLKA) = TMNLF  (IONID)
         ILVSTK (NPFLKA) = IEXLVL (IONID)
         LFRPHN (NPFLKA) = .FALSE.
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      ELSE IF ( IJBEAM .EQ. -2 ) THEN
         IJHION = IPROM  * 100000 + MOD ( IPROZ, 100 ) * 1000 + IPROA
         IJHION = IJHION * 100    + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
         EEXSTK (NPFLKA) = EXENRG (IONID)
         TMNSTK (NPFLKA) = TMNLF  (IONID)
         ILVSTK (NPFLKA) = IEXLVL (IONID)
         ILOFLK (NPFLKA) = IJHION
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
*  |  Parent radioactive isotope:
         IRDAZM (NPFLKA) = 0
*  |  Particle age (s)
         AGESTK (NPFLKA) = +ZERZER
*  |  Kinetic energy of the particle (GeV)
         TKEFLK (NPFLKA) = SQRT ( PBEAM**2 + AM (IONID)**2 )
     &                   - AM (IONID)
*  |  Particle momentum
         PMOFLK (NPFLKA) = PBEAM
*        PMOFLK (NPFLKA) = SQRT ( TKEFLK (NPFLKA) * ( TKEFLK (NPFLKA)
*    &                          + TWOTWO * AM (IONID) ) )
         LFRPHN (NPFLKA) = .FALSE.
*  |
*  +-------------------------------------------------------------------*
*  |  Normal hadron:
      ELSE
         IONID = IJBEAM
         EEXSTK (NPFLKA) = EXENRG (IJBEAM)
         TMNSTK (NPFLKA) = TMNLF  (IJBEAM)
         ILVSTK (NPFLKA) = IEXLVL (IJBEAM)
         ILOFLK (NPFLKA) = IJBEAM
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
*  |  Parent radioactive isotope:
         IRDAZM (NPFLKA) = 0
*  |  Particle age (s)
         AGESTK (NPFLKA) = +ZERZER
*  |  Kinetic energy of the particle (GeV)
* sample Kinetic energy of the particle or its  momentum (GeV)
 5       RNDSPE = FLRNDM (RNDSPE)
* find the correct energy bin 
         DO 2000 IG = 1, NOBINS
            IF ( RNDSPE .LE. CUMSPE (IG) ) THEN
               IGOK = IG
               GO TO 2500
            END IF
 2000    CONTINUE
         STOP ' STOP:SPECTRUM-IG ?????'
 2500    CONTINUE
* Sample linearly inside the bin
         RNDSPE = FLRNDM ( RNDSPE )
C      WRITE(*,*) "IGOK = ",IGOK
         IF ( IGOK .GT. 1 ) THEN
            SPEMAX = SPECTR(IGOK)
            SPEMIN = SPECTR(IGOK-1)
            CUMMAX = CUMSPE(IGOK)
            CUMMIN = CUMSPE(IGOK-1)
            EMAXIM = EMED(IGOK)
            EMINIM = EMED(IGOK-1)
         ELSE
            SPEMAX = SPECTR(IGOK+1)
            SPEMIN = SPECTR(IGOK)
            CUMMAX = CUMSPE(IGOK+1)
            CUMMIN = CUMSPE(IGOK)
            EMAXIM = EMED(IGOK+1)
            EMINIM = EMED(IGOK)
         ENDIF
C      write (*,*) " spemax, spemin, cummax, cummin, eminim, emaxim"
C      write (*,*)  spemax, spemin, cummax, cummin, eminim, emaxim
CGB      RNDSPE=(RNDSPE-CUMMIN)/(CUMMAX-CUMMIN)
         RNDSPE=FLRNDM(XDUMMY)
         IF(SPEMAX.EQ.SPEMIN) THEN
c   indeed, inside the bin I the spectrum is flat!
            VSPEC = RNDSPE
         ELSE
            VSPEC = RNDSPE*(SPEMAX-SPEMIN)*(SPEMAX+SPEMIN)
            VSPEC = SQRT(VSPEC+SPEMIN**2)-SPEMIN
            VSPEC = VSPEC/(SPEMAX-SPEMIN)
         ENDIF
         VSPEC = VSPEC*(EMAXIM-EMINIM)+EMINIM
*         WRITE (LUNOUT, *) 'VSPEC', VSPEC
*         CALL FLUSH(LUNOUT)
*         CALL FLNRR2(RGAUS1,RGAUS2)
*         VSPEC = VSPEC+RGAUS1*WHASOU(3)+WHASOU(4)
         
*         IF ((VSPEC .LT. 0) .OR. (VSPEC .GT. WHASOU(5))) THEN
*            GO TO 5
*         ENDIF
            
C      write(*,*) "vspec = ",vspec
         IF ( LEKINE ) THEN
            TKEFLK (NPFLKA) = VSPEC
            PMOFLK (NPFLKA) = SQRT ( TKEFLK (NPFLKA) * ( TKEFLK (NPFLKA)
     &           + TWOTWO * AM (IONID) ) )
         ELSE
            PMOFLK (NPFLKA) = VSPEC
            PBEAM = PMOFLK (NPFLKA)
            TKEFLK (NPFLKA) = SQRT ( PBEAM**2 + AM (IONID)**2 ) -             
     &           AM (IONID)
         ENDIF
*  |  +----------------------------------------------------------------*
*  |  |  Check if it is a neutrino, if so force the interaction
*  |  |  (unless the relevant flag has been disabled)
         IF ( LISNUT (IJBEAM) .AND. LNUFIN ) THEN
            LFRPHN (NPFLKA) = .TRUE.
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  Not a neutrino
         ELSE
            LFRPHN (NPFLKA) = .FALSE.
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
      END IF
*  |
*  +-------------------------------------------------------------------*
*  From this point .....
*  Particle generation (1 for primaries)
      LOFLK  (NPFLKA) = 1
*  User dependent flag:
      LOUSE  (NPFLKA) = 0
*  No channeling:
      KCHFLK (NPFLKA) = 0
      ECRFLK (NPFLKA) = ZERZER
*  Extra infos:
      INFSTK (NPFLKA) = 0
      LNFSTK (NPFLKA) = 0
      ANFSTK (NPFLKA) = ZERZER
*  Parent variables:
      IPRSTK (NPFLKA) = 0
      EKPSTK (NPFLKA) = ZERZER
*  User dependent spare variables:
      DO 100 ISPR = 1, MKBMX1
         SPAREK (ISPR,NPFLKA) = ZERZER
 100  CONTINUE
*  User dependent spare flags:
      DO 200 ISPR = 1, MKBMX2
         ISPARK (ISPR,NPFLKA) = 0
 200  CONTINUE
*  Save the track number of the stack particle:
      ISPARK (MKBMX2,NPFLKA) = NPFLKA
      NPARMA = NPARMA + 1
      NUMPAR (NPFLKA) = NPARMA
      NEVENT (NPFLKA) = 0
      DFNEAR (NPFLKA) = +ZERZER
*  ... to this point: don't change anything
      AKNSHR (NPFLKA) = -TWOTWO
*  Cosines (tx,ty,tz)
      TXFLK  (NPFLKA) = UBEAM
      TYFLK  (NPFLKA) = VBEAM
      TZFLK  (NPFLKA) = WBEAM
*     TZFLK  (NPFLKA) = SQRT ( ONEONE - TXFLK (NPFLKA)**2
*    &                       - TYFLK (NPFLKA)**2 )
*  Polarization cosines:
      TXPOL  (NPFLKA) = -TWOTWO
      TYPOL  (NPFLKA) = +ZERZER
      TZPOL  (NPFLKA) = +ZERZER
*  Particle coordinates
      XFLK   (NPFLKA) = XBEAM
      YFLK   (NPFLKA) = YBEAM
      ZFLK   (NPFLKA) = ZBEAM
*  Calculate the total kinetic energy of the primaries: don't change
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope:
      IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      ELSE IF ( ILOFLK (NPFLKA) .EQ. -2 .OR.
     &          ILOFLK (NPFLKA) .GT. 100000 ) THEN
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
*  |
*  +-------------------------------------------------------------------*
*  |  Standard particle:
      ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
         TKESUM = TKESUM + ( TKEFLK (NPFLKA) + AMDISC (ILOFLK(NPFLKA)) )
     &          * WTFLK (NPFLKA)
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      END IF
*  |
*  +-------------------------------------------------------------------*
      RADDLY (NPFLKA) = ZERZER
*  Here we ask for the region number of the hitting point.
*     NREG (NPFLKA) = ...
*  The following line makes the starting region search much more
*  robust if particles are starting very close to a boundary:
      CALL GEODRR ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
      CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &              NRGFLK(NPFLKA), IDISC )
*  Do not change these cards:
      CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
      NLATTC (NPFLKA) = MLATTC
      CMPATH (NPFLKA) = ZERZER
      CALL SOEVSV
      RETURN
*=== End of subroutine Source =========================================*
      END

*
*=== rdmysp ===========================================================*
*
      SUBROUTINE RDMYSP

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
      INCLUDE '(SOURCM)'
*     max number of energy bins
      PARAMETER ( MXEFIL = 1000 )
*
      DIMENSION EKNMIN (MXEFIL), EKNMAX (MXEFIL), CUMSPE(MXEFIL),
     &          EMED(MXEFIL), DELTAE (MXEFIL), SPECTR  (MXEFIL)
      LOGICAL LEKINE, LSNORM
      COMMON /MYSPEC/ EKNMIN, CUMSPE, SPECTR, EMED, DELTAE, TOTSPE, 
     &                ESUPSP, NOBINS, LEKINE, LSNORM

      CHARACTER *12 FILE

      IF (ABS(WHASOU(1)) .GT. ZERZER ) THEN
         IFORMT = ABS(WHASOU(1))
      ELSE 
         IFORMT = 1
      END IF
      INDF = LNNBLN ( SDUSOU )
      IF ( INDF .EQ. 0 )THEN
         FILE = '7MeV.dat'
         WRITE (LUNOUT, '(A)') FILE
      ELSE
         IF ( INDEX (SDUSOU,'.') .EQ. 0) THEN
            FILE=SDUSOU(1:INDF)//'.dat'
            WRITE (LUNOUT, '(A)') FILE
*            CALL FLUSH(LUNOUT)
         ELSE
            FILE = SDUSOU
         ENDIF
      ENDIF
      CALL OAUXFI ( FILE, 22, 'OLD', IERR )
      IF ( IERR .NE. 0 ) 
     &   CALL FLABRT ( 'SOURCE', 'IMPOSSIBLE TO OPEN FILE' )
*  +-------------------------------------------------------------------*
*  |  EMIN EMAX CONTENT:

*      WRITE (LUNOUT, *) 'find nbins emin emax'
*      CALL FLUSH(LUNOUT)
      IF ( IFORMT .EQ. 1 ) THEN
         DO IE = 1, MXEFIL
            READ (22,*, END=101, ERR=1000) EKNMIN (IE), EKNMAX (IE), 
     &                                     SPECTR (IE)
            DELTAE (IE) = EKNMAX (IE) - EKNMIN (IE)
            Write(*,*) EKNMIN(IE),EKNMAX(IE),DELTAE(IE),SPECTR(IE)
         ENDDO
 101     CONTINUE
         NOBINS = IE - 1

         
*  | 
*  +-------------------------------------------------------------------*
*  |  EMID  CONTENT , FIXED WIDTH:
      ELSE IF ( IFORMT .EQ. 2 ) THEN
         WRITE (LUNERR, '(A)') FILE
         DO IE = 1, MXEFIL
             READ (22,*, END=102, ERR=1000) EKNMAX (IE),  
     &           SPECTR (IE)
             EKNMAX(IE) = EKNMAX(IE)/1000.D+00
*             WRITE(LUNOUT,*) 'primo ciclo', EKNMAX(IE)
*             CALL FLUSH(LUNOUT)
             IF ( IE .EQ. 2 ) THEN
                DD = EKNMAX (IE) - EKNMAX (IE-1)
                DDM = HLFHLF * DD
                EKNMIN (1) = EKNMAX (1) - DDM
                DELTAE (1) = DD
             ENDIF
             EKNMIN (IE) = EKNMAX (IE) - DDM
             DELTAE (IE) = DD
          ENDDO
          DO IE = 1, MXEFIL
             Write(*,*) EKNMIN(IE),EKNMAX(IE),DELTAE(IE),SPECTR(IE)
*             WRITE(LUNOUT,*) 'secondo ciclo', EKNMAX(IE)
*             CALL FLUSH(LUNOUT)
          ENDDO
 102      CONTINUE
          NOBINS = IE -1
*          WRITE (LUNOUT, *) 'find nobinbs', NOBINS
*          CALL FLUSH(LUNOUT)
*  +-------------------------------------------------------------------*
*  |  NOBINS,  NOBINS+1 ENERGY LIMITS, NOBINS  CONTENTS,:
      ELSE IF ( IFORMT .EQ. 3 ) THEN
         READ (22,*, ERR=1000) NOBINS
         DO IE = 1, NOBINS 
            READ (22,*, ERR=1000) EKNMIN (IE)
            IF ( IE .NE. 1) THEN
               DELTAE (IE-1) = EKNMIN (IE) - EKNMIN (IE-1)
               EKNMAX (IE-1) = EKNMIN (IE)
            ENDIF
         ENDDO
         READ (22,*, ERR=1000, END=1000) EKNMAX (NOBINS)
         DELTAE (NOBINS) = EKNMAX (NOBINS) - EKNMIN (NOBINS)
         DO IE = 1, NOBINS 
            READ (22,*, ERR=1000, END=1000) SPECTR (IE)
         ENDDO
         DO IE = 1,NOBINS
            Write(*,*) EKNMIN(IE),EKNMAX(IE),DELTAE(IE),SPECTR(IE)
         ENDDO
      ENDIF
*  | 
*  +-------------------------------------------------------------------*
*  |
      ESUPSP = EKNMAX (NOBINS)
* NOW BUILD THE CUMULATIVE DISTR
      CUMSPE (1) = SPECTR (1)* DELTAE (1)
*     WRITE (LUNOUT, *) 'do cumulative'
*      CALL FLUSH(LUNOUT)
      DO 200 IG = 2, NOBINS
         CUMSPE(IG) = CUMSPE(IG-1) +  SPECTR (IG)* DELTAE (IG)
 200  CONTINUE
*      WRITE (LUNOUT, *) 'done cumulative'
*      CALL FLUSH(LUNOUT)
* NOW NORMALIZE IT
      TOTSPE = CUMSPE (NOBINS)
      IF ( TOTSPE .LT. AZRZRZ ) CALL FLABRT ( 'SOURCE' , ' 0 SPECTRUM') 
*      WRITE (LUNOUT, *) 'cumspe', TOTSPE
*      CALL FLUSH(LUNOUT)
      DO 201 IG = 1, NOBINS
         CUMSPE(IG) = CUMSPE(IG) / TOTSPE 
*        WRITE (LUNOUT, *) 'cumspe at ig',IG, '  ' ,CUMSPE(IG)
*        CALL FLUSH(LUNOUT)
*     write (67,*) emin(ig), emax(ig), SPECTR (IG), cumspe(ig)
 201  CONTINUE

      RETURN
 1000 CONTINUE
      WRITE (LUNOUT, *) 'source : error reading from file '
      WRITE (LUNOUT, '(A)') FILE
      WRITE (LUNERR, *) 'source : error reading from file '
      WRITE (LUNERR, '(A)') FILE
      CALL FLABRT ( 'SOURCE' , ' error reading file ')
      RETURN
*=== End of subroutine rdmysp  =========================================*
      END
