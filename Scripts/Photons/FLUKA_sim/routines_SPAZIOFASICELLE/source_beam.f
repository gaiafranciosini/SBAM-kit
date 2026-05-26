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
* *
* Copyright (C) 1990-2013      by   Alfredo Ferrari & Paola Sala  *
* All Rights Reserved.                                             *
* *
* *
* New source for FLUKA9x-FLUKA20xy:                                *
* *
* Created on 07 January 1990   by   Alfredo Ferrari & Paola Sala  *
* Infn - Milan       *
* *
* Last change on  19-May-13    by   Alfredo Ferrari               *
* *
* This is just an example of a possible user written source routine.  *
* note that the beam card still has some meaning - in the scoring the *
* maximum momentum used in deciding the binning is taken from the     *
* beam momentum.  Other beam card parameters are obsolete.            *
* *
* Output variables:                                              *
* *
* Nomore = if > 0 the run will be terminated              *
* *
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
      INCLUDE '(UNRTSF)'
      INCLUDE '(CASLIM)'

*
      LOGICAL LFIRST
*
      SAVE LFIRST
      DATA LFIRST / .TRUE. /
*======================================================================*
* *
* HIT eff VERSION MODIFIED                             *
* *
*======================================================================*
      character*30 FILE
      double precision x, y, z, ekin
      double precision vx, vy, vz, v_mod
      double precision cosx, cosy, cosz
      integer id_part
      character(100) line
            
      NOMORE = 0
* +-------------------------------------------------------------------*
* |  First call initializations:
      IF ( LFIRST ) THEN
* |  *** The following 3 cards are mandatory ***
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.
* |  *** User initialization ***
         INDF = LNNBLN ( SDUSOU )
         IF ( INDF .EQ. 0 )THEN
            FILE = 'beam.txt'
         ELSE
            IF ( INDEX (SDUSOU,'.') .EQ. 0) THEN
               FILE=SDUSOU(1:INDF)//'.txt'
            ELSE
               FILE = SDUSOU
            ENDIF
         ENDIF         
         CALL OAUXFI ( FILE, 22, 'OLD', IERR )
         IF ( IERR .NE. 0 ) then
            CALL FLABRT ( 'SOURCE', 'IMPOSSIBLE TO OPEN FILE' )
         endif
         Write(*,*)'source init: reading first line of input file:'
         READ (22,*, END=1000, ERR=1000) line
      END IF
* |
* +-------------------------------------------------------------------*
* Push one source particle to the stack. 
      NPFLKA = NPFLKA + 1
* Wt is the weight of the particle
      WTFLK  (NPFLKA) = ONEONE
      WEIPRI = WEIPRI + WTFLK (NPFLKA)

*-----------------------------------------------------------------------
* LETTURA DEI DATI DAL FILE DI INPUT
*-----------------------------------------------------------------------
      ekin = 0.D0
      do while(ekin.eq.0.D0)
         READ (22,*, END=1001, ERR=1000) id_part, x, y, z,
     &                                   vx, vy, vz, ekin
      end do
      goto 1002
      
 1000 continue
      write(*,*)'source: error in event number = ', ncase
      CALL FLABRT ( 'SOURCE', 'Error reading input file' )
      
 1001 CONTINUE
      nomore = 1
      RETURN
      
 1002 Continue

* Conversion MeV -> GeV
      ekin = ekin/1000.D0

*-----------------------------------------------------------------------
* ASSEGNAZIONE ID PARTICELLA
*-----------------------------------------------------------------------
      IF ( IJBEAM .NE. -2 ) THEN
         IONID = id_part
         ILOFLK (NPFLKA) = id_part
         LRADDC (NPFLKA) = .FALSE.
         IGROUP (NPFLKA) = 0
      ELSE
* | Heavy ion logica di FLUKA standard
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         ILOFLK (NPFLKA) = IJHION
         LRADDC (NPFLKA) = .FALSE.
         IGROUP (NPFLKA) = 0
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
      END IF

*-----------------------------------------------------------------------
* INIZIALIZZAZIONE VARIABILI DI STACK (Eseguita una sola volta qui!)
*-----------------------------------------------------------------------
      LOFLK  (NPFLKA) = 1
      LOUSE  (NPFLKA) = 0
      KCHFLK (NPFLKA) = 0
      ECRFLK (NPFLKA) = ZERZER
      INFSTK (NPFLKA) = 0
      LNFSTK (NPFLKA) = 0
      ANFSTK (NPFLKA) = ZERZER
      IPRSTK (NPFLKA) = 0
      EKPSTK (NPFLKA) = ZERZER
      DO 100 ISPR = 1, MKBMX1
         SPAREK (ISPR,NPFLKA) = ZERZER
 100  CONTINUE
      DO 200 ISPR = 1, MKBMX2
         ISPARK (ISPR,NPFLKA) = 0
 200  CONTINUE
      ISPARK (MKBMX2,NPFLKA) = NPFLKA
      NPARMA = NPARMA + 1
      NUMPAR (NPFLKA) = NPARMA
      NEVENT (NPFLKA) = 0
      DFNEAR (NPFLKA) = +ZERZER

*-----------------------------------------------------------------------
* CALCOLO COSENI DIRETTORI (Normalizzazione vettore)
*-----------------------------------------------------------------------
      v_mod = sqrt(vx**2 + vy**2 + vz**2)
      if (v_mod .le. 0.D0) then
         write(*,*) 'source: Errore modulo vettore direzione nullo!'
         CALL FLABRT ( 'SOURCE', 'Vettore direzione nullo nel file' )
      endif

      cosx = vx / v_mod
      cosy = vy / v_mod
      cosz = vz / v_mod

* Assegnazione dell'età e dell'energia cinetica
      AGESTK (NPFLKA) = +ZERZER
      AKNSHR (NPFLKA) = -TWOTWO
      TKEFLK (NPFLKA) = ekin

* Calcolo dinamico dell'impulso usando la massa automatica AM(IONID)
      PMOFLK (NPFLKA) = SQRT ( TKEFLK (NPFLKA) * ( TKEFLK (NPFLKA)
     &                       + TWOTWO * AM (IONID) ) )

* Coseni direttori nel vettore FLUKA
      TXFLK  (NPFLKA) = COSX
      TYFLK  (NPFLKA) = COSY
      TZFLK  (NPFLKA) = COSZ

* Polarizzazione cosines:
      TXPOL  (NPFLKA) = -TWOTWO
      TYPOL  (NPFLKA) = +ZERZER
      TZPOL  (NPFLKA) = +ZERZER

*-----------------------------------------------------------------------
* ASSEGNAZIONE COORDINATE SPAZIALI
*-----------------------------------------------------------------------
      if(whasou(2).eq.1.D0) then
         XFLK   (NPFLKA) = x
         YFLK   (NPFLKA) = y
         ZFLK   (NPFLKA) = z
      else
         XFLK   (NPFLKA) = -90.D0
         YFLK   (NPFLKA) = -90.D0
         ZFLK   (NPFLKA) = -90.D0
      endif

* Calcolo dell'energia cinetica totale dei primari
      IF ( ILOFLK (NPFLKA) .EQ. -2 .OR. ILOFLK (NPFLKA) .GT. 100000 )
     &   THEN
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
         TKESUM = TKESUM + ( TKEFLK (NPFLKA) + AMDISC (ILOFLK(NPFLKA)) )
     &            * WTFLK (NPFLKA)
      ELSE
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      END IF
      RADDLY (NPFLKA) = ZERZER
      
* Robustezza geometrica ai confini delle regioni
      CALL GEODRR ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
      CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &              NRGFLK(NPFLKA), IDISC )
      CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
      NLATTC (NPFLKA) = MLATTC
      CMPATH (NPFLKA) = ZERZER
      CALL SOEVSV
      RETURN

      END
