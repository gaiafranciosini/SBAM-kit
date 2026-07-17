*$ CREATE MGDRAW.FOR
*COPY MGDRAW
*                                                                      *
*=== mgdraw ===========================================================*
*                                                                      *
      SUBROUTINE MGDRAW ( ICODE, MREG )
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2006      by        Alfredo Ferrari           *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     MaGnetic field trajectory DRAWing: actually this entry manages   *
*                                        all trajectory dumping for    *
*                                        drawing                       *
*                                                                      *
*     Created on   01 march 1990   by        Alfredo Ferrari           *
*                                              INFN - Milan            *
*     Last change  05-may-06       by        Alfredo Ferrari           *
*                                              INFN - Milan            *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(CASLIM)'
      INCLUDE '(COMPUT)'
      INCLUDE '(SOURCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(GENSTK)'
      INCLUDE '(MGDDCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(QUEMGD)'
      INCLUDE '(SUMCOU)'
      INCLUDE '(TRACKR)'
      INCLUDE '(LTCLCM)'
c
      include 'mgdraw.inc'
*
      DIMENSION DTQUEN ( MXTRCK, MAXQMG )
      integer  ICPART, IBPART 
      character*8 NEWREGNAM,MREGNAM
      CHARACTER*8 MLATNAM
      logical ldead
      double precision erawfiber, equenchedfiber
      double precision erawsci, equenchsci
      RETURN
*----------------------------------------------------------------------*
*                                                                      *
*     Icode = 1: call from Kaskad                                      *
*     Icode = 2: call from Emfsco                                      *
*     Icode = 3: call from Kasneu                                      *
*     Icode = 4: call from Kashea                                      *
*     Icode = 5: call from Kasoph                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
c      write(*,*)"sono in mgdraw"
*  +-------------------------------------------------------------------*
      if(idbflg.gt.1) then
         call GEOR2N ( mreg, MREGNAM, IERR ) 
         write(*,*)' '
         write(*,*)'------------- mgdraw: Ev =',ncase,' -------------'
         write(*,*)'jtrack = ',jtrack,' regione = ',MREGNAM
         write(*,*)' '
      endif

c
c************* PARTICLES AND HEAVY IONS WITH JTRACK=-2
      IF((JTRACK.GE.1).OR.((JTRACK.LE.-2).AND.(JTRACK.GE.-6)))THEN
         ICPART = ICHRGE(JTRACK)
         IBPART = IBARCH(JTRACK)
         AMPART = AM(JTRACK)
c*************HEAVY IONS WITH JTRACK.LT.-6
      ELSEIF(JTRACK.LT.-6)THEN
         ICPART = ICHEAV(-JTRACK)
         IBPART = IBHEAV(-JTRACK)
         AMPART = AMNHEA(-JTRACK)
      ELSE
         ICPART=int(zfrttk)
         IBPART = IBARCH(JTRACK)
C         IBPART=0
         AMPART = AM(JTRACK)
      ENDIF

c**************************************************************
      CALL GEOR2N(MREG,MREGNAM,IERR)
      CALL GEON2R(MREGNAM,MREG,IERR)
*      CALL GEOL2N(MLATTC,MLATNAM,IRTLAT,IERR)
c      write(*,*) jtrack, mreg,' ',MREGNAM, MLATTC
c
      if(idbflg.gt.1) then
         write(*,*)'ptrack = ',ptrack,' ekin = ',etrack-ampart
         write(*,*)'zfrttk = ',zfrttk,' icpart = ',icpart,'  m =',ampart
         WRITE (*,*)'x,y,z = ', ( SNGL (XTRACK (I)), SNGL (YTRACK (I)),
     &        SNGL (ZTRACK (I)), I = 0, NTRACK )
         WRITE (*,*)'cx,cy,cz = ',SNGL(CXTRCK),SNGL(CYTRCK),SNGL(CZTRCK)
         WRITE (*,*)'px,py,pz = ',SNGL(ptrack*CXTRCK),
     &        SNGL(ptrack*CYTRCK),SNGL(ptrack*CZTRCK)
         write(*,*)'time = ',atrack,' step = ',ctrack,' path = ',cmtrck
      endif
c
c Find which particle of the part common we are following:
c let's assume it's the last one...
c

      call UpdateCurrPart(mreg,icpart,ibpart,ampart,icode,
     &     xtrack(0),ytrack(0),ztrack(0),JRUN)
c     
c      write(*,*)'Esco da Mgdraw'
      RETURN
*     
*================================================================*
*                                                                      *
*     Boundary-(X)crossing DRAWing:                                    *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             19: boundary crossing                                    *
*     Icode = 2x: call from Emfsco                                     *
*             29: boundary crossing                                    *
*     Icode = 3x: call from Kasneu                                     *
*             39: boundary crossing                                    *
*     Icode = 4x: call from Kashea                                     *
*             49: boundary crossing                                    *
*     Icode = 5x: call from Kasoph                                     *
*             59: boundary crossing                                    *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY BXDRAW ( ICODE, MREG, NEWREG, XSCO, YSCO, ZSCO )
c     scoring at the boundaries... select the region
c     
c         write(*,*)'Sono in Bxdraw'
c************* PARTICLES AND HEAVY IONS WITH JTRACK=-2
       IF((JTRACK.GE.1).OR.((JTRACK.LE.-2).AND.(JTRACK.GE.-6)))THEN
         ICPART = ICHRGE(JTRACK)
         IBPART = IBARCH(JTRACK)
         AMPART = AM(JTRACK)
c************* HEAVY IONS WITH JTRACK.LT.-6
       ELSEIF(JTRACK.LT.-6)THEN
         ICPART = ICHEAV(-JTRACK)
         IBPART = IBHEAV(-JTRACK)
         AMPART = AMNHEA(-JTRACK)
       ELSE
         ICPART= ICHRGE(JTRACK)
C         IBPART=0
         IBPART = IBARCH(JTRACK)
         AMPART = AM(JTRACK)
       ENDIF
c
c debug printing
c
      if(idbflg.gt.0) then
         call GEOR2N ( newreg, NEWREGNAM, IERR ) 
         call GEOR2N ( mreg, MREGNAM, IERR ) 
*         CALL GEOL2N ( MLATTC,MLATNAM,IRTLAT,IERR )
         write(*,*)' '
         write(*,*)'------------- bxdraw: Ev =',ncase,' -------------'
         write(*,*)'jtrack = ',jtrack,' icode = ',icode
         write(*,*)'da -> ',mregnam,'a -> = ',newregnam, ' lat_idx = ',
     &        MLATTC

         if(idbflg.gt.0) then
            write(*,*)'zfrttk = ',zfrttk,' icpart = ',icpart,
     &           ' m =',AMPART
            write(*,*)'ptrack = ',ptrack,' ekin = ',etrack-AMPART
            write(*,*)'x,y,z = ',XSCO,YSCO,ZSCO
            WRITE (*,*)'cx,cy,cz = ',SNGL(CXTRCK), SNGL(CYTRCK)
     &           ,SNGL(CZTRCK)
            write(*,*)'number of track segments: ntrack = ',NTRACK
            do ii =1,ntrack
               write(*,*)'track seg # ',ii,' track seg lenght = '
     &              ,ttrack(ii)
            end do
            write(*,*)'number of en deposition: mtrack = ',MTRACK
            do ii =1,mtrack
               write(*,*)'track dep # ',ii,' track dep val = ',
     &              dtrack(ii)
            end do
         endif
         write(*,*)' '
      endif
c
c scoring the crossing info.
c
      if((newreg.eq.nregtriumf).and.(mreg.eq.nregvuoto)) then
         call score_cross(
     &        ICPART,IBPART,AMPART,newreg,mreg,xsco,ysco,zsco)
      endif   
c      write(*,*)'Esco da Bxdraw'
c     
      RETURN
*     
* ======================================================================*
*                                                                      *
*     Event End DRAWing:                                               *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY EEDRAW ( ICODE )
      RETURN
      if(idbflg.gt.2) then
         call dump_common()
      endif
      RETURN
*
*======================================================================*
*                                                                      *
*     ENergy deposition DRAWing:                                       *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             10: elastic interaction recoil                           *
*             11: inelastic interaction recoil                         *
*             12: stopping particle                                    *
*             13: pseudo-neutron deposition                            *
*             14: escape                                               *
*             15: time kill                                            *
*     Icode = 2x: call from Emfsco                                     *
*             20: local energy deposition (i.e. photoelectric)         *
*             21: below threshold, iarg=1                              *
*             22: below threshold, iarg=2                              *
*             23: escape                                               *
*             24: time kill                                            *
*     Icode = 3x: call from Kasneu                                     *
*             30: target recoil                                        *
*             31: below threshold                                      *
*             32: escape                                               *
*             33: time kill                                            *
*     Icode = 4x: call from Kashea                                     *
*             40: escape                                               *
*             41: time kill                                            *
*             42: delta ray stack overflow                             *
*     Icode = 5x: call from Kasoph                                     *
*             50: optical photon absorption                            *
*             51: escape                                               *
*             52: time kill                                            *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY ENDRAW ( ICODE, MREG, RULL, XSCO, YSCO, ZSCO )
*  +-------------------------------------------------------------------*
      RETURN
c debug
c
      if(idbflg.gt.1) then
         write(*,*)'------------- endraw: Ev =',ncase,' -------------'
         write(*,*)'jtrack = ',jtrack,' mreg = ',mreg
         write(*,*)'rull = ',rull,' icode = ',icode,' idcurr = ',idcurr
         write(*,*)'idead(idcurr) = ',idead(idcurr)
         write(*,*)'x,y,x = ',sngl(XSCO), sngl(YSCO), sngl(ZSCO) 
         WRITE (*,*)'cx,cy,cz = ',SNGL(CXTRCK),SNGL(CYTRCK),SNGL(CZTRCK)
      endif
c
Cgb
Cgb   to preserve idead flag correctly
Cgb
      idead (idcurr) = icode
Cgb
      if(idbflg.gt.1) then
         write(*,*)'ENDRAW: Now idead(idcurr) = ',idead(idcurr)
      endif
Cgb
      if(jtrack.lt.62) then
c************* PARTICLES AND HEAVY IONS WITH JTRACK=-2
         IF((JTRACK.GE.1).OR.((JTRACK.LE.-2).AND.(JTRACK.GE.-6)))THEN
            ICPART = ICHRGE(JTRACK)
            IBPART = IBARCH(JTRACK)
            AMPART = AM(JTRACK)
c*************HEAVY IONS WITH JTRACK.LT.-6
         ELSEIF(JTRACK.LT.-6)THEN
            ICPART = ICHEAV(-JTRACK)
            IBPART = IBHEAV(-JTRACK)
            AMPART = AMNHEA(-JTRACK)
         ELSE
            ICPART= ICHRGE(JTRACK)
            IBPART = IBARCH(JTRACK)
C            IBPART=0
            AMPART = AM(JTRACK)
         ENDIF
c
c if not delta ray below threshold than update the current particle
c
         if( (JTRACK.ne.3).or.(ICODE.ne.22)) then
c            write(*,*)'***********************************************'
c            write(*,*)' NB UpdateCurrentParticle Called from ENDRAW!!!'
c            write(*,*)'***********************************************'
            call UpdateCurrPart(mreg,icpart,ibpart,ampart,icode,
     &           xsco,ysco,zsco,JEND)
            if (jtrack.eq.7) then
               px(idcurr) = pxf(idcurr)
               py(idcurr) = pyf(idcurr)
               pz(idcurr) = pzf(idcurr)
            endif
         endif
      endif 
c     
c     end debug
c     
      RETURN
*
*======================================================================*
*                                                                      *
*     SOurce particle DRAWing:                                         *
*                                                                      *
*======================================================================*
*
      ENTRY SODRAW
      RETURN
      if(idbflg.gt.1) then
         write(*,*)'------------- sodraw: Ev =',ncase,' -------------'
      endif
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope: it works only for 1 source particle on
*  |  the stack for the time being
      IF ( ILOFLK (NPFLKA) .GE. 100000 .AND. LRADDC (NPFLKA) ) THEN
         IARES  = MOD ( ILOFLK (NPFLKA), 100000  )  / 100
         IZRES  = MOD ( ILOFLK (NPFLKA), 10000000 ) / 100000
         IISRES = ILOFLK (NPFLKA) / 10000000
         IONID  = ILOFLK (NPFLKA)
c
*  |
*  +-------------------------------------------------------------------*
*  |  Patch for heavy ions: it works only for 1 source particle on
*  |  the stack for the time being
      ELSE IF ( ABS (ILOFLK (NPFLKA)) .GE. 10000 ) THEN
         IONID = ILOFLK (NPFLKA)
         CALL DCDION ( IONID )
      ELSE 
         IONID = ILOFLK (NPFLKA)
      END IF
c
c   first call to initialing the kinematic tracking 

c 
      call UpdateCurrPart(mreg,icpart,ibpart,ampart,icode,
     &     0.D+00,0.D+00,0.D+00,JSTART)
c
      if(idbflg.gt.0) then
         write(*,*)'reg = ',numreg(nump),' tprod = ',tempo(nump),
     &        ' jpa = ',jpa(nump)
         write(*,*)'vert = ',vxi(nump),vyi(nump),vzi(nump),
     &        ' p = ',px(nump),py(nump),pz(nump)
      endif
*  |
*  +-------------------------------------------------------------------*
      RETURN
*
*======================================================================*
*                                                                      *
*     USer dependent DRAWing:                                          *
*                                                                      *
*     Icode = 10x: call from Kaskad                                    *
*             100: elastic   interaction secondaries                   *
*             101: inelastic interaction secondaries                   *
*             102: particle decay  secondaries                         *
*             103: delta ray  generation secondaries                   *
*             104: pair production secondaries                         *
*             105: bremsstrahlung  secondaries                         *
*             110: decay products                                      *
*     Icode = 20x: call from Emfsco                                    *
*             208: bremsstrahlung secondaries                          *
*             210: Moller secondaries                                  *
*             212: Bhabha secondaries                                  *
*             214: in-flight annihilation secondaries                  *
*             215: annihilation at rest   secondaries                  *
*             217: pair production        secondaries                  *
*             219: Compton scattering     secondaries                  *
*             221: photoelectric          secondaries                  *
*             225: Rayleigh scattering    secondaries                  *
*     Icode = 30x: call from Kasneu                                    *
*             300: interaction secondaries                             *
*     Icode = 40x: call from Kashea                                    *
*             400: delta ray  generation secondaries                   *
*  For all interactions secondaries are put on GENSTK common (kp=1,np) *
*  but for KASHEA delta ray generation where only the secondary elec-  *
*  tron is present and stacked on FLKSTK common for kp=npflka          *
*                                                                      *
*======================================================================*


*
      ENTRY USDRAW ( ICODE, MREG, XSCO, YSCO, ZSCO )
      RETURN
Cgb     To mark correctly interaction code of current particle
Cgb
      intcode(idcurr) = icode
*     
c debug printing
c
      if(idbflg.gt.0) then
         call GEOR2N ( mreg, MREGNAM, IERR ) 
         write(*,*)' '
         write(*,*)'------------- usdraw: Ev =',ncase,' -------------'
         write(*,*)' Fluka index= ',ISPUSR(MKBMX2),' idcurr= ',idcurr,
     &        'jtrack = ',jtrack,' icode = ',icode
         write(*,*)'regione = ',mregnam, 'x,y,z = ',XSCO,YSCO,ZSCO
         if(idbflg.gt.2) then
            write(*,*)'zfrttk = ',zfrttk,' icpart = ',icpart,',m = ',
     &           AMPART
            write(*,*)'ptrack = ',ptrack,' ekin = ',etrack-AMPART,
     &           ' t = ',atrack
            WRITE (*,*)'cx,cy,cz = ',SNGL(CXTRCK), SNGL(CYTRCK)
     &           ,SNGL(CZTRCK)
         endif
         write(*,*)'  USDRAW secondary writing : NP=',NP
         do  ip = 1, NP
            WRITE(*, *) ip,' jpa sec= ',KPART(ip),' Ek sec= ',
     &           TKI(ip),' p sec= ',plr(ip)
            if(kpart(ip).lt.-1) then
               CALL USRDCI(IJ,IONA,IONZ,IONM)
               write(*,*)'A= ',iona,' Z= ',ionz
            endif
         end do
         write(*,*)'  USDRAW heavy secondary writing : NPHEAV=',NPHEAV
         do  ip = 1, NPHEAV
            WRITE(*, *) ip,' jpa sec= ',KHEAVY(ip),' Ek sec= ',
     &           tkheav(ip),' p sec= ',pheavy(ip)
         end do
         write(*,*)' '
      endif
*
*   Do nothing for elastic scattering
c      if ( icode .eq.  100) return
*  +-------------------------------------------------------------------*
*  |  Delta rays, Compton and Moller-Bhabba scattering:
      IF ( ICODE.EQ.103 .OR. ICODE.EQ.210 .OR. ICODE.EQ. 212 ) THEN
         ldead = .false.
*  |
*  +-------------------------------------------------------------------*
*  |  Bremsstrahlung:
      ELSE IF ( ICODE .EQ. 208 .OR. ICODE .EQ. 105 ) THEN
         NPA    = NP
         ldead = .false.
*  +-------------------------------------------------------------------*
*  |  Compton scattering:
      ELSE IF ( ICODE .EQ. 219 ) THEN
         idead (idcurr) = icode
         ldead = .true.
*  +-------------------------------------------------------------------*
*  |  PhotoElectric:
      ELSE IF ( ICODE .EQ. 221 ) THEN
         idead (idcurr) = icode
         ldead = .true.
*  |
*  +-------------------------------------------------------------------*
*  |  In flight or at rest annihilation:
      ELSE IF ( ICODE .EQ. 214 .OR. ICODE .EQ. 215 ) THEN
         idead (idcurr) = icode
         ldead = .true.
*  |
*  +-------------------------------------------------------------------*
*  |  Pair production:
      ELSE IF ( ICODE.EQ.217 .OR. ICODE.EQ.104 ) THEN
         idead (idcurr) = icode
         ldead = .true.
*  |
*  +-------------------------------------------------------------------*
*  |  Low energy neutron secondaries:
      ELSE IF ( ICODE.EQ.300 ) THEN
         idead (idcurr) = icode
         ldead =.true.
*  |
*  +-------------------------------------------------------------------*
*  |  Decays
      ELSE IF ( ICODE.EQ.102 .OR. icode.eq.110) THEN
         idead (idcurr) = icode
         ldead = .true.
*  |
*  +-------------------------------------------------------------------*
*  |  Inelastic interactions:
      ELSE IF ( ICODE.EQ.101 ) THEN
         idead (idcurr) = icode
         ldead = .true.
*  |
*  +-------------------------------------------------------------------*
*  |  Everything else:
      ELSE
         ldead = .false.
      END IF
*  |
*     +-------------------------------------------------------------------*
      if(numint.lt.maxint) then
         numint = numint +1
         xint(numint) = xsco
         yint(numint) = ysco
         zint(numint) = zsco
         intpa(numint) = idcurr
C     idead(numint) = idead(idcurr)
         if(idbflg.gt.1) then
            write(*,*)'USDRAW: Now numint= ', numint
         endif
         intcode(numint) = intcode(idcurr)
      else
         write(*,*)'WARNIG ::: USDRAW: ev',ncase,' numint exceeding= ',
     &        numint
      endif
Cgb
Cgb      idead (numint) = intcode(idcurr)
Cgb
      IF (icode.eq.219 .or. icode.eq.101 .or. icode.eq.102 .or. 
     &     icode.eq.110 .or. icode.eq.214 .or. icode.eq.215 .or.
     &     icode.eq.217 .or. icode.eq.221 ) then
Cgb
         call UpdateCurrPart(mreg,0,0,zerzer,icode,
     &        xsco,ysco,zsco,JEND)
Cgb
      endif
c
      if ( ldead ) then
         idead (idcurr) = icode
c         idead (idcurr) = intcode(idcurr)
         if(idbflg.gt.1) then
            write(*,*)'USDRAW: ldead Now idcurr, idead(idcurr) = ',
     &           idcurr,idead(idcurr),icode
         endif
         vxf (idcurr) = sngl(xsco)
         vyf (idcurr) = sngl(ysco)
         vzf (idcurr) = sngl(zsco)
      end if
c
      RETURN
*=== End of subrutine Mgdraw ==========================================*
      END

