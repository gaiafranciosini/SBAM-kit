!$ CREATE USBREA.FOR
!COPY USBREA
!
!=== Usbrea ===========================================================*
!
      PROGRAM USBREA_2_MHD
      

      CHARACTER FILE*1024,RUNTIT*80,RUNTIM*32,CHSTAT*10,CHCASE*15
      CHARACTER MYARG*1024

!
!----------------------------------------------------------------------*
!                                                                      *
!     Copyright (C) 1989-2008      by  Alberto Fasso`, Alfredo Ferrari *
!     All Rights Reserved.                    & Paola Sala             *
!                                                                      *
!----------------------------------------------------------------------*
!
!
!=== usrbin ===========================================================*
!
!----------------------------------------------------------------------*
!     Module USRBIN:                                                   *
!     A. Ferrari & A. Fasso': user defined binnings                    *
!                                                                      *
!          Last change A. Ferrari 12-feb-1997                          *
!                                                                      *
!                                                                      *
!     Up to MXUSBN user defined binnings are allowed                   *
!            lusbin = logical flag, .true. if at least 1 user defined  *
!                     binning is used                                  *
!            lusevt = logical flag, .true. if at least 1 user defined  *
!                     binning event by event is used                   *
!            lustkb = logical flag, .true. if at least 1 user defined  *
!                     track-length binning is used. Track length bin-  *
!                     ning are recognized as binnings where accurate   *
!                     deposition along the track is requested but      *
!                     for generalized particles other than 208/211     *
!            nusrbn = number of user defined binnings used             *
!            itusbn = type of binning: 0 = cartesian, .ne. 0 = RZ      *
!            idusbn = distribution to be scored: the usual values are  *
!                     allowed.                                         *
!            titusb = binning name                                     *
!            ipusbn = logical unit to print the results on: formatted  *
!                     if > 0, unformatted if < 0                       *
!            kbusbn = initial location in blank common of the consi-   *
!                     dered binning (for index 1, real*8 address for   *
!                     lsngbn = .false., real*4 address for lsngbn      *
!                            = .true.)                                 *
!            nxbin  = number of x (r for RZ) intervals                 *
!            nybin  = number of y (1 for RZ) intervals                 *
!            nzbin  = number of z intervals                            *
!         xlow/high = minimum and maximum x (r   for R-Phi-Z)          *
!         ylow/high = minimum and maximum y (phi for R-Phi-Z)          *
!         zlow/high = minimum and maximum z                            *
!            dxusbn = x (r) bin width                                  *
!            dyusbn = y bin width                                      *
!            dzusbn = z bin width                                      *
!            tcusbn = time cut-off (seconds) for this binning          *
!            bkusbn = 1st Birk's law parameter for this binning        *
!                     (meaningful only for energy scoring)             *
!            b2usbn = 2nd Birk's law parameter for this binning        *
!                     (meaningful only for energy scoring)             *
!            xaxusb = x-axis offset for R-Z binning (not possible for  *
!                     R-Phi-Z), just for back-compatibility, now it    *
!                     should be done through rototraslations           *
!            yaxusb = y-axis offset for R-Z binning (not possible for  *
!                     R-Phi-Z), just for back-compatibility, now it    *
!                     should be done through rototraslations           *
!            levtbn = logical flag for binning to be printed at the    *
!                     end of each event only                           *
!            lntzer = logical flag for printing only non zero cells    *
!            ltrkbn = logical flag for flagging track-length binnings  *
!            lsngbn = logical flag for single precision storing of     *
!                     data                                             *
!            krtnbn = flag of a possible pre-defined rotation (if >0)  *
!                     to be applied to the binnings before scoring     *
!                                                                      *
!----------------------------------------------------------------------*
!
      PARAMETER ( MXUSBN = 400 )
      LOGICAL LUSBIN, LEVTBN, LNTZER, LUSEVT, LUSTKB, LTRKBN, LSNGBN
      CHARACTER*10 TITUSB
      COMMON /UUSRBN/  XLOW  (MXUSBN), XHIGH (MXUSBN), YLOW  (MXUSBN),&
                      YHIGH (MXUSBN), ZLOW  (MXUSBN), ZHIGH (MXUSBN),&
                      DXUSBN(MXUSBN), DYUSBN(MXUSBN), DZUSBN(MXUSBN),&
                      TCUSBN(MXUSBN), BKUSBN(MXUSBN), B2USBN(MXUSBN),&
                      XAXUSB(MXUSBN), YAXUSB(MXUSBN),&
                      NXBIN (MXUSBN), NYBIN (MXUSBN), NZBIN (MXUSBN),&
                      ITUSBN(MXUSBN), IDUSBN(MXUSBN), KBUSBN(MXUSBN),&
                      IPUSBN(MXUSBN), KRTNBN(MXUSBN), LEVTBN(MXUSBN),&
                      LNTZER(MXUSBN), LSNGBN(MXUSBN), LTRKBN(MXUSBN),&
                      NUSRBN, LUSBIN, LUSEVT, LUSTKB
      COMMON /USRCH/  TITUSB(MXUSBN)
      double precision PIPIPI
      PARAMETER ( PIPIPI = 3.141592653589793D+00 )
      PARAMETER ( MXDUMM = 40000000 )
      PARAMETER ( AINFNT = 1.E+29 )
      PARAMETER ( ZERZER = 0.E+00 )
      DIMENSION GMSTOR (MXDUMM), GBSTOR (MXDUMM),&
                JB(MXUSBN), LIO(MXUSBN)
      CHARACTER CHFORM*20
      LOGICAL LOPEN, LSTATI, LREFLX, LREFLY, LREFLZ


    integer(4):: iargfile1,iargfile2
    logical dose2Gy

!	asch mhd format specifiers >>>
	integer(4):: nn(3)
	real(4):: off(3), hs(3)
!	asch mhd format specifiers <<<

	if(command_argument_count().lt.2 .or. command_argument_count().gt.3) then
	write(*,*) 'converter from fluka unformatted bnn file to mhd format'
    write(*,*) 'usage: bnn2mhd [-Gy] inputfname outputfname'
    write(*,*) 'option -Gy : convert fluka dose in Gy'
	stop
	endif

    dose2Gy = .false.
    iargfile1 = 1
    iargfile2 = 2

    do i = 1,command_argument_count()
        CALL get_command_argument(i,myarg)
        if(trim(myarg).eq.'-Gy') then
            dose2Gy = .true.
            select case(i)
            case(1)
                iargfile1 = 2
                iargfile2 = 3
            case(2)
                iargfile1 = 1
                iargfile2 = 3
            case(3)
                iargfile1 = 1
                iargfile2 = 2
            end select
        endif
    enddo
    if(dose2Gy .and. command_argument_count().ne.3) then
        write(*,*) 'missing file name'
        write(*,*) 'usage: bnn2mhd.x [-Gy] inputfname outputfname'
        stop
    endif




	
      
      DO 5555 I=1,MXDUMM
         GMSTOR(I)=0.E+00
         GBSTOR(I)=0.E+00
 5555 CONTINUE
      MIO   = 0
      NCTOT = 0
      MCTOT = 0
      WCTOT = 0.E+00
      LSTATI=.FALSE.
!      WRITE(*,'('' Type the input file: '',$)')
!      READ (*,'(A)')FILE
	CALL get_command_argument(iargfile1, file)
	write(*,*) 'Input file ',trim(file)
      LQ=LNNBLN(FILE)
      IF (LQ .LE. 0) GO TO 2000
      OPEN (UNIT=1,FILE=FILE,STATUS='OLD',FORM='UNFORMATTED')
      IRECRD = 0
      READ (1,ERR=100) RUNTIT,RUNTIM,WEIPRI,NCASE,MCASE,NBATCH
      IRECRD = IRECRD + 1
      WRITE(*,*)RUNTIT
      WRITE(*,*)RUNTIM
      WRITE(*,*)WEIPRI
      WRITE(*,*)NCASE
      WRITE(*,*)MCASE
      WRITE(*,*)NBATCH
      GO TO 200
 100  CONTINUE
      REWIND 1
      IRECRD = 0
      READ (1,ERR=110) RUNTIT,RUNTIM,WEIPRI,NCASE,NBATCH
      IRECRD = IRECRD + 1
      WRITE(*,*)RUNTIT
      WRITE(*,*)RUNTIM
      WRITE(*,*)WEIPRI
      WRITE(*,*)NCASE
      MCASE  = 0
      GO TO 200
 110  CONTINUE
      REWIND 1
      IRECRD = 0
      READ (1,ERR=120) RUNTIT,RUNTIM,WEIPRI,NCASE
      IRECRD = IRECRD + 1
      WRITE(*,*)RUNTIT
      WRITE(*,*)RUNTIM
      WRITE(*,*)WEIPRI
      WRITE(*,*)NCASE
      MCASE  = 0
      NBATCH = 1
      GO TO 200
 120  CONTINUE
      REWIND 1
      IRECRD = 0
      READ (1,ERR=150) RUNTIT,RUNTIM,NCASE
      IRECRD = IRECRD + 1
      WRITE(*,*)RUNTIT
      WRITE(*,*)RUNTIM
      WRITE(*,*)NCASE
      WEIPRI = REAL(NCASE)
      MCASE  = 0
      NBATCH = 1
      GO TO 200
 150  CONTINUE
      IRECRD = 0
      REWIND 1
      WRITE(*,*)' How many particle for this run?'
      READ(*,*)NCASE
      WEIPRI = REAL(NCASE)
      MCASE  = 0
      NBATCH = 1
 200  CONTINUE
      NCTOT  = NCTOT + NCASE
      MCTOT  = MCTOT + MCASE
      IF ( NCTOT .GE. 1000000000 ) THEN
         NCTOT = NCTOT - 1000000000
         MCTOT = MCTOT + 1
      END IF
      WCTOT  = WCTOT + WEIPRI
      KLAST  = 0
      KMAX   = 0
      DO 700 IB = 1, 1000
         NB = IB
         READ (1,ERR=400,END=1000) MB, TITUSB(NB), ITUSBN(NB),&
                         IDUSBN(NB),&
                         XLOW(NB),XHIGH(NB),NXBIN(NB),&
                         DXUSBN(NB),YLOW(NB),&
                         YHIGH(NB),NYBIN(NB),DYUSBN(NB),&
                         ZLOW(NB),ZHIGH(NB),NZBIN(NB),&
                         DZUSBN(NB),LNTZER(NB),BKUSBN(NB),&
                         B2USBN(NB),TCUSBN(NB)
         IRECRD = IRECRD + 1
         print *,MB,TITUSB(NB),ITUSBN(NB),IDUSBN(NB)
         print *,XLOW(NB),XHIGH(NB),NXBIN(NB),DXUSBN(NB)
         print *,YLOW(NB),YHIGH(NB),NYBIN(NB),DYUSBN(NB)
         print *,ZLOW(NB),ZHIGH(NB),NZBIN(NB),DZUSBN(NB)
         print *,LNTZER(NB),BKUSBN(NB)
         print *,B2USBN(NB),TCUSBN(NB)

         GO TO 500
  400    CONTINUE
!        BACKSPACE 1
         REWIND 1
         DO IR = 1, IRECRD
            READ (1)
         END DO
         READ (1,ERR=450,END=1000) MB, TITUSB(NB), ITUSBN(NB),&
                           IDUSBN(NB),&
                           XLOW(NB),XHIGH(NB),NXBIN(NB),&
                           DXUSBN(NB),YLOW(NB),&
                           YHIGH(NB),NYBIN(NB),DYUSBN(NB),&
                           ZLOW(NB),ZHIGH(NB),NZBIN(NB),&
                           DZUSBN(NB)
         IRECRD = IRECRD + 1
         LNTZER(NB)=.FALSE.
         TCUSBN(NB)=1.E+38
         BKUSBN(NB)=0.E+00
         B2USBN(NB)=0.E+00
         GO TO 500
  450    CONTINUE
!        BACKSPACE 1
         REWIND 1
         DO IR = 1, IRECRD
            READ (1)
         END DO
         READ (1,ERR=460) CHSTAT,IJX
         IRECRD = IRECRD + 1
         IF ( CHSTAT .EQ. 'STATISTICS' ) THEN
            LSTATI = .TRUE.
            GO TO 1000
         END IF
  460    CONTINUE
         STOP 'UNKNOWN RECORD TYPE'
  500    CONTINUE
         JB(IB)=NB
!        JB(IB)=MB
         IPUSBN(NB)=1
         IF ( IDUSBN(NB) .NE. 208 .AND. IDUSBN(NB)&
              .NE. 211 .AND. ITUSBN(NB) .GE. 10 ) THEN
            LTRKBN(NB) = .TRUE.
         ELSE
            LTRKBN(NB) = .FALSE.
         END IF
         K0 = KLAST + 1
         KBUSBN (NB) = K0
         K1 = NXBIN(NB) * NYBIN(NB) * NZBIN(NB) + K0 -1
         KLAST = K1
         KMAX  = MAX ( KMAX, KLAST )
         READ (1) (GMSTOR(J), J = K0, K1)
         IRECRD = IRECRD + 1
  700 CONTINUE
 1000 CONTINUE
      IB = IB-1
      NUSRBN = IB
      IF ( LSTATI ) THEN
         KLAST = 0
         DO JJ=1,IB
            NB = JB(JJ)
            K0 = KLAST + 1
            K1 = NXBIN(NB) * NYBIN(NB) * NZBIN(NB) + K0 -1
            KLAST = K1
            READ (1) (GBSTOR(J), J = K0, K1)
         END DO
      END IF
      CLOSE (UNIT=1)
 2000 CONTINUE
!      WRITE(*,'('' Type the output file name:'',$)')
!      READ (*,'(A)') FILE
	CALL get_command_argument(iargfile2, file)
! Start_VAX_seq
!     OPEN (UNIT=1,FILE=FILE,STATUS='NEW',FORM='FORMATTED',
!    &      CARRIAGECONTROL='LIST')
!     LUNOUT=6
! End_VAX_seq
! Start_UNIX_seq
!      OPEN (UNIT=1,FILE=FILE,STATUS='UNKNOWN',FORM='FORMATTED')
!      LUNOUT=11
! End_UNIX_seq
! write to m3d format
!     &                          XLOW(NB),XHIGH(NB),NXBIN(NB),DXUSBN(NB),
!     &                          YLOW(NB),YHIGH(NB),NYBIN(NB),DYUSBN(NB),
!     &                          ZLOW(NB),ZHIGH(NB),NZBIN(NB),DZUSBN(NB)
  NB=1
  nn=(/NXBIN(NB),NYBIN(NB),NZBIN(NB)/)
  hs=(/DXUSBN(NB),DYUSBN(NB),DZUSBN(NB)/)
  off=(/XLOW(NB),YLOW(NB),ZLOW(NB)/) + hs*0.5

  print *,'hs',hs
  ! print *,'jb',jb
    ! DO IB = 1,NUSRBN
    !  NB   = JB    (IB)
    !  print *,nb,KBUSBN(NB)
    ! enddo

  ! ITUHLP = MOD(ITUSBN(NB),10)        
  ! print *,'ITUHLP',ITUHLP 
  ! stop
  
  WRITE (*,*) 'saving map to ',TRIM(FILE)
  OPEN(UNIT=77,FILE=trim(FILE),STATUS='UNKNOWN',FORM='FORMATTED')
  write(77,*)
  write(77,*) 'ObjectType = Image'
  write(77,*) 'NDims = 3'
  write(77,*) 'DimSize =',nn
  write(77,*) 'BinaryData = True'
  write(77,*) 'BinaryDataByteOrderMSB = False'
  write(77,*) 'CompressedData = False'
  write(77,*) 'TransformMatrix = 1 0 0 0 1 0 0 0 1'
  write(77,*) 'Offset = ', off*10.
  write(77,*) 'ElementSpacing =',hs*10.
  write(77,*) 'ElementType = MET_FLOAT'
  write(77,*) 'ElementDataFile = LOCAL'
  close(77)



  OPEN(UNIT=77,FILE=trim(FILE),STATUS='OLD',FORM='UNFORMATTED',ACCESS='STREAM',POSITION='APPEND')
  K0 = KBUSBN(NB)
  K1 = NXBIN(NB) * NYBIN(NB) * NZBIN(NB) + K0 - 1
  if(dose2Gy) then
    write(*,*) 'dose converted from GeV/g to Gy'
    write(77) GMSTOR(K0:K1)*1.602176462E-7
  else
    write(77) GMSTOR(K0:K1)
  endif
  close(77) 
  stop 

!  +-------------------------------------------------------------------*
!  |                          Writes on formatted files the results!!!
      DO 6900 IB = 1,NUSRBN
         NB   = JB    (IB)
         LOUT = IPUSBN(NB)
!  |  +----------------------------------------------------------------*
!  |  |  Formatted output requested:
         IF ( LOUT .GT. 0 ) THEN
 6650       CONTINUE
!  |  |  +-------------------------------------------------------------*
!  |  |  |  Loop on already opened units:
            DO 6700 IO = 1,MIO
               IF (LIO(IO) .EQ. LOUT) GO TO 6800
 6700       CONTINUE
!  |  |  |
!  |  |  +-------------------------------------------------------------*
            INQUIRE ( UNIT = LOUT, OPENED = LOPEN, FORM = CHFORM )
!  |  |  +-------------------------------------------------------------*
!  |  |  |  This unit is already opened:
            IF ( LOPEN ) THEN
!  |  |  |  +----------------------------------------------------------*
!  |  |  |  |  Open for unformatted access,while formatted is requested:
               IF ( CHFORM .EQ. 'UNFORMATTED' ) THEN
                  WRITE (LUNOUT,8800)LOUT,NB,LOUT+1
                  LOUT = LOUT + 1
                  GO TO 6650
               END IF
!  |  |  |  |
!  |  |  |  +----------------------------------------------------------*
!  |  |  |
!  |  |  +-------------------------------------------------------------*
!  |  |  |  Unit LOUT is not yet opened:
            ELSE
               MIO = MIO + 1
               LIO(MIO) = LOUT
! Start_IBM_seq
!              WRITE(CHOUT,'(A3)') LOUT
!              LENGTH = MIN(8, INDEX(INPFIL,' ')-1, INDEX(INPFIL,'.')-1)
!              CALL FILEINF(IR,'RECFM','FBA','LRECL',133,'BLKSIZE',1330)
! End_IBM_seq
               OPEN ( UNIT = LOUT, STATUS = 'NEW', FILE='dummy',&
                      ERR = 6750 )
               WRITE ( CHCASE, '(F15.0)' ) DBLE (NCASE)&
                                + 1.D+09 * DBLE (MCASE)
               WRITE (LOUT,8950) RUNTIT, RUNTIM, CHCASE, WEIPRI
               GO TO 6800
 6750          CONTINUE
               WRITE (LUNOUT,8900)LOUT,NB,LOUT+1
               LOUT = LOUT + 1
               MIO  = MIO - 1
               GO TO 6650
            END IF
!  |  |  |
!  |  |  +-------------------------------------------------------------*
 6800       CONTINUE
            ITUHLP = MOD(ITUSBN(NB),10)
!  |  |  +-------------------------------------------------------------*
!  |  |  |  Plain region binning:
            IF ( ITUHLP .EQ. 2 ) THEN
               WRITE (LOUT,8500) NB,TITUSB(NB),IDUSBN(NB),&
                                MAX (NXBIN(NB),NYBIN(NB),NZBIN(NB)),&
                                NINT(XLOW(NB)),NINT(XHIGH(NB)),&
                                NINT(DXUSBN(NB)),NINT(YLOW(NB)),&
                                NINT(YHIGH(NB)),NINT(DYUSBN(NB)),&
                                NINT(ZLOW(NB)),NINT(ZHIGH(NB)),&
                                NINT(DZUSBN(NB))
               LREFLX = .FALSE.
               LREFLY = .FALSE.
               LREFLZ = .FALSE.
!  |  |  |
!  |  |  +-------------------------------------------------------------*
!  |  |  |  Region/lattice/user binning:
            ELSE IF ( ITUHLP .EQ. 8 ) THEN
               WRITE (LOUT,8520) NB,TITUSB(NB),IDUSBN(NB),NXBIN(NB),&
                                NINT(XLOW(NB)),NINT(XHIGH(NB)),&
                                NINT(DXUSBN(NB)),NYBIN(NB),&
                                NINT(YLOW(NB)),NINT(YHIGH(NB)),&
                                NINT(DYUSBN(NB)),&
                                ZLOW(NB),ZHIGH(NB),NZBIN(NB),&
                                DZUSBN(NB)
               LREFLX = .FALSE.
               LREFLY = .FALSE.
               LREFLZ = .FALSE.
!  |  |  |
!  |  |  +-------------------------------------------------------------*
!  |  |  |  R-Z only binning:
            ELSE IF ( ( ITUHLP .EQ. 1 .OR. ITUHLP .EQ. 7 ) .AND.&
                        NYBIN (NB) .LE. 1 ) THEN
               WRITE (LOUT,8550) NB, TITUSB(NB), IDUSBN(NB),&
                                XLOW(NB),XHIGH(NB),NXBIN(NB),DXUSBN(NB),&
                                ZLOW(NB),ZHIGH(NB),NZBIN(NB),DZUSBN(NB),&
                                XAXUSB(NB),YAXUSB(NB)
               LREFLX = .FALSE.
               LREFLY = .FALSE.
               LREFLZ = ITUHLP .EQ. 7
!  |  |  |
!  |  |  +-------------------------------------------------------------*
!  |  |  |  R-Phi-Z binning:
            ELSE IF ( ITUHLP .EQ. 1 .OR. ITUHLP .EQ. 7 ) THEN
               WRITE (LOUT,8570) NB, TITUSB(NB), IDUSBN(NB),&
                                XLOW(NB),XHIGH(NB),NXBIN(NB),DXUSBN(NB),&
                                YLOW(NB),YHIGH(NB),NYBIN(NB),DYUSBN(NB),&
                                ZLOW(NB),ZHIGH(NB),NZBIN(NB),DZUSBN(NB)
               LREFLX = .FALSE.
               LREFLY = .FALSE.
               LREFLZ = ITUHLP .EQ. 7
!  |  |  |
!  |  |  +-------------------------------------------------------------*
!  |  |  |  Cartesian binning:
            ELSE
               WRITE (LOUT,8600) NB, TITUSB(NB), IDUSBN(NB),&
                                XLOW(NB),XHIGH(NB),NXBIN(NB),DXUSBN(NB),&
                                YLOW(NB),YHIGH(NB),NYBIN(NB),DYUSBN(NB),&
                                ZLOW(NB),ZHIGH(NB),NZBIN(NB),DZUSBN(NB)
               LREFLX = ITUHLP .EQ. 3 .OR. ITUHLP .EQ. 6
               LREFLY = ITUHLP .EQ. 4 .OR. ITUHLP .EQ. 6
               LREFLZ = ITUHLP .EQ. 5 .OR. ITUHLP .EQ. 6
            END IF
!  |  |  |
!  |  |  +-------------------------------------------------------------*
            IF ( LREFLX ) WRITE (LOUT,8410)
            IF ( LREFLY ) WRITE (LOUT,8420)
            IF ( LREFLZ ) WRITE (LOUT,8430)
            IF ( ITUSBN(NB)/10 .GE. 1 ) WRITE (LOUT,8440)
            IF ( LTRKBN(NB) ) WRITE (LOUT,8465)
            IF ( LNTZER(NB) ) WRITE (LOUT,8470)
            IF ( BKUSBN(NB) .GT. ZERZER ) WRITE (LOUT,8480) BKUSBN(NB),&
                 B2USBN(NB)
            IF ( TCUSBN(NB) .LT. AINFNT ) WRITE (LOUT,8490) TCUSBN(NB)
            K0 = KBUSBN(NB)
            K1 = NXBIN(NB) * NYBIN(NB) * NZBIN(NB) + K0 - 1
            WRITE (LOUT,*) 'K0,K1',K0,K1
            
!            WRITE (LOUT,8700) (GMSTOR(J), J = K0, K1)
!  |  |  +-------------------------------------------------------------*
!  |  |  |  Statistics is available:
            IF ( LSTATI ) THEN
!  |  |  |  +----------------------------------------------------------*
!  |  |  |  |  Plain region binning:
               IF ( ITUHLP .EQ. 2 ) THEN
                  WRITE (LOUT,8710)
!  |  |  |  |
!  |  |  |  +----------------------------------------------------------*
!  |  |  |  |  Region/lattice/user binning:
               ELSE IF ( ITUHLP .EQ. 8 ) THEN
                  WRITE (LOUT,8720)
!  |  |  |  |
!  |  |  |  +----------------------------------------------------------*
!  |  |  |  |  R-Z only binning:
               ELSE IF ( ( ITUHLP .EQ. 1 .OR. ITUHLP .EQ. 7 ) .AND.&
                           NYBIN (NB) .LE. 1 ) THEN
                  WRITE (LOUT,8730)
!  |  |  |  |
!  |  |  |  +----------------------------------------------------------*
!  |  |  |  |  R-Phi-Z binning:
               ELSE IF ( ITUHLP .EQ. 1 .OR. ITUHLP .EQ. 7 ) THEN
                  WRITE (LOUT,8740)
!  |  |  |  |
!  |  |  |  +----------------------------------------------------------*
!  |  |  |  |  Cartesian binning:
               ELSE
                  WRITE (LOUT,8750)
               END IF
!  |  |  |  |
!  |  |  |  +----------------------------------------------------------*
               WRITE (LOUT,8700) (100.E+00*GBSTOR(J), J = K0, K1)
            END IF
!  |  |  |
!  |  |  +-------------------------------------------------------------*
         END IF
!  |  |
!  |  +----------------------------------------------------------------*
 6900 CONTINUE
!  |
!  +-------------------------------------------------------------------*
 8410 FORMAT (6X,'+/- X symmetry requested and implemented' )
 8420 FORMAT (6X,'+/- Y symmetry requested and implemented' )
 8430 FORMAT (6X,'+/- Z symmetry requested and implemented' )
 8440 FORMAT (6X,'accurate deposition along the tracks requested' )
 8450 FORMAT (6X,'unnormalized data will be printed event by event')
 8460 FORMAT (6X,'normalized (per unit volume) data will be printed',&
                 ' at the end of the run')
 8465 FORMAT (6X,'this is a track-length binning')
 8470 FORMAT (6X,'only non-zero cells will be printed')
 8480 FORMAT (6X,'Birk law quenching factors: ',1P,E9.2,' g/(GeV cmq)',&
              ',',E9.2,' [g/(GeV cmq)]^2')
 8490 FORMAT (6X,'Time cut off set at:',1P,E9.2,' s')
 8493 FORMAT (6X,'Coordinate transformation #',I4,&
                 ' applied to this binning')
 8496 FORMAT (6X,'Single precision storage requested for this binning')
 8500 FORMAT ('1',/,3X,'Region    binning n. ',I3,'  "',A10,&
              '" , generalized particle n. ',I4,&
              /,6X,I5,' bins corresponding to the region sets:'&
              /,6X,'from region ',I5,' to region ',I5,' in step of ',I5,&
              ' regions, or',&
              /,6X,'from region ',I5,' to region ',I5,' in step of ',I5,&
              ' regions, or',&
              /,6X,'from region ',I5,' to region ',I5,' in step of ',I5,&
              ' regions',&
              /,6X,'Data follow in an array A(ir), format ',&
              '(1(5x,1p,10(1x,e11.4)))',/)
 8520 FORMAT (/,3X,'Reg/Lat/U binning n. ',I3,'  "',A10,&
              '" , generalized particle n. ',I4,&
              /,6X,I5,' bins corresponding to the region-related set:'&
              /,6X,'from region ',I5,' to region ',I5,' in step of ',I5,&
              ' regions, and',&
              /,6X,I5,' bins corresponding to the lattice-related set:'&
              /,6X,'from lattice',I5,' to lattice',I5,' in step of ',I5,&
              ' lattices, and',&
              /,6X,'U coordinate: from ',1P,E11.4,&
              ' to ',E11.4,' ux, ',0P,I5,' bins (',1P,E11.4,' ux wide)',&
              /,6X,'Data follow in an array A(ir,il,iu), format ',&
              '(1(5x,1p,10(1x,e11.4)))',/)
 8550 FORMAT ('1',/,3X,'R - Z     binning n. ',I3,'  "',A10,&
              '" , generalized particle n. ',I4,&
              /,6X,'R coordinate: from ',1P,E11.4,&
              ' to ',E11.4,' cm, ',0P,I5,' bins (',1P,E11.4,' cm wide)',&
              /,6X,'Z coordinate: from ',1P,E11.4,&
              ' to ',E11.4,' cm, ',0P,I5,' bins (',1P,E11.4,' cm wide)',&
              /,6X,'axis coordinates: X =',1P,E11.4,', Y = ',1P,E11.4,&
              ' cm',&
              /,6X,'Data follow in a matrix A(ir,iz), format ',&
              '(1(5x,1p,10(1x,e11.4)))',/)
 8570 FORMAT ('1',/,3X,'R-Phi-Z   binning n. ',I3,'  "',A10,&
              '" , generalized particle n. ',I4,&
              /,6X,'R coordinate: from ',1P,E11.4,&
            ' to ',E11.4,'  cm, ',0P,I5,' bins (',1P,E11.4,' cm  wide)',&
              /,6X,'P coordinate: from ',1P,E11.4,&
            ' to ',E11.4,' rad, ',0P,I5,' bins (',1P,E11.4,' rad wide)',&
              /,6X,'Z coordinate: from ',1P,E11.4,&
            ' to ',E11.4,'  cm, ',0P,I5,' bins (',1P,E11.4,' cm  wide)',&
              /,6X,'Data follow in a matrix A(ir,ip,iz), format ',&
              '(1(5x,1p,10(1x,e11.4)))',/)
 8600 FORMAT ('1',/,3X,'Cartesian binning n. ',I3,'  "',A10,&
              '" , generalized particle n. ',I4,&
              /,6X,'X coordinate: from ',1P,E11.4,&
              ' to ',E11.4,' cm, ',0P,I5,' bins (',1P,E11.4,' cm wide)',&
              /,6X,'Y coordinate: from ',1P,E11.4,&
              ' to ',E11.4,' cm, ',0P,I5,' bins (',1P,E11.4,' cm wide)',&
              /,6X,'Z coordinate: from ',1P,E11.4,&
              ' to ',E11.4,' cm, ',0P,I5,' bins (',1P,E11.4,' cm wide)',&
              /,6X,'Data follow in a matrix A(ix,iy,iz), format ',&
              '(1(5x,1p,10(1x,e11.4)))',/)
 8650 FORMAT (/,3X,'Binning n:',I5,', "',A10,'",  Event #:',I6,&
              ', Primary(s) weight ',1P,E11.4)
 8700 FORMAT (1(5X,1P,10(1X,E11.4)))
 8710 FORMAT (/,6X,'Percentage errors follow in an array A(ir),',&
              ' format (1(5x,1p,10(1x,e11.4)))',/)
 8720 FORMAT (/,6X,'Percentage errors follow in an array A(ir,il,iu),',&
              ' format (1(5x,1p,10(1x,e11.4)))',/)
 8730 FORMAT (/,6X,'Percentage errors follow in a matrix A(ir,iz),',&
              ' format (1(5x,1p,10(1x,e11.4)))',/)
 8740 FORMAT (/,6X,'Percentage errors follow in a matrix A(ir,ip,iz),',&
              ' format (1(5x,1p,10(1x,e11.4)))',/)
 8750 FORMAT (/,6X,'Percentage errors follow in a matrix A(ix,iy,iz),',&
              ' format (1(5x,1p,10(1x,e11.4)))',/)
 8800 FORMAT (/,5X,'***** I/O unit ',I4,' cannot be opened to write',&
              ' formatted (unformatted) data for binning n. ',I3,&
              ' *****',/,5X,'***** it is ',&
              ' already opened in unformatted (formatted) I/O mode ',&
              ' *****',/,5X,'***** try with unit ',I4,' *****')
 8900 FORMAT (/,5X,'***** I/O unit ',I4,' cannot be opened to write',&
              ' formatted (unformatted) data for binning n. ',I3,&
              ' *****',/,5X,'***** because of an unrecognized error ',&
              ' *****',/,5X,'***** try with unit ',I4,' *****')
 8950 FORMAT (/,1X,'*****',2X,A80,2X,'*****',/,/,10X,A32,/,/,&
              10X,'Total number of particles followed ',A14,', for a ',&
              'total weight of ',1P,E11.4,/)
 8970 FORMAT (/,1X,'*****',2X,A80,2X,'*****',/,/,10X,A32,/,/,&
              10X,'Total number of particles to be followed ',I7,&
              ', event by event',/)
      STOP
      END

!$ CREATE LNNBLN.FOR
!COPY LNNBLN
!
!=== Lnnbln ===========================================================*
!
      INTEGER FUNCTION LNNBLN (CARD)

!
      CHARACTER CARD*(*)
!
      LENGTH = LEN (CARD)
      DO 100 LQ = LENGTH, 1, -1
         IF ( CARD (LQ:LQ).NE.' ') GO TO 200
  100 CONTINUE
      LQ     = 0
  200 CONTINUE
      LNNBLN = LQ
      RETURN
!=== End of function Lnnbln ===========================================*
      END

