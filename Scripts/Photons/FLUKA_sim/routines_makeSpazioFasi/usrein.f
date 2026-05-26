COPY USREIN
*
*=== Usrein ===========================================================*
*
      SUBROUTINE USREIN

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      INCLUDE '(CASLIM)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1991-2005      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     USeR Event INitialization: this routine is called before the     *
*     showering of an event is started, but after the source particles *
*     of that event have been already loaded on the stack              *
*                                                                      *
*     Created on 01 january 1991   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 09-apr-99     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
c
      include 'mgdraw.inc'
      integer ii  
c
      if(idbflg .gt. 0)then
         write(*,*)'***************************************************'
         write(*,*)'                   EVENT',ncase
         write(*,*)'***************************************************'
         write(*,*)''
      endif

      numev = ncase
c
      do ii = 1,min(nump,maxnump)
         idpa(ii) = 0
         igen(ii) = 0
         icha(ii) = 0
         numreg(ii) = 0
         iba(ii) = 0
         idead(ii) = 0
         jpa(ii) = 0
         vxi(ii) = 0.
         vyi(ii) = 0.
         vzi(ii) = 0.
         vxf(ii) = 0.
         vyf(ii) = 0. 
         vzf(ii) = 0.
         px(ii) = 0.
         py(ii) = 0.
         pz(ii) = 0.
         pxf(ii) = 0.
         pyf(ii) = 0.
         pzf(ii) = 0.
         amass(ii) = 0.
         tempo(ii) = 0.
         tof(ii) = 0.
         trlen(ii) = 0.
c
         idfluka(ii) = 0   ! aux variables for particle latching
c
      end do
      nump = 0
c
c      
      do ii = 1,min(numint,maxint)
         xint(ii) = 0
         yint(ii) = 0
         zint(ii) = 0
         intpa(ii) = 0
         intcode(ii) = 0
      end do
      numint = 0
c
      ncross = 0
      do ii= 1, nmaxcross
         idcross(ii) = 0
         nregcross(ii) = 0
         nregoldcross(ii) = 0
         nlattcross(ii) = 0
         nlattoldcross(ii) = 0
         pxcross(ii) = 0.
         pycross(ii) = 0.
         pzcross(ii) = 0.
         xcross(ii) = 0.
         ycross(ii) = 0.
         zcross(ii) = 0.
         tcross(ii) = 0.
         chcross(ii) = 0.
         amacross(ii) = 0.
      end do
      ncross = 0
c
c root event initialization
c
c      CALL myusrein(NCASE)
c

      RETURN
*=== End of subroutine Usrein =========================================*
      END


