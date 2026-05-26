*$ CREATE USREOU.FOR
*COPY USREOU
*
*=== Usreou ===========================================================*
*
      SUBROUTINE USREOU

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1991-2005      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     USeR Event OUtput: this routine is called at the end of each     *
*     event                                                            *
*                                                                      *
*     Created on 01 january 1991   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 09-apr-99     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
      include '(CASLIM)'
      include 'mgdraw.inc'
      integer ii
c
c filling the root tree for each block
c
c      write(*,*)'sono in usreou'
c
      if(idbflg.gt.0) then
         write(*,*)'ev= ',numev
      endif
c
      if( (ifragflag.eq.0) .or. ( (ifragflag.gt.0).and.(nump.gt.1)))then
         write(outunit,*) ncase,nump,ncross
         do ii = 1,nump
c            call fillpart(idpa(ii), igen(ii), icha(ii), 
c     &           numreg(ii), iba(ii), idead(ii), jpa(ii), vxi(ii), 
c     &           vyi(ii), vzi(ii), vxf(ii), vyf(ii), vzf(ii), px(ii), 
c     &           py(ii),pz(ii),pxf(ii),pyf(ii),pzf(ii),amass(ii), 
c     &           tempo(ii),tof(ii),trlen(ii))
c           

*            write(outunit,200)idpa(ii), igen(ii), icha(ii), 
             write(outunit,*)idpa(ii), igen(ii), icha(ii), 
     &            numreg(ii), iba(ii), idead(ii), jpa(ii), vxi(ii), 
     &            vyi(ii), vzi(ii), vxf(ii), vyf(ii), vzf(ii), px(ii), 
     &            py(ii),pz(ii),pxf(ii),pyf(ii),pzf(ii),amass(ii), 
     &            tempo(ii),tof(ii),trlen(ii)
         end do
C
c  boundary crossing info
c         
         do ii = 1,ncross
c            write(outunit,800) idcross(ii), nregcross(ii), xcross(ii), 
            write(outunit,*) idcross(ii), nregcross(ii),
     &           nregoldcross(ii),xcross(ii), 
     &           ycross(ii), zcross(ii),pxcross(ii), pycross(ii),
     &           pzcross(ii), amacross(ii), chcross(ii), tcross(ii)
         end do
c     
c     call myusreou()
      endif
c     
      if(idbflg .gt. 0)then
         call dump_common()
         write(*,*)'***************************************************'
         write(*,*)'***************** EVENT END ***********************'
         write(*,*)'***************************************************'
         write(*,*)''
      endif

c
c      write(*,*)'esco da  usreou'
      RETURN
*=== End of subroutine Usreou =========================================*
      END

