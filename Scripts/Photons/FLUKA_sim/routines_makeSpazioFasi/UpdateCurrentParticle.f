c-------------------------------------------------------------------------
      subroutine UpdateCurrPart(MREG,ICPART,IBPART,AMPART,ICODE,
     &     XX,YY,ZZ,JFLAG)
c--------------------------------------------------------------------------
c
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      INCLUDE '(TRACKR)'
      INCLUDE '(CASLIM)'
      INCLUDE '(BEAMCM)'
      INCLUDE '(RDCYCM)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(GENSTK)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(PAPROP)'

      include 'mgdraw.inc'
      integer JFLAG, JAPART, JZPART 
      logical OLDONE
c
c  JFLAG = 1 -> init tracking,  JFLAG = 0 -> tracking of the particle during the event
c
      integer idfluka_now, idpadre
c
c no optical photons
c
c      if(numev.eq.147) then
c         idbflg = 10
c         else
c            idbflg = 0
c      endif

      if(jtrack.eq.-1) then
         return       
      endif
c
      idfluka_now = ISPUSR(MKBMX2)
c      idbflg = 3
c
c    ______________     START  INIT   ________________________ 
c
c     Event initialinzing of the kinematics tracking
c
      if(jflag.eq.JSTART) then
         IONID  = ILOFLK (NPFLKA)
         idfluka(1) = 1
         nump = 1
         idcurr = 1
         idpa(idcurr) = 0
         tempo(idcurr) = sngl( agestk(npflka))
         tof(idcurr) = 0.
         trlen(idcurr) = 0.
         numreg(idcurr) = nrgflk(npflka)
         igen(idcurr) = 1
         idead(idcurr) = 0
         idflukamax = 1

C         write(*,*) 'ionid 0 = ',ionid
C         if( ionid .GE. 100000 .AND. LRADDC (NPFLKA) ) then
       if( ionid .GE. 100000 ) then
C            write(*,*) 'ionid 1 = ',ionid
            JAPART = MOD ( ILOFLK (NPFLKA), 100000  )  / 100
            JZPART = MOD ( ILOFLK (NPFLKA), 10000000 ) / 100000
            AAA = JAPART
            ZZZ = JZPART
            AMPART = JAPART * AMUC12
     &           + EMVGEV * EXMSAZ (AAA,ZZZ,.TRUE.,IZDUMM) 
            jpa (idcurr) = -2
         else if ( abs(ionid) .ge. 10000 ) then
            call dcdion ( ionid )
C            write(*,*) 'ionid 3 = ',ionid
            JZPART = ichrge (ionid)
            JAPART =  ibarch (ionid)
            AMPART  = am (ionid)
            jpa (idcurr) = -2
         else if ( ionid .lt. -6) then
C            write(*,*) 'ionid 2 = ',ionid
            JZPART = icheav(-ionid)
            JAPART = ibheav (-ionid)
            AMPART  = amnhea (-ionid)
            jpa (idcurr) = ionid
         else
            JZPART = ichrge (ionid)
            JAPART =  ibarch (ionid)
            AMPART  = am (ionid)            
            jpa (idcurr) = ionid
         endif
c         
         vxi (idcurr) = sngl(xflk (npflka))
         vyi (idcurr) = sngl(yflk (npflka))
         vzi (idcurr) = sngl(zflk (npflka))
         vxf (idcurr) = vxi(idcurr)
         vyf (idcurr) = vyi(idcurr)
         vzf (idcurr) = vzi(idcurr)
c
c when decay the momentum is defined <0 ?
c
         ppp          = abs(pmoflk (npflka))
         px  (idcurr) = sngl(txflk (npflka) * ppp)
         py  (idcurr) = sngl(tyflk (npflka) * ppp)
         pz  (idcurr) = sngl(tzflk (npflka) * ppp)
         pxf(idcurr) = px(idcurr)
         pyf(idcurr) = py(idcurr)
         pzf(idcurr) = pz(idcurr)
         icha(idcurr) = JZPART
         iba(idcurr) = JAPART
         amass(idcurr) = sngl(ampart)
c
         if(idbflg.gt.1) then
            write(*,*)' '
            write(*,*)'===== UpdateCurrPart : Ev =',numev,' == INIT = '
            write(*,*)' '
            write(*,*)' ionid= ',ionid,' npflka= ',npflka,
     &           ' Z= ',jzpart,' A= ',jApart            
            write(*,*)'age= ',tempo(idcurr),' reg= ',numreg(idcurr)
            write(*,*)'vert= ',xflk (npflka),yflk (npflka),zflk (npflka)
            write(*,*)'pmod= ',pmoflk (npflka)
            ppp= abs(pmoflk (npflka))
            write(*,*)px(idcurr),py(idcurr),pz(idcurr)
            write(*,*)'mass= ',ampart,' iloflk= ',iloflk(npflka),
     &           ' jtrack= ',jtrack
            write(*,*)'idcurr= ',idcurr,'idfluka(idcurr)= ',
     &           idfluka(idcurr),' idfluka_now= ',idfluka_now
            write(*,*)' '
            write(*,*)'===== UpdateCurrPart : END INIT = '
            write(*,*)' '
         endif
         return
      else
c    _____________________ END  INIT   ________________________ 
c
c    normal Kinematic tracking
c            
         if(idbflg.gt.1) then
            write(*,*)' '
            write(*,*)
     &           '--------- UpdateCurrPart : Ev =',numev,' ------- '
            write(*,*)' '
            write(*,*)'JFLAG= ',JFLAG,' Icode = ',icode," nump= ",nump
            write(*,*)'idcurr= ',idcurr,'idfluka(icurr)= ',
     &           idfluka(idcurr),'idfluka_now= ',idfluka_now
            write(*,*)'jtrack = ',jtrack,' jpa(idcurr)= ',jpa(idcurr)
            write(*,*)'PART_CHG = ',icpart,' PART_MASS = ',ampart
            write(*,*)'ptrack = ',ptrack,'px,y,z = ',
     &           ptrack*cxtrck,ptrack*cytrck,ptrack*cztrck
            write(*,*)'x(0) = ',xtrack(0),ytrack(0),ztrack(0)
            write(*,*)'x(1) = ',xtrack(1),ytrack(1),ztrack(1)
            write(*,*)'atrack = ',atrack,' time = ',tempo(idcurr),
     %           ' tof = ',tof(idcurr)
            write(*,*)'totlen = ',trlen(idcurr),' ctrack = ',ctrack,
     %           ' cmtrack = ',cmtrck
         endif
c
c  =========:  FOLLOWING THE SAME PARTICLE or an OLD ONE ?: ===============   
c
         OLDONE = .false.

         if (idfluka(idcurr).eq.idfluka_now ) then
            if (icode.ne.219) then
               oldone = .true.
            else
               if (idbflg.gt.1)  
     &         write(*,*) '--------- UpdateCurrPart : Icode 1 =',icode
            endif
c
c  new charged particle below threshold, giving the energy deposition to father
c
         else if( (ptrack.eq.0).and. (JFLAG.eq.JEND) )  then
            if (icode.ne.219) then
               oldone = .true.
            else
               if (idbflg.gt.1)  
     &         write(*,*) '--------- UpdateCurrPart : Icode 2 =',icode
            endif
            if(idbflg.gt.0) then
       write(*,*)'UpdateCurrentParticle ptrack = 0 searching for father'
            endif
            do ii = 1,numint
               if(idbflg.gt.0) then
                  write(*,*)'int= ',ii, 'x int= ',xint(ii),yint(ii),
     &                 zint(ii),' intpa= ',intpa(ii)
               endif
               if( (xx.eq.xint(ii)) .and. (yy.eq.yint(ii)) .and. 
     &              (zz.eq.zint(ii)) ) then
                  if(idgflg.gt.0) write(*,*)'father= ',ii
                  idcurr = intpa(ii)
                  return
               endif
            end do
         else
            do kk = 1,nump
               if (idfluka(kk).eq.idfluka_now ) then
                  if (icode.ne.219) then
                     oldone = .true.
                  else
               if (idbflg.gt.1)  
     &             write(*,*) 
     &             '--------- UpdateCurrPart : Icode 3 =',icode
                  endif
                  idcurr = kk
               endif
            end do
         endif
c
         if(oldone) then

            if(JFLAG.eq.JRUN) then
               vxf(idcurr) = sngl(xtrack(1))
               vyf(idcurr) = sngl(ytrack(1))
               vzf(idcurr) = sngl(ztrack(1))
            elseif(JFLAG.eq.JEND) then   
               vxf(idcurr) = sngl(xx)
               vyf(idcurr) = sngl(yy)
               vzf(idcurr) = sngl(zz)
C               idead(idcurr) = icode
            else                ! Non dovremmo essere qua
            write(*,*)'UpdateCurrentParticle: ERROR -> wrong JFLAG= ',
     &              JFLAG
            endif
            pxf(idcurr) = sngl(ptrack*cxtrck)
            pyf(idcurr) = sngl(ptrack*cytrck)
            pzf(idcurr) = sngl(ptrack*cztrck)
            tof(idcurr) = sngl(atrack) - tempo(idcurr)
            trlen(idcurr) = sngl(cmtrck)
         else
c     
c     =========:  FOLLOWING A NEW PARTICLE  : ===============   
c     
            if(idbflg.gt.0) then
               write(*,*)' '
               write(*,*)
     &     'UPdateCurrentParticle= --------- New part found ----------'
               write(*,*)' '
               write(*,*)
               write(*,*)'jtrack= ',jtrack,' cha= ',icpart,
     &              ' mass= ',ampart
               write(*,*)'nump= ',nump,' idfluka_now',idfluka_now
               write(*,*)'JRUN= ',JRUN,'numint= ',numint
               write(*,*) 'intcode =',icode
               write(*,*)'ispusr(mkbmx2)= ',ISPUSR(MKBMX2),
     &              'ispusr(1)= ',ISPUSR(1)
               write(*,*)'idcurr= ',idcurr,'idfluka(idcurr)= ',
     &              idfluka(idcurr),' idpa(idcurr) = ',idpa(idcurr)
               write(*,*)'vi= ',xx,yy,zz
               write(*,*)'pmod= ',ptrack,' ekin= ',etrack-ampart
            endif
c
            if(nump.eq.maxnump) then
               write(*,*)'Error: UpdateCurrPart:'
               write(*,*)'Max particle number exceeded : ev= ',ncase
               return
            endif
c
            nump = nump + 1
c
            idcurr = nump
            idfluka(idcurr) = idfluka_now 
            tempo(idcurr) = sngl(ATRACK)
            tof(idcurr) = 0.
            numreg(idcurr) = mreg
            icha(idcurr) = icpart
            iba(idcurr) = ibpart
            amass(idcurr) = sngl(ampart)
            jpa(idcurr) = jtrack
            
            if(JFLAG.eq.JRUN) then
               vxi(idcurr) = sngl(xtrack(0))
               vyi(idcurr) = sngl(ytrack(0))
               vzi(idcurr) = sngl(ztrack(0))
               vxf(idcurr) = sngl(xtrack(1))
               vyf(idcurr) = sngl(ytrack(1))
               vzf(idcurr) = sngl(ztrack(1))
            elseif(JFLAG.eq.JEND) then
               vxi(idcurr) = sngl(xx)
               vyi(idcurr) = sngl(yy) 
               vzi(idcurr) = sngl(zz)
               vxf(idcurr) = sngl(xx)
               vyf(idcurr) = sngl(yy)
               vzf(idcurr) = sngl(zz)
C               idead(idcurr) = icode
            endif
            px(idcurr) = sngl(ptrack*cxtrck)
            py(idcurr) = sngl(ptrack*cytrck)
            pz(idcurr) = sngl(ptrack*cztrck)
            pxf(idcurr) = px(idcurr)
            pyf(idcurr) = py(idcurr)
            pzf(idcurr) = pz(idcurr)
            trlen(idcurr) = sngl(cmtrck)
c
c particle latching: father, igen
c
            do ii = 1,numint
               if(idbflg.gt.0) then
                  write(*,*)'int= ',ii, 'x int= ',xint(ii),yint(ii),
     &                 zint(ii),' intpa= ',intpa(ii)
               endif
               if( (xx.eq.xint(ii)) .and. (yy.eq.yint(ii)) .and. 
     &              (zz.eq.zint(ii)) ) then
                  if(idgflg.gt.0) write(*,*)'father= ',ii
                  idpa(idcurr) = intpa(ii)
                  igen(idcurr) = igen(intpa(ii)) + 1
               endif
            end do
            
         endif  ! end of normal tracking
c
      endif
c
      return
      end

