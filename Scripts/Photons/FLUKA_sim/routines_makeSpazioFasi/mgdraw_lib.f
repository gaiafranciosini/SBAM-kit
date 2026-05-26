c----------------------------------------------------
      subroutine dump_common_part() 
c-----------------------------------------------------
c
      integer ii
      include 'mgdraw.inc'
      double precision ptot
c
      write(*,*)' '
      write(*,*)'------- dump_common_part ---------'
      write(*,*)' '
c
      write(*,*)'DUMP common PART'
      write(*,*)'numero particelle = ',nump
      do ii = 1, nump
         write(*,*)'particle = ', ii
         write(*,*)'idparent = ', idpa(ii),' igen = ',igen(ii), 
     &        ' charge = ',icha(ii),' reg =',numreg(ii)
         write(*,*)'idead = ',idead(ii),' jpa =',jpa(ii),
     &        ' mass = ',amass(ii),' time = ', tempo(ii)
         write(*,*)'vert in = ',vxi(ii),vyi(ii),vzi(ii)
         write(*,*)'vert fin = ',vxf(ii),vyf(ii),vzf(ii)
         write(*,*)'p ini = ',px(ii),py(ii),pz(ii)
         ptot = sqrt(px(ii)**2+py(ii)**2+pz(ii)**2)
         write(*,*)'p tot ini = ',ptot,' ekin ini= ',sqrt(ptot*ptot +
     &        amass(ii)*amass(ii))-amass(ii)
         write(*,*)'p fin = ',pxf(ii),pyf(ii),pzf(ii)
         write(*,*)'p tot fin = ',sqrt(pxf(ii)**2+pyf(ii)**2+pzf(ii)**2)
      end do
      return 
      end
      


c----------------------------------------------------
      subroutine dump_common_cross() 
c-----------------------------------------------------
c
      integer ii
      include 'mgdraw.inc'
c
      write(*,*)' '
      write(*,*)'------- dump_common_cross ---------'
      write(*,*)' '
c
      write(*,*)'numero crossing: ',ncross
      do ii = 1,ncross
         write(*,*)'crossing = ',ii,' part id = ',idcross(ii)
         write(*,*)'x,y,z = ',xcross(ii),ycross(ii),zcross(ii)
         write(*,*)'px,py,pz  in = ',pxcross(ii),pycross(ii),
     &        pzcross(ii)
         write(*,*)'reg = ',nregcross(ii),' m = ',amacross(ii),
     &        ' time = ',tcross(ii),' cha= ',chcross(ii)
         write(*,*)
      end do      
c     
      return
      end


c----------------------------------------------------
      subroutine dump_common() 
c-----------------------------------------------------
c
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      INCLUDE '(TRACKR)'
      include 'mgdraw.inc'
c
      if(idbflg .gt. 0) then
         call dump_common_part()
c
         call dump_common_cross()
      endif
c
      return
      end

            
c-------------------------------------------------------------------------
      SUBROUTINE score_cross(
     &        icharge,numbar,ampart,newreg,mreg,xsco,ysco,zsco)
c--------------------------------------------------------------------------
c
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      INCLUDE '(TRACKR)'
      include 'mgdraw.inc'
c
      integer  icharge, numbar, newreg, mreg 
      double precision ampart, xsco, ysco, zsco
c
c      write(*,*)'sono in score_cross'
      if(ncross.ge.maxcross) then
         write(*,*)"SCORE_CROSS: max number of crossing exceeded"
         return
      endif
      ncross = ncross + 1
      idcross(ncross) = idcurr
      xcross(ncross)  = sngl(xsco)
      ycross(ncross)  = sngl(ysco)
      zcross(ncross)  = sngl(zsco)
      pxcross(ncross) = sngl(ptrack*cxtrck)
      pycross(ncross) = sngl(ptrack*cytrck)
      pzcross(ncross) = sngl(ptrack*cztrck)
      tcross(ncross)  = sngl(atrack)
      chcross(ncross) = icharge
      amacross(ncross) = sngl(ampart)
      nregcross(ncross) = newreg
      nregoldcross(ncross) = mreg
c
      if(idbflg.gt.1) then
         write(*,*)' '
         write(*,*)'--------------- Score_cross -----------------'
         write(*,*)'idcurr = ',idcurr,' ncross= ',ncross
         write(*,*)'reg= ',newreg,' pxy,z= ',pxcross(ncross),
     &        pycross(ncross),pzcross(ncross),' mass= ',ampart
         write(*,*)'x,y,zcross= ',xsco,ysco,zsco,' t= ',atrack,
     &        ' cha= ',icharge
      endif             
c
c      write(*,*)'esco da score_cross'
      return
      end


