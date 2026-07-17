c
c     idcurr = position of the current particle in the particle common block
c           if idcurr<=0  -> new particle to be followed
c
c
      integer idbflg, outunit
      integer ifragflag
      integer numev
      integer idcurr
      integer JSTART, JRUN, JEND
      parameter (JSTART=1, JRUN=0, JEND=-1)
      common /comunico/idbflg, numev, idcurr, outunit, ifragflag

c *************************************************************
c
c particle common block
c
c   nump     = number of tracks
c 
c  for each of the produced track we register:
c
c   idpa     = index in the part common of the track parent
c   icha     = charge
c   iba      = barionic number
c   idead    = with which interaction the particle dies
c   jpa      = kind of particle
c   igen     = kind of interaction
c   vxi, vyi, vzi = production position of the track
c   vxf, vyf, vzf = final position of the track
c   px,py,pz, = production momentum of the track
c   amass = track mass
c   tempo = production time of the track
c   tof   = tof of the particle
c   trlen = track lenght of the particle
c
      integer nump, maxnump
      parameter(maxnump=1000)

      integer idpa(maxnump), igen(maxnump)
      integer icha(maxnump), numreg(maxnump), iba(maxnump)
      integer idead(maxnump), jpa(maxnump)
      real vxi(maxnump),vyi(maxnump),vzi(maxnump)
      real vxf(maxnump),vyf(maxnump),vzf(maxnump)
      real px(maxnump),py(maxnump),pz(maxnump)
      real pxf(maxnump),pyf(maxnump),pzf(maxnump)
      real amass(maxnump), tempo(maxnump), tof(maxnump)
      real trlen(maxnump)
      common /particle_common/ vxi, vyi, vzi, 
     &  vxf, vyf, vzf,px, py, pz, pxf, pyf, pzf,
     &  amass, tempo, tof, trlen, nump, idpa, igen, 
     &  icha, numreg, iba, idead, jpa
c *************************************************************
c common needed to score the boundary crossing info 
c and put the right pointer to the particle common of the produced electrons
c
      integer ncross, maxcross
      parameter (maxcross=10000)
      integer idcross(maxcross), nregcross(maxcross) 
      integer nregoldcross(maxcross), nlattcross(maxcross)
      integer nlattoldcross(maxcross)
      real pxcross(maxcross), pycross(maxcross), pzcross(maxcross)
      real xcross(maxcross), ycross(maxcross), zcross(maxcross)
      real tcross(maxcross), amacross(maxcross), chcross(maxcross)
      common /crssing_common/ 
     & pxcross, pycross, pzcross, 
     & xcross, ycross, zcross, tcross, chcross, amacross, ncross,  
     & idcross, nregcross, nregoldcross, nlattcross, nlattoldcross
c
c Numbers of fluka region to be scored
c
      integer nregvuoto,nregfanto
      common/region_common/nregvuoto,nregfanto
c
c
c  interaction stack: common block for the parent particle
c  recording with the info on the interaction point
c
      integer idfluka(maxnump)
      integer numint, maxint
      parameter (MAXINT=2000)
      double precision xint(maxint),yint(maxint),zint(maxint)
      integer intpa(maxint),intcode(maxint)
      common /particle_stack/xint,yint,zint,intpa,intcode,
     &       idfluka,numint
