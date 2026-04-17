c   Modified in Barcelona on December 2023
c   Modified by Dan, Volodymyr and Angel to run event-by-event simulations. 
c   kslab parameter must be set to 8 (non-symmetric inItial sate)
c   tnc parameter is now read from inhydro file: line 5, colunm 4
c   krandm parameter is now read from inhydro file: line 5, colunm 3
      program hydro
      implicit real*8 (a-h,o-z)
      parameter (maxfsiz=100)
      character infile*(maxfsiz), outfile*(maxfsiz), instat*(maxfsiz),
     1 instaf
      character kstrng, iss, isu, isr, ibl
c
      common /files/ outfile, infile, linfile, loutfile
      common /files1/ instaf(100)
      common /myl/ kstrng(100)
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
      common /pscmz/ zz
      common /countnevent/ nevent
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
      common /thestf/ rap00, Tfin, dxv, dyv, dzv, epsilon00, b0, rho00,
     1 e00
      common /gss/ sigx, sigy, sigz, radx, rady, radz, kspher
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
      common /other/ eusual, epn0, t0, epion, ppion, edelt, delta, pdelt
     1 , temp, tempkn, eps, xth, epsc, dnsmax, ncycpt, ncycrp, ncycaf
      common /numcmm/ numbk, numcel, numcst, numext, numfz1, numlv,
     1 numout, numouu, numou2, jdim1, kfreez, kfrzcy, press, freez,
     2 compmy, koast, kext, cpntot, compmx, compe,iundee(58)
      common /newcrp/ entot, dentot, tempav, numcll
      common /idnt/ ident(8)
      common /pscm2/ cycen, istop, isume
      common /corneliuss/ totenercornelius,totpxcornelius,
     1 totpycornelius,totpzcornelius,totbarcornelius

      dimension ibl(80), iss(20), isr(12), isu(13)
      data isr /'I','N','P','U','T',' ','f','i','l','e',':',' '/
      data iss /'I','N','I','T','I','A','L','_','S','T','A','T','E',' ',
     1'f','i','l','e',':',' '/
      data isu /'O','U','T','P','U','T',' ','f','i','l','e',':',' '/
c
c  icysum is the total number of cylcles
c  bcoef  is the coef. of I (thermal energy) in the eq. of state.
c  eta    defines freeze-out criteria
c  keqost specifies the equation of state  -x := QGP
c  comprs is the nuclear compressibility coef.
c  ______
c

      totenercornelius=0.
      totpxcornelius=0.
      totpycornelius=0.
      totpzcornelius=0.
      totbarcornelius=0.
      nevent=-1
      kresrt=0
      do 10 i=1,80
      ibl(i)=' '
   10 continue
c                  READING FROM INHYDRO FILE
      infile='../input/inhydro'
      linfile=mylenc(infile)
      open (155,file=infile(1:linfile),iostat=ios,status='old',err=190)
      rewind 155
      if (ios.ne.0) go to 190
      innf=1
      read (155,*,err=290) ioptos1,ioptos2,ioptos3,ioptos4,kdebug,kslab
      ioptosq=kslab
      hsizez=0.
      call opnfil
      innf=2
      read (155,500,err=290) id
      innf=3
      read (155,370,err=300) outfile
      loutfile=mylenc(outfile)
      open (156,file=outfile(1:loutfile),status='unknown',err=200)
c       write (156,*) ioptos1,ioptos2,ioptos3,ioptos4,kdebug,kslab
c       write (156,320) id
c       write (156,530)
      linfile=mylenc(infile)
      lg=linfile+12
      mg=(82-lg)/2
      mg1=82-lg-mg
c       write (156,350) (ibl(i),i=1,mg),isr,(kstrng(i),i=1,linfile),(
c      1 ibl(i),i=1,mg1)
      loutfile=mylenc(outfile)
      lg=loutfile+13
      mg=(82-lg)/2
      mg1=82-lg-mg
c       write (156,350) (ibl(i),i=1,mg),isu,(kstrng(i),i=1,loutfile),(
c      1 ibl(i),i=1,mg1)
      read (155,370,err=310) instat
      lstfile=mylenc(instat)
      lg=lstfile+20
      mg=(82-lg)/2
      mg1=82-lg-mg
      if (kslab.eq.2.or.kslab.eq.6) then
c         write (156,350) (ibl(i),i=1,mg),iss,(kstrng(i),i=1,lstfile),(
c      1  ibl(i),i=1,mg1)
      endif
      innf=6
c       write (156,380)
c
c    ioptos sets the printing flags, a 0 -> no output, a 1 -> output
c    ioptosq-> if non-zero, input from unit 2: Magas' code
c    kw -> if 0, plots to unit 6, otherwise to 33
c
c        tapex -> position of particles for every 10 cycles
c        tape2 -> (4, conditions as output from Magas' code
c                              _____________________________INPUT
c
c   Note that kpmax is three times the number of markers because of
c    three dimensions
      nq=3
      lcmb=91627233
      ijkmax=lcmb/nq
      if (ijkmax.gt.30542452) stop 248
      kpmax=85877055
      kpm=kpmax/3
      call dmill (ijkmax,ij0,ij1a,ij2)
      call dmill (kpmax,kp0,kp1,kp2)
c       write (156,360) ij0,ij1a,ij2,kp0,kp1,kp2
c
c    kdedug.NE.0 write to fort.11 cells of the type 1,
c                to fort.12 cells of types 2,3, if =0 - don't write
c
      pi=3.1415926535897932384626
      rho0=.145299917857864
      r0=(3/(4*pi*rho0))**(1./3.)
      convpu=((3./(4.*pi))**.3333333333333333)/1.18
c   convpu converts from fermi to program units
c
      kw=33
c                       ID of the run
c  10-> move: all cells
c  11-> move: type 1
c  12-> move: types 2 & 3
c
      if (ioptos4.ne.0) write (10,420) id
      if (kdebug.ne.0) write (11,420) id
      if (kdebug.ne.0) write (12,420) id
      if (kdebug.ne.0) write (13,420) id
      if (kdebug.ne.0) write (14,420) id
      time=0.0001
      innf=5
c
c                 READING LINE 5 FROM INHYDRO FILE
c
      if (kresrt.eq.0) read (155,*,err=180) keqost,icysum,krandm,tnc
      write(*,*) 'krandm= ',krandm,'tnc=',tnc
      if (icysum.gt.0) go to 40
c       write (156,*) ' ICYSUM < 0 - natural stop '
      stop 2
   40 bcoef=1./3.
      innf=6
      interp=1
      interv=1
      nonsym=0
c                 READING LINE 6 FROM INHYDRO FILE
      if (kresrt.eq.0) read (155,*,err=220) jdim1
c
c      kslab: if 0, ellipsoids - default
c             if 1, slabs
c             if 2, uniform streaks from Vladimir
c             if 3, uniform streaks, same E, N
c             if 4, Landau - removed
c             if 5, Spheroidal
c             if 6, torus
c             if 7, arbitrary - file must be read in via inarbt
c             if 8, non-symmetric initia state
c
      if (kslab .eq. 8) go to 120
      write(*,*) 'error: kslab must be set to 8'
      stop 189
c      go to (70,50,50,160,90,100,110,120), kslab
  120 continue
      nonsym=1
      ap=197.
      at=197.
      b0=0.5
      rho00=0.145
      epsilon00=100.*1000.  
      call cnucl
      go to 140
  130 continue
  140 continue
      call read_chiral_eos
      call picrel
  150 continue
c
c	to continue coombined calculations
c	istop=1
c       write (156,*) ' at 80'
      write (*,*) ' at 80 natural stop '
      write (*,*) 'ncycaf=',ncycaf
 
      return
c      stop 1
  160 continue
c       write (156,*) 'Error: kslab=',kslab,' option not implemented'
      write (*,*) 'Error: kslab=',kslab,' option not implemented'
      stop 4
  170 continue
c       write (156,*) ' error in read, ap=',ap,at,gevp,gevt
      write (*,*) ' error in read, ap=',ap,at,gevp,gevt
      stop 5
  180 continue
c       write (156,*) 'error in reading line5 keqost',keqost,icysum,
c      1 'ioptosq=',ioptosq
      write (*,*) 'error in reading line 5 keqost',keqost,icysum,'ioptos
     1=',ioptosq
      write (*,340)
      write (*,500) id
      write (*,*) 'Original input file was: ',infile
      write (*,*) 'Original output file was: ',outfile
      write (*,500) id
      stop 6
  190 continue
c       write (156,*) ' MISSING INPUT FILE ',infile
      write (*,*) ' MISSING INPUT FILE ',infile
      write (*,340)
      write (*,500) id
      write (*,*) 'Original input file was: ',infile
      write (*,*) 'Original output file was: ',outfile
      write (*,500) id
      stop 469
  200 continue
c       write (156,*) ' Problem with output file',outfile
      write (*,*) ' Problem with output file',outfile
      write (*,340)
      write (*,500) id
      write (*,*) 'Original input file was: ',infile
      write (*,*) 'Original output file was: ',outfile
      write (*,500) id
      stop 21
  210 continue
      write (*,*) ' '
      write (*,*) ' ERROR'
c       write (156,*) ' MISSING  FILE for INITIAL STATE: ',instat
      write (*,*) ' MISSING  FILE for INITIAL STATE: ',instat
      write (*,340)
      write (*,500) id
      write (*,*) 'Original input file was: ',infile
      write (*,*) 'Original output file was: ',outfile
      stop 488
  220 write (*,*) 'Error in input line 6, jdim1,...'
      stop 5
  230 continue
      write (*,*) ' error in reading Tfin'
      stop 493
  240 continue
      write (*,*) ' error in reading second line of instat'
      stop 511
  250 write (*,*) ' kslab = 0 no longer supported'
c       write (156,*) 'kslab = 0 no longer supported'
      stop 498
  260 write (*,*) 'kslab too large',kslab
c       write (156,*) 'kslab too large',kslab
      stop 502
  270 write (*,*) 'Error in input, line with npart,...'
      stop 504
  280 write (*,*) ' Error in opening file Restart_old'
      stop 506
  290 write (*,*) 'Error reading line 3',innf,' of input'
      stop 523
  300 write (*,*) 'Error reading line 4',innf,' of input'
      stop 525
  310 write (*,*) 'Error reading line 5',innf,' of input'
      stop 527
c

  320 format ('|',82x,'|',/,'|',1x,20a4,1x,'|')
  330 format ('|',29x,'Begin reading input file',29x,'|')
  340 format (//,'input already read:',/)
  350 format ('|',82a1,'|')
  360 format ('|',3x,'Maximum # of cells =',i2,',',i3,',',i3,'  Array si
     1ze =',i2,',',i3,',',i3,' = 3 times the # pts.'4x,'|',/,'|',82x,'|'
     2 )
  370 format (a)
  380 format ('|',82x,'|')
  390 format (/,20x,'________ RESTART _______',/)
  400 format (2i5,i3)
  410 format (3i3,2f8.3)
  420 format ('#',3x,'"Rel.HI-PIC-HYDRO" Code  vers.: hydro8.3  ','FREEZ
     1E OUT data'/'#',3x,' ID of this run: ',20a4)
  430 format (2x,'i_x ',2x,'l_y ',2x,'indic',1x,'z_min     ',2x,'z_max
     1   ',2x,'n [1/fm^3]  ',2x,'e [GeV/fm^3]',2x,'th(y)=v_z   ',2x,'  R
     2_plus    ',2x,'  R_minus   ',2x,' T_max [fm] ')
  440 format (7x,E9.3,5x,10a4)
  450 format (18x,E8.2,5x,e8.2,11x,i3,6x,i3,6x,i2)
  460 format (5x,e11.5,8x,e11.5,6x,e11.5)
  470 format (5x,e8.3,6x,e8.3,6x,e8.3)
  480 format (21x,i4,2x,i4,5x,i4,5x,i4)
  490 format (4i1,4i2)
  500 format (20a4)
  510 format ('new dxv of mini cell=',f10.6,'fm  r_0=',f7.4,' R=',f6.4,'
     1fm  cells=',f5.0,'   ktnc=',i3)
  520 format (10a4)
  530 format ('|',82x,'|',/,'|',27('___'),'_|',/,'|',4x,' One fluid hydr
     1odynamic calculation of relativistic heavy ion collisions'6x,'|',/
     2 ,'|',34x,'Version 8f.4 ',35x,'|',/,'|',82x,'|')
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine xiter (i,j,k,x,y)
c
c  A subroutine that, for a given xi, iterates to calculate x from
c  which one can then obtain y
c
      implicit real*8 (a-h,o-z)

      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
      common /prstbl/ prs(99,102), trs(99,102), pres, temp, chem, amurs(
     1 99,102), entrs(99,102), entro, ipres
      common /rot/ a, b, c
      common /hxi33/ hxi3a, xi, hxi
c
      bsq=b*b
      hxi=eqost(xi)
      hxi3a=a*hxi
      hxi3a=hxi3a/xi
      numit=0
      if (c.gt.bsq) go to 50
   10 xnew=c*((1.+bcoef*x)/((1.+bcoef)*b+hxi3a*dsqrt(1.-x)))**2
      delx=abs(xnew-x)
      x=xnew
      numit=numit+1
      if (delx.lt.1.e-6) go to 40
      if (numit.eq.300) go to 20
      if (x.lt.1.0) go to 10
   20 x=1.0
      numit=0
   30 xnew=1.-(((1.+bcoef*x)*dsqrt(c/x)-(1.+bcoef)*b)/hxi3a)**2
      delx=abs(xnew-x)
      x=xnew
      numit=numit+1
      if (delx.lt.1.e-5) go to 40
      if (numit.eq.500) go to 50
      if (x.gt.0.) go to 30
      go to 50
c
c  here for all cases
c
   40 y=(xi*b/a*dsqrt(1.-x)-x*hxi)/(1.+bcoef*x)
      pt1=xi*b/a*dsqrt(1.-x)
      pt2=x*hxi
      pt3=pt1-pt2
      return
   50 continue
c       write (156,60) a,b,c,numit
      stop 60
c
   60 format (' no solution a,b,c,numit',3f12.7,i5)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine wtnwt (jdim1,npower)
c
c  npower determines the polynomial
c
      implicit real*8 (a-h,o-z)
c
c  routine to calculate n'th degree weighting for velocity transport
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
      common /scaled/ hsclz, dzmod
      common /transp/ xte, yte, zte, hpx, hmx, hpy, hmy, hpz, hmz, i, j,
     1 k, l, ueff, veff, weff, cneff, xmeff, ymeff, zmeff, eneff, koko
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233)
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c
      ktrak=0
      sumu=0.
      sumv=0.
      sumw=0.
      sume=0.
      sumvol=0.
      ijkb1=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c       if (capn(ijkb1).eq.0.) write (156,*) ' N zero wtnwt',ijkb1,i,j,k
      hx=(aint(xte)+.5-xte)
      hy=(aint(yte)+.5-yte)
      hz=(aint(zte)+.5-zte)
      if (npower.eq.2) go to 40
      if (npower.ne.1) stop 4321
      if (hx.gt.0.) go to 10
c x >1/2
      sxl=0.
      sx0=1.+hx
      sxr=-hx
c x <1/2
      go to 20
   10 sxl=hx
      sx0=1.-hx
      sxr=0.
   20 if (hz.gt.0.) go to 30
c z >1/2
c    Changes made 5 Jan 08; the scaling requires more thought when dzmod
c    is small
      szl=0.
      sz0=1.+hz
      szr=-hz
      go to 50
c z <1/2
   30 szl=hz
      sz0=1.-hz
      szr=0.
      zps=zte-aint(zte)
c    Changes made 5 Jan 08; the scaling requires more thought when dzmod
c    is small
      go to 50
   40 sxl=0.5*(0.5+hx)*(0.5+hx)
      sx0=0.75-hx*hx
      sxr=0.5*(0.5-hx)*(0.5-hx)
      sx=sxl+sx0+sxr
      syl=0.5*(0.5+hy)*(0.5+hy)
      sy0=0.75-hy*hy
      syr=0.5*(0.5-hy)*(0.5-hy)
      sy=syl+sy0+syr
      szl=0.5*(0.5+hz)*(0.5+hz)
      sz0=0.75-hz*hz
      szr=0.5*(0.5-hz)*(0.5-hz)
      sz=szl+sz0+szr
   50 continue
c   Look at the j-1, j, j+1 planes
      do 170 jj=1,3
      if (jdim1.eq.1.and.jj.ne.2) go to 170
      if (jj.eq.3.and.yte.lt.0.5) go to 170
      if (jj.eq.1.and.yte.gt.0.5) go to 170
      ijkb=ijkb1+(jj-2)*nqib
      if (yte.le.0.5.and.jj.eq.1) ijkb=ijkb1
      if (hy.gt.0.) go to 60
      syl=0.
      syr=-hy
      syz=1.+hy
      if (jj.eq.1) syz=syl
      if (jj.eq.3) syz=syr
      go to 70
   60 continue
c yte<.5
      syz=0.
      if (jj.eq.2) syz=1.-hy
      if (jj.eq.1) syz=hy
   70 continue
      if (jdim1.ne.1.and.j.eq.1.and.jj.eq.3) syz=syl+syr
c  evaluate quantities for the level j'=j
c  center cell- 1
      ijk=ijkb
      if (capn(ijk).eq.0.) go to 80
      sq=sx0*syz*sz0
      sumu=sumu+u(ijk)*sq
      sumv=sumv+v(ijk)*sq
      sumw=sumw+w(ijk)*sq
      sume=sume+e(ijk)*sq
      sumvol=sumvol+sx0*syz*sz0
c  nearest neighbors
c  (i-1,j,k) - 2
   80 ijk=ijkb-nq
      if (i.eq.1) go to 90
      if (ijk.le.0) go to 90
      if (capn(ijk).eq.0.) go to 90
      sq=sxl*syz*sz0
      sumu=sumu+u(ijk)*sq
      sumv=sumv+v(ijk)*sq
      sumw=sumw+w(ijk)*sq
      sume=sume+e(ijk)*sq
      sumvol=sumvol+sq
c  (i+1,j,k) - 3
   90 ijk=ijkb+nq
      if (i.eq.ibar) go to 100
      if (capn(ijk).eq.0.) go to 100
      sq=sxr*syz*sz0
      sumu=sumu+u(ijk)*sq
      sumv=sumv+v(ijk)*sq
      sumw=sumw+w(ijk)*sq
      sume=sume+e(ijk)*sq
      sumvol=sumvol+sq
c  (i,j,k-1)
 100   ijk=ijkb-nqijb
      if (k.eq.1) go to 110
      if (ijk.le.0) go to 110
      if (capn(ijk).eq.0.) go to 110
      sq=sx0*syz*szl
      sumu=sumu+u(ijk)*sq
      sumv=sumv+v(ijk)*sq
      sumw=sumw+w(ijk)*sq
      sume=sume+e(ijk)*sq
      sumvol=sumvol+sq
c  (i,j,k+1)
  110 ijk=ijkb+nqijb
      if (k.eq.kbar) go to 120
      if (capn(ijk).eq.0.) go to 120
      sq=sx0*syz*szr
      sumu=sumu+u(ijk)*sq
      sumv=sumv+v(ijk)*sq
      sumw=sumw+w(ijk)*sq
      sume=sume+e(ijk)*sq
      sumvol=sumvol+sq
c  corners
c  (i-1,j,k-1)
  120 ijk=ijkb-nq-nqijb
      if (i.eq.1) go to 130
      if (k.eq.1) go to 130
      if (ijk.le.0) write (*,*) 'PROBLEM in WTNWT',i,j,k,ijk,ijkb
      if (ijk.le.0) go to 130
      if (capn(ijk).eq.0.) go to 130
      sq=sxl*syz*szl
      sumu=sumu+u(ijk)*sq
      sumv=sumv+v(ijk)*sq
      sumw=sumw+w(ijk)*sq
      sume=sume+e(ijk)*sq
      sumvol=sumvol+sq
c  (i+1,j,k+1)
  130 ijk=ijkb+nq+nqijb
      if (i.eq.ibar) go to 150
      if (k.eq.kbar) go to 140
      if (capn(ijk).eq.0.) go to 140
      sq=sxr*syz*szr
      sumu=sumu+u(ijk)*sq
      sumv=sumv+v(ijk)*sq
      sumw=sumw+w(ijk)*sq
      sume=sume+e(ijk)*sq
      sumvol=sumvol+sq
c  (i+1,j,k-1)
  140 ijk=ijkb+nq-nqijb
      if (i.eq.ibar) go to 150
      if (ijk.le.0) go to 150
      if (k.eq.1) go to 150
      if (capn(ijk).eq.0.) go to 150
      sq=sxr*syz*szl
      sumu=sumu+u(ijk)*sq
      sumv=sumv+v(ijk)*sq
      sumw=sumw+w(ijk)*sq
      sume=sume+e(ijk)*sq
      sumvol=sumvol+sq
c  (i-1,j,k+1)
  150 ijk=ijkb-nq+nqijb
      if (i.eq.1) go to 160
      if (k.eq.kbar) go to 160
      if (capn(ijk).eq.0.) go to 160
      sq=sxl*syz*szr
      sumu=sumu+u(ijk)*sq
      sumv=sumv+v(ijk)*sq
      sumw=sumw+w(ijk)*sq
      sume=sume+e(ijk)*sq
      sumvol=sumvol+sq
  160 continue
      if (jj.eq.1.and.yte.lt.0.5) sumv=-sumv
  170 continue
      ueff=sumu/sumvol
      veff=sumv/sumvol
      weff=sumw/sumvol
      eneff=sume/sumvol
      ktrak=0
      return
      end
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tablt
      implicit real*8 (a-h,o-z)
c
      common /rtrt/ iqwe
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
      common /prstbl/ prs(99,102), trs(99,102), pres, temp, chem, amurs(
     1 99,102), entrs(99,102), entro, ipres
      common /rot/ a, b, c
      common /eqqcd/ tq, eq1, pq, sq, anq1, amuq
      common /transp/ xte, yte, zte, hpx, hmx, hpy, hmy, hpz, hmz, i, j,
     1 k, l, ueff, veff, weff, cneff, xmeff, ymeff, zmeff, eneff, koko
      common /cnuclr/ hvc, hv2, pi0, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
c     
      del=0.01*a
      iqwe=1
      do 10 i=0,100
      an=del*(i)
      if (i.eq.0) an=0.0001
      x=1.-an*an/(a*a)
      y=b-dsqrt(c*x) ! Energy density in the rest frame in MeV/fm3
c      z=0.1453*939.5731*y
c      az=0.1453*an
c      z=1000.*y
      z=y ! is already in Mev/fm3
      az=an

c Angel, 2024, not needed with new EoS
c
c      presx=y/3.-4.*bbag/3. !pressure in MeV/fm3
c      if (presx .lt. 0.) presx=0.
c      if (presx .gt. 0.) z=y
c      if (presx .eq. 0.) z=3.*bbag*((y+bbag)/(4.*bbag))**(4./3.)
      call eqsqcd3 (az,z)
      pijk=pq !/(939.5731*0.1453)
      pres=pq !/(939.5731*0.1453)
c       write (156,*) 'line 1762'
      temp=tq
      chem=amuq
      entro=sq
c     ADD ON FEBRUARY 19, 2024
      if (pijk .ge. z) then
      write (*,*) 'pijk>e (line640):', 'pijk= ',pijk,'eijk= ',z
c       pause
      endif



c       if (az.eq.0.) write (156,*) ' az=0-i',i,j,k,az,z
      fun=c-(b+pres)*(b+pres)*x
      gsq=1./(1.-x)
c       write (156,20) an,fun,a,b,c,pres,x,y,gsq,c/(gsq*gsq*x)
   10 continue
      iqwe=0
      return
c
   20 format ('an-',f12.7,e14.4,'  a,b,c',3f10.4,' p',3f11.6,3x,f11.7,2f
     1 16.6)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine strmomb (ijk,i,c)
      implicit real*8 (a-h,o-z)
c
      common /bccket/ lablbk(1527120), cn(1527120)
      data modnum /304139/
c
      k=mod(ijk,modnum)+1
      do 10 i=k,1527120
      if (lablbk(i).eq.ijk) go to 30
      if (lablbk(i).ne.-1) go to 10
      lablbk(i)=ijk
      cn(i)=c
      go to 40
   10 continue
      do 20 i=1,k
      if (lablbk(i).eq.ijk) go to 30
      if (lablbk(i).ne.-1) go to 20
      lablbk(i)=ijk
      cn(i)=c
      go to 40
   20 continue
c       write (156,50) ijk,c
      write (*,50) ijk,c
      stop 47
   30 continue
      cn(i)=cn(i)+c
   40 continue
      return
c
   50 format (' out of space in strmomb ',i9,5f10.6)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine strmoma (ijk,i,c)
c
      implicit real*8 (a-h,o-z)
c
      common /backet/ lablbk(22), cn(22), convol(22)
      data modnum /17/
c
      k=mod(ijk,modnum)+1
c       if (c.eq.0.) write (156,*) ' c is zero  help'
      do 10 i=k,21
      if (lablbk(i).eq.ijk) go to 30
      if (lablbk(i).ne.-1) go to 10
      lablbk(i)=ijk
      cn(i)=c
      go to 40
   10 continue
c
      do 20 i=1,k
      if (lablbk(i).eq.ijk) go to 30
      if (lablbk(i).ne.-1) go to 20
      lablbk(i)=ijk
      cn(i)=c
      go to 40
   20 continue
c       write (156,50) ijk,c
      write (*,50) ijk,c
      stop 46
   30 continue
      cn(i)=cn(i)+c
   40 continue
      return
c
   50 format (' out of space in strmom ',i9,5f10.6)
      end


c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine strmom (ijk,i,x,y,z,e,c)
c
      implicit real*8 (a-h,o-z)
c
      common /bucket/ lablbk(1527120), xm(1527120), ym(1527120), zm(1527120)
     1 , en(1527120), cn(1527120), cnzvl(1527120)
      data modnum /304139/
c
      k=mod(ijk,modnum)+1
      if (k.gt.1527120) go to 50
c       if (c.eq.0.) write (156,*) ' c is zero in strmom  help',ijk,x,y,
c      1 z,e
      do 10 i=k,1527120
      if (lablbk(i).eq.ijk) go to 30
      if (lablbk(i).ne.-1) go to 10
      lablbk(i)=ijk
      xm(i)=x
      ym(i)=y
      zm(i)=z
      en(i)=e
      cn(i)=c
      go to 40
   10 continue
      do 20 i=1,k
      if (lablbk(i).eq.ijk) go to 30
      if (lablbk(i).ne.-1) go to 20
      lablbk(i)=ijk
      xm(i)=x
      ym(i)=y
      zm(i)=z
      en(i)=e
      cn(i)=c
      go to 40
   20 continue
      write (*,*) ' problem in STRMOM',i,lablbk(i)
c       write (156,*) ' problem in STRMOM',i,lablbk(i)
      stop 46
   30 continue
      xm(i)=xm(i)+x
      ym(i)=ym(i)+y
      zm(i)=zm(i)+z
      en(i)=en(i)+e
      cn(i)=cn(i)+c
   40 continue
      return
   50 continue
c       write (156,*) 'Error in STRMOM  ijk=',ijk,k,x,y,z,e
      write (*,*) 'Error in STRMOM  ijk=',ijk,k,x,y,z,e
      stop 1456
c
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rtmj (x,f,fct,xli,xri,eps,iend,ier)
c
      implicit real*8 (a-h,o-z)
c
c  routine identical to rtmi, but name is changed
c  to prevent recursive call
c
c    prepare iteration
      ier=0
      xl=xli
      xr=xri
      x=xl
      tol=x
      f=fct(tol)
      if (f.eq.0.) go to 70
      fl=f
      x=xr
      tol=x
      f=fct(tol)
      if (f.eq.0.) go to 70
      fr=f
      if (sign(1.d+00,fl)+sign(1.d+00,fr).ne.0.) go to 120
c    basic assumption fl*fr less than 0 is satisfied.
c    generate tolerance for function values.
      i=0
      tolf=100.*eps
c    start iteration loop
   10 i=i+1
c    start bisection loop
      do 50 k=1,iend
      x=.5*(xl+xr)
      tol=x
      f=fct(tol)
      if (f.eq.0.) go to 70
      if (sign(1.d+00,f)+sign(1.d+00,fr).ne.0.) go to 20
c    interchange xl and xr in order to get the same sign in f and fr
      tol=xl
      xl=xr
      xr=tol
      tol=fl
      fl=fr
      fr=tol
   20 tol=f-fl
      a=f*tol
      a=a+a
      if (a-fr*(fr-fl).ge.0.) go to 30
      if (i-iend.le.0) go to 80
   30 xr=x
      fr=f
c    test on satisfactory accuracy in bisection loop
      tol=eps
      a=abs(xr)
      if (a-1.0.le.0.) go to 40
      tol=tol*a
   40 if (dabs(xr-xl)-tol.gt.0.) go to 50
      if (dabs(fr-fl)-tolf.le.0.) go to 60
   50 continue
c    end of bisection loop
c    no convergence after iend iteration steps followed by iend
c    successive steps of bisection or steadily increasing function
c    values at right bounds. error return.
      ier=1
   60 if (dabs(fr)-abs(fl).le.0.) go to 70
      x=xl
      f=fl
   70 return
c    computation of iterated x-value by inverse parabolic interpolation
   80 a=fr-f
      dx=(x-xl)*fl*(1.+f*(a-tol)/(a*(fr-fl)))/tol
      xm=x
      fm=f
      x=xl-dx
      tol=x
      f=fct(tol)
      if (f.eq.0.) go to 70
c    test on satisfactory accuracy in iteration loop
      tol=eps
      a=abs(x)
      if (a-1..le.0.) go to 90
      tol=tol*a
   90 if (dabs(dx)-tol.gt.0.) go to 100
      if (dabs(f)-tolf.le.0.) go to 70
c    preparation of next bisection loop
  100 if (sign(1.d+00,f)+sign(1.d+00,fl).ne.0.) go to 110
      xr=x
      fr=f
      go to 10
  110 xl=x
      fl=f
      xr=xm
      fr=fm
      go to 10
c    end of iteration loop
c    error return in case of wrong input data
  120 ier=2
      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rtmi (x,f,fct,xli,xri,eps,iend,ier)
c
      implicit real*8 (a-h,o-z)
c
c    prepare iteration
      ier=0
      xl=xli
      xr=xri
      x=xl
      tol=x
      f=fct(tol)
      if (f.eq.0.) go to 70
      fl=f
      x=xr
      tol=x
      f=fct(tol)
c      if (f) 20,160,20
      if (f.eq.0.) go to 70
      fr=f
      if (sign(1.d+00,fl)+sign(1.d+00,fr).ne.0.) go to 120
c    basic assumption fl*fr less than 0 is satisfied.
c    generate tolerance for function values.
      i=0
      tolf=100.*eps
c    start iteration loop
   10 i=i+1
c    start bisection loop
      do 50 k=1,iend
      x=.5*(xl+xr)
      tol=x
      f=fct(tol)
      if (f.eq.0.) go to 70
      if (sign(1.d+00,f)+sign(1.d+00,fr).ne.0.) go to 20
c    interchange xl and xr in order to get the same sign in f and fr
      tol=xl
      xl=xr
      xr=tol
      tol=fl
      fl=fr
      fr=tol
   20 tol=f-fl
      a=f*tol
      a=a+a
      if (a-fr*(fr-fl).ge.0.) go to 30
      if (i-iend.le.0) go to 80
   30 xr=x
      fr=f
c    test on satisfactory accuracy in bisection loop
      tol=eps
      a=abs(xr)
      if (a-1.0.le.0.) go to 40
      tol=tol*a
   40 if (dabs(xr-xl)-tol.gt.0.) go to 50
      if (dabs(fr-fl)-tolf.le.0.) go to 60
   50 continue
c    end of bisection loop
c
c    no convergence after iend iteration steps followed by iend
c    successive steps of bisection or steadily increasing function
c    values at right bounds. error return.
      ier=1
   60 if (dabs(fr)-abs(fl).le.0.) go to 70
      x=xl
      f=fl
   70 return
c    computation of iterated x-value by inverse parabolic interpolation
   80 a=fr-f
      dx=(x-xl)*fl*(1.+f*(a-tol)/(a*(fr-fl)))/tol
      xm=x
      fm=f
      x=xl-dx
      tol=x
      f=fct(tol)
c      if (f) 180,160,180
      if (f.eq.0.) go to 70
c
c    test on satisfactory accuracy in iteration loop
      tol=eps
      a=abs(x)
      if (a-1..le.0.) go to 90
      tol=tol*a
   90 if (dabs(dx)-tol.gt.0.) go to 100
      if (dabs(f)-tolf.le.0.) go to 70
c    preparation of next bisection loop
  100 if (sign(1.d+00,f)+sign(1.d+00,fl).ne.0.) go to 110
      xr=x
      fr=f
      go to 10
  110 xl=x
      fl=f
      xr=xm
      fr=fm
      go to 10
c    end of iteration loop
c    error return in case of wrong input data
  120 ier=2
      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function rootr (s)
      implicit real*8 (a-h,o-z)
      common /roots/ rmaj, rmin, x0, y0
      t=rmin*dsqrt(1.-s**2/rmaj**2)
      rootr=rmaj*dsqrt(rmaj**2-s**2)*(s-x0)-(rmin*s)*(t-y0)
      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function retrev (ijk)
      implicit real*8 (a-h,o-z)
c
      common /bucket/ lablbk(1527120), ensort(1527120), eza(1527120), ezb(1527120), ezc(1527120), ezd(1527120), cnzzz(1527120)
c
      data modnum /304139/
c
      k=mod(ijk,modnum)+1
      retrev=0.
      do 10 i=k,1527120
      if (ijk.eq.lablbk(i)) go to 30
      if (lablbk(i).eq.-1) go to 40
   10 continue
      do 20 i=1,k
      if (ijk.eq.lablbk(i)) go to 30
      if (lablbk(i).eq.-1) go to 40
   20 continue
c       write (156,50) ijk
      write (*,50) ijk
      stop 44
   30 retrev=ensort(i)
   40 continue
      return
   50 format (' failure in retrev ',i5)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine readinput
c   Modified February 2023 by Volodymyr & Angel: we've changed the units of
c   the energy density 
c
c   Subroutine to set up initial state for QGP calculation using input
c   Vladimir's code.
c       Updated on 26/02/2003 for the expanding streaks initial state
c     Feb 2005 - changed to read input from fort.100
      implicit real*8 (a-h,o-z)
      common /mssss/ outn, porton, npxz, npt, kout
      common /averps/ xaveps, yaveps, zaveps
      common /other/ eusual, epn0, t0, epion, ppion, edelt, delta, pdelt
     1 , temp, tempkn, eps, xth, epsc, dnsmax, ncycpt, ncycrp, ncycaf
      common /rot/ a, b, c
      common /plottt/ kbarx, kbar1, kfrst, imidpt, kmidpt, klmin
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
      common /nosym/ h1, h2, h3
      common /numcmm/ numbk, numcel, numcst, numext, numfz1, numlv,
     1 numout, numouu, numou2, jdim1, kfreez, kfrzcy, press, freez,
     2 compmy, koast, kext, cpntot, compmx, compe,iundee(58)
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
      common /matric/ zmin(300,300), zmax(300,300), energ(300,300), velz
     1 (300,300), dens(300,300), zsmz(200), adens(300,700), kmove(700),
     2 Rpl(300,300), Rmi(300,300), chyz(300,300), dummy(1)
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
      common /extrl/ zkmin(300,300), zkmax(300,300)
      Integer imid, kmid, indicator(300,300)
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
      common /thestf/ rap00, Tfin, dxv, dyv, dzv, epsilon00, b0, rho00,
     1 e00
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(9162
     2 7233), yncont(700), yecont(700), yzcont(700), xncont(700), xecont
     3 (700), xzcont(700), zncont(700), zecont(700), zzcont(700)
      dimension pressure(91627233)
      dimension px(85877055), py(85877055), pz(85877055), alpha(18)
      real*8 eijk, anijk, pressijk, wijk, uijk, vijk, zmomijk, xmomijk,
     1 ymomijk, pNtot, Etot
      character*256 hline1, hline2, hline3, hline4, hline5, hline6,
     1 hline7, hline8
      integer ios_hdr, ipos, ieq, ieq2
      logical hash_header
      dimension nypt(700)
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c     
      if (nredux.eq.0) nredux=1
c     nredux=1
      econg=.001*econ
      nreduy=nredux !nreduy=1
      nreduz=nredux ! nreduz=1
      npart=2
      npartz=3
      read (155,*) npart,npartz ! npart=3, npartz=3      
      if (npart.le.0) stop 99
      nparp=npart ! nparp=3
      nparpz=npartz ! nparpz=3
      nzero=0
      mwr=0
      if (kslab.eq.8) then
c   Define ibar, jbar, kbar
c       ijkmax -> # of fluid cells
        zmxx=ijkmax/3 ! ijkmax=30542411, zmxx=2036160
c   extraz is the ratio of ibar to kbar and allows more expansion in the
c     z direction.  Can change if needed.
        extraz=3.0
        ip=int((zmxx/extraz)**.3333333) 
        if (mod(ip,2).eq.1) ip=ip-1 ! ip=86
        ibar=ip ! ibar=86
        jbar=ibar ! jbar=86
        kbar=int(zmxx/float(ibar*jbar))
        if (mod(kbar,2).eq.1) kbar=kbar-1 ! kbar=274
        ijbar=ibar*jbar
        ijkbar=ijbar*kbar
        ijkba3=nq*ijkbar
        ijk=ibar*jbar*kbar
        nqib=nq*ibar
        nqijb=nq*ijbar
      endif
c       write (156,1200) ibar,jbar,kbar,ijk
      fibar=ibar
      fjbar=jbar
      fkbar=kbar
c       write (156,1190)
c       write (156,1210)
c       write (156,1190)
c       write (156,1050) nredux
c       write (156,1190)
c       if (krandm.eq.0) write (156,1060)
c       if (krandm.ne.0) write (156,1070)
c       write (156,1190)
c       write (156,1080) npart,npartz,nparp,nparpz
c       write (156,1190)
      if (nparpz.ne.npartz) nparpz=npartz
      open(unit=100,file='../input/fort.100')
      rewind 100
      read (100,'(A)',iostat=ios_hdr) hline1
      if (ios_hdr.ne.0) then
      write(*,*) 'Error reading header from ../input/fort.100'
      stop 9011
      endif
      ipos=1
 9110 if (ipos.le.256.and.hline1(ipos:ipos).eq.' ') then
      ipos=ipos+1
      go to 9110
      endif
      hash_header=.false.
      if (ipos.le.256) then
      if (hline1(ipos:ipos).eq.'#') hash_header=.true.
      endif
      if (.not.hash_header) then
      backspace 100
      read (100,1520) epsilon00
      read (100,1520) b0
      read (100,1520) ap
      read (100,1520) at
      read (100,1240) rap00
      read (100,1250) dxv
      read (100,1250) dyv
      read (100,1250) dzv
      else
      read (100,'(A)',iostat=ios_hdr) hline2
      read (100,'(A)',iostat=ios_hdr) hline3
      read (100,'(A)',iostat=ios_hdr) hline4
      read (100,'(A)',iostat=ios_hdr) hline5
      read (100,'(A)',iostat=ios_hdr) hline6
      read (100,'(A)',iostat=ios_hdr) hline7
      read (100,'(A)',iostat=ios_hdr) hline8
      if (ios_hdr.ne.0) then
      write(*,*) 'Error reading comment header from ../input/fort.100'
      stop 9012
      endif
      ieq=index(hline3,'=')
      if (ieq.eq.0) then
      write(*,*) 'Missing epsilon0 in fort.100 header'
      stop 9013
      endif
      read (hline3(ieq+1:),*,iostat=ios_hdr) epsilon00
      if (ios_hdr.ne.0) then
      write(*,*) 'Bad epsilon0 value in fort.100 header'
      stop 9014
      endif
      ieq=index(hline3,'b=')
      if (ieq.eq.0) then
      write(*,*) 'Missing b value in fort.100 header'
      stop 9015
      endif
      read (hline3(ieq+2:),*,iostat=ios_hdr) b0
      if (ios_hdr.ne.0) then
      write(*,*) 'Bad b value in fort.100 header'
      stop 9016
      endif
      ieq=index(hline3,'Ap, At')
      if (ieq.eq.0) ieq=index(hline3,'Ap,At')
      if (ieq.eq.0) then
      write(*,*) 'Missing Ap,At values in fort.100 header'
      stop 9017
      endif
      ieq2=index(hline3(ieq:),'=')
      if (ieq2.eq.0) then
      write(*,*) 'Malformed Ap,At field in fort.100 header'
      stop 9018
      endif
      read (hline3(ieq+ieq2:),*,iostat=ios_hdr) ap,at
      if (ios_hdr.ne.0) then
      write(*,*) 'Bad Ap,At values in fort.100 header'
      stop 9019
      endif
      ieq=index(hline4,'=')
      if (ieq.eq.0) then
      write(*,*) 'Missing beam y0 value in fort.100 header'
      stop 9020
      endif
      read (hline4(ieq+1:),*,iostat=ios_hdr) rap00
      if (ios_hdr.ne.0) then
      write(*,*) 'Bad beam y0 value in fort.100 header'
      stop 9021
      endif
      ieq=index(hline5,'dx=')
      if (ieq.eq.0) then
      write(*,*) 'Missing dx value in fort.100 header'
      stop 9022
      endif
      read (hline5(ieq+3:),*,iostat=ios_hdr) dxv
      if (ios_hdr.ne.0) then
      write(*,*) 'Bad dx value in fort.100 header'
      stop 9023
      endif
      ieq=index(hline5,'dy=')
      if (ieq.eq.0) then
      write(*,*) 'Missing dy value in fort.100 header'
      stop 9024
      endif
      read (hline5(ieq+3:),*,iostat=ios_hdr) dyv
      if (ios_hdr.ne.0) then
      write(*,*) 'Bad dy value in fort.100 header'
      stop 9025
      endif
      ieq=index(hline5,'dz=')
      if (ieq.eq.0) then
      write(*,*) 'Missing dz value in fort.100 header'
      stop 9026
      endif
      read (hline5(ieq+3:),*,iostat=ios_hdr) dzv
      if (ios_hdr.ne.0) then
      write(*,*) 'Bad dz value in fort.100 header'
      stop 9027
      endif
      endif

      nparp2=nparp**2 ! nparp2=9
c     To get correct cell sizes ???????
      dxdydv=dxv*dyv*dzv ! dxdydv=0.2596940720 fm3
      dx=dxv ! dx=0.638 fm
      dy=dyv ! dy=0.638 fm
      dz=dzv ! dz=0.638 fm
      dxdydz=dx*dy*dz ! dxdydz=0.2596940720 fm3
      rdxdydz=1./dxdydz ! rdxdydz=3.8506847395423027 fm-3
      dxdy=dxv*dyv
      dzvfm=2.*float(nreduz)*dzv !
      volxyz=dxdydz
      vlxyz=volxyz*float(nredux*nreduy)
c       write (156,1430) vlxyz,volxyz,dxdy,dxdydz
c       write (156,1190)
      rvolxyz=1./volxyz
      do 20 i=1,300
      do 10 j=1,700
      adens(i,j)=0.
   10 continue
   20 continue
      do 30 i=1,lcmb-2
      aa1(i)=0. ! In this vector is stored the three components of the 
c     flow vwlocity
      aa2(i)=0. ! In this vector is stored N0 (capn),pressure (p) and T00 (e)
c     N0 [fm-3], pressure [MeV/fm3], T00 [MeV/fm3]
      aa3(i)=0. ! In this vector is stored the three components of the 
c     momentum
   30 continue
      do 40 i=1,700
      yncont(i)=0.
      yecont(i)=0.
      yzcont(i)=0.
      nypt(i)=0.
   40 continue
      do 50 i=1,700
      xncont(i)=0.
      xecont(i)=0.
      xzcont(i)=0.
      zncont(i)=0.
      zecont(i)=0.
      zzcont(i)=0.
   50 continue
c       Just to fill in since in common
      do 70 i=1,300
      do 60 j=1,70
      dens(i,j)=0.
      energ(i,j)=0.
      indicator(i,j)=0
      velz(i,j)=0.
      chyz(i,j)=0.
      Rpl(i,j)=0.
      Rmi(i,j)=0.
      zmin(i,j)=0.
      zmax(i,j)=0.
   60 continue
   70 continue
      if (krandm.eq.0) go to 80 ! krandm is set to 1
      call random_seed
   80 continue
      imid=(ibar+1)/2 ! imid=43
      jmid=(jbar+1)/2 ! jmid=43
      kmid=(kbar+1)/2 ! kmid=137
      klmin=100000
      krmax=0
      ipmin=100000
      ipmax=-1
c      some part is taken away for reading of initial state
c        Instead of ESRM calc. we have now
      do 110 k=1,kbar
      do 100 j=1,jbar
      do 90 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      capn(ijk)=0.
      e(ijk)=0.
      p(ijk)=0.
      pressure(ijk)=0.
      u(ijk)=0.
      v(ijk)=0.
      w(ijk)=0.
      zmom(ijk)=0.
      xmom(ijk)=0.
      ymom(ijk)=0.
   90 continue
  100 continue
  110 continue
      i1min=9999
      j1min=9999
      k1min=9999
      i1max=0
      j1max=0
      k1max=0
      i7=1000000000
      zmomt=0
      umaxmx=0.
      umaxmy=0.
      umaxmz=0.
      umaxw=0.
      aan=0.
      een=0.
      zmomps=0.
      zmomng=0
      asmerr=0.
      presin=0.
      dpt=dxdydz/float(npart*npart*npartz) ! npart=npartz=3, dpt=9.6182989629629637E-003 (fm3)
      npoito=0
      anmax=0.
      egsum=0.
      agsum=0.
      i6=0
      egsum_prueba=0. ! used to check energy conservation
      agsum_prueba=0. ! used to check charge baryon  conservation
      sum_zmomijk=0.
      sum_press=0.
      
      sum_wijk=0.
c                      READING DATA FROM FORT.100 FILE   
  120 read (100,*,end=130) i1,j1,k1,eijk,anijk,pressijk,wijk,uijk,vijk
     1 ,zmomijk,xmomijk,ymomijk 
c     eijk=T00 in GeV/fm3
c     anijk=N0 in fm-3
c     pressijk=pressure in GeV/fm3
c     wijk,uijk,vijk are given in units c
c     zmomijk,xmomijk,ymomijk are given in GeV/fm3  
c     To get MeV
      eijk=eijk*1000. ! T00 [MeV/fm3]
      pressijk=pressijk*1000. ! pressure [MeV/fm3]
      zmomijk=zmomijk*1000. ! T0Z [MeV/fm3]
      xmomijk=xmomijk*1000. ! T0X [MeV/fm3]
      ymomijk=ymomijk*1000. ! T0Y [MeV/fm3]
      ijk=nq*(ijbar*(k1-1)+ibar*(j1-1)+i1-1)+1
      if (uijk.ne.0..or.vijk.ne.0..or.xmomijk.ne.0..or.ymomijk.ne.0.) go
     1 to 890
      if (anijk.eq.0.) anijk=.00004
      egsum=egsum+eijk
      agsum=agsum+anijk
      ansum=ansum+anijkreadinput
      egsum_prueba=egsum_prueba+eijk
      agsum_prueba=agsum_prueba+anijk
      sum_zmomijk=sum_zmomijk+zmomijk
      sum_press=sum_press+pressijk
      sum_wijk=sum_wijk+abs(wijk)
      anijk=anijk/0.145 ! N0/n0
      npoin=int(dxdydz*anijk/dpt+.5000001)
      erras=anijk-float(npoin)*dpt/dxdydz
      asmerr=asmerr+erras
      anijk=anijk*0.145 ! N0 [fm-3]
      presin=presin+pressijk ! MeV/fm3
      i6=i6+1
      e(ijk)=eijk ! T00 given in MeV/fm3
      p(ijk)=anijk ! baryon density !!!!!
c   Note that the anijk is stored temporarily in array p
c   Pressure is stored in the array pressure, but later moved to p(ijk)
c   It could be stored in u() if there is no transverse initial velocity
      if (anijk.gt.anmax) anmax=anijk
      if (abs(xmomijk).gt.umaxmx) umaxmx=abs(xmomijk)
      if (abs(ymomijk).gt.umaxmy) umaxmy=abs(ymomijk)
      if (abs(zmomijk).gt.umaxmz) umaxmz=abs(zmomijk)
      if (abs(wijk).gt.umaxw) umaxw=abs(wijk)
      aan=aan+anijk
      een=een+eijk
      npoito=npoito+npoin
      pressure(ijk)=pressijk
      zmomt=zmomt+zmomijk
      u(ijk)=uijk
      v(ijk)=vijk
      w(ijk)=wijk ! z-component of flow velocity in units c
      if (u(ijk)**2+v(ijk)**2+w(ijk)**2 .gt. 1.) then
      write (*,*) 'v^2>1 in line 5554',i,j,k,u(ijk),v(ijk),w(ijk)
      wijk=0.
      uijk=0.
      vijk=0.
      w(ijk)=0.
      u(ijk)=0.
      v(ijk)=0.
c       pause 1314
      endif
c     ------------------------------------------------------------------
c                Only used in this subroutine
      xncont(i1)=xncont(i1)+anijk*dxdydz ! baryon charge
      xecont(i1)=xecont(i1)+eijk*dxdydz ! T00 [MeV]
      xzcont(i1)=xzcont(i1)+zmomijk ! T0Z [MeV/fm3]
      yncont(j1)=yncont(j1)+anijk*dxdydz ! baryon charge
      yecont(j1)=yecont(j1)+eijk*dxdydz ! T00 [MeV]
      yzcont(j1)=yzcont(j1)+zmomijk ! T0Z [MeV/fm3]
      zncont(k1)=zncont(k1)+anijk*dxdydz ! baryon charge
      zecont(k1)=zecont(k1)+eijk*dxdydz ! T00 [MeV]
      zzcont(k1)=zzcont(k1)+zmomijk ! T0Z [MeV/fm3]
c     ------------------------------------------------------------------
      zmom(ijk)=zmomijk
      xmom(ijk)=xmomijk
      ymom(ijk)=ymomijk
      if (zmom(ijk).gt.0.) zmomps=zmomps+zmom(ijk)
      if (zmom(ijk).lt.0.) zmomng=zmomng+zmom(ijk)
      if (i1min.gt.i1) i1min=i1
      if (i1max.lt.i1) i1max=i1
      if (j1min.gt.j1) j1min=j1
      if (j1max.lt.j1) j1max=j1
      if (k1min.gt.k1) k1min=k1
      if (k1max.lt.k1) k1max=k1
      if (k1.lt.klmin) klmin=k1
      if (k1.gt.krmax) krmax=k1
      if (ipmin.gt.i1) ipmin=i1
      if (ipmax.lt.i1) ipmax=i1
      go to 120
  130 continue
c                    FINISH READING FORT.100 FILE 
c       write (*,*) 'CHECKING ENERGY, MOMENTO AND BARYON CHARGE READ FROM'
c       write (*,*) '                 FORT.100 FILE'
c       write (*,*) 'T00 [MeV*fm-3], T0Z [MeV*fm-3], T0Y [MeV*fm-3]'
c       write (*,*) 'T0X [MeV*fm-3],N0 [fm-3], press [MeV*fm-3]'
c       write (*,*) 'sumT00*V=',egsum_prueba*dxdydz,'sumT0Z*V',
c      1 sum_zmomijk*dxdydz,'sumpress=',sum_press*dxdydz,'sumN0*V=',
c      2 agsum_prueba*dxdydz,'sum(abs(wijk))',sum_wijk
c       write (156,1190)
      numps=kpmax/3 ! kpmax=85877055, numps=28625685
      numdes=porton*dfloat(numps) ! porton=0.9, numdes=5152623
      dpt=agsum*dxdydz/dfloat(numdes) ! Now dpt is dimensionless dpt=3.3330075366957684E-005
c       write (156,1260) i6
c       write (156,940) anmax
c       write (156,1190)
c       write (156,1090) umaxmx,umaxmy,umaxmz,umaxw
c       write (156,960) zmomng*dxdydz,zmomps*dxdydz,(zmomng+zmomps)*dxdydz
c       write (156,950) aan*dxdydz,een*dxdydz,een/aan
c       write (156,1190)
      i1mid=(i1max+i1min+1)/2
      j1mid=(j1max+j1min+1)/2
      k1mid=(k1max+k1min+1)/2
      kwidth=k1max-k1min
c       write (156,1270) i1min,i1mid,i1max,j1min,j1mid,j1max,k1min,k1mid
c      1 ,k1max,kwidth
      iadjst=imid-i1mid
      jadjst=jmid-j1mid
      kadjst=kmid-k1mid
c       write (156,1280) iadjst,jadjst,kadjst
c       write (156,1190)
      write(*,*) k1min,k1max, k1mid,kmin,kmax,kmid,kadjst,kwidth/2,jbar
      if (k1mid+kadjst+kwidth/2+2.ge.k1min) go to 920
c     I am not sure that xmom(ijk), ymom(ijk) will have any effect
c     even if these are given non-zero
c       Coming back to the units used in the code
c   Change units and center initial state
c


c     SHITTING INITIAL STATE GRID TO COINCIDE WITH THE CENTER OF THE HYDRO GRID
      do 160 i1=i1min,ibar
      do 150 j1=j1min,jbar
      do 140 k1=k1min,kbar
      i=i1+iadjst
      j=j1+jadjst
      k=k1+kadjst
      ijk1=nq*(ijbar*(k1-1)+ibar*(j1-1)+i1-1)+1
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (ijk1.le.0) write (*,*) i,j,k,i1,j1,k1
      if (ijk1.le.0) stop 392
      if (ijk.le.0) write (*,*) i,j,k,i1,j1,k1,iadjst,jadjst,kadjst
      if (ijk.le.0) stop 4132
      if (ijk.gt.ijkmax) go to 140
      if (ijk1.gt.ijkmax) go to 140
      if (p(ijk1).eq.0.) go to 140
      if (p(ijk).ne.0.) go to 930
      e(ijk)=e(ijk1) ! T00 in MeV/fm3
      p(ijk)=p(ijk1) ! N0 (baryon density)
      u(ijk)=u(ijk1) ! ADDED ON APRIL 30, 2023
      v(ijk)=v(ijk1) ! ADDED ON APRIL 30, 2023
      w(ijk)=w(ijk1)
      pressure(ijk)=pressure(ijk1) ! in MeV/fm3
      zmom(ijk)=zmom(ijk1) ! in MeV/fm3
      xmom(ijk)=xmom(ijk1) ! in MeV/fm3
      ymom(ijk)=ymom(ijk1) ! in MeV/fm3
      if (xmom(ijk).ne.0.) write (*,*) 'xmom',i,j,k,xmom(ijk)
      if (ymom(ijk).ne.0.) write (*,*) 'ymom',i,j,k,ymom(ijk)
      p(ijk1)=0.
      e(ijk1)=0.
      zmom(ijk1)=0.
      u(ijk1)=0.
      v(ijk1)=0.
      w(ijk1)=0.
  140 continue
  150 continue
  160 continue


      cmax=0.
      cmin=1110.
      klmin=klmin+1
      krmax=krmax+1
      ptimid=(ipmin+ipmax+1)/2
      ptkmid=(klmin+krmax+1)/2
      xmid=ptimid
      zmid=ptkmid
      kmidpt=int(ptkmid+.1)
      imidpt=int(ptimid)
      imidii=(ipmin+ipmax)/2
      asum=0.
      asumv=0.
      esum=0.
      pnu(1)=dpt ! dpt is defined in line 4067 and is dimensionless
      pnu(2)=dpt
      koun=0
      pmax=0.
      rewind 8
      i2min=i1min+iadjst
      i2max=i1max+iadjst
      j2min=j1min+jadjst
      j2max=j1max+jadjst
      k2min=k1min+kadjst
      k2max=k1max+kadjst
      mchk=0
      do 190 j2=j2min,j2max
      do 180 i2=i2min,i2max
      do 170 k2=k2min,k2max
      ijk=nq*(ijbar*(k2-1)+ibar*(j2-1)+i2-1)+1
      if (p(ijk).eq.0..and.e(ijk).ne.0.) write (*,1300) i2,j2,k2,p(ijk)
     1 ,e(ijk)
      if (p(ijk).eq.0.) go to 170
  170 continue
  180 continue
  190 continue
c    
      axtot=0
      extot=0.
      aytot=0
      eytot=0.
      aztot=0
      eztot=0.
      zero=0.
      i1mean=i1max+i1min
      j1mean=j1max+j1min
      k1mean=k1max+k1min
      ireflc=i1max-i1min
      jreflc=j1max-j1min
      kreflc=k1max-k1min
      do 200 i1=1,700
      if (xecont(i1).eq.0.) go to 200
      ireflc=i1mean-i1
      axtot=axtot+xncont(i1) ! baryon charge
      extot=extot+xecont(i1) ! energy in MeV
      epern=0. ! energy per nucleon
      if (xncont(i1).ne.0.) epern=xecont(i1)/xncont(i1) ! MeV/nucl
      zpern=0. ! z-component of momentum per nucleon
      if (xncont(i1).ne.0.) zpern=xzcont(i1)/xncont(i1) ! MeV/fm3/nucl
  200 continue
      do 210 j1=1,700
      if (yecont(j1).eq.0.) go to 210
      jreflc=j1mean-j1
      epern=0.
      aytot=aytot+yncont(j1)
      eytot=eytot+yecont(j1)
      if (yncont(j1).ne.0.) epern=yecont(j1)/yncont(j1)
      zpern=0.
      if (yncont(j1).ne.0.) zpern=yzcont(j1)/yncont(j1)
  210 continue
      do 220 j1=1,700
      if (zecont(j1).eq.0.) go to 220
      kreflc=k1mean-j1
      epern=0.
      aztot=aztot+zncont(j1)
      eztot=eztot+zecont(j1)
      if (zncont(j1).ne.0.) epern=zecont(j1)/zncont(j1)
      zpern=0.
      if (zncont(j1).ne.0.) zpern=zzcont(j1)/zncont(j1)
      if (k1.eq.k1mean) write (*,*) ' '
  220 continue
c       write (156,1320) axtot,aytot,aztot,extot,eytot,eztot
c       write (156,1190)
      do 280 i=1,ibar
      do 270 j=1,jbar
      do 260 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (ijk.gt.ijkmax) stop 21
      if (p(ijk).eq.0.) go to 260
      if (cmax.lt.p(ijk)) cmax=p(ijk)
      if (cmin.gt.p(ijk)) cmin=p(ijk)
      asum=asum+p(ijk) ! baryon density
      esum=esum+e(ijk) ! Mev/fm3
      koun=koun+1
  260 continue
  270 continue
  280 continue
c
c       write (156,1190)
c       write (156,1450) volxyz,rdxdydz,vlxyz
c       write (156,1190)
c       write (156,1410) cmax,cmin,cmax/dpt,cmin/dpt
c       write (156,1190)
      i9min=9999
      j9min=9999
      k9min=9999
      i9max=0
      j9max=0
      k9max=0
      do 310 i=1,ibar
      do 300 j=1,jbar
      do 290 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (ijk.gt.ijkmax) go to 290
      if (p(ijk).eq.0.) go to 290
      if (i9min.gt.i) i9min=i
      if (i9max.lt.i) i9max=i
      if (j9min.gt.j) j9min=j
      if (j9max.lt.j) j9max=j
      if (k9min.gt.k) k9min=k
      if (k9max.lt.k) k9max=k
      if (kslab.ne.8) write (105,1380) i,j,k,e(ijk)*rdxdydz
     1 ,p(ijk)*rdxdydz/float(nredux*nreduz),e(ijk)/3.,w
     2 (ijk),0.e0,0.e0,zmom(ijk)*rdxdydz*1e3/939.5731,0.e0,0.e0
  290 continue
  300 continue
  310 continue
      i9mid=(i9max+i9min+1)/2
      j9mid=(j9max+j9min+1)/2
      k9mid=(k9max+k9min+1)/2
      k9widt=k9max-k9min
c       write (156,1340) i9min,i9mid,i9max,j9min,j9mid,j9max,k9min,k9mid
c      1 ,k9max,k9widt
c     each marker particle will carry dpt
c       Please vary if there are too many or too few particles
c       write (156,1190)
      if (kslab.lt.8) then
c         write (156,1480) dpt,cmax*volxyz/dpt,cmin*volxyz/dpt,asumv/dpt,
c      1  3.*asumv/dpt
        write (*,1320) dpt,cmax*volxyz/dpt,cmin*volxyz/dpt,asumv/dpt,3.
     1  *asumv/dpt
      endif
      if (kslab.eq.8) then
        asumv=asum
c         write (156,1480) dpt,cmax/dpt,cmin/dpt,asumv/dpt,3.*asumv/dpt
      endif
      pnu(1)=dpt
      pnu(2)=dpt
c       write (156,1190)
c       write (156,1510)
c       write (156,1500) gam(1),pnu(1),e0(1),zmom0(1),w0(1)
c       if (kslab.lt.8) write (156,1500) gam(2),pnu(2),e0(2),zmom0(2),
c      1 w0(2)
c       write (156,1190)
c
      npxz=0
      esumin=0.
      asumot=0.
      asumin=0.
      npcmin=0
      esumot=0.
      anerr=0
      summz=0.
      summzo=0.
      ko=-2
      mzmax=0
      jom=0
      iom=0
      kout=0
      parpi=1./float(nparp)
      parpzi=1./float(nparpz)
      ypmin=1./float(nparp)
      ypmin=0.45
      if (nparp.eq.1) parpi=0.
      if (nparpz.eq.1) parpzi=0.
c       if (nparp.gt.6) write (156,*) ' nparp > 6'
c       if (nparpz.gt.6) write (156,*) ' npaprz >6'
      if (nparp.gt.7) stop 1236
      if (nparpz.gt.7) stop 1238
c
c_______________________________________________________________________
c
c    Assign marker particles to cells
c
c_______________________________________________________________________
c
      ibar2=(ibar)/2
      kbar2=(kbar+3)/2
      ibarx=ibar
      if (b0.eq.0.) ibarx=ibar2
      write (*,*) 'b0=',b0
      if (kslab.eq.8) then
        ibarx=ibar
        kbar2=kbar
      endif
      rrmax=0.
      egsum=0.
      agsum=0.
      nmpc=0
      agsumpru=0;
      do 550 i=1,ibarx
      do 540 j=1,jbar
      do 530 k=1,kbar2
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (ijk.gt.ijkmax) go to 910
      if (p(ijk).eq.0.) go to 530
      egsum=egsum+e(ijk)*dxdydz ! sum(T00)*V in MeV
      agsum=agsum+p(ijk)*dxdydz ! sum(N0)*V
      if (ijk.lt.0) write (105,1350) i,j,k,p(ijk),e(ijk)/1000. 
c	Number of marker particles per cell, mz >=1 
       mz=int(p(ijk)*dxdydz/(float(nredux*nreduz)*dpt)+.5)
       nmpc=nmpc+mz
      if (mz.ne.0) go to 330
c   insufficient matter to warrant a marker
      nzero=nzero+1
c       if (asumot.eq.0.) write (156,1490) dpt
      asumot=asumot+p(ijk)*dxdydz ! sum(N0)*V
      summzo=summzo+zmom(ijk)*dxdydz ! T0Z in MeV
      esumot=esumot+e(ijk) ! sum(T00)*V in MeV
      e(ijk)=0.
      p(ijk)=0.
      xmom(ijk)=0.
      ymom(ijk)=0.
      zmom(ijk)=0.
      capn(ijk)=0.
      go to 520
  330 continue
      if (mz.gt.mzmax) mzmax=mz
      lji=0
      laccpt=0
      laccp2=0
  340 lji=lji+1
c   Distribute marker particles
      if (krandm.eq.0) go to 360
      call random_number (x)
      call random_number (y)
      call random_number (z)
      if (krandm.eq.1) go to 500
      if (krandm.eq.2) go to 500
      f=1./float(krandm-1)
      x=x*f
      y=y*f
      z=z*f
      go to 500
  360 continue
c
c   no randomization ___________________________________________________
c   Keep only for backward compatibility
c
      x=.5
      y=.5
      z=.5
      if (nparp.eq.1) go to 370
      if (nparp.eq.2) go to 380
      if (nparp.eq.3) go to 390
      if (nparp.eq.4) go to 400
      if (nparp.eq.5) go to 410
      if (nparp.eq.7) go to 420
      stop 4542
  370 xdf=0
      ydf=0.
      go to 430
  380 xdf=parpi*mod((lji-1)/2,2)-.25
      ydf=parpi*mod((lji-1)/(nparp*nparpz),2)-.25
      go to 430
  390 xdf=parpi*(mod((lji-1)/nparpz,3)-1)
      ydf=parpi*(mod((lji-1)/(nparp*nparpz),3)-1)
      go to 430
  400 ydf=.25*mod((lji-1)/(nparp*nparpz),4)-.375
      xdf=.25*mod((lji-1)/nparpz,4)-.375
      go to 430
  410 ydf=.2*mod((lji-1)/(nparp*nparpz),5)-.4
      xdf=.2*mod((lji-1)/nparpz,5)-.4
      go to 430
  420 ydf=(mod((lji-1)/(nparp*nparpz),7)-3)/7.
      xdf=(mod((lji-1)/nparpz,7)-3)/7.
c
  430 if (i.lt.imid) xdf=-xdf
      if (nparpz.eq.1) go to 440
      if (nparpz.eq.2) go to 450
      if (nparpz.eq.3) go to 460
      if (nparpz.eq.4) go to 470
      if (nparpz.eq.5) go to 480
      stop 4567
  440 zdf=0.
      go to 490
  450 zdf=0.5*mod(lji,2)-.25
      go to 490
  460 zdf=.3333333333*(mod(lji,3)-1)
      go to 490
  470 zdf=.25*mod((lji-1),4)-.375
      go to 490
  480 zdf=.2*mod((lji-1),5)-.4
  490 x=xdf+x
      y=ydf+y
      z=zdf+z
c
c   ___________________________________________________
c
  500 continue
      xc=i-1.+x
      yc=j-1.+y
      zc=k-1.+z
c
c   To which of the original cells does this correspond?
c
      if (kslab.lt.8.and.j.eq.1) npxz=npxz+1
      if (kslab.eq.8.and.j.eq.jmid) npxz=npxz+1
      jo=int(yc*float(nreduy)+0.5)
      if (jo.gt.jom) jom=jo
      io=int(dfloat(nredux)*(xc+float(imid)-0.1)+.5)
      if (io.gt.iom) iom=io
      zko=(zc-float(kmid)-1.)*dzvfm
      if (io.le.0) go to 340
      if (ymin.gt.yc) ymin=yc
      if (ymax.lt.yc) ymax=yc
      kkk=int(zc+1.)
      if (kslab.ne.8) then
        dist=yc*yc+(xc-imid)**2
        sdist=sqrt(dist)
        write (103,*) xc,yc,zc,ko,sdist
        if (abs(zc-113.1).lt.0.05) write (103,*) xc,yc,zc,ko
        if (k.eq.kmid.and.dist.gt.12.5) write (102,*) xc,yc
        if (b0.lt.0.01.and.sdist.gt.7..and.ktnc.eq.12) write (102,*) 'om
     1 it',xc,yc,zc,sdist
        if (b0.lt.0.01.and.sdist.gt.13..and.ktnc.eq.16) go to 340
        if (b0.lt.0.01.and.sdist.gt.7..and.ktnc.eq.12) go to 340
      endif
c
      ko=ko+3
      capn(ijk)=capn(ijk)+dpt
      asumin=asumin+dpt
      npcmin=npcmin+1
      if (ko.gt.kpmax) go to 900
      laccpt=laccpt+1
      px(ko)=xc
      py(ko)=yc
      pz(ko)=zc
      jps=int(yc+1.)
      nypt(jps)=nypt(jps)+1
c       if (py(ko).le.0.) write (156,*) ' error in y',y
      if (kslab.eq.8) go to 510
  510 if (laccpt.lt.mz) go to 340
  520 continue
  530 continue
  540 continue
  550 continue
c       write (156,1190)
      egsum=egsum/econg
c
c__________________Finished with assignments____________________________
c______________Do checks and print out summaries________________________
c
      kp=ko
      npt=ko
c    
c       write (156,1190)
      if (kslab.ne.8) then
c         write (156,1420) kout,kout*dpt
      endif
      sumnut=dpt*float(ko)/6.
      sumnup=dpt*float(ko)/6.
      nyptmn=0
      nyptmx=0
      do 560 j=1,95
      if (nypt(j).eq.0) go to 560
      if (nyptmn.eq.0) nyptmn=j
      nyptmx=j
  560 continue
c       write (156,1120) nyptmn,nyptmx
      do 580 i=1,300
      do 570 j=1,700
      adens(i,j)=0.
  570 continue
  580 continue
c   Following only for graphics
      k=0
c       write (156,1190)
c       write (156,1130)
      do 600 j=nyptmn-1,nyptmx+1
      k=k+1
c       write (156,1110) j,nypt(j),3.0*float(nypt(j)-nypt(nyptmx-k+2))/
c      1 float(ko)
      do 590 koo=1,npt,3
      jpss=py(koo)+1.
      if (jpss.ne.j) go to 590
  590 continue
  600 continue
      do 610 koo=1,npt,3
      jpss=py(koo)+1.
  610 continue
c       write (156,1190)
      do 620 koo=1,npt,3
      jpss=py(koo)+1.
  620 continue
c
      asumpt=asumin
      esum=0.
      asumin=0.
      bc0tot=0.
      asumer=0.
      summz=0.
      presou=0.
      kcells=0
      cnmax=0.
 
      do 650 i=1,ibar
      do 640 j=1,jbar
      do 630 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (ijk.gt.ijkmax) go to 630
      if (capn(ijk).eq.0.) go to 630
      kcell=kcell+1
      asumer=asumer+capn(ijk)-p(ijk)
      gamm=1./sqrt(1.-w(ijk)**2) 
       if (j.eq.1) adens(i,k)=capn(ijk)/gamm ! baryon density in the rest frame
c     ------------------------------------------------------------------
c               From here p(ijk)=pressure given in MeV
c     ------------------------------------------------------------------
c   NOTE change
      if (kslab.ne.8) p(ijk)=e(ijk)/3. ! e(ijk)=T00  in MeV !!!!!!!
      if (kslab.eq.8) p(ijk)=pressure(ijk) ! in MeV/fm3
      p(ijk)=pressure(ijk) ! in MeV/fm3
      esum=esum+e(ijk)*dxdydz ! sum(T00)*V in MeV
      asumin=asumin+capn(ijk) ! sum(N0)
      summz=summz+zmom(ijk)*dxdydz ! T0Z in MeV
      presou=presou+p(ijk)  ! pressure in MeV/fm3
c     To get MeV  + averaging
      e(ijk)=e(ijk) ! in MeV/fm3
      zmom(ijk)=zmom(ijk) ! in MeV/fm3
c
      n=capn(ijk)/pnu(1)+.001
c       if (abs(capn(ijk)-float(n)*pnu(1)).gt.0.0001) write (156,*) 'ERROR
c      1 ',i,j,k,capn(ijk),n,capn(ijk)/pnu(1)
      if(cnmax.lt.capn(ijk))then
         cnmax=capn(ijk)
         icmax=i
         jcmax=j
         kcmax=k
         ncmax=n
      endif
  101 format(' N max',4i4,5f10.5)
      capn(ijk)=capn(ijk)*rdxdydz
      bc0tot=bc0tot+capn(ijk)
      e(ijk)=e(ijk) ! T00 in MeV/fm3
      p(ijk)=p(ijk) ! pressure in MeV/fm3
      xmom(ijk)=xmom(ijk) ! in MeV/fm3
      ymom(ijk)=ymom(ijk) ! in MeV/fm3
      zmom(ijk)=zmom(ijk) ! in MeV/fm3
      if (e(ijk).eq.0.) go to 630
c       check for fort.100
  630 continue
  640 continue
  650 continue
c       write (156,1190)
c       write (156,1460) mzmax,nzero
c       write (156,990) asumin,asumot,-asumer
      errr=asum-(asumin+asumot-asumer)
c       if (abs(errr).gt.0.0000000001) write (*,*) 'error sum matter'
      esum=esum
      esumot=esumot
c       write (156,1020) esum/econg,esumot/econg,(esum+esumot)/econg
c       write (156,1010) summz,summzo,summz+summzo
c       write (156,1000) presou,presin-presou
c       write (156,1190)
c       write (156,1030)
c       write (156,1440) (npt+2)/3,npxz,kpmax/3,100.*float(npt)/
c      1 float(kpmax),asum,esum,esum/asum
c       write (156,1040) kcell
c       write (156,1190)
      smz=0.
      kplmin=0
      kplmax=0

c___Following statements only for symmetric case________________________
c
c   Find center of points
      xave=0
      yave=0.
      zave=0.
      anumkp=float((npt+2))/3.
      firstx=1000.
      firsty=1000.
      firstz=1000.
      alastx=0
      alasty=0.
      alastz=0.
      do 850 ko=1,npt+2,3
      xave=xave+px(ko)
      yave=yave+py(ko)
      zave=zave+pz(ko)
      if (px(ko).lt.firstx) firstx=px(ko)
      if (py(ko).lt.firsty) firsty=py(ko)
      if (pz(ko).lt.firstz) firstz=pz(ko)
      if (px(ko).gt.alastx) alastx=px(ko)
      if (py(ko).gt.alasty) alasty=py(ko)
      if (pz(ko).gt.alastz) alastz=pz(ko)
  850 continue
c       write (156,1150)
      xaveps=xave/anumkp
      yaveps=yave/anumkp
      zaveps=zave/anumkp
c       write (156,1160) firstx,xave/anumkp,alastx
c       write (156,1170) firsty,yave/anumkp,alasty
c       write (156,1180) firstz,zave/anumkp,alastz
c       write (156,1190)
      do 880 i=1,ibar
      do 870 j=1,jbar
      do 860 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (xmom(ijk).ne.0.) write (*,*) 'xmom',i,j,k,xmom(ijk)
      if (ymom(ijk).ne.0.) write (*,*) 'ymom',i,j,k,ymom(ijk)
  860 continue
  870 continue
  880 continue
      call plane
      return
  890 write (*,*) 'x or y momenta or velocity .ne.0. - reprogram'
      stop 99
  900 continue!write (156,*) ' too many particles',i,j,k,ko,'>',kpmax
      stop 98
  910 continue !write (156,*) ' number of cells too large ',ijk
      stop 5113
  920 continue
  996 format (3i4,1x,f11.4,1x,f11.4,1x,f11.4,f11.4,2x,f11.4,1x,
     1 f11.4,f11.4,2x,f11.4,1x,f11.4)
  921 format('cells in shift will overlap',i6,' k1min',i4,' ik1mid',i4,
     1 ' kadjst',i5,'  kwidth/2',i4)
c       write (*,921) k1mid+kadjst+kwidth/2+2,k1min,k1mid,
c      1 kadjst,kwidth/2
c       write (156,921) k1mid+kadjst+kwidth/2+2,k1min,k1mid,
c      1 kadjst,kwidth/2

      stop 999
  930 continue
c       write (*,*) 'ijk not empty - need to rethink reassignments',i,j,k
c      1 ,i1,j1,k1,p(ijk),p(ijk1)
      stop 'assignment error'
c
  251 format('|',26x,' max p ',3i4,f11.7,26x,'|')
  102 format( '|',4x,'Max density at',2i3,i4,'  # markers=',i5, 
     1'  max density=',f7.4,'  n/n0=',f5.3,4x,'|')
  940 format ('|',23x,'Maximum density=',f8.4,' baryons/cell',22x,'|')
  950 format ('|',15x,'Initial N=',f7.2,3x,'E=',f9.2,' GeV   E/N=',f7.2,
     1' GeV',14x,'|')
  960 format ('|',12x,'Total Z Momentum  Z-',f9.2,'  Z+',f8.2,'  Sum',f9
     1 .4,' GeV',11x,'|')
  970 format (i6,3f11.6)
  980 format (3i4,' asumin',2f10.5,2f14.4,9f10.5)
  990 format ('|',6x,'sum capn',f9.4,'  omitted matter',f7.4,'  discreti
     1zation error',f8.4,6x,'|')
 1000 format ('|',4x,'Summed pressure=',f9.2,'  pressure lost from exclu
     1ded cells=',f8.2,' GeV',5x,'|')
 1010 format ('|',5x,'Momentum included=',f8.3,'  Momentum excluded=',f7
     1 .3,' Sum',f9.4,' GeV/c',5x,'|')
 1020 format ('|',7x,'Energy included=',f10.3,'  Energy excluded=',f7.3,
     1' Sum',f10.3,' GeV',6x,'|')
 1030 format ('|',22x,'Final Description of the Initial State:',21x,'|')
 1040 format ('|',26x,'Number of occupied cells=',i6,25x,'|')
 1050 format ('|',36x,' nredux=',i2,36x,'|')
 1060 format ('|',20x' Marker particles at centers of grids',20x,'|')
 1070 format ('|',24x' Randomization of marker particles',24x,'|')
 1080 format ('|',23x,'npart=',i1,2x,'npartz=',i1,2x,'nparp=',i1,2x,'npa
     1rpz=',i1,23x,'|')
 1090 format ('|',15x,'Max Momenta x:',f6.3,' y:',f6.3,' z:',f8.3,' Max
     1w',f6.3,15x,'|')
 1100 format (i3,i10,f10.4)
 1110 format ('|',30x,i2,i9,f11.7,30x,'|')
 1120 format ('|',25x,'min, max j, non-zero cells',2i3,25x,'|')
 1130 format ('|',26x,'Amount of asymmetry in y plane',26x,'|',/,'|',31x
     1 ,'j','    # pts','  asym/total number',22x,'|')
 1140 format (3f8.3,4i9)
 1150 format ('|',16x,'Particle Distribution   Lowest   Average   Highes
     1t',16x,'|')
 1160 format ('|',29x,'x plane',3f10.3,16x,'|')
 1170 format ('|',29x,'y plane',3f10.3,16x,'|')
 1180 format ('|',29x,'z plane',3f10.3,16x,'|')
 1190 format ('|',82x,'|')
 1200 format ('|',18x,'ibar=',i3,'  jbar=',i3,'  kbar=',i3,'  product='
     1 ,i8,18x,'|')
 1210 format ('|',16x,'_____Begin setup NON-SYMMETRIC STATE __________',
     1 19x,'|')
 1220 format ('|',2x,'INPUT:',1x,18a4,1x,'|')
 1230 format ('|',2x,'INPUT:',22x,'dxv=',f5.3,' dyv=',f5.3,' dzv=',f5.3s
     1 ,23x,'|')
 1520 format (f11.4)
 1521 format (f9.4)
 1522 format (f5.0)
 1523 format (f5.0)

 1240 format (10x,f12.5)
 1250 format (f13.6)
 1251 format (f13.6)
 1252 format (f13.6)
 1260 format ('|',19x,'Number of lines of input file fort.100=',i5,19x,'
     1|')
 1270 format ('|',12x,'Original imin=',i2,' imid=',i2,' imax=',i2,'  jmi
     1n=',i2,' jmid=',i2,' jmax=',i2,13x,'|',/,'|',19x,'  kmin=',i3,' km
     2id=',i3,' kmax=',i3,' kidth=',i2,26x,'|')
 1280 format ('|',14x,'Shifts to move to center of grid: x:',i3,' y:',i3
     1 ,'  z:',i4,15x,'|')
 1290 format (2i4,3f11.5)
 1300 format ('N=0',3i4,4f11.5)
 1310 format (2i4,8f11.4)
 1320 format ('|',12x,'Calculated N in profiles:',3f11.3,12x,'|',/'|',12
     1 x,'Calculated E in profiles:',3f11.1,12x,'|')
 1330 format ('|',24x,'pNtot=',f10.4,'   Etot=',f10.2,24x,'|')
 1340 format ('|',17x,'imin=',i2,' imid=',i2,' imax=',i2,'  jmin=',i2,'
     1jmid=',i2,' jmax=',i2,17x,'|',/,'|',22x,'  kmin=',i3,' kmid=',i3,'
     2 kmax=',i3,' kwidth=',i2,22x,'|')
 1350 format (3i4,2f10.4,i2)
 1360 format (i6,3i4,i10,3f10.3,f10.4)
 1370 format (3i4,3f10.3)
 1380 format (3i4,1x,f11.4,1x,f11.4,1x,f11.4,f11.4,2x,f11.4,1x,f11.4,f11
     1 .4,2x,f11.4,1x,f11.4)
 1390 format (20a4)
c 1400 format ('randm',3i4,3f10.4,3x,3f10.4)
 1410 format ('|',3x,' Density extremes of matter in cells: cmax cmin'
     1 ,f8.5,f11.8,' nucl/fm^3',3x,'|',/,'|',27x,'-> number of markers '
     2 ,f9.4,f12.8,13x,'|')
 1420 format (' # markers outside boundaries',i6,' amt of matter',f6.2)
 1430 format ('|',2x,'INPUT: vlxyz',f8.5,' volxyz',f8.5,' dxdy',f9.6,'
     1dzvfm'f9.5,15x,'|')
 1440 format ('|',5x,'number markers=',i8,' xz plane=',i6,'   kpmax=',i8
     1 ,' % used=',f7.3,'%',5x,'|',/,'|',9x,'A_calc=',f7.2,'  Esum=',f11
     2 .2,' GeV  E/n=',f9.3,' GeV/nucleon',10x,'|')
 1450 format ('|',14x,'volxyz=',f8.5,' fm^3  rdxdydz=',f8.5,2x,'vlxyz='
     1 ,f8.5,14x,'|')
 1460 format ('|',22x,'Maximum number of markers in a cell',i6,21x,'|',/
     1 ,'|',4x,'Number of cells for which amount of matter was too small
     2 for a marker=',i6,4x,'|')
 1470 format (i5,' omitted ',3i3,2f10.5,2f6.2)
 1480 format ('|',10x,'each marker will represent',f14.6,' amount ofbary
     1onic matter',11x,'|',/,'|',10x,'projected max. number of markers/c
     2ell=',f12.1,' min number=',f9.1,10x,'|',/,'|',17x,'predicted # mar
     3ker particles=',2f14.1,16x,'|')
 1490 format ('|',5x,'In the setup, the cells omitted because the amount
     1 of matter in the cell',5x,'|',/,'|',12x,'was smaller than one-hal
     2l of',f8.5,' are listed in fort.8',13x,'|')
 1500 format ('|',6x,5f14.7,6x,'|')
 1510 format ('|',21x,' Parameters of the projectile and target',21x,'|'
     1 ,/,'|',15x,'gamma      baryon/MP  gam(mn-B)/mn','  moment. mn/c',
     2 7x,'v/c',9x,'|')
      end


c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine prxx (anprz,enprz)
c
c    fast continuous eos but less accurate than press
c
c    input anpr    baryon density  [ 1./fm3]
c          enpr    energy density  [mev/fm3]
c         jtemp (in common /yestemp/ )
c             =1 temperature is calculated
c             =0 temperature is not calculated in qgp (:= faster)
c
c    calculates pressure at density anpr and energy density aepr
c    the result is in pres,temp and
c    ipres is an integer =
c
c        +99    :  hadronic phase
c        +99 or 77: mixed phase  last digit is qgp abundance
c        -22 - -26: qgp phase
c        -51 - -54: energy too low (t<0) (or T<1mev for mixture)
c        -61 - -63: density too low at T=Tcr no data in table
c
      implicit real*8 (a-h,o-z)
c
      common /prstbl/ prs(99,102), trs(99,102), pres, temp, chem, amurs(
     1 99,102), entrs(99,102), entro, ipres
      common /pne/ boumq(1001,3), bouhl(1001,3), vbmq(1000), vbhl(1000),
     1 boumh(1001,3), vbmh(1000)
      common /eqqcd/ tq, eq, pq, sq, anq, amuq
      common /cnuclr/ hvc, hv2, pi, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
      common /eqmix/ ancrh, ancrq, ecrh, ecrq, em, pm, sm, amum, imix,
     1 iph, ch, cq, cph, cpq, cmu, ceh, ceq, csh, csq, antq, anth, entq,
     2 enth, itsr, ier, tcr1, tcr2, tcr3, anth2, anth3, antq2, antq3
      common /pin/ anmin, emin, dn, de, nn, ne, opdne, dnbq, dnbh
      common /n0spec/ i1, j1, an1, e1, xi1, xi2, yj1, yj2, sprs(5), strs
     1 (5)
      common /yestemp/ jtemp
c
      anpr=anprz
      enpr=enprz
      kntz=0
   10 continue
      t=-0.0001
      temp=-0.0001
      pres=-0.001
      ipres=-99
c    initial density branching
      anin1=anpr
      enin1=enpr
      aepr=enpr
      if (anpr-ancrq.gt.0.) go to 20
      if (anpr-ancrh.gt.0.) go to 50
      if (anpr-antq.ge.0.) go to 90
      if (anpr-anth.ge.0.) go to 110
c                                    n>nq(t=1)
   20 continue
      if (aepr-ecrq.ge.0.) go to 30
      ipres=-51
      go to 270
   30 eqt0=hvc*1.5/pi**2*(1.5*pi**2*anpr)**(4./3.)+bbag
      if (aepr-eqt0.ge.0.) go to 40
      ipres=-52
      go to 270
   40 continue
      ipres=-22
      go to 190
c                                 nh(t=1)<n<nq(t=1)
   50 continue
      alamq=(anpr-ancrh)/(ancrq-ancrh)
      elimt1=alamq*ecrq+(1.-alamq)*ecrh
      if (aepr-elimt1.ge.0.) go to 60
      ipres=-53
      go to 270
   60 temp=((anpr-antq)+0.5*dnbq)/dnbq+1
      in=temp
      if (in-999.le.0.) go to 70
      in=999
   70 dd=anpr-antq-(in-1)*dnbq
      elim=((boumq(in,3)*dd+boumq(in,2))*dd+boumq(in,1))*dd+vbmq(in)
      if (aepr-elim.ge.0.) go to 80
      ipres=10
      go to 210
   80 ipres=-23
      go to 190
c                            nq(tc)<n<nh(t=1)
   90 continue
c
      if (ancrh.eq.ancrq) go to 140
c                            detour for pure hadronic m
      tempin=((anpr-antq)+0.5*dnbq)/dnbq+1
      in=tempin
      if (in-999.le.0.) go to 100
      in=999
  100 dd=anpr-antq-(in-1)*dnbq
      elim=((boumq(in,3)*dd+boumq(in,2))*dd+boumq(in,1))*dd+vbmq(in)
      ipres=10
      if (aepr-elim.lt.0.) go to 160
      ipres=-25
      go to 190
c                                     nq(tc)>n>nh(tc)
  110 continue
      if (aepr-entq.lt.0.) go to 120
      ipres=-26
      go to 190
  120 alamq=(anpr-anth)/(antq-anth)
      elitcr=alamq*entq+(1.-alamq)*enth
      if (enpr-elitcr.lt.0.) go to 130
      ipres=-63
      go to 270
  130 ipres=20
      go to 160
c                              n<nh(tc)
  140 continue
      if (aepr-entq.lt.0.) go to 150
      ipres=-24
      go to 190
  150 if (aepr-enth.le.0.) go to 170
      ipres=-61
      go to 270
c    calculation of boundary of hadronic phase at T=0
  160 ipres=ipres+10
  170 tempin=((anpr-0.0001)+0.5*dnbh)/dnbh+1
      in=tempin
      if (in-999.le.0.) go to 180
      in=999
  180 dd=anpr-0.0001-(in-1)*dnbh
      eht0=((bouhl(in,3)*dd+bouhl(in,2))*dd+bouhl(in,1))*dd+vbhl(in)
      if (aepr-eht0.gt.0.) go to 200
      ipres=-54
      kntz=kntz+1
c  CHANGES TO ENSURE REASONABLE PRESSURE  9 MAY 2000
      if (kntz.ge.11) go to 270
      enpr=1.01*enpr
      go to 10
c                        Q G P
 190  continue
c Angel, 2024, not needed with new EoS
c 190 pres=(enpr-4.*bbag)/3.
c      if (pres .lt. 0.) pres=0.
c      if (pres .gt. 0.) enpr=enpr
c      if (pres .eq. 0.) enpr=3.*bbag*((enpr+bbag)/(4.*bbag))**(4./3.)
      if (jtemp.eq.0) return
      call eqsqcd3 (anpr,enpr)
      pres=pq
      temp=tq
      entrop=sq
c     ADD ON FEBRUARY 19, 2024
      if (pres .ge. enpr) then
c       write (*,*) 'pijk>e (line2206):', 'pijk= ',pres,'eijk= ',enpr
c       pause
      endif
      return
c                      mixed phase
c                     hadronic phase
  200 continue
  210 continue
c                      standard interpolation
      an=anpr
      e=enpr
      ipres=99
      tempi=(an-anmin)/dn+1.
      i=int(tempi)
      tempe=(e-emin)/de+1.
      j=int(tempe)
      if (i-nn.lt.0.) go to 220
      i=nn-1
  220 if (j-ne.lt.0.) go to 230
      j=ne-1
  230 xi=anmin+(i-1)*dn
      yj=emin+(j-1)*de
      if (i.eq.i1.and.j.eq.j1) go to 250
      dxm=an-xi
      dxp=dxm-dn
      dym=e-yj
      dyp=dym-de
c       if (i.gt.98) write (156,*) ' PRS ',i,j
c       if (j.gt.101) write (156,*) ' PRS ',i,j
c       if (i.lt.1) write (156,*) ' PRS ',i,j
c       if (j.lt.1) write (156,*) ' PRS ',i,j,an,anmin,e,emin,dn,de
      z=opdne*(prs(i,j)*dxp*dyp+prs(i+1,j+1)*dxm*dym-prs(i,j+1)*dxp*dym
     1 -prs(i+1,j)*dxm*dyp)
      zt=opdne*(trs(i,j)*dxp*dyp+trs(i+1,j+1)*dxm*dym-trs(i,j+1)*dxp*dym
     1 -trs(i+1,j)*dxm*dyp)
      zmu=opdne*(amurs(i,j)*dxp*dyp+amurs(i+1,j+1)*dxm*dym-amurs(i,j+1)
     1 *dxp*dym-amurs(i+1,j)*dxm*dyp)
      temp=zt
      chem=zmu
      if (ancrh.eq.ancrq) go to 240
c                        detour for pure hadronic m.
      eup=elim
c      if (xi-antq) 440,440,420
      if (xi-antq.le.0.) go to 240
c                                spec interpolation
c                                at the edge of qgp,.
c                                to ensure continuous
c                                transition to qgp.
      diff=e-eup+2.*de
c      if (diff) 440,440,430
      if (diff.le.0.) go to 240
      pup=(eup-4.*bbag)/3.
      ipres=77
      corr=(pup-z)*diff/(2.*de)
      z=z+corr
  240 continue
      go to 260
  250 continue
c
c      spec iteration around n0, e0
c
      jel=1
      if (an.lt.an0) jel=jel+2
      if (e.lt.e1) jel=jel+1
      if (jel.eq.1) then
        dxm=an-an1
        dxp=an-xi2
        dym=e-e1
        dyp=e-yj2
        pij=0.
        tij=0.
        pipjp=prs(i1+1,j1+1)
        tipjp=trs(i1+1,j1+1)
        pijp=sprs(2)
        tijp=strs(2)
        pipj=sprs(3)
        tipj=strs(3)
        opdnes1=1./((xi2-an1)*(yj2-e1))
      endif
      if (jel.eq.2) then
        dxm=an-an1
        dxp=an-xi2
        dym=e-yj1
        dyp=e-e1
        pij=sprs(4)
        tij=strs(4)
        pipjp=sprs(3)
        tipjp=strs(3)
        pijp=0.
        tijp=0.
        pipj=prs(i1+1,j)
        tipj=trs(i1+1,j)
        opdnes1=1./((xi2-an1)*(e1-yj1))
      endif
      if (jel.eq.3) then
        dxm=an-xi1
        dxp=an-an1
        dym=e-e1
        dyp=e-yj2
        pij=sprs(5)
        tij=strs(5)
        pipjp=sprs(2)
        tipjp=strs(2)
        pijp=prs(i1,j1+1)
        tijp=trs(i1,j1+1)
        pipj=0.
        tipj=0.
        opdnes1=1./((an1-xi1)*(yj2-e1))
      endif
      if (jel.eq.4) then
        dxm=an-xi1
        dxp=an-an1
        dym=e-yj1
        dyp=e-e1
        pij=prs(i1,j1)
        tij=trs(i1,j1)
        pipjp=0.
        tipjp=0.
        pijp=sprs(5)
        tijp=strs(5)
        pipj=sprs(4)
        tipj=strs(4)
        opdnes1=1./((an1-xi1)*(e1-yj1))
      endif
      zs1=opdnes1*(pij*dxp*dyp+pipjp*dxm*dym-pijp*dxp*dym-pipj*dxm*dyp)
      zts1=opdnes1*(tij*dxp*dyp+tipjp*dxm*dym-tijp*dxp*dym-tipj*dxm*dyp)
      temp=zts1
      z=zs1
  260 continue
      pres=z
      return
c              no solution
  270 pres=-0.001
      return
      end
      block data quasw
c
ccc This block data gives the weights and points for calculating the
ccc Guassian type integral @int_0^@infinit f(x)*x*dexp(-x)*dx
c
      dimension qx(100), qw(100)
      common /imsl/ qx, qw
      data qw /0.221730E-02,0.133951E-01,0.407054E-01,0.913457E-01,0.172
     1 547E+00,0.291542E+00,0.455661E+00,0.672259E+00,0.948720E+00,0.129
     2 269E+01,0.171163E+01,0.221321E+01,0.280523E+01,0.349560E+01,0.429
     3 247E+01,0.520364E+01,0.623796E+01,0.740324E+01,0.870891E+01,0.101
     4 633E+02,0.117749E+02,0.135543E+02,0.155097E+02,0.176518E+02,0.199
     5 900E+02,0.225348E+02,0.252968E+02,0.282871E+02,0.315160E+02,0.349
     6 968E+02,0.387399E+02,0.427594E+02,0.470665E+02,0.516735E+02,0.566
     7 005E+02,0.618548E+02,0.674539E+02,0.734146E+02,0.797517E+02,0.864
     8 822E+02,0.936259E+02,0.101203E+03,0.109227E+03,0.117723E+03,0.126
     9 708E+03,0.136215E+03,0.146260E+03,0.156872E+03,0.168063E+03,0.179
     $ 884E+03,0.192351E+03,0.205492E+03,0.219361E+03,0.233968E+03,0.249
     $ 358E+03,0.265578E+03,0.282681E+03,0.300696E+03,0.319682E+03,41*0.
     $ 00/
c
      data qx /.36330E-1,.121828,.256231,.439564,.671871,0.953205E+00,0.
     1 128361E+01,0.166320E+01,0.209203E+01,0.257024E+01,0.309795E+01,0.
     2 367529E+01,0.430238E+01,0.497938E+01,0.570648E+01,0.648384E+01,0.
     3 731168E+01,0.819019E+01,0.911961E+01,0.101002E+02,0.111321E+02,0.
     4 122158E+02,0.133514E+02,0.145392E+02,0.157796E+02,0.170730E+02,0.
     5 184195E+02,0.198198E+02,0.212741E+02,0.227828E+02,0.243464E+02,0.
     6 259653E+02,0.276401E+02,0.293711E+02,0.311590E+02,0.330042E+02,0.
     7 349073E+02,0.368690E+02,0.388899E+02,0.409706E+02,0.431118E+02,0.
     8 453143E+02,0.475787E+02,0.499059E+02,0.522967E+02,0.547520E+02,0.
     9 572727E+02,0.598597E+02,0.625141E+02,0.652368E+02,0.680290E+02,0.
     $ 708919E+02,0.738267E+02,0.768346E+02,0.799170E+02,0.830754E+02,0.
     $ 863112E+02,0.896260E+02,0.930215E+02,0.964995E+02,0.100062E+03,0.
     $ 103711E+03,0.107448E+03,0.111275E+03,0.115196E+03,0.119213E+03,0.
     $ 123328E+03,0.127544E+03,0.131864E+03,0.136292E+03,0.140832E+03,0.
     $ 145486E+03,0.150260E+03,0.155158E+03,0.160185E+03,0.165346E+03,0.
     $ 170647E+03,0.176094E+03,0.181695E+03,0.187458E+03,0.193390E+03,0.
     $ 199503E+03,0.205806E+03,0.212312E+03,0.219036E+03,0.225993E+03,0.
     $ 233201E+03,0.240684E+03,0.248465E+03,0.256577E+03,0.265056E+03,0.
     $ 273950E+03,0.283315E+03,0.293231E+03,0.303799E+03,0.315165E+03,0.
     $ 327552E+03,0.341322E+03,0.357178E+03,0.376936E+03/
      end


c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine plane
      implicit real*8 (a-h,o-z)
c
      integer*2 koost
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
      common /plcm7/ kextrn(30542411), koost(30542411)
      common /averps/ xaveps, yaveps, zaveps
      common /gamsol/ enzz, bbagzz, rrr, amom, bbagde
      common /gss/ sigx, sigy, sigz, radx, rady, radz, kspher
      common /prstbl/ prs(99,102), trs(99,102), pres, tempzz, chem,
     1 amurs(99,102), entrs(99,102), entro, ipres
      common /mincl/ imincl, imaxcl, jmincl, jmaxcl, kmincl, kmaxcl
      common /conser/ xmoint, ymoint, zmoint
      common /cnuclr/ hvc, hv2, pi0, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
      common /rot/ a, b, c
      common /numcmm/ numbk, numcel, numcst, numext, numcli, numlv,
     1 numout, numouu, numou2, jdim1, kfreez, kfrzcy, press, freez,
     2 compmy, koast, kext, cpntot, compmx, compe,iundee(58)
      common /maxinf/ icmax, jcmax, kcmax
      common /leaves/ xmomlv(300), ymomlv(300), zmomlv(300)
      common /epslon/ eps1, epsss
      common /pscm2/ cycen, istop, isume
      common /scaled/ hsclz, dzmod
      common /pscmz/ zz
      common /extrmes/ xextrn, xextrx, yextrx, zextrn, zextrx, kox
      common /hxi33/ hxi3a, xi, hxi
      common /newcrp/ entot, dentot, tempav, numcll
      common /bucket/ lablbk(1527120), xmomcl(1527120), ymomcl(1527120),
     1 zmomcl(1527120), encell(1527120), capncl(1527120), cnvvl(1527120)
      common /backet/ lbblbk(22), cpnncl(22), convol(22)
      common /bccket/ lablbc(1527120), celvol(1527120)
      common /midle/ imid, kmid, kend, imidpt, kmidpt, klmin, krmax,
     1 kdisp, ipmin, ipmax, imid1
      common /spcl/ delx, delz
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
      common /other/ eusual, epn0, t0, epion, ppion, edelt, delta, pdelt
     1 , temp, tempkn, eps, xth, epsc, dnsmax, ncycpt, ncycrp, ncycaf
      common /mssss/ outn, porton, npxz, npt, kout
      common /thestf/ rap00, Tfin, dxv, dyv, dzv, epsilon00, b0, rho00,
     1 e00
      common /pbarz/ ksuppr
      common /eqqcd/ tq, eq1, pq, sq, anq1, amuq
      common /misca/ kmove(700), numlvv(52)
      common /dmisca/ adens(300,700), eedns(300,700), aplot(300,700),
     1 bplot(300,700), cplot(300,700)
      common /miscnw/ cnfst, numfst
c
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233)
c
      dimension px(85877055), py(85877055), pz(85877055)
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
c
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c
c   Calculate average transverse momentum
      avept=0.
      avepz=0.
      numclz=0
      xmompl=0
      xmommn=0.
      ymompl=0
      ymommn=0.
      zmompl=0
      zmommn=0.
      numxpl=0
      numxmn=0
      numypl=0
      numymn=0
      numzpl=0
      numzmn=0
      anumzp=0
      anumzm=0
      totcel=0
      totpar=0.
      umaxmx=0.
      umaxmy=0.
      umaxmz=0.
      imidcl=(imincl+imaxcl)/2
      jmidcl=(jmincl+jmaxcl)/2
      kmidcl=(kmincl+kmaxcl)/2
      sumx=0.
      sumy=0.
      sumz=0.
      dpt=pnu(1)
      factop=dxdydz/dpt
      do 10 ko=1,npt,3
      i=px(ko)+1.
      j=py(ko)+1.
      k=pz(ko)+1.
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      anumq=capn(ijk)*factop
      totcel=totcel+1./anumq
      totzm=totzm+zmom(ijk)/anumq
c   X plane
      if (px(ko).gt.xaveps) then
        xmompl=xmompl+0.1*xmom(ijk)/capn(ijk)
        numxpl=numxpl+1
      else
        xmommn=xmommn+xmom(ijk)
        numxmn=numxmn+1
      endif
c   Y plane
      if (py(ko).gt.yaveps) then
        ymompl=ymompl+anumq*ymom(ijk)
        anumyp=anumyp+capn(ijk)/anumq
        numypl=numypl+1
      else
        ymommn=ymommn+anumq*ymom(ijk)
        anumym=anumym+capn(ijk)/anumq
        numymn=numymn+1
      endif
c   Z plane
      if (pz(ko).gt.zaveps) then
        zmompl=zmompl+anumq*zmom(ijk)
        anumzp=anumzp+capn(ijk)/anumq
        numzpl=numzpl+1
      else
        zmommn=zmommn+anumq*zmom(ijk)
        anumzm=anumzm+capn(ijk)/anumq
        numzmn=numzmn+1
      endif
   10 continue
      totxm=xmommn+xmompl
      totym=ymommn+ymompl
      totzm=zmommn+zmompl
      kwrit=0
      if (kwrit.ne.0) then
        write (*,80) numxpl,numxmn,xmommn,xmompl,totxm
        write (*,90) numypl,numymn,ymommn,ymompl,totym
        write (*,100) numzpl,numzmn,zmommn,zmompl,totzm
      endif
c       write (156,110) numxpl,numxmn,numypl,numymn,numzpl,numzmn
      sumz=0.
      sumx=0.
      sumy=0.
      do 40 i=1,ibar
      do 30 j=1,jbar
      do 20 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 20
      avept=avept+xmom(ijk)**2+ymom(ijk)**2
      avepz=avepz+zmom(ijk)**2
      numclz=numclz+1
      sumx=sumx+xmom(ijk)
      sumy=sumy+ymom(ijk)
      sumz=sumz+zmom(ijk)
      if (abs(xmom(ijk)).gt.umaxmx) umaxmx=abs(xmom(ijk))
      if (abs(ymom(ijk)).gt.umaxmy) umaxmy=abs(ymom(ijk))
      if (abs(zmom(ijk)).gt.umaxmz) umaxmz=abs(zmom(ijk))
   20 continue
   30 continue
   40 continue
      if (ncyc.ne.0.) go to 50
c   store initial momenta to check conservation laws
      xmoint=sumx*dxdydz
      ymoint=sumy*dxdydz
      zmoint=sumz*dxdydz
c       write (156,140) xmoint,ymoint,zmoint
      go to 60
   50 continue
      errxm=abs(sumx*dxdydz-xmoint)
      errym=abs(sumy*dxdydz-ymoint)
      errzm=abs(sumz*dxdydz-zmoint)
   60 continue
      if (mod(ncyc,5).eq.0) then
c         write (156,120) sumx*dxdydz,sumy*dxdydz,sumz*dxdydz,umaxmx*
c      1  dxdydz,umaxmy*dxdydz,umaxmz*dxdydz
c         write (156,130) sqrt(avept)/float(numclz),sqrt(avepz)/
c      1  float(numclz),(avept+avepz)/float(numclz)
      endif
      return
c
   70 format (i7,f7.1,f8.1,' -',f8.3,f8.3,f5.1,2f10.4,f7.1,2f8.3)
   80 format ('|',5x,'Numpl/mn',2i6,5x,'  Mx',3f13.4)
   90 format ('|',5x,'Numpl/mn',2i6,5x,'  My',3f13.4)
  100 format ('|',5x,'Numpl/mn',2i6,5x,'  Mz',3f13.4,/)
  110 format ('|',4x,'Numb. markers in +/- x direction:',2i6,'  y:',2i6,
     1'  z:',2i6,1x,'|')
  120 format ('|',5x,'Total momenta',3e10.2,2x,'Max. ',2f7.3,f8.3,5x,'|'
     1 )
  130 format ('|',14x,'Ave. pT=',f8.4,'  Ave. pZ=',f8.4,'   ave. p^2='
     1 ,f8.3,14x,'|')
  140 format ('|',12x,'Initial Momenta:',3e14.6,12x,'|')
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine picrel
c
c   Main subroutine; it controls the calculation and calls all 
c   subroutines that perform the calculations.
c
      implicit real*8 (a-h,o-z)
      real*4 secnds, base, start
c
      integer*2 koost
      INTEGER :: Nambi, Ndisc,Nsurf
      INTEGER :: NCornCallCycle, NCornSurfCycle, NCornSkipCycle
      INTEGER :: NCornInvalidCycle
      INTEGER :: IProgLen, IProgFill, ip, ITermCols
      INTEGER :: iEmax, jEmax, kEmax, iPmax, jPmax, kPmax
      INTEGER :: IAdaptiveFO, NInactiveFO, NInactiveReq, NMinFOCheck
      INTEGER :: ICycleExtendStep, ICycleMaxHard, IClosedFO
      LOGICAL :: has_cross_pre, has_cross_post, cube_valid
      REAL*8 ProgFrac, ProgPerc, maxEcyc, maxPcyc
      REAL*8 EPSFOE, EPSFOB
      CHARACTER*122 ProgBar
      CHARACTER*5 EscGreen
      CHARACTER*4 EscReset
      CHARACTER*16 SCols
      common /chktsm/ sumchem,sumtemp,sumentro,nlpc
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
      common /plcm7/ kextrn(30542411), koost(30542411)
      common /averps/ xaveps, yaveps, zaveps
      common /gamsol/ enzz, bbagzz, rrr, dumm1, dumm2
      common /gss/ sigx, sigy, sigz, radx, rady, radz, kspher
      common /prstbl/ prs(99,102), trs(99,102), pres, tempzz, chem,
     1 amurs(99,102), entrs(99,102), entro, ipres
      common /cnuclr/ hvc, hv2, pi0, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
      common /rot/ a, b, c
      common /numcmm/ numbk, numcel, numcst, numext, numcli, numlv,
     1 numout, numouu, numou2, jdim1, kfreez, kfrzcy, press, freez,
     2 compmy, koast, kext, cpntot, compmx, compe,iundee(58)
      common /maxinf/ icmax, jcmax, kcmax
      common /leaves/ xmomlv(300), ymomlv(300), zmomlv(300)
      common /epslon/ eps1, epsss
      common /pscm2/ cycen, istop, isume
      common /scaled/ hsclz, dzmod
      common /pscmz/ zz
      common /extrmes/ xextrn, xextrx, yextrx, zextrn, zextrx, kox
      common /hxi33/ hxi3a, xi, hxi
      common /newcrp/ entot, dentot, tempav, numcll
      common /bucket/ lablbk(1527120), xmomcl(1527120), ymomcl(1527120),
     1 zmomcl(1527120), encell(1527120), capncl(1527120), cnvvl(1527120)
      common /backet/ lbblbk(22), cpnncl(22), convol(22)
      common /bccket/ lablbc(1527120), celvol(1527120)
      common /volbuk/ lblvol(1527120), retvol(1527120)
      common /momen/ dumdum(61), npion, kmin, kmax, i1, i2, k1, k2
      common /yestemp/ jtemp
      common /juet/ dxyz1, pwq, amnW, beta2, gammaf, betat, cosfi, gjut,
     1 y0, hj, sumths
      common /timing/ diffp1, diffp2, diffp3, diffp4, diffp5, start
      common /plasma/ sumplast, iplasma, kplasma
      common /matric/ storev(929801), idensz(500), edensz(500)
      common /transp/ xte, yte, zte, hpx, hmx, hpy, hmy, hpz, hmz, i, j,
     1 k, l, ueff, veff, weff, cneff, xmeff, ymeff, zmeff, eneff, ko
      common /plottt/ kbarx, kbar1, kfrstt, imidptz, kmidptz, klminz
      common /midle/ imid, kmid, kend, imidpt, kmidpt, klmin, krmax,
     1 kdisp, ipmin, ipmax, imid1
      common /spcl/ delx, delz
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
      common /other/ eusual, epn0, t0, epion, ppion, edelt, delta, pdelt
     1 , temp, tempkn, eps, xth, epsc, dnsmax, ncycpt, ncycrp, ncycaf
      common /mssss/ outn, porton, npxz, npt, kout
      common /thestf/ rap00, Tfin, dxv, dyv, dzv, epsilon00, b0, rho00,
     1 e00
      common /countnevent/ nevent
      common /pbarz/ ksuppr
      common /eqqcd/ tq, eq1, pq, sq, anq1, amuq
      dimension ensort(1527120)
      common /misca/ kmove(700), numlvv(52)
      common /dmisca/ adens(300,700), eedns(300,700), aplot(300,700),
     1 bplot(300,700), cplot(300,700)
      common /miscnw/ cnfst, numfst
      common /corneliuss/ totenercornelius,totpxcornelius,
     1 totpycornelius,totpzcornelius,totbarcornelius
      common /foclose/ EPSFOE, EPSFOB, IAdaptiveFO, NInactiveFO,
     1 NInactiveReq, NMinFOCheck, ICycleExtendStep, ICycleMaxHard,
     2 IClosedFO
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233)
     
      common /timeconv/ dtconv
      
      
      common /posevent/ enerdensrfpos(1:300,1:300,1:700),enerdensrfpost(
     & 1:300,1:300,1:700),bardensrfpos(1:300,1:300,1:700),bardensrfpost(
     & 1:300,1:300,1:700),enerdenslfpos(1:300,1:300,1:700),bardenslfpos(
     & 1:300,1:300,1:700),vxpos(1:300,1:300,1:700),vxpost(1:300,1:300,1:
     & 700),vypos(1:300,1:300,1:700),vypost(1:300,1:300,1:700),vzpos(1:3
     & 00,1:300,1:700),vzpost(1:300,1:300,1:700),presspos(1:300,1:300,1:
     & 700),xmomdenslfpos(1:300,1:300,1:700),ymomdenslfpos(1:300,1:300,1
     & :700),zmomdenslfpos(1:300,1:300,1:700),posN0(1:300,1:300,1:700),p
     & osNx(1:300,1:300,1:700),posNy(1:300,1:300,1:700),posNz(1:300,1:30
     & 0,1:700),xmubpos(1:300,1:300,1:700),temppos(1:300,1:300,1:700)
     
      common /prevent/ enerdensrfpre(1:300,1:300,1:700),enerdensrfprev(1
     & :300,1:300,1:700),bardensrfpre(1:300,1:300,1:700),bardensrfprev(1
     & :300,1:300,1:700),enerdenslfpre(1:300,1:300,1:700),bardenslfpre(1
     & :300,1:300,1:700),vxpre(1:300,1:300,1:700),vxprev(1:300,1:300,1:7
     & 00),vypre(1:300,1:300,1:700),vyprev(1:300,1:300,1:700),vzpre(1:30
     & 0,1:300,1:700),vzprev(1:300,1:300,1:700),xmomdenslfpre(1:300,1:30
     & 0,1:700),ymomdenslfpre(1:300,1:300,1:700),zmomdenslfpre(1:300,1:3
     & 00,1:700),presspre(1:300,1:300,1:700),preNx(1:300,1:300,1:700),pr
     & eNy(1:300,1:300,1:700),preNz(1:300,1:300,1:700),preN0(1:300,1:300
     & ,1:700),xmubpre(1:300,1:300,1:700),temppre(1:300,1:300,1:700)
     
      dimension hypercube(0:1,0:1,0:1,0:1),dSigma(0:3,8)    
      dimension vmid(0:3,8)
     
      dimension px(85877055), py(85877055), pz(85877055)
      dimension numkt1(2), numkt2(2), icnt(20)
      
      DIMENSION Hypervx(0:1,0:1,0:1,0:1),Hypervy(0:1,0:1,0:1,0:1)
      DIMENSION Hypervz(0:1,0:1,0:1,0:1),Hypern(0:1,0:1,0:1,0:1)
      DIMENSION Hypertemp(0:1,0:1,0:1,0:1),Hypermub(0:1,0:1,0:1,0:1)
      DIMENSION Hyperpress(0:1,0:1,0:1,0:1)

     
      equivalence (ensort(1),xmomcl(1))
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
      equivalence (nx(1),nxt), (nx(2),nxp), (ny(1),nyt), (ny(2),nyp),
     1 (nz(1),nzt), (nz(2),nzp), (nplot(1),nplott), (nplot(2),nplotp),
     2 (gam(1),gamt), (gam(2),gamp), (smaj(1),smajt), (smaj(2),smajp),
     3 (smin(1),smint), (smin(2),sminp), (fx(1),fxt), (fx(2),fxp), (fz(1
     4 ),fzt), (fz(2),fzp), (w0(1),wt), (w0(2),wp), (sumnu(1),sumnut),
     5 (sumnu(2),sumnup)
      external fct1
c
c  Following are variables that are used as dummy constants;
c  they are zeroed here to simply avoid warning messages
c  from the compiler
c
      data en, merge, implict, marknm, radann, xannih, zannih, radan2 /0
     1 .,0,0,0,0.,0.,0.,0./
c
c  if a two-dimensional calculation is desired, set jdim1 to 1
c  units - UNITS - Units:
c  Pressure in code: p pijk p(ijk) = Press.in[MeV/fm3]/(939.5731*0.1453)
c  Energy dens. E in code e(ijk) = E in[MeV/fm3]/(0.1453*939.5731)
c  Momentum density M in code zmom = M[MeV/c/fm3]/(939.5731*0.1453)
c  Baryon density N in code capn(ijk)= N [fm-3]/0.1453
c  Rest energy dens. in code eps = e[MeV/fm3]/(0.1453*939.5731)
c  Rest baryon dens in code xi = n/n_0
c  Time in code dt t            = Time in[E-23sec]*(0.1453*2.99792458
c  Time in code dt t            = Time in[fm/c]*(0.1453)
c  Cell size in code is 1./0.1453(.145^1/3); Cell volume is 1./(.1453)
c
      Nambi=0
      Ndisc=0
c      xn0 = 0.15891d0 ! fm-3
c      ene0 = 146.51751415742d0 ! MeV/fm3
      xn0 = 0.15d0 ! fm-3
      ene0 = 146d0 ! MeV/fm3
      ecrit = 500d0!*ene0 !600.d0 ! energy density for the hypersurface (MeV/fm3)
      start=secnds(base)
      econ=931.5731
      speedc=2.99792458
      radius=1.18*dfloat(197)**(1./3.)
      kdisp=0
      kbarmx=700
      kskipp=1
      kbarex=11
      iwr=0
      kap=ap+.00000001
      kat=at+.00000001
c       write (*,*) kap,kat
c moved forward by 18 lines, to act if ioptosq.ne.0 - 29/7/2008, LC
      if (xdisp.eq.0.) kdisp=0
      if (kslab.eq.1) jdim1=1
c       if (kslab.eq.0) write (156,1250) ktnc
c       if (ioptosq.ne.0.and.kslab.eq.2) write (156,1240) ktnc
c       if (ioptosq.ne.0.and.kslab.eq.3) write (156,1240) ktnc
c       if (kslab.eq.3) write (156,850)
c       if (jdim1.le.2) write (156,1480)
c       if (jdim1.gt.2) write (156,1490)
c       if (kslab.eq.1) write (156,860)
c       if (kslab.eq.4.or.kslab.eq.5) write (156,1500)
c       if (kslab.eq.6) write (156,1510)
      kpbar=0
      if (ioptosq.ne.0) go to 10
c       if (ap.le.0.) write (156,1530)
c       write (156,1550) id
c       write (156,1560) ap,at,gevp,gevt,tnc,xdisp,keqost,icysum,eta,
c      1 comprs,bcoef,ioptos1,ioptos2,kslab,kw,kdebug
      if (ap.le.0) kpbar=1
      if (ioptosq.ne.0) kpbar=0
c                              _____________________________INPUT
      if (kpbar.eq.1) read (155,*) ksuppr
      go to 20
   10 continue
c     Input here is from (2): Initial state table
c       if (kslab.lt.4) write (156,1260) ap,at,icysum,epsilon00,b0,ktnc
      dy1=rap00
      beta=(dexp(dy1)-dexp(-dy1))/(dexp(dy1)+dexp(-dy1))
      e1=1./dsqrt(1.-beta*beta)
      e1=939.5731*(e1-1.)
      e1gev=.001*e1
      xdisp=b0
c       if (kslab.lt.4) write (156,1270) rap00,rho00,e00,beta,e1gev
   20 amnW=939.5731
      amnuc=939.5731
      wp=1.18*ap**(.3333333)
      wt=1.18*at**(.3333333)
      endn=.01
      ffac=2./3./pi**2
      radius=wt
      delx=radius/tnc
      volcub=(radius/tnc)**3
      em0=wp+wt
      en0=xdisp*em0
      vol1=4.*pi*wt**3/3.
      vol2=4.*pi*wp**3/3.
c       if (kslab.lt.4) write (156,1600) wp,wt,em0,en0,vol1,vol2
c       if (kslab.eq.4.or.kslab.eq.5) write (156,1610) wp,vol1
c       if (interp.eq.1) write (156,1340)
c       if (interp.eq.2) write (156,1350)
c       if (interp.eq.3) write (156,1360)
c       if (interv.eq.1) write (156,1320)
c       if (interv.eq.2) write (156,1330)
      delypl=tnc/wt
c      delypl is the number of cells to equal one fm.
      pwq=pi
      kfreez=0
      kfrzcy=20
c
c      kfreez = 0  no freeze out
c             = 1 freeze out base on density  (see 1900 lines on)
c
      twopi=pi+pi
      econ=931.5731
      if (keqost.le.0) econ=939.5731
      econg=0.9315731
      amass=(1.18/tnc)**3*at
      enmaxc=0.
      cpntot=0.
      leave=0
      mleave=0
      epn0=-8.
      zmovemx=0.
      xmovemx=0.
      if (keqost.lt.0) epn0=0.
      eusual=1.+epn0/econ
      eps1=eusual
      epsss=eps1
      emevp=1000.*gevp
      emevt=1000.*gevt
      wp=-dsqrt(1.-(econ/(emevp+econ))**2)
      wt=dsqrt(1.-(econ/(emevt+econ))**2)
      if (ioptosq.eq.0) go to 30
      gevp=e1gev
      gevt=e1gev
      wp=-beta
      wt=beta
   30 wtp=wp
      if (ap.le.0.) wt=0.
      if (kslab.ge.5) wp=.9
      amome=econ*wtp
      kleave=0
      em0=1.0
      en0=1.0
      em0n0=em0*en0
      ren0=1./en0
      rem0n0=1./em0n0
      remnsq=rem0n0*rem0n0
      eosgm1=bcoef
      compmx=0.
      compe=0.
      compmy=0.
      hxieq1=-2.*(931.5731/939.5731)/3.
      if (keqost.gt.0) hxieq1=eqost(dfloat(1))
c
c   Change  18/2/2006
c      hxieq1=1.
c       if (kslab.lt.3.and.kdebug.ne.0) write (156,1520) keqost,hxieq1
c      1 ,eusual
      eps0=-em0n0*hxieq1/eosgm1
      if (kresrt.ne.0) go to 40
c
c__7_____LINE 7  DEPENDS ON VALUE OF KSLAB
c                             _____________________________INPUT 7
      if (kslab.eq.0) then
        if (kresrt.eq.0) read (155,*,err=740) npart,npartz,nparp,nparpz
     1  ,krandm,dtmod,ncycpt,ncycrp,ncycaf
        go to 50
      endif
      if (kslab.le.3) read (155,*,err=740) nparp,nparpz,dtmod,ncycpt
     1 ,ncycrp,ncycaf
      if (kslab.gt.3) read (155,*,err=740) porton,dtmod,ncycpt,ncycrp
     1 ,ncycaf
      if (ncycaf.eq.0) ncycaf=100
      nparp=max0(nparp,1)
      nparpz=max0(nparpz,1)
      npart=nparp
      npartz=nparpz
      dzmod=1.
c      if (kslab.gt.0) krandm=1
c      krandm=0
c   Change above if one wants a lattice arrangement
      go to 50
   40 continue
      if (kslab.eq.0) go to 750
      read (3,*) dx,dy,dz
      dxdydz=dx*dy*dz
      read (3,*) ibar,jbar,kbar,ijbar
   50 continue
      if (ncycpt.eq.0) ncycpt=10
      if (ncycrp.eq.0) ncycrp=10
      dzmod=1.
c
c   Input variables pmin and fztmod removed, 28 Jan '08
c
c      apmin=1.+(pmin/939.5731)**2
c
c      if (ioptosq.lt.0.and.kdebug.ne.0) write (156,1230) pmin,apmin
c
c   Define a scale factor for use in phase
      hsclz=0.5/dzmod
c
c  ktemp ne 0 -> temperatures will written on tape 3
c

      if (dtmod.le.0..or.dtmod.gt.1.) dtmod=1.
c      Following are because of changes not yet implemented in
c      prelude to call to inter
      if (interp.eq.2.and.npart.ne.nparp) go to 60
      if (interp.eq.2.and.npartz.ne.nparpz) go to 60
      go to 70
   60 write (*,*) 'interp.eq.2.and.npart.ne.nparp or npartz.ne.nparpz'
      stop 51
   70 continue
      volnrm=1./dfloat(npart*npart*npartz)
      if (jdim1.eq.1) volnrm=volnrm*dfloat(npart)
      dzmod=1.
c
c  ncycpt: particle positions written every ncycpt cycles except for
c          for those to turb and turbf
c  ncycrp: write to 10 called every ncycrp cycles
c
c       if (kslab.ge.4) write (156,1430) porton,dtmod,ncycpt,ncycrp,ncycaf
c      1 ,ioptosq,ioptos1,ioptos2,ioptos3,ioptos4,kdeb
c       if (kslab.lt.5.and.kslab.ne.0) write (156,1440) nparp,nparpz,
c      1 krandm,dzmod,dtmod,ncycpt,ncycrp,ncycaf,kslab,ioptos1,ioptos2,
c      2 ioptos3,ioptos4,kdebug
c       if (kslab.eq.0) write (156,1420) npart,npartz,nparp,nparpz,krandm
c      1 ,dtmod,ncycpt,ncycrp,ncycaf,ioptosq,ioptos1,ioptos2,ioptos3
c      2 ,ioptos4,kdebug
      if (kresrt.eq.1) go to 80
c
c  0 < dtmod < 2  if < 1, a smaller time step is used
      sbase=3./(4.*pi*en0)
c  with this choice of sbase, the volume is at/en0
      smajt=(sbase*at)**(1./3.)
      smajp=(sbase*ap)**(1./3.)
      gamt=1./dsqrt(1.-wt**2)
      gamp=1./dsqrt(1.-wp**2)
      smint=smajt/gamt
      sminp=smajp/gamp
      dx1=dx
c  dx, dy, dz are cell dimensions in program units [p.u.]
      dx=smajt/tnc
      dy=dx
      dz=dx
      delz=delx
c
      dxdydz=dx*dy*dz
      rdxdydz=1./dxdydz
      dxdydz=dx*dy*dz
      rdx=1./dx
      rdy=1./dy
      rdz=1./dz
c       write (*,*) ' smajt ', smajt,'    tnc ', tnc
c       write (*,*) ' dx ', dx,'    dy ', dy,'    dz ', dz
c       write (*,*) 'rdx ',rdx,'   rdy ',rdy,'   rdz ',rdz
c       pause
      rdxdydz=rdx*rdy*rdz
      hsizex=1.
      hsizez=1.
      if (hsizex.le.0.) hsizex=1.
      if (hsizez.le.0.) hsizez=1.
c
      if (mod(ibar,2).eq.1) ibar=ibar-1
      possbl=(dfloat(ijkmax)/3./2./dfloat(kbarex))**.3333333333333333
      if (kslab.eq.8) possbl=0.5*possbl
      ip=possbl+.0001
      jbar=ip
      ibar=2*ip
      if (kslab.eq.8) jbar=2*ip
      kmx=kbarex*ip
c
      if (mod(ibar,2).eq.1) ibar=ibar-1
      ibar=min0(ibar,300)
      if (mod(kmx,2).eq.1) kmx=kmx+1
      kbar=min0(kmx,kbarmx)
c
   80 continue
      ijbar=ibar*jbar
      ijkbar=ijbar*kbar
      ijkba3=nq*ijkbar
      ijk=nq*(ijbar*(kbar-1)+ibar*(jbar-1)+ibar-1)+1
      if (ijkba3.gt.ijkmax) go to 730
      fibar=dfloat(ibar)
      fjbar=dfloat(jbar)
      fkbar=dfloat(kbar)
      imid=(ibar+1)/2
      kmid=(kbar+1)/2
      nqib=nq*ibar
      nqijb=nq*ijbar
      nxt=npart
      nyt=npart
      nzt=npartz
      nxp=nparp
      nyp=nparp
      nzp=nparpz
      if (jdim1.ne.1) go to 90
      nyt=1
      nyp=1
   90 continue
      nplott=1
      nplotp=3
      nsc=8
      kplmin=0
      kplmax=0
      kplminx=0
      kplmaxx=0
c
      nlc=ijkbar*3
      if (ap.eq.0.) go to 100
      realt=0.5*radius/(tnc*abs(wp)*2.997924581)
      go to 110
  100 realt=0.5*radius/(tnc*2.997924581)
  110 continue
      do 120 i=1,lcmb-2
      aa1(i)=0.
      aa2(i)=0.
      aa3(i)=0.
  120 continue
c       if (ioptosq.eq.0) write (156,1470)
      if (ioptosq.ne.0) go to 130
      go to 240
  130 continue
      if (kresrt.eq.1) go to 140
      diff1=secnds(start)
c     ------------------------------------------------------------------
       if (kslab.eq.8.and.kresrt.eq.0) call readinput
c     ------------------------------------------------------------------
c     From here dxdydz=0.259694072 fm3 (0.638^3)
c      e(ijk)! T00 in MeV/fm3
c      capn(ijk) N0 in fm-3
c      p(ijk) ! pressure in MeV/fm3
c      xmom(ijk) ! in MeV/fm3
c      ymom(ijk) ! in MeV/fm3
c      zmom(ijk) ! in MeV/fm3

      rdx=1./dx
      rdy=1./dy
      rdz=1./dz
      rdxdydz = rdz * rdy * rdz
c       write (*,*) ' dx ', dx,'    dy ', dy,'    dz ', dz
c       write (*,*) 'rdx ',rdx,'   rdy ',rdy,'   rdz ',rdz
c       pause




c                  write header for OSCAR output
c        call oschydro_header(thydro_start,itmax,number)

c      rewind 334
c       write (*,*) 'timepre= ',ncyc*dtconv
c      write (*,*) 'ncyc iteration loop= ',ncyc
      do i=1,ibar
      do j=1,jbar
      do k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      bardensrfpre(i,j,k)=0.
      bardensrfpre(i,j,k)=sqrt(1-w(ijk)**2-u(ijk)**2-v(ijk)**2)*
     1 capn(ijk)
      bardenslfpre(i,j,k)=0.
      bardenslfpre(i,j,k)=capn(ijk)
      enerdenslfpre(i,j,k)=0.
      enerdenslfpre(i,j,k)=e(ijk)
      enerdensrfpre(i,j,k)=0.
      enerdensrfpre(i,j,k)=(1.-w(ijk)**2-u(ijk)**2-v(ijk)**2)*(e(ijk)+
     1 p(ijk))-p(ijk)
      xmomdenslfpre(i,j,k)=0.
      xmomdenslfpre(i,j,k)=xmom(ijk)
      ymomdenslfpre(i,j,k)=0.
      ymomdenslfpre(i,j,k)=ymom(ijk)
      zmomdenslfpre(i,j,k)=0.
      zmomdenslfpre(i,j,k)=zmom(ijk)
      vxpre(i,j,k)    = 0.
      vxpre(i,j,k)    = u(ijk)
      vypre(i,j,k)    = 0.
      vypre(i,j,k)    = v(ijk)
      vzpre(i,j,k)    = 0.
      vzpre(i,j,k)    = w(ijk)
      presspre(i,j,k) = 0.
      presspre(i,j,k) = p(ijk)
      temppre(i,j,k)  = 0.
      temppre(i,j,k)  = temptab(enerdensrfpre(i,j,k),bardensrfpre(i
     & ,j,k))
      xmubpre(i,j,k)   = 0.
      xmubpre(i,j,k)   = chemtab(enerdensrfpre(i,j,k),bardensrfpre(
     & i,j,k))
      if (e(ijk) .ne. 0. .or. capn(ijk) .ne. 0.) then 
      write (333,*) i,j,k,e(ijk),capn(ijk),p(ijk)
      endif
c.......................................................................
c       collect cells below switching criterion before evolution
      if(e(ijk) .gt. 1d-8 .and. enerdensrfpre(i,j,k) .le. ecrit)then
      tf=ncyc*dtconv
      xifr=(float(i)-float(ibar)*0.5d0-0.5d0)*sngl(dx)
      yjfr=(float(j)-float(jbar)*0.5d0-0.5d0)*sngl(dy)
      zkfr=(float(k)-float(kbar)*0.5d0-0.5d0)*sngl(dz)
      dstfr=dx*dy*dz
      dsxfr=0.0d0
      dsyfr=0.0d0
      dszfr=0.0d0
      tfr=temptab(enerdensrfpre(i,j,k),bardensrfpre(i,j,k))
      chemqfr=chemtab(enerdensrfpre(i,j,k),bardensrfpre(i,j,k))
      chemsfr=0.
      vxfr= u(ijk)
      vyfr= v(ijk)
      vzfr= w(ijk)

      call oschyper(221,tf,xifr,yjfr,zkfr,enerdensrfpre(i,j,k),
     1 presstab(enerdensrfpre(i,j,k),bardensrfpre(i,j,k)),
     2 bardensrfpre(i,j,k),vxfr,vyfr,vzfr,tfr,chemqfr,chemsfr,
     3 dstfr,dsxfr,dsyfr,dszfr)

      end if
      enddo
      enddo
      enddo


      nlc=ijkbar*3
c       write (156,760)
      kdisp=1
      dpt=pnu(1) ! dimensionless
  140 continue
      if (kresrt.eq.0) go to 190
  190 centrx=(ibar+1)/2
      centrz=(kbar+1)/2
c       write (156,870)
      etot=0.
      cntot=0.
      zmptot=0.
      wptot=0.
      xmptot=0.
      do 220 k=1,kbar
      do 210 j=1,jbar
      do 200 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 200
      coverp=capn(ijk)/pnu(1)
      cntot=cntot+capn(ijk) ! baryon density in fm-3
      etot=etot+e(ijk) ! energy density in MeV/fm3
      if (zmom(ijk).lt.0.) go to 200
      zmptot=zmptot+zmom(ijk) ! in MeV/fm3
      wptot=wptot+w(ijk) 
      xmptot=xmptot+xmom(ijk) ! in MeV/fm3
  200 continue
  210 continue
  220 continue
c       write (156,770) cntot*dxdydz,etot*dxdydz,zmptot*dxdydz,wptot
      dxv=nredux*dxv
      dyv=nredux*dyv
      dzv=nreduz*dzv
      xcen=0.
      zcen=0.
      kp=(npt+2)/3
      do 230 ko=1,npt,3
      zcen=zcen+pz(ko)
      xcen=xcen+px(ko)
  230 continue
c       if (kdebug.ne.0) write (156,*) 'cen',3.*xcen/dfloat(kp),3.*zcen
c      1 /dfloat(kp)
      if (kslab.gt.1) then
        delx=dzv
        delz=dzv
        volcub=delz**3
      endif
c       write (156,840)
c c       write (156,880)
c       write (156,1130) dx,dy,dz,dxv,dyv,dzv
c       write (156,1140) dx/dxv,dy/dyv,dz/dzv,dx*dy*dz/(dxv*dyv*dzv),
c      1 rdxdydz,rdx,dxdydz
      if (ioptos4.ne.0) write (10,780) alam,Bbag
      if (ioptos4.ne.0) write (10,1280) dx,dy,dz,dxv,dyv,dzv
  240 continue
      jkskp=1
      ko=-2
      ymax=0.
      kp=npt
      if (kresrt.eq.1) go to 330
c  write to files the initial particle positions
  250 ko=ko+3
      if (ko.ge.kp) go to 260
      xte=px(ko)
      yte=py(ko)
      zte=pz(ko)
      ix=px(ko)+1.
      jx=py(ko)+1.
      kx=pz(ko)+1.
      ijk=nq*(ijbar*(kx-1)+ibar*(jx-1)+ix-1)+1
c       if (capn(ijk).eq.0.) write (156,970) ko,ix,jx,kx,px(ko),py(ko),
c      1 pz(ko),ijk,ibar,ijbar,ko,kp,kp-ko
      if (capn(ijk).eq.0.) write (*,970) ko,ix,jx,kx,px(ko),py(ko),pz(ko
     1 ),ijk,ibar,ijbar
      if (capn(ijk).eq.0.) write (*,*) 'CAPN = 0 in picrel  '
      if (capn(ijk).eq.0.) stop 1186
      if (yte.gt.ymax) ymax=yte
      go to 250
  260 continue
      if (jkskp.eq.0) go to 330
      nsumcl=0
      sumcl=0.
      sumcld=0.
      jmxx=min0(jbar,10)
      scapn=0.
      sen=0.
      do 320 j=1,jmxx
      if (jdim1.lt.3.and.j.gt.1) go to 320
      if (ioptos2.ne.0) write (kw,*)
      if (ioptos2.ne.0) write (kw,*) '_____ j=',j,'________'
      nsum=0
      kijend=0
      do 310 i=1,ibar
      do 300 k=1,kbar
      eedns(i,k)=0.
      adens(i,k)=0.
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (ijk.gt.lcmb) write (*,*) 'PICREL: ijk too large',i,j,k,ijbar
     1 ,ibar,ijk,ijkmax,lcmb
c       if (ijk.gt.lcmb) write (156,*) 'PICREL: ijk too large',i,j,k,ijbar
c      1 ,ibar,ijk,ijkmax,lcmb
      if (ijk.gt.lcmb) stop 1222
      u(ijk)=0.
      if (capn(ijk).eq.0.) go to 300
c       if (kdebug.gt.0) write (156,1150) i,j,k,ijk,capn(ijk),e(ijk),
c      1 p(ijk),u(ijk),v(ijk),w(ijk),xmom(ijk),ymom(ijk),zmom(ijk)
      scapn=scapn+capn(ijk)
      sen=sen+e(ijk)
      if (k.gt.kijend) kijend=k
      gamm=1./dsqrt(1.-w(ijk)*w(ijk))
      nsumcl=nsumcl+1
      sumcl=sumcl+capn(ijk)
      sumcld=sumcld+capn(ijk)*dxdydz
      capn(ijk)=dxdydz*capn(ijk)
      adens(i,k)=capn(ijk) ! N0 [fm-3]
      n=capn(ijk)/pnu(1)+.0001
      nsum=nsum+n
      eedns(i,k)=n
      if (dabs(capn(ijk)-dfloat(n)*pnu(1)).gt.0.0001.and.interp.eq.1)
     1 write (16,1230) ncyc,i,j,k,ijk,capn(ijk),n,capn(ijk)/pnu(1),pnu(1
     2 ),rdxdydz*capn(ijk),rdxdydz
      capn(ijk)=rdxdydz*capn(ijk)
      if (ioptosq.eq.0) go to 300
      u(ijk)=0.
  300 continue
  310 continue
      if (nsum.eq.0) go to 330
      if (ioptos2.ne.0) then
        write (kw,*) ' Nsum=',nsum
        write (kw,*) ' CAPN ',ncyc
        write (kw,*)
        write (kw,*) 'Marker particles/cell   j=',j
      endif
  320 continue
  330 continue
c       write (156,870)
      kend=kijend
c
      do 340 i=1,1527120
      lblvol(i)=lablbk(i)
      retvol(i)=cnvvl(i)
      lablbk(i)=0.
      cnvvl(i)=0.
  340 continue
c
      icentr=0
      if (dabs(gevp-gevt).lt.0.001) icentr=1
      if (gevp.ne.0..and.gevt.ne.0.) icentr=1
c  icentr=1 -> centre of mass
      iequal=0
      if (dabs(ap-at).lt..01.and.icentr.eq.1) iequal=1
      pmom=dsqrt(gevp*gevp+2*gevp*econg)
      xrapid=(gevp+econg+pmom)/(gevp+econg-pmom)
      xrapid=-0.5*dlog(xrapid)
      pmom=dsqrt(gevt*gevt+2*gevt*econg)
      trapid=(gevt+econg+pmom)/(gevt+econg-pmom)
      trapid=0.5*dlog(trapid)
      frapid=2.*dmax1(xrapid,trapid)
      i=5.*frapid+.5
      frapid=dfloat(i)*.2
      rapmin=-4.6
      if (frapid.lt.2.3) rapmin=-2.3
      rapdel=0.1
      if (frapid.lt.2.3) rapdel=0.05
      if (ioptosq.ne.0) xrapid=-rap00
      if (ioptosq.ne.0) trapid=rap00
c       if (kslab.le.3) write (156,1410) xrapid,trapid,frapid
      y1=abs(xrapid-trapid)
      deltarap=y1
      dy1=y1
      beta=(dexp(dy1)-dexp(-dy1))/(dexp(dy1)+dexp(-dy1))
      e1=1./dsqrt(1.-beta*beta)
      e1=econ*(e1-1.)
      e1gev=0.001*e1
      e1tev=0.001*e1gev
c       if (kslab.le.3.and.e1gev.le.1000000.) write (156,1390) e1,e1gev,
     1 beta
      if (kslab.le.3.and.e1gev.gt.1000000.) write (156,1400) e1gev,e1tev
c      1 ,beta
      kbr=kbar/2+1
      kbrx=3
      if (2*(kbar/2).ne.kbar) kbrx=2
      kx=(kbar+1)/2
      volcubz=volcub*dzmod
      if (kslab.gt.3) volcubz=delx**3
c       write (156,1310) radius,delx,delz,volcubz,delz/speedc,10.*delz
c      1 /speedc/speedc
      numkt1(1)=1
      numkt2(1)=3*npxz/2-2
      numkt1(2)=3*npxz+1
      numkt2(2)=3*(npt-npxz)/2+numkt2(1)-3
c       if (kdebug.ne.0) write (156,1380) numkt1,numkt2
      fimid=imid
      fkmid=kmid
      if (ap.eq.0.) go to 350
      dtote=2.5*dz/(dabs(wp)+dabs(wt))
      if (kslab.gt.0) dtote=1.25*dz
c  with this choice of dt, a marker particle requires at least four time
c  steps  (usually four to five) to travel the length of a cell.
      go to 360
  350 dtote=2.5*dz
  360 dtote=dtote*dtmod ! dtmod=0.5 if need to change this value do it from inhydro file
      dto=dtote
      dto=aint(10.*dtote)*0.1
      write(*,*)'dto',dto,' dtote',dtote

      if (dtote.lt.0.1) dto=aint(100.*dtote)*0.01
      if (dtote.lt.0.01) dto=aint(1000.*dtote)*0.001
      if (dtote.lt.0.001) dto=aint(10000.*dtote)*0.0001
      if (dtote.lt.0.0001) dto=aint(100000.*dtote)*0.00001
      if (kresrt.eq.0) t=0.
      if (kresrt.eq.0) ncyc=0


      dt=0.2*dto !dto=0.3

c       write(*,*)'dt ',dt,'dto',dto,' dtote',dtote
      tout=dto
      iprt=1
      numtd=0
      zpminz=-1212
      zpmaxz=1212
      iexpnd=0

      dtodx=dt*rdx
      dtody=dt*rdy
      dtodz=dt*rdz
c       write(*,*)'rdx ',rdx,'   rdy ',rdy,'   rdz ',rdz
c       write(*,*)'TEST dt',dt,'derivadas',dtodx,dtody,dtodz
c       pause
c
      xmomot=0.
      zmomot=0.
      freez=eta
      if (freez.eq.0.) kfreez=0
      if (freez.gt.0.01.and.freez.lt.0.99) kfreez=1
      if (kpbar.eq.0) go to 370
      if (ksuppr.eq.0) go to 370
c      call removs
  370 continue
      sumpar=sumnut+sumnup
      if (kslab.eq.7) convpu=1.
c
c   CHANGE  8 DEC 2010
c       write (156,1570) ibar,jbar,kbar,npxz,(npt+2)/3,dfloat(ibar)*dx
c      1 /convpu,float(jbar)*dy/convpu,dfloat(kbar)*dz/convpu
c       write (156,1290) ijkbar,ijkmax,dfloat(100*nq*ijkbar)/
c      1 dfloat(ijkmax)/3.,(npt+2)/3,kpmax/3,100.*dfloat(npt+2)/
c      2 dfloat(kpmax)
c  hsize is the smoothing distance - not implemented yet
c  if zero, it is set to the cell size



C   time problem 

c       rzero=1.18
c       dens0=3./(4.*pi*rzero**3)
c       dens0=dens0**.33333333
        dens0=1. !not needed anymore 02062024

      time=1./(dens0*speedc)
      timedt=time*dt        !time=  0.33356409672594067 en 10^(-23) s

c      timfmc=1./dens0
c      timdfmc=timfmc*dt
c      dtodz=dtodx
c       write (156,1190)
      diffs=secnds(start)
c       if(ioptos4.gt.10)write (156,800) diffs

c       dtconv=10.*dt*time/speedc !02062024
c      dtconv=5.*dt*time/speedc !02062024
c      dt = dtconv
c      write (*,*) 'dtconv ',dtconv
c       pause
 
       dtconv=1d0*dt !en fm ! 02062024
      
c       write (156,1170) dtote,dtodx,dtodz
c       write (156,1181) dt,timedt,dtconv
     
c   Old old      
 1180 format ('|',2x,'Each time step=',f7.5,' in program units; ',f8.6,'
     1 X 10(-23) s = ',f7.5,' fm/c',4x,'|'/,'|',82x,'|') 

c   New new
 1181 format ('|',2x,'Each time step=',f7.5,' in fm; ',f8.6,'
     1 X 10(-23) s = ',f7.5,' fm',4x,'|'/,'|',82x,'|')
      
C   time problem solved      
      
      
      if (ioptos4.ne.0) write (10,1450) timedt,dtconv,ibar,jbar,kbar,pnu
     1 (1),pnu(2)
      if (kdebug.ne.0) then
        write (13,1450) timedt
        write (11,1450) timedt
        write (12,1450) timedt
        write (14,1450) timedt
      endif
      timel=gamt*timedt
      if (gamp.gt.1.01.and.gamt.eq.0.) write (156,1460) timel
      IAdaptiveFO=1
      NInactiveFO=0
      IClosedFO=0
      NInactiveReq=20
      NMinFOCheck=40
      ICycleExtendStep=50
      fdt=0.5d0/dtmod
      ICycleMaxHard=max0(int(3.d0*icysum*fdt*1.5d0),
     1 icysum+int(200.d0*fdt*1.5d0))
      if (ICycleMaxHard.lt.icysum) ICycleMaxHard=icysum
      EPSFOE=1d-6
      EPSFOB=1d-8
      if (icysum.eq.0) go to 530
c
c _____________________________________________________________________
c|                                                                     |
c|     Begining of iteration loop        (ends at 3795)                |
c|_____________________________________________________________________|
c


  380 if (ncyc.gt.0) iprt=0

c     
      iplasma=0
      sumplast=0.
      timou1=secnds(start)
      modcyc=mod(ncyc,ncycaf)
c       write(*,*)'istop=',istop,ncyc,icysum
      kphs=0
      smx=0.
      smz=0.
      leave=0
      sme=0.
      sme1=0.
      capnto=0.
      freez=eta
c     _________________________________________________________________
c    |                  FREEZE OUT condition                           |
      if (compmy.gt.compmx) compmx=compmy
      numc=0
      nwrit=0
      sumsz=0.
      afin=0.
c    |_________________________________________________________________|
c
c    calculate total x, z momentum and energy
c
      ncell=0
      angx=0.
      angy=0.
      angz=0.
      angxp=0.
      angyp=0.
      angzp=0.
      ncellps=0
      do 480 k=1,kbar
      do 470 j=1,jbar
      do 460 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 460
      angx=-dfloat(j)*zmom(ijk)+dfloat(k-kmid)*xmom(ijk)+angx
      angy=dfloat(i-imid)*zmom(ijk)-dfloat(k-kmid)*xmom(ijk)+angy
      angz=-dfloat(j)*xmom(ijk)+(dfloat(i-imid)-.5)*ymom(ijk)+angz
      if (i-imid.le.0) go to 390
      if (k-kmid.le.0) go to 390
      ncellps=ncellps+1
      angxp=-dfloat(j)*zmom(ijk)+dfloat(k-kmid)*xmom(ijk)+angxp
      angyp=dfloat(i-imid)*zmom(ijk)-dfloat(k-kmid)*xmom(ijk)+angyp
      angzp=-dfloat(j)*xmom(ijk)+dfloat(i-imid)*ymom(ijk)+angzp
  390 ncell=ncell+1
      capnto=capnto+capn(ijk)*dxdydz
      smx=smx+xmom(ijk)*dxdydz
      smz=smz+zmom(ijk)*dxdydz
      sme=sme+e(ijk)*dxdydz
      if (j.ne.1) go to 430
      if (i.ne.imid+10) go to 400
      timlpt=dfloat(ncyc)*timedt
      timfmc=dfloat(ncyc)*dtconv
c   Contours
      zcl=(dfloat(k-kmid)-0.5)*dz/convpu
      zcl=(dfloat(k-kmid)-1.5)*dz/convpu
  400 continue
      if (i.ne.imid) go to 410
      timlpt=float(ncyc)*timedt
      timfmc=float(ncyc)*dtconv
      zcl=(dfloat(k-kmid)-0.5)*dz/convpu
      zcl=(dfloat(k-kmid)-1.5)*dz/convpu
      i1=i+1
      ijk1=nq*(ijbar*(k-1)+ibar*(j-1)+i1-1)+1
  410 continue
      if (k.ne.kmid) go to 420
      k1=k+1
      ijk1=nq*(ijbar*(k1-1)+ibar*(j-1)+i-1)+1
      k2=k+2
      ijk2=nq*(ijbar*(k2-1)+ibar*(j-1)+i-1)+1
      timlpt=float(ncyc)*timedt
      timfmc=float(ncyc)*dtconv
      xgp=(float(i-imid)-0.5)*dx/convpu
      xcl=(dfloat(i-imid)-0.5)*dx/convpu
  420 if (k.ne.kmid+20) go to 430
      timlpt=float(ncyc)*timedt
      timfmc=float(ncyc)*dtconv
      xgp=(float(i-imid)-0.5)*dx/convpu
      xcl=(dfloat(i-imid)-0.5)*dx/convpu
  430 continue
      if (i.ne.imid) go to 440
      if (k.ne.kmid) go to 440
      k1=k+1
      ijk1=nq*(ijbar*(k1-1)+ibar*(j-1)+i-1)+1
      k2=k+2
      k3=k+3
      ijk2=nq*(ijbar*(k2-1)+ibar*(j-1)+i-1)+1
      i1=i+1
      ijk3=nq*(ijbar*(k1-1)+ibar*(j-1)+i1-1)+1
      ijk4=nq*(ijbar*(k2-1)+ibar*(j-1)+i1-1)+1
      im=i-1
      ijk5=nq*(ijbar*(k1-1)+ibar*(j-1)+im-1)+1
      ijk6=nq*(ijbar*(k2-1)+ibar*(j-1)+im-1)+1
      ygp=(float(j)-0.5)*dx/convpu
      eave=0.25*(e(ijk1)+e(ijk2)+e(ijk3)+e(ijk4))
      cave=0.25*(capn(ijk1)+capn(ijk2)+capn(ijk3)+capn(ijk4))
      if (j.gt.3) go to 460
      if (ncyc.ge.0) go to 460
c       write (156,*) 'j=',j
      ijk7=nq*(ijbar*(k-1)+ibar*(j-1)+im-1)+1
      ijk8=nq*(ijbar*(k-1)+ibar*(j-1)+i1-1)+1
      ijks=nq*(ijbar*(k3-1)+ibar*(j-1)+im-1)+1
      ijkt=nq*(ijbar*(k3-1)+ibar*(j-1)+i-1)+1
      ijku=nq*(ijbar*(k3-1)+ibar*(j-1)+i1-1)+1
      i2=i+2
      ijkv=nq*(ijbar*(k3-1)+ibar*(j-1)+i2-1)+1
      ijkb=nq*(ijbar*(k-1)+ibar*(j-1)+i2-1)+1
      ijk9=nq*(ijbar*(k1-1)+ibar*(j-1)+i2-1)+1
      ijka=nq*(ijbar*(k2-1)+ibar*(j-1)+i2-1)+1
  440 continue
      if (j.gt.2) go to 460
      if (k.ne.kmid) go to 460
      if (i.ne.imid) go to 460
c       write (156,*) ' '
      do 450 indx=-4,5
      i1=imid+indx
      k1=kmid-3
      ijk1=nq*(ijbar*(k1-1)+ibar*(j-1)+i1-1)+1
      k2=k1+1
      k3=k1+2
      k4=k1+3
      k5=k1+4
      k6=k1+5
      k7=k1+6
      k8=k1+7
      k9=k1+8
      ijk2=nq*(ijbar*(k2-1)+ibar*(j-1)+i1-1)+1
      ijk3=nq*(ijbar*(k3-1)+ibar*(j-1)+i1-1)+1
      ijk4=nq*(ijbar*(k4-1)+ibar*(j-1)+i1-1)+1
      ijk5=nq*(ijbar*(k5-1)+ibar*(j-1)+i1-1)+1
      ijk6=nq*(ijbar*(k6-1)+ibar*(j-1)+i1-1)+1
      ijk7=nq*(ijbar*(k7-1)+ibar*(j-1)+i1-1)+1
      ijk8=nq*(ijbar*(k8-1)+ibar*(j-1)+i1-1)+1
      ijk9=nq*(ijbar*(k9-1)+ibar*(j-1)+i1-1)+1
c       if (indx.eq.0) write (156,*) ' '
c       if (indx.eq.2) write (156,*) ' '
c       write (156,1110) i1,e(ijk2),e(ijk3),e(ijk4),e(ijk5),e(ijk6),
c      1 e(ijk7),e(ijk8),e(ijk9)
      if (kdebug.eq.0) go to 450
      i1=imid-indx
      k1=kmid-indx
      ijk1=nq*(ijbar*(k1-1)+ibar*(j-1)+i1-1)+1
      k2=k1-1
      k3=k1-2
      ijk2=nq*(ijbar*(k2-1)+ibar*(j-1)+i1-1)+1
      ijk3=nq*(ijbar*(k3-1)+ibar*(j-1)+i1-1)+1
c       write (156,1110) i1,e(ijk1),e(ijk2),e(ijk3)
      i2=i1-1
      ijk4=nq*(ijbar*(k1-1)+ibar*(j-1)+i2-1)+1
      ijk5=nq*(ijbar*(k2-1)+ibar*(j-1)+i2-1)+1
      ijk6=nq*(ijbar*(k3-1)+ibar*(j-1)+i2-1)+1
c       write (156,1110) i1,e(ijk4),e(ijk5),e(ijk6)
  450 continue
  460 continue
  470 continue
  480 continue
c
      smx=smx+xmomot
      smz=smz+zmomot
      eeerr=origen-sme
      if (ncyc.eq.0) eeerr=0.
c       if (dabs(eeerr).gt.0.0000001) write (156,*) ' error in E _______'
c      1 ,origen,sme,origen-sme
      num=0
      do 490 i=1,1527120
      if (lablbk(i).eq.-1) go to 490
      num=num+1
  490 continue
c
      if (dtconv.ge.0.005) then
c         write (156,1580) t,ncyc*dtconv,ncyc,smx,smz,sme,eeerr,compmy
c      1  ,compmx,icmax,jcmax,kcmax
c        write (*,1580) t,ncyc*dtconv,ncyc,smx,smz,sme,eeerr,compmy
c      1  ,compmx,icmax,jcmax,kcmax
c      porcomp = 100. - float(icysum-ncyc)/float(icysum)*100.
c      write (*,*) porcomp,' % completed'
c         if (kdebug.gt.0) write (156,1210) sme,sme1,sme-sme1,capnto
        if (cnfst/cntot.gt.21.0) write (*,810) numfst,cnfst,100.*cnfst
     1  /cntot
        if (cnfst.gt.cnfstx) cnfstx=cnfst
      else
c         write (156,1590) t,ncyc*dtconv,ncyc,smx,smz,sme,eeerr,compmy
c      1  ,compmx,numcst,num,numext,numbk,icmax,jcmax,kcmax
c         if (kdebug.gt.0) write (156,1210) sme,sme1,sme-sme1,capnto
c         write (156,990) ncyc,ncell,ncellps,angx,angy,angz,angxp,angyp
c      1  ,angzp
        write (11,990) ncyc,ncell,ncellps,angx,angy,angz,angxp,angyp
     1  ,angzp
      endif
c       if (kdebug.ne.0) write (156,*) origen,sme
c
      if (ncyc.eq.0) origen=sme
      gamm=0.
      if (sumpar.ne.0.) gamm=capnto/sumpar
c
      cpntot=cpntot*dxdydz
c         write for table printed at the end of the run to (1)
      isume=0
      if (numcel+numext.le.1.and.ncyc.gt.3) isume=1
      compmy=0.
      numcst=0
      numbk=0
      numcel=0
      numext=0
      numlv=0
      call outfin (timfmc,centrz,ktemp)
c       write (*,*) 'istop value in subroutine picrel:', istop
      if (istop.eq.0) go to 620
c      stop 3562 ! added on December, 2023 If hydro pass through here we
c     have a problem with units
c                           GO TO 620
C
c   End of calculation; write quantities needed for restart
  530 continue
      if (icysum.gt.0) go to 570
      sme=0.
      do 560 i=1,ibar
      do 550 j=1,jbar
      do 540 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 540
      sme=sme+e(ijk)*dxdydz
  540 continue
  550 continue
  560 continue
      origen=sme
  570 continue
      anucsum=0.
      enucsum=0.
      do 600 k=1,kbar
      do 590 j=1,jbar
      do 580 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c       if (capn(ijk).lt.0.2) capn(ijk)=0.2
c       if (e(ijk).lt.13.0) e(ijk)=13.0
      if (capn(ijk).eq.0.) go to 580
      capn(ijk)=capn(ijk)
      e(ijk)=e(ijk)
      p(ijk)=p(ijk)
      xmom(ijk)=xmom(ijk)
      ymom(ijk)=ymom(ijk)
      zmom(ijk)=zmom(ijk)
      anucsum=anucsum+capn(ijk)
      enucsum=enucsum+e(ijk)
  580 continue
  590 continue
  600 continue
      i=-1
      j=-1
      k=-1
c       write (156,890) anucsum,enucsum,enucsum*.939573,enucsum*.939573
c      1 /anucsum
      if (timtt.eq.0.) timtt=timph1+timph2+timph3
      timt0=0.01*timtt
c       if(ioptos4.gt.10)write (156,820) timph1,timph2,timph3,timtt,timtt+
c      1 diffs,diff6
c       write (156,830) timph1/timt0,timph2/timt0,timph3/timt0,timtt/(
c      1 diff6+diffs)
      write (*,*) ' finished'
      return
c                            CONTINUE HERE
  620 numfil=0
      kstrz=0
c     
c     __________________________________________________________________
c                             STARTING PHASE I
c
c     NOTE THAT ENERGY, BARYON CHARGE AND MOMENTA ARE GIVEN IN:
c     e(ijk) T00 [MeV/fm3]
c     xmom(ijk) T0X [MeV/fm3]
c     ymom(ijk) T0Y [MeV/fm3]
c     zmom(ijk) T0Z [MeV/fm3]
c     capn(ijk) N0 [fm-3]
c     pressijk(ijk) [MeV/fm3]
c
      timou2=secnds(start)
      timout=timou2-timou1
      diff2=secnds(start)
C     ------------------------------------------------------------------
C       CHECKING FOR ENERGY, MOMEMTUM AND BARYON CHARGE CONSERVATION
      sumN0=0
      sumT00=0
      sumT0Z=0
      sumpress=0
      sumwijk=0
      do i=1,ibar
      do j=1,jbar
      do k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      sumN0=sumN0+capn(ijk)
      sumT00=sumT00+e(ijk)
      sumT0Z=sumT0Z+zmom(ijk)
      sumpress=sumpress+p(ijk)
      sumwijk=sumwijk+abs(w(ijk))
      enddo
      enddo
      enddo
c     verbose per-phase print disabled

c                STORING THE DATA OF THE LAST CYCLE  (ncyc-1)
c      rewind 334
      write (333,*) 'ncyc =',99999999
      write (334,*) 'ncyc =',99999999
      do k=1,kbar
      do j=1,jbar
      do i=1,ibar 
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      bardensrfpre(i,j,k)=0.
      bardensrfpre(i,j,k)=sqrt(1-w(ijk)**2-u(ijk)**2-v(ijk)**2)*
     1 capn(ijk)
      bardenslfpre(i,j,k)=0.
      if (ncyc.le.2) then
      bardenslfpre(i,j,k) = capn(ijk)
      else
      bardenslfpre(i,j,k) = bardenslfpos(i,j,k)
      endif
      enerdenslfpre(i,j,k)=0.
      if (ncyc.le.2) then
      enerdenslfpre(i,j,k)=e(ijk)
      else
      enerdenslfpre(i,j,k)=enerdenslfpos(i,j,k)
      endif
      enerdensrfpre(i,j,k)=0.
      enerdensrfpre(i,j,k)=(1.-w(ijk)**2-u(ijk)**2-v(ijk)**2)*(e(ijk)+
     1 p(ijk))-p(ijk)
      xmomdenslfpre(i,j,k)=0.
      xmomdenslfpre(i,j,k)=xmom(ijk)
      ymomdenslfpre(i,j,k)=0.
      ymomdenslfpre(i,j,k)=ymom(ijk)
      zmomdenslfpre(i,j,k)=0.
      zmomdenslfpre(i,j,k)=zmom(ijk)
      vxpre(i,j,k)    = 0.
      vxpre(i,j,k)    = u(ijk)
      vypre(i,j,k)    = 0.
      vypre(i,j,k)    = v(ijk)
      vzpre(i,j,k)    = 0.
      vzpre(i,j,k)    = w(ijk)
      presspre(i,j,k) = 0.
      presspre(i,j,k) = p(ijk)
      temppre(i,j,k)  = 0.
      temppre(i,j,k)  = temptab(enerdensrfpre(i,j,k),bardensrfpre(i
     & ,j,k))
      xmubpre(i,j,k)   = 0.
      xmubpre(i,j,k)   = chemtab(enerdensrfpre(i,j,k),bardensrfpre(
     & i,j,k))

      
c       if (vxpre(i,j,k) .ne. 0. .or. vypre(i,j,k) .ne. 0. .or.
c      1 vzpre(i,j,k) .ne. 0.) then
c       write (334,*) vxpre(i,j,k),vypre(i,j,k),vzpre(i,j,k)
c       endif
c      if (e(ijk) .ne. 0. .or. capn(ijk) .ne. 0. .or. p(ijk) .ne. 0.
c     1 .or. w(ijk) .ne. 0. .or. u(ijk) .ne. 0. .or. v(ijk) .ne. 0.
c     2 .or. zmom(ijk) .ne. 0. .or. ymom(ijk) .ne. 0. .or. xmom(ijk)
c     3 .ne. 0.) then
c      write (10,999)i,j,k,e(ijk),capn(ijk),
c     1 p(ijk),wijk,uijk,vijk,zmom(ijk),ymom(ijk),xmom(ijk)
c      endif
      enddo
      enddo
      enddo

      nege=0
      do i=1,ibar
      do j=1,jbar
      do k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (e(ijk).lt. 0.) then
      nege=nege+1
      e(ijk)=0.
      endif
      enddo
      enddo
      enddo

C     ------------------------------------------------------------------
      call phase1 (mleave,kpbar)

      nege2=0
      do i=1,ibar
      do j=1,jbar
      do k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (e(ijk).lt. 0.) then
      nege2=neg2+1
      e(ijk)=0.
      endif
      enddo
      enddo
      enddo
c
c     Note that after calling phase1 e(ijk),zmom(ijk),ymom(ijk) and 
c     xmom(ijk) are given in MeV. capn(ijk) is multiplied by dxdydz 
c     and pressure is given in MeV/fm3 b is given in MeV*fm-3, cmom and 
c     c in MeV2*fm-6
C     ------------------------------------------------------------------
C       CHECKING FOR ENERGY, MOMEMTUM AND BARYON CHARGE CONSERVATION
      sumN0=0
      sumT00=0
      sumT0Z=0
      sumpress=0
      sumwijk=0
      do i=1,ibar
      do j=1,jbar
      do k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      sumN0=sumN0+capn(ijk)
      sumT00=sumT00+e(ijk)
      sumT0Z=sumT0Z+zmom(ijk)
      sumpress=sumpress+p(ijk)
      sumwijk=sumwijk+abs(w(ijk))
      enddo
      enddo
      enddo
c     verbose per-phase print disabled

c     ------------------------------------------------------------------
c
      diff3=secnds(start)
      timph1=timph1+diff3-diff2
      numcll=0
      entot=0.
      dentot=0.
      tempav=0.
c     __________________________________________________________________
c                            STARTING PHASE II 
      if (iwr.eq.0) go to 670
      do 660 jv=1,3
c       write (156,*) 'Content for j=',jv
      do 650 iv=90,115
      do 630 jvx=1,15
      icnt(jvx)=0
  630 continue
      iwc=0
      do 640 kv=1,kbar
      ijk=nq*(ijbar*(kv-1)+ibar*(jv-1)+iv-1)+1
      if (capn(ijk).eq.0.) go to 640
      iwc=iwc+1
      if (iwc.gt.19) write (*,*) kv,iwc,ijk
      icnt(iwc)=kv
  640 continue
      if (iwc.eq.0) go to 650
c       if (icnt(1).eq.93) write (156,1060) iv,(icnt(kv),kv=1,iwc)
c       if (icnt(1).eq.94) write (156,1070) iv,(icnt(kv),kv=1,iwc)
c       if (icnt(1).eq.95) write (156,1080) iv,(icnt(kv),kv=1,iwc)
c       if (icnt(1).eq.96) write (156,1090) iv,(icnt(kv),kv=1,iwc)
c       if (icnt(1).eq.97) write (156,1100) iv,(icnt(kv),kv=1,iwc)
c       if (icnt(1).eq.98) write (156,1060) iv,(icnt(kv),kv=1,iwc)
  650 continue
  660 continue
  670 continue
c       if (kdebug.ne.0) write (156,*) 'phase 2'
      diff3=secnds(start)

      call phase2 (mleave,zpmin,zpmax,enmaxc)

      nege3=0
      do i=1,ibar
      do j=1,jbar
      do k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (e(ijk).lt. 0.) then
      nege3=nege3+1
      e(ijk)=0.
      endif
      enddo
      enddo
      enddo

C     ------------------------------------------------------------------
C       CHECKING FOR ENERGY, MOMEMTUM AND BARYON CHARGE CONSERVATION
      sumN0=0
      sumT00=0
      sumT0Z=0
      sumpress=0
      sumwijk=0
      do i=1,ibar
      do j=1,jbar
      do k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      sumN0=sumN0+capn(ijk)
      sumT00=sumT00+e(ijk)
      sumT0Z=sumT0Z+zmom(ijk)
      sumpress=sumpress+p(ijk)
      sumwijk=sumwijk+abs(w(ijk))
      enddo
      enddo
      enddo
c     verbose per-phase print disabled

c     ------------------------------------------------------------------
      diff4=secnds(start)
      timph2=timph2+diff4-diff3
      if (kslab.lt.8) call chksym
      diff5=secnds(start)

      do i=1,ibar
      do j=1,jbar
      do k=1,kbar
      bardensrfpost(i,j,k)=0.
      enerdensrfpost(i,j,k)=0.
      bardensrfprev(i,j,k)=0.
      enerdensrfprev(i,j,k)=0.
      bardensrfpos(i,j,k)=0.
      enerdensrfpos(i,j,k)=0.
      xmomdenslfpos(i,j,k)=0.
      ymomdenslfpos(i,j,k)=0.
      zmomdenslfpos(i,j,k)=0.
      vxpos(i,j,k)=0.
      vypos(i,j,k)=0.
      vzpos(i,j,k)=0.
      vxpost(i,j,k)=0.
      vypost(i,j,k)=0.
      vzpost(i,j,k)=0.
      vxprev(i,j,k)=0.
      vyprev(i,j,k)=0.
      vzprev(i,j,k)   =0.
      presspos(i,j,k) =0.
      temppos(i,j,k)  = 0.
      xmubpos(i,j,k)   = 0.

      enddo
      enddo
      enddo

c          ________________________ phase 3 _________________________
      call phase3 (kdisp,zpmin,zpmax,zpminz,zpmaxz,iexpnd)

      nege4=0
      do i=1,ibar
      do j=1,jbar
      do k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (e(ijk).lt. 0.) then
      nege4=nege4+1
      e(ijk)=0.
      endif
      enddo
      enddo
      enddo
      write (869,*) nege,nege2,nege3,nege4
c
C     ------------------------------------------------------------------
C       CHECKING FOR ENERGY, MOMEMTUM AND BARYON CHARGE CONSERVATION
      sumN0=0
      sumT00=0
      sumT0Z=0
      sumpress=0
      sumwijk=0
      maxEcyc=0.d0
      maxPcyc=0.d0
      iEmax=0
      jEmax=0
      kEmax=0
      iPmax=0
      jPmax=0
      kPmax=0
      do i=1,ibar
      do j=1,jbar
      do k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (isnan(e(ijk)).or.isnan(p(ijk)).or.isnan(capn(ijk))) then
      write(*,*) 'NaN bulk field at cyc,i,j,k=',ncyc,i,j,k,
     1 ' e=',e(ijk),' p=',p(ijk),' n=',capn(ijk)
      stop 9151
      endif
      if (dabs(e(ijk)).gt.huge(1.d0).or.dabs(p(ijk)).gt.huge(1.d0)) then
      write(*,*) 'Inf bulk field at cyc,i,j,k=',ncyc,i,j,k,
     1 ' e=',e(ijk),' p=',p(ijk)
      stop 9152
      endif
      sumN0=sumN0+capn(ijk)
      sumT00=sumT00+e(ijk)
      sumT0Z=sumT0Z+zmom(ijk)
      sumpress=sumpress+p(ijk)
      sumwijk=sumwijk+abs(w(ijk))
      if (dabs(e(ijk)).gt.maxEcyc) then
      maxEcyc=dabs(e(ijk))
      iEmax=i
      jEmax=j
      kEmax=k
      endif
      if (dabs(p(ijk)).gt.maxPcyc) then
      maxPcyc=dabs(p(ijk))
      iPmax=i
      jPmax=j
      kPmax=k
      endif
      enddo
      enddo
      enddo
      if (maxEcyc.gt.1.d12.or.maxPcyc.gt.1.d12) then
      write(*,*) 'WARNING unphysical bulk scale at cyc=',ncyc,
     1 ' |e|max=',maxEcyc,' at ',iEmax,jEmax,kEmax,
     2 ' |p|max=',maxPcyc,' at ',iPmax,jPmax,kPmax
      endif

      ProgFrac=0.d0
      if (icysum.gt.0) ProgFrac=dfloat(ncyc)/dfloat(icysum)
      if (ProgFrac.lt.0.d0) ProgFrac=0.d0
      if (ProgFrac.gt.1.d0) ProgFrac=1.d0
      SCols=' '
      call getenv('COLUMNS',SCols)
      ITermCols=0
      read(SCols,*,err=9112,end=9112) ITermCols
 9112 continue
      if (ITermCols.le.0) ITermCols=120
      IProgLen=ITermCols-32
      if (IProgLen.lt.10) IProgLen=10
      if (IProgLen.gt.80) IProgLen=80
      IProgFill=int(ProgFrac*dfloat(IProgLen)+0.5d0)
      if (IProgFill.lt.0) IProgFill=0
      if (IProgFill.gt.IProgLen) IProgFill=IProgLen
      ProgBar=' '
      ProgBar(1:1)='['
      do ip=1,IProgLen
      ProgBar(ip+1:ip+1)='.'
      enddo
      ProgBar(IProgLen+2:IProgLen+2)=']'
      do ip=1,IProgFill
      ProgBar(ip+1:ip+1)='#'
      enddo
      ProgPerc=100.d0*ProgFrac
      EscGreen=char(27)//'[92m'
      EscReset=char(27)//'[0m'

      write(*,1002) EscGreen,ProgBar(1:IProgLen+2),EscReset,ProgPerc,
     1 ncyc,icysum
      write(*,1003) sumT00*dxdydz,sumT0Z*dxdydz,sumpress*dxdydz,
     1 sumN0*dxdydz
c     ------------------------------------------------------------------
c
      call plane
      diff6=secnds(start)
      timph3=timph3+diff6-diff5
      timtt=timph1+timph2+timph3
      if (mod(ncyc,5).ne.0) go to 680
      timt0=0.01*timtt
c       if(ioptos4.gt.10)write (156,820) timph1,diffp3,timph2,timph3,timtt
c      1 ,diff6+diffs,diff6
c       write (156,830) timph1/timt0,diffp3/timt0,timph2/timt0,timph3/
c      1 timt0
  680 continue
      if (kslab.lt.8) call chksym
      kphs=7
      i=-1
      if (iplasma.gt.0) then
        sumplast=sumplast/float(iplasma)
c         write (156,1370) iplasma,sumplast,ncyc
      endif
c                       IMPLEMENTING CORNELIUS ALGORITHM


c      rewind 333
c      rewind 333
      rewind 336
      rewind 337
c      rewind 338
c     file number for OSCAR output
      number=21
c
      hypercube   = 0.
      dsigma      = 0.
      enetot      = 0.
      pxtot       = 0.
      pytot       = 0.
      pztot       = 0.
      bardenstot  = 0.
      deltaenergy = 0.
      rnegbflow   = 0.
      posbflow    = 0.
      rnegene     = 0.
      posene      = 0.
      t00prev     = 0.
      t00post     = 0.
      xn0prev      = 0.
      xn0post      = 0.
      NCornCallCycle = 0
      NCornSurfCycle = 0
      NCornSkipCycle = 0
      NCornInvalidCycle = 0
      
      
      do i=1,ibar
      do j=1,jbar
      do k=1,kbar


      if (isnan(vxpos(i,j,k))) stop 3997
      if (isnan(vxpre(i,j,k))) stop 3998
      if (isnan(vypos(i,j,k))) stop 3999
      if (isnan(vypre(i,j,k))) stop 4000
      if (isnan(vzpos(i,j,k))) stop 4001
      if (isnan(vzpre(i,j,k))) stop 4002
      if (isnan(vxpost(i,j,k))) stop 4003
      if (isnan(vxprev(i,j,k))) stop 4004
      if (isnan(vypost(i,j,k))) stop 4005
      if (isnan(vyprev(i,j,k))) stop 4006
      if (isnan(vzpost(i,j,k))) stop 4007
      if (isnan(vzprev(i,j,k))) stop 4008
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1

c     CHECK ENERGY AND BARYON CHARGE CONSERVATION

      gammaprev = 1./(1. - vxprev(i,j,k) - vyprev(i,j,k) 
     &                   - vzprev(i,j,k))
      gammapost = 1./(1. - vxpost(i,j,k) - vypost(i,j,k) 
     &                   - vzpost(i,j,k))
c      t00prev = t00prev + 
c     & ((enerdensrfprev(i,j,k) + presspre(i,j,k))*gammaprev**2 -
c     & presspre(i,j,k))*dx*dy*dz

      t00prev = t00prev + enerdenslfpre(i,j,k)*dx*dy*dz
c      t00post = t00post + 
c     & ((enerdensrfpost(i,j,k) + presspos(i,j,k))*gammapost**2 -
c     & presspos(i,j,k))*dx*dy*dz

      t00post = t00post + enerdenslfpos(i,j,k)*dx*dy*dz
     
c      xn0prev = xn0prev + gammaprev*bardensrfprev(i,j,k)*dx*dy*dz
c      xn0post = xn0post + gammapost*bardensrfpost(i,j,k)*dx*dy*dz


      xn0prev = xn0prev + bardenslfpre(i,j,k)*dx*dy*dz
      xn0post = xn0post + bardenslfpos(i,j,k)*dx*dy*dz
      
      
      if (i.lt.ibar.and.j.lt.jbar.and.k.lt.kbar) then

      premin=enerdensrfpre(i,j,k)
      premax=premin
      postmin=enerdensrfpos(i,j,k)
      postmax=postmin

      premin=dmin1(premin,enerdensrfpre(i,j,k+1))
      premax=dmax1(premax,enerdensrfpre(i,j,k+1))
      premin=dmin1(premin,enerdensrfpre(i,j+1,k))
      premax=dmax1(premax,enerdensrfpre(i,j+1,k))
      premin=dmin1(premin,enerdensrfpre(i,j+1,k+1))
      premax=dmax1(premax,enerdensrfpre(i,j+1,k+1))
      premin=dmin1(premin,enerdensrfpre(i+1,j,k))
      premax=dmax1(premax,enerdensrfpre(i+1,j,k))
      premin=dmin1(premin,enerdensrfpre(i+1,j,k+1))
      premax=dmax1(premax,enerdensrfpre(i+1,j,k+1))
      premin=dmin1(premin,enerdensrfpre(i+1,j+1,k))
      premax=dmax1(premax,enerdensrfpre(i+1,j+1,k))
      premin=dmin1(premin,enerdensrfpre(i+1,j+1,k+1))
      premax=dmax1(premax,enerdensrfpre(i+1,j+1,k+1))

      postmin=dmin1(postmin,enerdensrfpos(i,j,k+1))
      postmax=dmax1(postmax,enerdensrfpos(i,j,k+1))
      postmin=dmin1(postmin,enerdensrfpos(i,j+1,k))
      postmax=dmax1(postmax,enerdensrfpos(i,j+1,k))
      postmin=dmin1(postmin,enerdensrfpos(i,j+1,k+1))
      postmax=dmax1(postmax,enerdensrfpos(i,j+1,k+1))
      postmin=dmin1(postmin,enerdensrfpos(i+1,j,k))
      postmax=dmax1(postmax,enerdensrfpos(i+1,j,k))
      postmin=dmin1(postmin,enerdensrfpos(i+1,j,k+1))
      postmax=dmax1(postmax,enerdensrfpos(i+1,j,k+1))
      postmin=dmin1(postmin,enerdensrfpos(i+1,j+1,k))
      postmax=dmax1(postmax,enerdensrfpos(i+1,j+1,k))
      postmin=dmin1(postmin,enerdensrfpos(i+1,j+1,k+1))
      postmax=dmax1(postmax,enerdensrfpos(i+1,j+1,k+1))

      has_cross_pre=(premin.le.ecrit.and.premax.ge.ecrit.and.
     1 (premin.lt.ecrit.or.premax.gt.ecrit))
      has_cross_post=(postmin.le.ecrit.and.postmax.ge.ecrit.and.
     1 (postmin.lt.ecrit.or.postmax.gt.ecrit))

      if (has_cross_pre.or.has_cross_post) then
      call fillcube(i,j,k,1,enerdensrfpre,enerdensrfpos,hypercube)

      cube_valid=.true.
      do it=0,1
      do ix=0,1
      do iy=0,1
      do iz=0,1
      if (isnan(HyperCube(it,ix,iy,iz))) cube_valid=.false.
      if (dabs(HyperCube(it,ix,iy,iz)).ge.huge(1d0)) cube_valid=.false.
      enddo
      enddo
      enddo
      enddo

      if (cube_valid) then
      NCornCallCycle=NCornCallCycle+1
c                CALLING CORNELIUS SUBROUTINE
      CALL Cornelius(ecrit,hypercube,dSigma,Nsurf,Vmid,dt,dx,
     1 dy,dz,Nambi,Ndisc)
      NCornSurfCycle=NCornSurfCycle+Nsurf
c      if (Nsurf .gt. 1) then
c      write (*,*) 'Nsurf= ',Nsurf
c     endif
c INTERPOLATE VALUES ON HYPERSURFACE  and write them to hypersurface array
      call fillcube(i,j,k,1,bardensrfpre,bardensrfpos,Hypern) ! n is de baryon density in the rest frame
      call fillcube(i,j,k,1,vxpre,vxpos,Hypervx) ! x-component of velocity flow
      call fillcube(i,j,k,1,vypre,vypos,Hypervy) ! y-component of velocity FLOW
      call fillcube(i,j,k,1,vzpre,vzpos,Hypervz) ! z-component of velocity FLOW
      call fillcube(i,j,k,1,presspre,presspos,Hyperpress) ! pressure density
      call fillcube(i,j,k,1,temppre,temppos,Hypertemp) ! Temperature [MeV]
      call fillcube(i,j,k,1,xmubpre,xmubpos,Hypermub) ! Baryon chemical potential [MeV]


      t0=ncyc*dtconv ! time step evolution
c       x0=(float(i)-float(ibar)*.5-.5)*sngl(dx) !SNGL(A) converts the double precision real A to a default real value
c       y0=(float(j)-float(jbar)*.5-.5)*sngl(dy)
c       z0=(float(k)-float(kbar)*.5-.5)*sngl(dz)
c ibar=jbar=86 kbar=274
c     coordinates of vertices
      x0=(float(i)-float(ibar)*.5-.5)*sngl(dx)
      y0=(float(j)-float(jbar)*.5-.5)*sngl(dy)
      z0=(float(k)-float(kbar)*.5-.5)*sngl(dz)

c      if (enerdensrfpos(i,j,k).ne.0.) then
c      if (abs(enerdensrfpos(i,j,k)-ecrit)/enerdensrfpos(i,j,k)
c     1 .le.0.1) then
c      write (943,*) ncyc,t0,x0,y0,z0
c      endif
c      endif



      do ii=1,Nsurf

      Hypert = t0 + Vmid(0,ii)
      Hyperx = x0 + Vmid(1,ii)
      Hypery = y0 + Vmid(2,ii)
      Hyperz = z0 + Vmid(3,ii)

      emid       = Quadrilinear(Vmid(:,ii),HyperCube,dt,dx,dy,dz)
      barcharmid = Quadrilinear(Vmid(:,ii),Hypern,dt,dx,dy,dz)
      pressmid2  = Quadrilinear(Vmid(:,ii),Hyperpress,dt,dx,dy,dz)
      Vxmid      = Quadrilinear(Vmid(:,ii),Hypervx,dt,dx,dy,dz)
      Vymid      = Quadrilinear(Vmid(:,ii),Hypervy,dt,dx,dy,dz)
      Vzmid      = Quadrilinear(Vmid(:,ii),Hypervz,dt,dx,dy,dz)
      tempmid    = Quadrilinear(Vmid(:,ii),Hypertemp,dt,dx,dy,dz)
      mubmid     = Quadrilinear(Vmid(:,ii),Hypermub,dt,dx,dy,dz)
      if (emid .lt. 0.) stop
      if (barcharmid .lt. 0.) stop

      vsquare = Vxmid**2 + Vymid**2 + Vzmid**2

      if (vsquare >= 1.) then

          Vxmid = Vxmid*sqrt(0.99 / vsquare);
          Vymid = Vymid*sqrt(0.99 / vsquare);
          Vzmid = Vzmid*sqrt(0.99 / vsquare);
          vsquare = 0.99;

      endif

      gammid=1.d0/dsqrt(1.d0-Vxmid**2 - Vymid**2 - Vzmid**2)

      umudsigmu=gammid*(dSigma(0,ii)+Vxmid*dSigma(1,ii)
     1 +Vymid*dSigma(2,ii)+Vzmid*dSigma(3,ii))

      bflow=barcharmid*umudsigmu
      if((bflow).lt.0)then
      rnegbflow=rnegbflow+bflow
      else
      posbflow=posbflow+bflow
      endif

      tnn=(emid+pressmid2)*gammid**2-pressmid2
      tnx=(emid+pressmid2)*Vxmid*gammid**2
      tny=(emid+pressmid2)*Vymid*gammid**2
      tnz=(emid+pressmid2)*Vzmid*gammid**2
      ene=(dSigma(0,ii)*tnn+dSigma(1,ii)*tnx+dSigma(2,ii)*tny+
     1     dSigma(3,ii)*tnz)

      txx=(emid+pressmid2)*Vxmid*Vxmid*gammid**2+pressmid2
      txy=(emid+pressmid2)*Vxmid*Vymid*gammid**2
      txz=(emid+pressmid2)*Vxmid*Vzmid*gammid**2
      tyy=(emid+pressmid2)*Vymid*Vymid*gammid**2+pressmid2
      tyz=(emid+pressmid2)*Vymid*Vzmid*gammid**2
      tzz=(emid+pressmid2)*Vzmid*Vzmid*gammid**2+pressmid2

      pxflow=(dSigma(0,ii)*tnx+dSigma(1,ii)*txx+dSigma(2,ii)*txy+
     1     dSigma(3,ii)*txz)
      pyflow=(dSigma(0,ii)*tny+dSigma(1,ii)*txy+dSigma(2,ii)*tyy+
     1     dSigma(3,ii)*tyz)
      pzflow=(dSigma(0,ii)*tnz+dSigma(1,ii)*txz+dSigma(2,ii)*tyz+
     1     dSigma(3,ii)*tzz)

      if((ene).lt.0)then
      rnegene=rnegene+ene
      else
      posene=posene+ene
      endif

      enetot     = enetot     + ene
      pxtot      = pxtot      + pxflow
      pytot      = pytot      + pyflow
      pztot      = pztot      + pzflow
      bardenstot = bardenstot + bflow

      tau = sqrt(1d0-Vxmid**2-Vymid**2-Vzmid**2)*Hypert
      eta = 0.5d0*log((1d0+Vzmid)/(1d0-Vzmid))

c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c                   OUTPUT FOR SMASH HADRON SAMPLER
      write (949,*) Hypert,Hyperx,Hypery,Hyperz,dSigma(0,ii),
     & dSigma(1,ii),dSigma(2,ii),dSigma(3,ii),gammid,Vxmid*gammid,
     & Vymid*gammid,Vzmid*gammid,tempmid/1000d0,mubmid/1000d0,0d0,0d0,
     & 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  293 continue
      enddo

      else
      NCornInvalidCycle=NCornInvalidCycle+1
      endif
      else
      NCornSkipCycle=NCornSkipCycle+1
      endif
      endif
      enddo
      enddo
      enddo
      
      totenercornelius = totenercornelius + enetot
      totpxcornelius   = totpxcornelius   + pxtot
      totpycornelius   = totpycornelius   + pytot
      totpzcornelius   = totpzcornelius   + pztot
      totbarcornelius  = totbarcornelius  + bardenstot
      write(*,*) 'FO flux cyc=',ncyc,' (P0,Px,Py,Pz,B)=',
     1 enetot,pxtot,pytot,pztot,bardenstot
      write(*,*) 'FO flux accum=',ncyc,' (P0,Px,Py,Pz,B)=',
     1 totenercornelius,totpxcornelius,totpycornelius,totpzcornelius,
     2 totbarcornelius


      write (734,*) ncyc*dtconv,rnegene,posene,enetot,
     1 totenercornelius,rnegbflow,posbflow,bardenstot,
     2 totbarcornelius

      if (IAdaptiveFO.ne.0.and.ncyc.ge.NMinFOCheck) then
      if (dabs(enetot).gt.EPSFOE.or.dabs(bardenstot).gt.EPSFOB.or.
     1 NCornSurfCycle.gt.0) then
      NInactiveFO=0
      IClosedFO=0
      else
      NInactiveFO=NInactiveFO+1
      if (NInactiveFO.eq.NInactiveReq) then
      IClosedFO=1
      write (*,*) 'FO close candidate confirmed at cyc=',ncyc,
     1 ' inactive_cycles=',NInactiveFO,' enetot=',enetot,
     2 ' bardenstot=',bardenstot,' Nsurf=',NCornSurfCycle
      endif
      endif
      endif

      if (ncyc.gt.0.and.modcyc.eq.0.and.IAdaptiveFO.ne.0.and.
     1 (mod(ncyc,10).eq.0.or.ncyc.le.NMinFOCheck.or.IClosedFO.eq.1))
     2 then
      write (*,*) 'FO close monitor cyc=',ncyc,' inactive=',
     1 NInactiveFO,'/',NInactiveReq,' closed=',IClosedFO
      endif
C      write(339,*) ncyc,enetot,rnegene,posene,bardenstot,rnegbflow,
c     1 posbflow
      
      
c  995 format (2x,16f16.2,20x,32f16.2,20x,i6)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      ikkk=9999
      sumtemp=0.
      sumchem=0.
      sumentro=0.
      nlpc=0.
      if (ioptos4.ne.0.and.modcyc.eq.0) write (10,1300) ikkk,j,k,e(ijk)
     1 ,capn(ijk),p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk)
     2 ,xmom(ijk),tq,sq,amuq
      if (ioptos4.ne.0.and.modcyc.eq.0) then
      write (118,1300) ikkk,j,k,e(ijk)
     1 ,capn(ijk),p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk)
     2 ,xmom(ijk),tq,sq,amuq
      endif
      
      
      nlpc=nlpc+1
      ipmaxz=ipmax
      numfz1=0
      if (numfz1.eq.0) go to 700
      do 690 j=1,numfz1
      i=idensz(j)/100
      k=idensz(j)-100*i
      adens(i,k)=edensz(j)
  690 continue
  700 continue
      if (numlv.eq.0) go to 720
      do 710 i=1,numlv
      ijk=numlvv(i)
      capn(ijk)=0.
      e(ijk)=0.
      xmom(ijk)=0.
      ymom(ijk)=0.
      zmom(ijk)=0.
      u(ijk)=0.
      v(ijk)=0.
      w(ijk)=0.
  710 continue
  720 continue
      go to 380
  730 continue
c       write (156,1540) ibar,jbar,kbar,ijkba3,ijkmax
c       write (156,*) 'fixed, ibar, jbar, kbarmax=',ijkmax/(3*ibar*jbar),
c      1 ' fixed, jbar, kbar, ibarmax=',ijkmax/(3*jbar*kbar)
      stop 22
      return
  740 write (*,*) 'Error in input line 7, subroutine picrel'
      stop 7
  750 write (*,*) 'Cannot use restart with kslab=0'
c       write (156,*) 'Cannot use restart with kslab=0'
      stop 97
c
  760 format ('|',34x,'Back in Picrel',34x,'|')
  770 format ('|',82x,'|',/,'|',5x,'Initial values  N=',f9.4,' E=',f16.
     1 3,'  Zm+=',f16.4,'  Wp+=',f9.4,5x,'|',/,'|',82x,'|')
  780 format ('# alam= ',f8.4,'    Bbag= ',f10.5)
  790 format (3i4,' capn',f10.4,i5,'  E=',2f11.4,2f11.3)
  800 format ('|',29x,'start-up time',f6.2,' secs',29x,'|')
  810 format (20x,'cells with c>b^2 =',i3,'  amount of matter='f9.5,' ='
     1 ,f5.2,'%')
  820 format (/,'Time phase1:',f6.2,'s wt:',f8.2,'s phase2:',f8.2,'s  ph
     1ase3:',f7.2,'s  phase tot:',6f8.2,'s')
  830 format (12x,f16.2,'%',6x,f16.2,'%',9x,f17.2,'%',9x,f17.2,'%',13x,
     1 3f16.2)
  840 format ('|',33x,'__Begin hydro__',34x,'|',/,'|',82x,'|')
  850 format (14x,' NOTE: THIS HAS THE TEST SETUP OF UNIFORM STREAKS',20
     1 x,'|',/,83x,'|')
  860 format (19x,'Initial state consists of two uniform slabs',21x,'|',
     1 /,83x,'|')
  870 format ('|',82x,'|')
  880 format ('|',4x,22('_'),' HYDRO CELL SIZE DEFINITION: ',22('_'),5x,
     1'|')
  890 format ('|',4x,'at end N=',f13.7,'  E=',2f18.6,f18.7' MeV/nuc. |')
  900 format (10a4)
  910 format (5i3,2f10.5,8i3)
  920 format (2f10.5,8i3)
  930 format (20('__'),/,f10.4)
  940 format (2f10.5,3i4)
  950 format (5x,i3,3x,f10.6,8x,f24.15,11x,i7,8x,f9.6,9x,f17.14,10x,3f18
     1 .14)
  960 format (3i4,i8,5f14.7)
  970 format ('zero CAPN when it should not be',i8,3i4,3f9.4,6i10)
  980 format (10f12.5)
  990 format (3i5,' ang=',3f13.3,2x,3f13.3)
  999 format (3i4,' ',9(E12.5,' '))
 1000 format ('ncyc=',i3,' t=',f10.6,'  E_tot=',f24.15,' [GeV]  kp=',i7,
     1' compmx=',f9.6,'  dxdydz=',f17.14,'  rdxdydz=',4f18.14)
 1002 format ('Progress ',a,a,a,1x,f6.2,'%',2x,'(',i7,'/',i7,')')
 1003 format ('sumT00*V=',e16.8,' sumT0Z*V=',e16.8,' sumpress*V=',
     1 e16.8,' sumN0*V=',e16.8)
 1010 format (' Following are quantities in each cell; to convert to fm^
     1 -3, multiply by',f10.6)
 1020 format (3f21.16)
 1030 format (3i4,i9,f19.15,5f21.15,3f18.15)
 1040 format ('   i   j   k                  N                E[GeV]
     1           P[GeV]              xmom                 ymom
     2      zmom',16x,'u',16x,'v',16x,'w')
 1050 format (' cell coords. ncyc=',i3,3f8.5,'  tnc=',i2,4i3,' b=',f4.1)
 1060 format (i3,'|',20i4)
 1070 format (i3,'|',4x,20i4)
 1080 format (i3,'|',8x,20i4)
 1090 format (i3,'|',12x,20i4)
 1100 format (i3,'|',16x,20i4)
 1110 format ('  energy im=',i3,4f13.6,' |',5f13.6)
 1120 format (i4,2f11.5,f12.7,3f13.5)
 1130 format ('|',15x,'cell sizes [p.u.]: ',3f11.7,15x,'|'/,'|',17x,'cel
     1l sizes [fm]: ',3f11.7,15x,'|')
 1140 format ('|',15x,'ratios pu/fm ',2x,3f11.3,19x,'|'/,'|',21x,'volume
     1 ratio ',f11.7,37x,'|'/,'|',13x,'rdxdydz  ',f11.7,'  rdx=',f11.7,'
     2  dxdydz=',f11.7,12x,'|'/,'|',82x'|')
 1150 format (3i3,i8,9f11.7)
 1160 format (i4,f10.3,3f10.2,3f10.2)
 1170 format ('|',82x,'|'/,'|',19x,'  dtote=',f7.5,' dtodx=',f7.5,' dtot
     1z=',f7.5,20x,'|')
 1190 format ('|',82('_'),'|',/,'|',82x,'|')
 1200 format (i4,2(f10.5,1x),f12.5,9f14.5)
 1210 format (' sme',3f12.3,' captot',f12.5)
 1220 format (8f12.6)
 1230 format (' ERROR in PICREL',4i4,i8,f12.7,i5,1x,5f17.10)
 1240 format ('|',22x,'Input from streak program of V. Magas',23x,'|'/'|
     1',82x,'|',/,'|',29x,'ktnc (# of cells) = ',i3,30x,'|',/'|',82x,'|'
     2 )
 1250 format (27x,'Original Input - Two Nuclei',29x,'|',/,29x,'ktnc (# o
     1f cells) = ',i3,31x,'|',/,83x,'|')
 1260 format ('|',2x,'Ap=',f4.0,' At=',f4.0,' icysum=',i4,' epsilon0='
     1 ,f7.2,2x,'b=',f4.2,' (R1+R2)  ktnc=',i2,' cells ',4x,'|')
 1270 format ('|',19x,'y0=',f9.5,3x,'rho0=',f8.5,5x,'e0=',f8.5,19x,'|'/,
     1'|',17x,'beta=',f10.8,'  Equiv. lab energy =',f8.3,' GeV',17x,'|'/
     2 ,'|',82x,'|')
 1280 format ('# Cell size ','dx= ',f10.8,3x,'dy= ',f10.8,3x,'dz= ',
     1 f10.8,'[fm]',3f10.8)
 1290 format ('|',13x,'Number of cells=',i7,'  Maximum number=',i7,f6.1,
     1'%',15x,'|',/,'|',8x,' Number of particles=',i7,'  Maximum number=
     2',i7,f6.1,'%',15x,'|'/,'|',82x,'|')
 1300 format (3i4,' ',12(E12.5,' '))
 1310 format ('|',5x,'Nuclear radius =',f6.3,'fm  cell: x dimension=',f6
     1 .3,'fm','  z dimen. =',f6.3'fm',5x,'|'/,'|',30x,'cell volume=',f6
     2 .3,'fm^3 ',29x,'|'/,'|',17x,'time for light to travel dz is',f7.4
     3 ,' 10^{-23} s',17x,'|',/,'|',54x,'=',f7.4,' fm/c',15x,'|')
 1320 format ('|',25x,' Velocity interpolation by area',26x,'|')
 1330 format (21x,' Velocity interpolation by M3 (quadratic)',21x,'|')
 1340 format ('|',82x,'|',/,'|',25x,' Quantity interpolation by area',26
     1 x,'|')
 1350 format (/,' Quantity interpolation by M3 (quadratic)')
 1360 format (/,' Quantity interpolation by M4 (cubic)')
 1370 format (' QGP in ',i4,' cells, of avg. Temp.= ',f6.1,'MeV,  at','
     1Cyc=',i3)
 1380 format (' numkt(i)',i2,i6,2i8,' pointers to marker particle storag
     1e indices')
 1390 format ('|',2x,'Equiv. lab energy=',f10.2,' MeV/nucl.=',f9.3,' GeV
     1/nucl. beta=',f12.9,4x,'|'/,83x,'|')
 1400 format ('|',2x,'Equiv. lab energy=',f11.1,' GeV/nucl.=',f8.1,' TeV
     1/nucl.  beta=',f12.10,3x,'|',/,'|',82x,'|')
 1410 format ('|',11x,'Projectile rapidity=',f6.3,2x,'Target rapidity='
     1 ,f6.3,f10.3,11x,'|')
 1420 format ('|',82x,'|'/'|',10x,'No. of marker particles per fluid cel
     1l in different directions',10x,'|',/,'|',22x,'npart=',i2,' npartz=
     2',i2,' nparp=',i2,' nparpz=',i2,23x,'|'/,'|',82x,'|'/18x,'Distribu
     3ting MP-s randomly y/n=(1,0), krandm=',i2,18x,'|',/,'|',23x,'Times
     4tep modification, dtmod=',f6.3,24x,'|'/,'|',82x,'|',/,'|',8x'Outpu
     5t regulating pars., ncycpt=',i3,' ncycrp=',i3,' ncycaf=',i3,'  ksl
     6ab=',i2,7x,'|'/,'|',32x,' ioptos=',5i2,32x,'|'/,'|',82x,'|')
 1430 format ('|',82x,'|',/,'|',10x,'No. of marker particles per fluid c
     1ell determined by density',12x,'|',/,'|',82x,'|'/,'|',15x,'Distrib
     2uting MP-s randomly y/n=(1,0),   porton=',f5.3,15x,'|'/,'|',23x,'T
     3imestep modification dtmod=',f6.3,25x,'|'/,'|',82x,'|',/,'|',8x,'O
     4utput regulating pars., ncycpt=',i3,' ncycrp=',i3,' ncycaf=',i3,'
     5kslab=',i2,8x,'|'/,'|',33x,' ioptos=',5i2,31x,'|'/'|',82x,'|')
 1440 format ('|',82x,'|'/'|',82x,'|'/,'|',5x,'No. of marker particles p
     1er fluid cell in different directions',15x,'|',/,'|',28x,' nparp='
     2 ,i2,' nparpz=',i2,35x,'|'/,'|',82x,'|'/,'|',18x,'Distributing MP-
     3s randomly y/n=(1,0),krandm=',i2,18x,'|'/,'|',5x,'Cell distortion
     4(dz/dx), dzmod=',f6.3,' Timestep modification, dtmod=',f6.3,4x,'|'
     5 /,'|',82x,'|',/,'|',7x,'Output regulating pars., ncycpt=',i3,' nc
     6ycrp=',i3,'  ncycaf=',i3,'  kslab=',i2,7x,'|'/,'|',33x,'ioptos=',5
     7 i2,32x,'|'/'|',82x,'|')
 1450 format ('# Time increment dt= ',2f12.9,' [fm/c] ',3i5,2f13.10)
 1460 format (' Laboratory Increment is ',E9.3,' 10(-23) sec'/)
 1470 format ('|',82x,'|'/'|',21x,'Parameters of the projectile and targ
     1et',22x,'|'/'| gamma baryon/MP  gam(mn-B)/mn   mom./c   v/c    nuc
     2l/cell   capn   #cells   pts. |')
 1480 format ('|',27x,'Two dimensional calculation',28x,'|')
 1490 format ('|',26x,'Three dimensional calculation',27x,'|')
 1500 format ('|',21x,'Initial conditions from Spheroidal Model',21x,'|'
     1 ,/,'|',82x,'|')
 1510 format ('|',26x,'Initial conditions of a Torus',27x,'|',/,83x,'|')
 1520 format (/,' EOS parameters: keqost=',i2,'  hxieq1=',f12.9,'
     1eusual=',f12.8)
 1530 format (/,47x,'pbar annihilation',/)
 1540 format (' maximum no. of cells exceeded. ibar=',i3,', jbar=',i3,'
     1, kbar=',i3,'  ijkba3=',2i8)
 1550 format (20a4,3x,'|',/,83x,'|')
 1560 format (' Ap=',f8.3,'  At=',f8.3,'  Ep=',f10.6,' GeV/nuc.   Et='
     1 ,f10.6,' GeV/nuc.',9x,'|'/' tnc (target semi-major axis, in cells
     2)=',f4.0,'  impact parameter, (0-1.)=',f4.2,8x,'|'/,' Eq. of state
     3 (EOS),   keqost=',i3,50x,'|'/,' Number of cycles (time steps) icy
     4sum=',i4,41x,'|'/,' eta=      ',f9.3,'  kappa=',f11.1,2x,' comprs=
     5',f9.4,'  bcoef=',f12.9,5x,'|',/,10x,'options:',4i3,53x,'|'/,83x,'
     6|')
 1570 format ('|',82x,'|',/,'|',26x,'ibar=',i3,'  jbar=',i3,'  kbar=',i3
     1 ,28x,'|',/,'|',10x,'no. of x-z particles=',i6,'  total no. of par
     2ticles=',i8,12x,'|'/,'|',15x,' xgrid=',f6.2,' fm','   ygrid=',f6.2
     3 ,' fm   z grid=',f6.2,'fm',15x,'|',/,'|',82x,'|')
 1580 format (/'t=',f6.3,'[T] or',f7.3,'[fm/c]  cyc=',i4,/' x=',e11.5
     1 ,' z=',e11.5,' E=',e11.5,' err=',e10.4,1x,'cmprs',2f5.2,
     2 ' (',i2,i3,i4')')
 1590 format (' t=',f6.4,'[T] or',f7.4,'[fm/c]  cyc=',i4,/'   x=',e12.6,
     1' z=',e11.5,' E=',e13.7,' err=',e10.4,1x,'cmprs',2f6.2,' num',i2
     2 ,i6,2i5,' (',2i3,i4,')')
 1600 format ('|',19x,'R_p=',f5.2,' R_t=',f5.2,' R_tot=',f5.2,'  b=',f5.
     1 2,'[fm]',19x,'|',/'|',21x,'Vol_t=',f9.3,'  Vol_p=',f9.3,' [fm^3]'
     2 ,22x,'|')
 1610 format ('|',22x,'  R_p=',f5.2,' [fm]   Vol_p=',f9.3,'fm^3',22x,'|'
     1 )
      end
c c.......................................................................
c
c        subroutine oschydro_header(thydro_start,itmax,number)
c c
c        implicit none
c
c       integer Ap, At, Zp, Zt, npart, nbar, nmes, ctag
c       integer nsteps,ranseed,event,eos,dectag,uid_cnt
c       integer NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
c       real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm
c c 7 integer
c c
c c       common /sys/ npart, nbar, nmes, ctag,nsteps,uid_cnt,
c c      +             ranseed,event,Ap,At,Zp,Zt,eos,dectag,
c c      +             NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
c c       common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm
c c
c c
c       character*12 d12_1,d12_2,d12_3
c       character*80 outputline
c       real*8 e0,n0,Bag,hc3,pi2
c       real*8 dx,dt,t,vol,cut
c       real*8 thydro_start
c       real*8 ecrit
c       integer itmax
c       integer number
c c
c c       common /grstate/ e0,n0,Bag,hc3,pi2
c c       common /gitter/ dx,dt,t,vol,cut
c c
c c chp header for OSCAR hydro output
c  920  format(3a12)
c  921  format(a80)
c  922  format(i6,6i4)
c  923  format(8f10.3)
c  924  format(a,f5.1,a)
c  925  format(a,f9.4,a,f9.4,a,f9.4)
c  926  format(a,f9.4,a,f9.4)
c  927  format(a,f9.4,a,f9.4,a,f9.4,a,f9.4)
c  928  format(a,2i5,a,2i5,a,f8.2,a)
c c
c       d12_1='OSCAR2008H '
c       d12_2='ideal      '
c       d12_3='final_hs'
c       write(number,920) d12_1,d12_2,d12_3
c
c c  INIT
c
c       write(outputline,928) 'INIT: UrQMD,',Ap,Zp,' on',At,Zt,
c      1              ' at E_CM= ',ecm,' AGEV'
c       write(number,921) outputline
c       write(outputline,925)'INIT: b= ',bimp,' n_B,0= ',n0,' e_0= ',
c      1 e0*1.d-3
c       write(number,921) outputline
c       write(outputline,926)'INIT: start time',thydro_start,
c     1     ' iso-energy density FO crit.',CTParam(64)*e0*1.d-3
c       write(number,921) outputline
c
c c  EOS
c
c       outputline='EOS: Lattice Eos from chiral model'
c       write(number,921) outputline
c
c c  CHARGES
c
c       outputline='CHARGES: baryon'
c       write(number,921) outputline
c       write(outputline,926)'HYPER: iso-energy density FO crit.',
c      1 ecrit*1.d-3
c       write(number,921) outputline
c
c c  GEOM
c
c       outputline='GEOM: 3d-cart'
c       write(number,921) outputline
c
c c  GRID
c
c c          outputline='GRID: Euler'
c c          write(number,921) outputline
c c          write(number,922) itmax,int(ngr),int(ngr),
c c      &    int(ngr),1,0,0
c c          write(number,923) thydro_start+t,thydro_start+t+itmax*dt,
c c      &                 (float(0)-float(ngr)*.5-.5)*sngl(dx),
c c      &                 (float(ngr)-float(ngr)*.5-.5)*sngl(dx),
c c      &                 (float(0)-float(ngr)*.5-.5)*sngl(dx),
c c      &                 (float(ngr)-float(ngr)*.5-.5)*sngl(dx),
c c      &                 (float(0)-float(ngr)*.5-.5)*sngl(dx),
c c      &                 (float(ngr)-float(ngr)*.5-.5)*sngl(dx)
c c
c c
c c  VISCOSITY
c
c       outputline='VISCOSITY: none'
c       write(number,921) outputline
c
c c c COMM
c c          outputline='COMM: SHASTA by D. Rischke'
c c          write(number,921) outputline
c c          outputline='COMM: modified by H. Petersen & M. Bleicher'
c c          write(number,921) outputline
c c          outputline='COMM: axis order, z-y-x'
c c          write(number,921) outputline
c
c c  END
c
c          outputline='END_OF_HEADER'
c          write(number,921) outputline
c c
c c
c        return
c        end

c..................................................................


       subroutine oschyper(number,t,x,y,z,e,p,n,vx,vy,vz,
     1                       temp,chem,schem,dst,dsx,dsy,dsz)

       integer number
       real*8 t,x,y,z,temp,chem,schem
       real*8 dst,dsx,dsy,dsz
       real*8 vx,vy,vz,e,p,n

       real*8 e0,n0,Bag,hc3,pi2
       common /grstate/ e0,n0,Bag,hc3,pi2

 997   format(6F11.4,2F11.6,3E24.16,3F11.6,4F11.6)


        write(number,997) t,x,y,z,e*e0*1.d-3,p*e0*1.d-3,
     1       temp*1.d-3,0.0d0,vx,vy,vz,n*n0,chem*1.d-3,schem*1.d-3,
     2       dst,dsx,dsy,dsz

      return
      end
c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function temptab(e,n)
      implicit none

      real*8 temptab
      real*8 p_chiral_eos_big  (0:2000,0:400),
     &       t_chiral_eos_big  (0:2000,0:400),
     &       mub_chiral_eos_big(0:2000,0:400),
     &	     s_chiral_eos_big  (0:2000,0:400),
     &       mus_chiral_eos_big(0:2000,0:400),
     &	     e_chiral_eos_big  (0:2000),
     &	     n_chiral_eos_big  (0:400)
      real*8 p_chiral_eos_small  (0:200,0:200),
     &       t_chiral_eos_small  (0:200,0:200),
     &       mub_chiral_eos_small(0:200,0:200),
     &       s_chiral_eos_small  (0:200,0:200),
     &       mus_chiral_eos_small(0:200,0:200),
     &	     e_chiral_eos_small  (0:200),
     &	     n_chiral_eos_small  (0:200)
      real*8 wes(0:1),wns(0:1),we(0:1),wn(0:1)
      real*8 e,n,e0,n0
      real*8  emin_big ,emax_big ,nmin_big ,nmax_big,
     &        emin_small,emax_small,nmin_small,nmax_small
      integer inn0,iee0
      integer Nes,Nns,Ne,Nn,ies,ins,ieb,inb,je,jn
      real*8  des,dns,de,dn,ems,nms,em,nm
      integer counting1,counting2
      common /chiraleos/ p_chiral_eos_big    ,t_chiral_eos_big  ,
     & 			 mub_chiral_eos_big  ,s_chiral_eos_big  ,
     &			 mus_chiral_eos_big  ,e_chiral_eos_big  ,
     &			 n_chiral_eos_big    ,
     &                   p_chiral_eos_small  ,t_chiral_eos_small,
     &			 mub_chiral_eos_small,s_chiral_eos_small,
     &                   mus_chiral_eos_small,e_chiral_eos_small,
     &                   n_chiral_eos_small  ,
     &			 emin_big ,emax_big ,nmin_big ,nmax_big,
     &        	         emin_small,emax_small,nmin_small,nmax_small
      if (e < 0.) then
      	  temptab = 0
      	  goto 10
      endif
      Ne  = 2001
      Nn  = 401
      Nes = 201
      Nns = 201
      if (e < 1460 .and. n < 0.3)then
      	 des = (emax_small - emin_small) / dfloat(Nes - 1)
         dns = (nmax_small - nmin_small) / dfloat(Nns - 1)
         ies = ((e - emin_small) / des)
         ins = ((n - nmin_small) / dns)
         if (ies < 0) then
             ies = 0
         endif
         if (ins < 0) then
             ins = 0
         endif
         if (ies > Nes - 2) then
             ies = Nes - 2
         endif
         if (ins > Nns - 2) then
             ins = Nns - 2
         endif
         ems = e - emin_small - dfloat(ies) * des
         nms = n - nmin_small - dfloat(ins) * dns
         wes(0) = 1. - ems/des
         wes(1) =      ems/des
         wns(0) = 1. - nms/dns
         wns(1) =      nms/dns
         temptab = 0.
         counting1=0
         do je = 0,1,1
            do jn = 0,1,1
                counting1 = counting1 + 1
                temptab = temptab + wes(je) * wns(jn) *
     1                    t_chiral_eos_small(ies + je,ins + jn)
            enddo
        enddo
        if (temptab < 0.)then
            temptab = 0.
        endif
        goto 10
      else if (e < 146000 .and. n < 6) then
        de  = (emax_big - emin_big) / dfloat(Ne - 1 )
        dn  = (nmax_big - nmin_big) / dfloat(Nn - 1 )
        ieb = ((e - emin_big) / de);
        inb = ((n - nmin_big) / dn);
        if (ieb  < 0) then
            ieb  = 0
        endif
        if (inb  < 0) then
            inb  = 0
        endif
        if (ieb  > Ne - 2 ) then
            ieb  = Ne - 2
        endif
        if (inb  > Nn - 2 ) then
            inb  = Nn - 2
        endif
        em = e - emin_big - dfloat(ieb) * de
        nm = n - nmin_big - dfloat(inb) * dn
        we(0) = 1. - em/de
        we(1) =      em/de
        wn(0) = 1. - nm/dn
        wn(1) =      nm/dn
        temptab = 0.
        do je = 0,1,1
            do jn = 0,1,1
                temptab = temptab + we(je) * wn(jn) *
     &                    t_chiral_eos_big(ieb + je,inb + jn)
            enddo
        enddo
        if (temptab < 0.)then
            temptab = 0.
        endif
        goto 10
      else
        temptab = 0.15120476935 * (e**0.25)
        goto 10
      endif

   10 continue
      return
      end
c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function chemtab(e,n)
      implicit none
c
c     This function-subprogram determines the baryon chemical potential
c     of the underlying EoS.
c     It is used in subroutines fileo and entro.
c
c     Type declarations for variables in common-blocks.
c
      real*8 chemtab
      real*8 p_chiral_eos_big  (0:2000,0:400),
     &       t_chiral_eos_big  (0:2000,0:400),
     &       mub_chiral_eos_big(0:2000,0:400),
     &	     s_chiral_eos_big  (0:2000,0:400),
     &       mus_chiral_eos_big(0:2000,0:400),
     &	     e_chiral_eos_big  (0:2000),
     &	     n_chiral_eos_big  (0:400)
      real*8 p_chiral_eos_small  (0:200,0:200),
     &       t_chiral_eos_small  (0:200,0:200),
     &       mub_chiral_eos_small(0:200,0:200),
     &       s_chiral_eos_small  (0:200,0:200),
     &       mus_chiral_eos_small(0:200,0:200),
     &	     e_chiral_eos_small  (0:200),
     &	     n_chiral_eos_small  (0:200)
      real*8 wes(0:1),wns(0:1),we(0:1),wn(0:1)
      real*8 e,n,e0,n0
      real*8  emin_big ,emax_big ,nmin_big ,nmax_big,
     &        emin_small,emax_small,nmin_small,nmax_small
      integer inn0,iee0
      integer Nes,Nns,Ne,Nn,ies,ins,ieb,inb,je,jn
      real*8  des,dns,de,dn,ems,nms,em,nm
      integer counting1,counting2
      common /chiraleos/ p_chiral_eos_big    ,t_chiral_eos_big  ,
     & 			 mub_chiral_eos_big  ,s_chiral_eos_big  ,
     &			 mus_chiral_eos_big  ,e_chiral_eos_big  ,
     &			 n_chiral_eos_big    ,
     &                   p_chiral_eos_small  ,t_chiral_eos_small,
     &			 mub_chiral_eos_small,s_chiral_eos_small,
     &                   mus_chiral_eos_small,e_chiral_eos_small,
     &                   n_chiral_eos_small  ,
     &			 emin_big ,emax_big ,nmin_big ,nmax_big,
     &        	         emin_small,emax_small,nmin_small,nmax_small
      n0 = 0.15d0 ! fm-3
      e0 = 146d0 ! MeV/fm3
c       e = e * e0
c       n = n * n0
      if (e < 0.) then
      	  chemtab = 0
      	  goto 10
      endif
      Ne  = 2001
      Nn  = 401
      Nes = 201
      Nns = 201
      if (e < 1460 .and. n < 0.3)then
      	 des = (emax_small - emin_small) / dfloat(Nes - 1)
         dns = (nmax_small - nmin_small) / dfloat(Nns - 1)
         ies = ((e - emin_small) / des)
         ins = ((n - nmin_small) / dns)
         if (ies < 0) then
             ies = 0
         endif
         if (ins < 0) then
             ins = 0
         endif
         if (ies > Nes - 2) then
             ies = Nes - 2
         endif
         if (ins > Nns - 2) then
             ins = Nns - 2
         endif
         ems = e - emin_small - dfloat(ies) * des
         nms = n - nmin_small - dfloat(ins) * dns
         wes(0) = 1. - ems/des
         wes(1) =      ems/des
         wns(0) = 1. - nms/dns
         wns(1) =      nms/dns
         chemtab = 0.
         do je = 0,1,1
            do jn = 0,1,1
                chemtab = chemtab + wes(je) * wns(jn) *
     1                     mub_chiral_eos_small(ies + je,ins + jn)
            enddo
        enddo
        goto 10
      else if (e < 146000 .and. n < 6) then
        de  = (emax_big - emin_big) / dfloat(Ne - 1 )
        dn  = (nmax_big - nmin_big) / dfloat(Nn - 1 )
        ieb = ((e - emin_big) / de);
        inb = ((n - nmin_big) / dn);
        if (ieb  < 0) then
            ieb  = 0
        endif
        if (inb  < 0) then
            inb  = 0
        endif
        if (ieb  > Ne - 2 ) then
            ieb  = Ne - 2
        endif
        if (inb  > Nn - 2 ) then
            inb  = Nn - 2
        endif
        em = e - emin_big - dfloat(ieb) * de
        nm = n - nmin_big - dfloat(inb) * dn
        we(0) = 1. - em/de
        we(1) =      em/de
        wn(0) = 1. - nm/dn
        wn(1) =      nm/dn
        chemtab = 0.
        do je = 0,1,1
            do jn = 0,1,1
                chemtab = chemtab + we(je) * wn(jn) *
     &                     mub_chiral_eos_big(ieb + je,inb + jn)
            enddo
        enddo
        goto 10
      else
        chemtab = 0.
        goto 10
      endif

   10 continue
      return
      end

c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function presstab(e,n)
      implicit none

      real*8 presstab
      real*8 p_chiral_eos_big  (0:2000,0:400),
     &       t_chiral_eos_big  (0:2000,0:400),
     &       mub_chiral_eos_big(0:2000,0:400),
     &	     s_chiral_eos_big  (0:2000,0:400),
     &       mus_chiral_eos_big(0:2000,0:400),
     &	     e_chiral_eos_big  (0:2000),
     &	     n_chiral_eos_big  (0:400)
      real*8 p_chiral_eos_small  (0:200,0:200),
     &       t_chiral_eos_small  (0:200,0:200),
     &       mub_chiral_eos_small(0:200,0:200),
     &       s_chiral_eos_small  (0:200,0:200),
     &       mus_chiral_eos_small(0:200,0:200),
     &	     e_chiral_eos_small  (0:200),
     &	     n_chiral_eos_small  (0:200)
      real*8 wes(0:1),wns(0:1),we(0:1),wn(0:1)
      real*8 e,n,e0,n0
      real*8  emin_big ,emax_big ,nmin_big ,nmax_big,
     &        emin_small,emax_small,nmin_small,nmax_small
      integer inn0,iee0
      integer Nes,Nns,Ne,Nn,ies,ins,ieb,inb,je,jn
      real*8  des,dns,de,dn,ems,nms,em,nm
      integer counting1,counting2
      common /chiraleos/ p_chiral_eos_big    ,t_chiral_eos_big  ,
     & 			 mub_chiral_eos_big  ,s_chiral_eos_big  ,
     &			 mus_chiral_eos_big  ,e_chiral_eos_big  ,
     &			 n_chiral_eos_big    ,
     &                   p_chiral_eos_small  ,t_chiral_eos_small,
     &			 mub_chiral_eos_small,s_chiral_eos_small,
     &                   mus_chiral_eos_small,e_chiral_eos_small,
     &                   n_chiral_eos_small  ,
     &			 emin_big ,emax_big ,nmin_big ,nmax_big,
     &        	         emin_small,emax_small,nmin_small,nmax_small
c      write(*,*) 'e[e0]= ',e,' n[n0]= ',n
      n0 = 0.15d0 ! fm-3
      e0 = 146d0 ! MeV/fm3
c       e = e * e0
c       n = n * n0
c      write(*,*) 'e[MeV/fm3]= ',e,' n[fm-3]= ',n
c      do inn0 = 0,200,1
c         do iee0 = 0,200,1
c            write(*,*) p_chiral_eos_big  (iee0,inn0)/e0,' ',
c     &       	       p_chiral_eos_small(iee0,inn0)/e0
c         enddo
c      enddo
c      stop 5033
c      write(*,*) emax_big  ,' ',emin_big  ,' ',nmax_big  ,' ',
c     &           nmin_big  ,' ',emax_small,' ',emin_small,' ',
c     &           nmax_small,' ',nmin_small
c      stop 5037
      if (e < 0.) then
      	  presstab = 0
      	  goto 10
      endif
      Ne  = 2001
      Nn  = 401
      Nes = 201
      Nns = 201
      if (e < 1460 .and. n < 0.3)then
c         write(*,*) 'first if statement'
      	 des = (emax_small - emin_small) / dfloat(Nes - 1)
         dns = (nmax_small - nmin_small) / dfloat(Nns - 1)
c         write(*,*) 'des= ',des,' dns= ',dns
         ies = ((e - emin_small) / des)
         ins = ((n - nmin_small) / dns)
c         write(*,*) 'ies= ',ies,' ins= ',ins
         if (ies < 0) then
             ies = 0
         endif
         if (ins < 0) then
             ins = 0
         endif
         if (ies > Nes - 2) then
             ies = Nes - 2
         endif
         if (ins > Nns - 2) then
             ins = Nns - 2
         endif
         ems = e - emin_small - dfloat(ies) * des
         nms = n - nmin_small - dfloat(ins) * dns
c         write(*,*) 'ems= ',ems,' nms= ',nms
         wes(0) = 1. - ems/des
         wes(1) =      ems/des
         wns(0) = 1. - nms/dns
         wns(1) =      nms/dns
c         write(*,*) 'wes(0)= ',wes(0),' wes(1)= ',wes(1)
c         write(*,*) 'wns(0)= ',wns(0),' wns(1)= ',wns(1)
         presstab = 0.
c         counting1=0
         do je = 0,1,1
            do jn = 0,1,1
c                counting1 = counting1 + 1
                presstab = presstab + wes(je) * wns(jn) *
     1                     p_chiral_eos_small(ies + je,ins + jn)
            enddo
        enddo
        if (presstab < 0.) then
            presstab = 0.
        endif
c        write(*,*) 'counting1= ',counting1,'presstab= ', presstab
        goto 10
      else if (e < 146000 .and. n < 6) then
c        write(*,*) 'e[MeV/fm3]= ',e,' n[fm-3]= ',n
c        write(*,*) 'second if statement'
        de  = (emax_big - emin_big) / dfloat(Ne - 1 )
        dn  = (nmax_big - nmin_big) / dfloat(Nn - 1 )
c        write(*,*) 'de= ',de,' dn= ',dn
        ieb = ((e - emin_big) / de);
        inb = ((n - nmin_big) / dn);
c        write(*,*) 'ieb= ',ieb,' inb= ',inb
        if (ieb  < 0) then
            ieb  = 0
        endif
        if (inb  < 0) then
            inb  = 0
        endif
        if (ieb  > Ne - 2 ) then
            ieb  = Ne - 2
        endif
        if (inb  > Nn - 2 ) then
            inb  = Nn - 2
        endif
        em = e - emin_big - dfloat(ieb) * de
        nm = n - nmin_big - dfloat(inb) * dn
c        write(*,*) 'em= ',em,' nm= ',nm
        we(0) = 1. - em/de
        we(1) =      em/de
        wn(0) = 1. - nm/dn
        wn(1) =      nm/dn
c        write(*,*) 'we(0)= ',we(0),' we(1)= ',we(1)
c        write(*,*) 'wn(0)= ',wn(0),' wn(1)= ',wn(1)
        presstab = 0.
c        counting2 = 0
        do je = 0,1,1
            do jn = 0,1,1
c                counting2 = counting2 + 1
                presstab = presstab + we(je) * wn(jn) *
     &                     p_chiral_eos_big(ieb + je,inb + jn)
            enddo
        enddo
        if (presstab < 0.) then
            presstab = 0.
        endif
c        write(*,*) 'counting2= ',counting2,'presstab= ', presstab
c        stop 5134
        goto 10
      else
        presstab = 0.2964 * e
c        write(*,*) 'e= ',e,'p= ',presstab
c        stop 5140
        goto 10
      endif

   10 continue
c      write(*,*) 'presstab= ', presstab
c      write(*,*) 'e= ',e,'n= ',n,'presstab= ', presstab
c      stop 5122
      return
      end
c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine read tabulated eos from file:
c     Files:
c     References:
      subroutine read_chiral_eos()
      implicit none 
      real*8 p_chiral_eos_big  (0:2000,0:400),
     &       t_chiral_eos_big  (0:2000,0:400),
     &       mub_chiral_eos_big(0:2000,0:400),
     &	     s_chiral_eos_big  (0:2000,0:400),
     &       mus_chiral_eos_big(0:2000,0:400),
     &	     e_chiral_eos_big  (0:2000),
     &	     n_chiral_eos_big  (0:400)
      real*8 p_chiral_eos_small  (0:200,0:200),
     &       t_chiral_eos_small  (0:200,0:200),
     &       mub_chiral_eos_small(0:200,0:200),
     &       s_chiral_eos_small  (0:200,0:200),
     &       mus_chiral_eos_small(0:200,0:200),
     &	     e_chiral_eos_small  (0:200),
     &	     n_chiral_eos_small  (0:200)
      integer inn0,iee0,j,Ne,Nes,Nn,Nns
      real*8  t,mub,e,p,n,s,mus,dummy
      real*8  e0,n0
      real*8  emin_big ,emax_big ,nmin_big ,nmax_big,
     &        emin_small,emax_small,nmin_small,nmax_small

      common /chiraleos/ p_chiral_eos_big    ,t_chiral_eos_big  ,
     & 			 mub_chiral_eos_big  ,s_chiral_eos_big  ,
     &			 mus_chiral_eos_big  ,e_chiral_eos_big  ,
     &			 n_chiral_eos_big    ,
     &                   p_chiral_eos_small  ,t_chiral_eos_small,
     &			 mub_chiral_eos_small,s_chiral_eos_small,
     &                   mus_chiral_eos_small,e_chiral_eos_small,
     &                   n_chiral_eos_small  ,
     &			 emin_big ,emax_big ,nmin_big ,nmax_big,
     &        	         emin_small,emax_small,nmin_small,nmax_small

      open(unit=51,file='../eos/chiraleos.dat')
      open(unit=52,file='../eos/chiralsmall.dat')

      n0 = 0.15d0 ! fm-3
      e0 = 146d0 ! MeV/fm3

C     T [MeV], mu_q [MeV], energy density [e_0], pressure [e_0], baryon
C     density [n_0], entropy density [n_0], mu_S [MeV], placeholder (not
c     important).
c
      do 1109 inn0 = 0,400,1
         do 1110 iee0 = 0,2000,1
            read(51,*) t,mub,e,p,n,s,mus,dummy
            p_chiral_eos_big  (iee0,inn0) = p * e0
            t_chiral_eos_big  (iee0,inn0) = t
            mub_chiral_eos_big(iee0,inn0) = mub
            s_chiral_eos_big  (iee0,inn0) = s * n0
            mus_chiral_eos_big(iee0,inn0) = mus
            e_chiral_eos_big  (iee0)      = e
            n_chiral_eos_big  (inn0)      = n
c            write(*,*) t_chiral_eos_big  (iee0,inn0),' ',
c     &       	       mub_chiral_eos_big(iee0,inn0),' ',
c     &       	       e_chiral_eos_big  (iee0)     ,' ',
c     &       	       p_chiral_eos_big  (iee0,inn0),' ',
c     &       	       n_chiral_eos_big  (inn0)     ,' ',
c     &                 s_chiral_eos_big  (iee0,inn0),' ',
c     &                 mus_chiral_eos_big(iee0,inn0),' ',
c     &                 dummy
 1110    continue
 1109 continue
      do 1209 inn0 = 0,200,1
         do 1210 iee0 = 0,200,1
            read(52,*) t,mub,e,p,n,s,mus,dummy
            p_chiral_eos_small  (iee0,inn0) = p * e0
            t_chiral_eos_small  (iee0,inn0) = t
            mub_chiral_eos_small(iee0,inn0) = mub
            s_chiral_eos_small  (iee0,inn0) = s * n0
            mus_chiral_eos_small(iee0,inn0) = mus
            e_chiral_eos_small  (iee0)      = e
            n_chiral_eos_small  (inn0)      = n
c            write(*,*) t_chiral_eos_small  (iee0,inn0),' ',
c     &       	       mub_chiral_eos_small(iee0,inn0),' ',
c     &       	       e_chiral_eos_small  (iee0)     ,' ',
c     &       	       p_chiral_eos_small  (iee0,inn0),' ',
c     &       	       n_chiral_eos_small  (inn0)     ,' ',
c     &                 s_chiral_eos_small  (iee0,inn0),' ',
c     &                 mus_chiral_eos_small(iee0,inn0)
 1210    continue
 1209 continue
      close(51)
      close(52)
      Ne  = 2001
      Nn  = 401
      Nes = 201
      Nns = 201
      emin_big   = e_chiral_eos_big  (0)       * e0;
      emax_big   = e_chiral_eos_big  (Ne - 1)  * e0;
      nmin_big   = n_chiral_eos_big  (0)       * n0;
      nmax_big   = n_chiral_eos_big  (Nn - 1)  * n0;
      emin_small = e_chiral_eos_small(0)       * e0;
      emax_small = e_chiral_eos_small(Nes - 1) * e0;
      nmin_small = n_chiral_eos_small(0)       * n0;
      nmax_small = n_chiral_eos_small(Nns - 1) * n0;
c      write(*,*) emax_big  ,' ',emin_big  ,' ',nmax_big  ,' ',
c     &           nmin_big  ,' ',emax_small,' ',emin_small,' ',
c     &           nmax_small,' ',nmin_small
      return
      end


      subroutine fillcube(i,j,k,nhyperstep,elast,e,HyperCube)
c      include 'defs.f'      
      
      real*8 e(1:300,1:300,1:700),elast(1:300,1:300,1:700)
      REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1,0:1) :: HyperCube
      integer nhyperstep  
      

      HyperCube(0,0,0,0)=elast(i,j,k)
      HyperCube(0,0,0,1)=elast(i,j,k+nhyperstep)
      HyperCube(0,0,1,0)=elast(i,j+nhyperstep,k) 
      HyperCube(0,0,1,1)=elast(i,j+nhyperstep,k+nhyperstep) 
      HyperCube(0,1,0,0)=elast(i+nhyperstep,j,k) 
      HyperCube(0,1,0,1)=elast(i+nhyperstep,j,k+nhyperstep) 
      HyperCube(0,1,1,0)=elast(i+nhyperstep,j+nhyperstep,k) 
      HyperCube(0,1,1,1)=elast(i+nhyperstep,j+nhyperstep,k+nhyperstep) 
      HyperCube(1,0,0,0)=e(i,j,k)
      HyperCube(1,0,0,1)=e(i,j,k+nhyperstep)
      HyperCube(1,0,1,0)=e(i,j+nhyperstep,k) 
      HyperCube(1,0,1,1)=e(i,j+nhyperstep,k+nhyperstep) 
      HyperCube(1,1,0,0)=e(i+nhyperstep,j,k) 
      HyperCube(1,1,0,1)=e(i+nhyperstep,j,k+nhyperstep) 
      HyperCube(1,1,1,0)=e(i+nhyperstep,j+nhyperstep,k) 
      HyperCube(1,1,1,1)=e(i+nhyperstep,j+nhyperstep,k+nhyperstep) 




      return 
      end
c    ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      FUNCTION Quadrilinear(Vmid,HyperCube,dt,dx,dy,dz)

      IMPLICIT NONE
      REAL(KIND(0D0)) :: Quadrilinear
      REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1,0:1) :: HyperCube
      REAL(KIND(0D0)),DIMENSION(0:3) :: Vmid
      REAL(KIND(0D0)) :: dt, dx, dy, dz, F
      REAL(KIND(0D0)) :: t, x, y, z
      INTEGER :: n,i,j,k

      Quadrilinear = 0D0
      t = Vmid(0)/dt
      x = Vmid(1)/dx
      y = Vmid(2)/dy
      z = Vmid(3)/dz

      DO k = 0,1
      DO j = 0,1
      DO i = 0,1
      DO n = 0,1
      Quadrilinear = Quadrilinear + (n*t+(1-n)*(1-t))*(i*x+(1-i)*(1-x)) 
     1 *(j*y+(1-j)*(1-y))*(k*z+(1-k)*(1-z))*HyperCube(n,i,j,k)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      END FUNCTION Quadrilinear









c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phase3 (kdisp,zpmin,zpmax,zpminz,zpmaxz,iexpnd)
c
c   Calculates rest-frame quantities
c
      implicit real*8 (a-h,o-z)
c       real*4 aa1, u, v, w
      integer*2 koost
      parameter (EFLOOR=1.d0, BFLOOR=1.d-4)
c     
      common /chktsm/ sumchem,sumtemp,sumentro,nlpc
      common /hxi33/ hxi3a, xi, hxi
      common /redund/ a, b, c
      common /other/ eusual, epn0, t0, epion, ppion, edelt, delta, pdelt
     1 , temp, tempkn, eps, xth, epsc, dnsmax, ncycpt, ncycrp, ncycaf
      common /prstbl/ prs(99,102), trs(99,102), pres, tempzz, chem,
     1 amurs(99,102), entrs(99,102), entro, ipres
      common /cnuclr/ hvc, hv2, pi0, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
      common /gamsol/ enzz, bbagzz, rrr, amom, bbagde
      common /eqqcd/ tq, eq1, pq, sq, anq1, amuq
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
      common /midle/ imid, kmid, kend, imidpt, kmidpt, klmin, krmax,
     1 mkdisp, ipmin, ipmax, imid1
      common /pscm2/ cycen, istop, isume
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
      common /dmisca/ adens(300,700), eedns(300,700), aplot(300,700),
     1 bplot(300,700), cplot(300,700)
      common /misca/ kmove(700), numlvv(52)
c
c   numfz1 changed to numcli, 4/7/22
c
      common /numcmm/ numbk, numcel, numcst, numext, numcli, numlv,
     1 numout, numouu, numou2, jdim1, kfreez, kfrzcy, press, freez,
     2 compmy, koast, kext, cpntot, compmx, compe,iundee(58)
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c      
      
      common /posevent/ enerdensrfpos(1:300,1:300,1:700),enerdensrfpost(
     & 1:300,1:300,1:700),bardensrfpos(1:300,1:300,1:700),bardensrfpost(
     & 1:300,1:300,1:700),enerdenslfpos(1:300,1:300,1:700),bardenslfpos(
     & 1:300,1:300,1:700),vxpos(1:300,1:300,1:700),vxpost(1:300,1:300,1:
     & 700),vypos(1:300,1:300,1:700),vypost(1:300,1:300,1:700),vzpos(1:3
     & 00,1:300,1:700),vzpost(1:300,1:300,1:700),presspos(1:300,1:300,1:
     & 700),xmomdenslfpos(1:300,1:300,1:700),ymomdenslfpos(1:300,1:300,1
     & :700),zmomdenslfpos(1:300,1:300,1:700),posN0(1:300,1:300,1:700),p
     & osNx(1:300,1:300,1:700),posNy(1:300,1:300,1:700),posNz(1:300,1:30
     & 0,1:700),xmubpos(1:300,1:300,1:700),temppos(1:300,1:300,1:700)
     
      common /prevent/ enerdensrfpre(1:300,1:300,1:700),enerdensrfprev(1
     & :300,1:300,1:700),bardensrfpre(1:300,1:300,1:700),bardensrfprev(1
     & :300,1:300,1:700),enerdenslfpre(1:300,1:300,1:700),bardenslfpre(1
     & :300,1:300,1:700),vxpre(1:300,1:300,1:700),vxprev(1:300,1:300,1:7
     & 00),vypre(1:300,1:300,1:700),vyprev(1:300,1:300,1:700),vzpre(1:30
     & 0,1:300,1:700),vzprev(1:300,1:300,1:700),xmomdenslfpre(1:300,1:30
     & 0,1:700),ymomdenslfpre(1:300,1:300,1:700),zmomdenslfpre(1:300,1:3
     & 00,1:700),presspre(1:300,1:300,1:700),preNx(1:300,1:300,1:700),pr
     & eNy(1:300,1:300,1:700),preNz(1:300,1:300,1:700),preN0(1:300,1:300
     & ,1:700),xmubpre(1:300,1:300,1:700),temppre(1:300,1:300,1:700)
     
     
      common /timeconv/ dtconv
     
      dimension u(91627233), v(91627233), w(91627233), capn(91627233),p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233), px(85877055), py(85877055), pz(85877055)
c
      common /bucket/ lablbk(1527120), xmomcl(1527120), ymomcl(1527120),
     1 zmomcl(1527120), encell(1527120), capncl(1527120), cnvvl(1527120)
      common /backet/ lbblbk(22), cpnncl(22), convol(22)
      common /bccket/ lablbc(1527120), celvol(1527120)
      common /volbuk/ lblvol(1527120), retvol(1527120)
      common /matric/ storev(929801), idensz(500), edensz(500)
      common /plcm7/ kextrn(30542411), koost(30542411)
      common /newcrp/ entot, dentot, tempav, numcll
c
      dimension ensort(1527120), largi(200)
c
      equivalence (ensort(1),xmomcl(1))
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
c
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c
      external gammf, gammfp, gammh, gamnix, gamnix1, funalt
c

	  write(156,*)' in the phase 3'
      xn0 = 0.15d0 ! fm-3
      ene0 = 146d0 ! MeV/fm3
      fast=0.
      if (ioptos3.ne.0) write (26,*) 'ncyc=',ncyc,'capn,capn*dxdydz/pnu(
     11),xni,  gamma*xni ,   capn/xni/gamma,   capn/xni/gamma'
      dentot=0.
      eq1=0.
      tempav=0.
      xmomot=0.
      numcll=1
      volnrm=1./dfloat(npart*npart*npartz)
      erat=1.! GeV ! OLD 939.5731*0.1453
      bbagzz=bbag
c      bbagde=bbag/1000.  !GeV/fm**3
      bbagde=bbag  ! in MeV/fm**3
      if (kslab.ne.8) call checks
c       if (kdebug.ne.0) write (156,*) ' Begin Phase 3, ncyc=',ncyc,
c      1 ' ___________________'
c
c      sumN0=0
c      sumT00=0
c      sumT0Z=0
c      sumpress=0
c      sumwijk=0
c      do i=1,ibar
c      do j=1,jbar
c      do k=1,kbar
c      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c      sumN0=sumN0+capn(ijk)
c      sumT00=sumT00+e(ijk)
c      sumT0Z=sumT0Z+zmom(ijk)
c      sumpress=sumpress+p(ijk)
c      sumwijk=sumwijk+abs(w(ijk))
c      enddo
c      enddo
c      enddo
c      write(*,*) 'N0phas3=',sumN0,'T00=',sumT00,'sumT0Z=',sumT0Z,
c     1 sumpress*dxdydz,sumwijk

      do 20 i=1,ibar
      do 10 k=1,kbar
      eedns(i,k)=0.
      adens(i,k)=0.
      aplot(i,k)=0.
      cplot(i,k)=0.
   10 continue
   20 continue
      do 50 i=1,ibar
      do 40 j=1,jbar
      do 30 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
   30 continue
   40 continue
   50 continue
      suume=0.
      suumec=0.
      suumx=0.
      suumz=0.
      numn=0
      ees=0.
      eet=0.
      if(nq.eq.0)stop 5068
      do 100 i=1,1527120
      ees=ees+encell(i) ! encell(i)=0 => ees=0
      if (lablbk(i).eq.-1) go to 100
      eet=eet+encell(i) ! eet=0
      numn=numn+1
      j=lablbk(i)
      ijk=(j-1)/nq
      k=ijk/ijbar+1
      ii=ijk-ijbar*(k-1)+1
      jj=ii/ibar+1
      iii=ii-ibar*(jj-1)+1
      ijkn=nq*(ijbar*(k-1)+ibar*(jj-1)+iii-1)+1
      enold=capn(j) ! []
      eeold=e(j) ! [MeV]
      xmom(j)=xmom(j)+xmomcl(i) ! [MeV]
      ymom(j)=ymom(j)+ymomcl(i) ! [MeV]
      zmom(j)=zmom(j)+zmomcl(i) ! [MeV]
      e(j)=e(j)+encell(i) ! encell(i)=0
      capn(j)=capn(j)+capncl(i)
      if (capn(j).gt.0.) go to 90
      if (capn(j).gt.-0.00001*pnu(1)) go to 80
      if (capn(j).lt.0.) write (*,630) i,j,iii,jj,k,ijkn,capn(j),capn(j)
     1 /pnu(1),enold,capncl(i),e(j),eeold,encell(i),pnu(1)
c       if (capn(j).lt.0.) write (156,630) i,j,iii,jj,k,ijkn,capn(j),
c      1 capn(j)/pnu(1),enold,capncl(i),e(j),eeold,encell(i),pnu(1)
      stop 121
   80 capn(j)=0.
   90 b=e(j) ! [MeV]
      bsq=b*b
      c=xmom(j)**2+ymom(j)**2+zmom(j)**2 ! [MeV2]a
      if (capn(ijkn).gt.cpmax) cpmax=capn(ijk)
      suumec=suumec+encell(i)
      suume=suume+capncl(i)
      suumx=suumx+xmomcl(i)
      suumz=suumz+zmomcl(i)
      if (capn(j).gt.BFLOOR) go to 100
      if (dabs(e(j)).gt.EFLOOR) write (*,640) j,e(j),capn(j),i
     1 ,capncl(i),encell(i)
c
      e(j)=0.
      capn(j)=0.
      xmom(j)=0.
      ymom(j)=0.
      zmom(j)=0.
      u(j)=0.
      v(j)=0.
      w(j)=0.
  100 continue
c      sumN0=0
c      sumT00=0
c      sumT0Z=0
c      sumpress=0
c      sumwijk=0
c      do i=1,ibar
c      do j=1,jbar
c      do k=1,kbar
c      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c      sumN0=sumN0+capn(ijk)
c      sumT00=sumT00+e(ijk)
c      sumT0Z=sumT0Z+zmom(ijk)
c      sumpress=sumpress+p(ijk)
c      sumwijk=sumwijk+abs(w(ijk))
c      enddo
c      enddo
c      enddo
c      write(*,*) 'N0phas3=',sumN0,'T00=',sumT00,'sumT0Z=',sumT0Z,
c     1 sumpress*dxdydz,sumwijk
      do 110 i=1,1527120
      lblvol(i)=lablbc(i)
      retvol(i)=celvol(i)
      sumv=sumv+retvol(i)
      lablbc(i)=-1
      celvol(i)=0.
      j=lblvol(i)
c  New
      if (j.le.0) go to 110
      if (capn(j).eq.0.) go to 110
c  Ensure velocities are non-zero if momenta are not
      if (w(j).ne.0..or.u(j).ne.0..or.v(j).ne.0.) go to 110
      if(e(j).eq.0.)stop 5135
      u(j)=xmom(j)/e(j)
      v(j)=ymom(j)/e(j)
      w(j)=zmom(j)/e(j)
      cmom=u(j)**2+v(j)**2+w(j)**2
c       if (cmom.gt.1.) write (156,790) j,u(j),v(j),w(j),cmom
      ensort(i)=0.
  110 continue
c   section removed, 25 June 2013 - see earlier versions
      if (zpmin.ge.zpminz) go to 120
      if (iexpnd.eq.1) go to 120
      if (zpmaxz.gt.zpmax) go to 120
c       write (156,860)
      iexpnd=1
  120 continue
      zpmaxz=zpmax
      zpminz=zpmin
      jtemp=1
      ktemp=1
      if (ktemp.ne.0.and.mod(ncyc,ncycpt).eq.0.and.keqost.le.0) jtemp=1
c  jtemp=1 -> temperatures are calculated in qgp
      ko=0
      kphs=3
      j=1
      smx=0.
      smz=0.
      numtt=0
c      sumN0=0
c      sumT00=0
c      sumT0Z=0
c      sumpress=0
c      sumwijk=0
c      do i=1,ibar
c      do j=1,jbar
c      do k=1,kbar
c      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c      sumN0=sumN0+capn(ijk)
c      sumT00=sumT00+e(ijk)
c      sumT0Z=sumT0Z+zmom(ijk)
c      sumpress=sumpress+p(ijk)
c      sumwijk=sumwijk+abs(w(ijk))
c      enddo
c      enddo
c      enddo
c      write(*,*) 'N0phas3=',sumN0,'T00=',sumT00,'sumT0Z=',sumT0Z,
c     1 sumpress*dxdydz,sumwijk
c
c   For k=kmid, locate for each j, the largest i for which a cell
c      is non-empty
c
      do 130 j=1,jbar
      largi(j)=0
  130 continue
      do 170 k=1,kbar
      do 160 j=1,jbar
      do 150 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (ijk0.gt.ijkmax) write (*,650) i,j,k,ijk0,ijkmax
      if (ijk.gt.3*ijkmax) write (*,*) i,j,k,ijk,ijkmax
      koost(ijk0)=0
      if (capn(ijk).eq.0.) go to 150
      if (k.ne.kmid) go to 140
      if (largi(j).gt.i) go to 140
      largi(j)=i
  140 continue
      num=capn(ijk)/pnu(1)+.1 ! pnu is dpt that is dimensionless
      numtt=numtt+num
      smx=smx+xmom(ijk)
      smz=smz+zmom(ijk)
      capn(ijk)=rdxdydz*capn(ijk) ! NOTE THAT FROM HERE N0 [fm-3] AGAIN
  150 continue
  160 continue
  170 continue
      smx=smx+xmomot
c       if (kdebug.gt.0) write (156,*) ' largi',(largi(j),j=1,25)
      numcli=0
      numcle=0
      numclt=0
      totn=0.
      sume=0.
      numtot=0
      do 200 k=1,kbar
      do 190 j=1,jbar
      do 180 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      kextrn(ijk0)=0
      if (capn(ijk).eq.0.) go to 180
      numtot=numtot+1
c      at this point capn (the density) assumes the value
c       of 1.0000 for normal nuclear density
      e(ijk)=rdxdydz*e(ijk) ! NOTE THAT FROM HERE T00 [MeV*fm-3] AGAIN
      xmom(ijk)=rdxdydz*xmom(ijk) ! NOTE THAT FROM HERE T0X [MeV*fm-3] AGAIN
      ymom(ijk)=rdxdydz*ymom(ijk) ! NOTE THAT FROM HERE T0Y [MeV*fm-3] AGAIN
      zmom(ijk)=rdxdydz*zmom(ijk) ! NOTE THAT FROM HERE T0Z [MeV*fm-3] AGAIN
      sume=sume+e(ijk)
  180 continue
  190 continue
  200 continue
c
c     FROM HERE T00, T0Z, T0Y AND T0X ARE GIVEN IN MeV*fm-3 AND N0 
c                          IN fm-3 AGAIN 
c
c      write (*,*) 'INSIDE SUBROUTIN PHASE3 T00,T0Z,T0Y,T0X AND N0 ARE'
c      write (*,*) 'DIVIDED BY DXDYDZ AGAIN SO T00 [Mev*fm-3],'
c      write (*,*) 'T00 [Mev*fm-3],T00 [Mev*fm-3],T00 [Mev*fm-3],'
c      write (*,*) 'T00 [Mev*fm-3] AND press [Mev*fm-3]' 
c      write(*,*) 'TNC4253',ktnc
c      sumN0=0
c      sumT00=0
c      sumT0Z=0
c      sumpress=0
c      sumwijk=0
c      do i=1,ibar
c      do j=1,jbar
c      do k=1,kbar
c      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c      sumN0=sumN0+capn(ijk)
c      sumT00=sumT00+e(ijk)
c      sumT0Z=sumT0Z+zmom(ijk)
c      sumpress=sumpress+p(ijk)
c      sumwijk=sumwijk+abs(w(ijk))
c      enddo
c      enddo
c      enddo     
c      write(*,*) 'N0phas3=',sumN0*dxdydz,'T00=',sumT00*dxdydz,
c     1 'sumT0Z=',sumT0Z*dxdydz,sumpress*dxdydz,sumwijk
c      if (kdebug.ne.0) then
c        write (15,810) ncyc
c        write (16,810) ncyc
c      endif
      delvmx=0.
      captot=0.
      eptot=0.
      captoe=0.
      eptote=0.
      do 420 iphs3=1,2
      kstrz=0
      ko=0
c
c  iphs3 = 1  -> interior cell
c  iphs3 = 2  -> exterior cell
c
c     Modified Magas 01.11.2005
c   Note: we use ncyc-1 because at this stage all arrays, etc. contain
c    the results for the previous cycle.
      if (ioptos4.ne.0.and.modcyc.eq.0) then
        if (iphs3.eq.1) write (10,810) ncyc-1,iphs3
        if (iphs3.eq.1) write (10,850) totn,keqost,ibar,jbar,kbar
        if (iphs3.eq.2.and.kdebug.ne.0) write (12,850) totn,keqost,ibar
     1  ,jbar,kbar
        if (iphs3.eq.2.and.kdebug.ne.0) write (12,810) ncyc,iphs3
      endif
      if (modcyc.eq.0.and.iphs3.eq.1.and.kdebug.ne.0) then
        write (13,850) totn,keqost,ibar,jbar,kbar
        write (14,850) totn,keqost,ibar,jbar,kbar
        if (iphs3.eq.1) write (11,850) totn,keqost,ibar,jbar,kbar
c         if (kslab.lt.4) write (156,850) totn,keqost,ibar,jbar,kbar
        write (11,810) ncyc,iphs3
        write (13,810) ncyc,iphs3
        write (14,810) ncyc,iphs3
      endif
      sume=0.
      do 370 k=1,kbar
      do 360 j=1,jbar
      do 350 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 350
      sume=sume+e(ijk) ! sumT00 [MeV*fm-3]
      a=capn(ijk)*ren0 ! N0 [fm-3], ren0=1
      b=e(ijk)*rem0n0 ! rem0n0=1, e is T00 [MeV*fm-3]
      cmom=xmom(ijk)**2+ymom(ijk)**2+zmom(ijk)**2 ! [MeV2*fm-6]
      c=cmom*remnsq ! remnsq=1, 
      if (iphs3.eq.1) numclt=numclt+1
      inter=0
      koast=0
c      at this point capn (the density) assumes the value
c       of 1.0000 for normal nuclear density
c
      call koccupy (i,j,k,iow,io6w,jdim1)
      if (iphs3.eq.1.and.iow.eq.0) go to 350
      if (iphs3.eq.2.and.iow.ne.0) go to 350
      if (iphs3.eq.1) numcli=numcli+1
      if (iphs3.eq.2) numcle=numcle+1
      temp=0.
      if (interp.eq.2) go to 300 ! interp=1                                                            |
c                                                                      |
c    Determine if a surface cell
c
      ipjk=ijk+nq
      imjk=ijk-nq
      ijpk=ijk+nqib
      ijmk=ijk-nqib
      ijkp=ijk+nqijb
      ijkm=ijk-nqijb
      dellx=1.
      delly=1.
      dellz=1.
      iversn=0
      if (iversn.eq.0) go to 250
      if (capn(ipjk).eq.0..and.capn(imjk).eq.0.) go to 220
      if (capn(ipjk).ne.0..and.capn(imjk).ne.0.) go to 220
      if (capn(ipjk).ne.0.) go to 210
      call counms (i,j,k,nm0,nm1,0)
      dellx=xextrx-float(i)+1.+delpt
      if (dellx.gt.1.) dellx=1.
      go to 220
c
  210 if (capn(imjk).ne.0.) go to 220
      call counms (i,j,k,nm0,nm1,0)
      dellx=float(i)-xextrn+delpt
      if (dellx.gt.1.) dellx=1.
c   at this point, know i-1 and i+1 are occupied, or both unoccupied
  220 continue
      if (capn(ijkp).eq.0..and.capn(ijkm).eq.0.) go to 240
      if (capn(ijkp).ne.0..and.capn(ijkm).ne.0.) go to 240
      if (capn(ijkp).ne.0.) go to 230
      call counms (i,j,k,nm0,nm1,-1)
      go to 240
c
  230 if (capn(ijkm).ne.0.) go to 240
      call counms (i,j,k,nm0,nm1,-1)
  240 continue
      if (capn(ijpk).ne.0.) go to 250
      call counms (i,j,k,nm0,nm1,-2)
      delly=yextrx-float(j)+1.+delpt
      if (delly.gt.1.) delly=1.
  250 continue
      delv=dellx*delly*dellz
      rdelv=0.
      if(delv.ne.0.)rdelv=1./delv
      trdelv=1000.*rdelv
      itrd=trdelv+.001
      rdelv=0.001*float(itrd)
      if (rdelv.gt.10.) rdelv=10.
      if (rdelv.gt.delvmx) delvmx=rdelv
c     rdelv is equal to 1
      call intern (i,j,k,ijk,iphs3,adens,ib4,numlvv,uijk,vijk,wijk,pijk
     1 ,rdelv)
c      if (iphs3.eq.2) pijk=eps/3.-4.*bbag/3.
c      if (kdebug.ne.0) write (156,690) i,j,k,pijk,(eps-4.*bbag)/3.,bbag,
c     1 pq
c
      if (pijk.lt.0.) pijk=0.
      if (e(ijk).gt.0..and.pijk.gt.e(ijk)) then
      if (pijk.gt.1.d12) write(*,*) 'WARNING p capped in phase3 cyc=',
     1 ncyc,' i,j,k=',i,j,k,' p=',pijk,' e_lab=',e(ijk)
      pijk=e(ijk)
      endif
      if (vijk.ne.0..and.jdim1.eq.1) write (*,*) 'PROBLEM in intern-vijk
     1',i,j,k,vijk,ymom(ijk)
      if (iphs3.eq.1) go to 290 !iphs3=1
      if (iphs3.eq.2) go to 290
      write (*,*) 'intern',iphs3
      stop 99
      usq=uijk**2+vijk**2+wijk**2
      if(usq.le.1.)gg=1./(1.-usq)
c      if (j.eq.1) write (156,760) iphs3,i,j,k,temp,entro,chem,gg,eps,
c     1 pres
      if(c.ne.0.)dd=16.*(b-bbagde)**2/c ! dimensionless
      ds=dsqrt(dd*dd-12.*dd)
      usq=0.5*(dd-6.-ds)
      cc=usq*16.*(b-bbagde)**2/(3.+usq)**2
c       write (156,*) 'dd=',dd,ds,'  usq=',usq,b,c
c
      epsx=b-dsqrt(c*usq) ! energy density in the rest frame in MeV/fm3
      if(usq.le.1.)gxx=1./(1.-usq)
      amom=xmom(ijk)**2+ymom(ijk)**2+zmom(ijk)**2
      if(gxx.ne.0.)xia=a/gxx
      repx=0.
      if(amom.ne.0.)repx=dsqrt(usq/amom)
      uijk=repx*xmom(ijk)
      vijk=repx*ymom(ijk)
      wijk=repx*zmom(ijk)
      vel2=uijk**2+vijk**2+wijk**2
      if (isnan(vel2).or.vel2.ge.1.d0) then
      write (*,*) 'v^2>=1 corrected (5292) cyc/i/j/k=',ncyc,i,j,k,
     1 ' vel2=',vel2
      if (isnan(vel2).or.vel2.le.0.d0) then
      uijk=0.d0
      vijk=0.d0
      wijk=0.d0
      else
      vfac=dsqrt(0.9999d0/vel2)
      uijk=uijk*vfac
      vijk=vijk*vfac
      wijk=wijk*vfac
      endif
      endif

c      write (156,750) usq,gxx,epsx,presx,uijk,vijk,wijk,xia
c Angel, 2024, not needed with new EoS

c      if (presx .lt. 0.) presx=0.
c      presx=epsx/3.-4.*bbagde/3. !pressure in MeV/fm3

c      if (presx .gt. 0.) epsx=epsx
c      if (presx .eq. 0.) epsx=3.*bbagde*((epsx+bbagde)/(4.*bbagde))**(4.
c     1 /3.)
      call eqsqcd3 (xi,epsx) ! epsx is already in  MeV/fm3, xi????
      if(erat.eq.0.) stop 5374

            presx=pq

c     ADD ON FEBRUARY 19, 2024
      if (presx .ge. epsx) then
      write (*,*) 'pijk>e (line5297):', 'pijk= ',presx,'eijk= ',epsx
c       pause
      endif


      Ecalc=(gxx*epsx+(gxx-1.)*presx)/erat  ! MeV/fm**3
      cmaybe=16.*usq*gxx*gxx*(epsx-bbagde)**2/9.
      if(a.ne.0.)xib=a/dsqrt(a)
c       write (156,*) 'xia=',xia,'  xib=',xib
c       if (dabs(Ecalc-b).gt.0.0000000001) write (156,660) b,Ecalc,
c      1 b-Ecalc,a,c
      if (dabs(Ecalc-b).gt.0.0000000001) write (*,660) b,Ecalc,b-Ecalc,a
     1 ,c
c       if (dabs(cmaybe-c).gt.0.0000001) write (156,670) c,cmaybe,c-cmaybe
      if (dabs(cmaybe-c).gt.0.0000001) write (*,670) c,cmaybe,c-cmaybe
c       write (156,*) ' '
      write (*,*)'iphs3=',iphs3
      if (iphs3.ne.3) go to 260
      xlower=0.00000000000000001
      deps=.00000000001
      call rtmi (egnn,froot,gamnix,xlower,a,deps,40,ier)
      if (ier.ne.0) then
        xupper=9.99999999999
        call rtmi (egnn,froot,gamnix,xlower,xupper,deps,40,ier)
      endif
      if(a.ne.0.)xv=1.-egnn*egnn/(a*a)
      gnn=1./(1.-xv)
      usq=xv
c       write (156,730) egnn,xv,gnn,gnn*gnn,gnn-gg
      if (dabs(gnn*gnn-gg).gt.0.001) write (*,*) egnn,gnn,gnn*gnn,gnn
     1 *gnn-gg
c       if (dabs(gnn*gnn-gg).gt.0.00001) write (156,*) '***********'
      gxx=gnn

      epsx=(3.*erat*e(ijk)+4.*bbagde*(gxx-1.))/(4.*gxx-1.)   !MeV/fm**3


      call eqsqcd3 (anq,epsx)      !presx=epsx/3.-4.*bbagde/3.    !MeV/fm**3

      presx=pq

c     ADD ON FEBRUARY 19, 2024
      if (presx .ge. epsx) then
      write (*,*) 'pijk>e (line5340):', 'pijk= ',presx,'eijk= ',epsx
c       pause
      endif

      Ecalc=(gxx*epsx+(gxx-1.)*presx)/erat    !MeV/fm**3
      cmaybe=usq*(epsx*(4./3.)-4.*bbag/3.)**2
      c1=usq*gxx*gxx*(epsx+presx)**2
c       write (156,*) ' c1',c1,c1/erat/erat
c       write (156,740) gxx,b,Ecalc,eps,epsx/erat,pres,presx
c       write (156,*) 'Ecalc',c,(cmaybe*gxx*gxx)/erat**2,c-(
c      1 cmaybe*gxx*gxx)/erat**2,b,Ecalc,b-Ecalc


  260 continue


c Angel, 2024, not needed with new EoS

c       presx=epsx/3.-4.*bbagzz/3.     !GeV/fm**3
c       if (presx.gt.0.) go to 270
c       if (presx.lt.0.) write (156,*) 'qqqqqqqqqqqq'
c       bbagzz=0.25*bbag
c
c
c       call rtmi (gxx,froot,gammfp,dfloat(1),dfloat(2),deps,40,ier)
c       epsx=(3.*erat*e(ijk)+4.*bbag*(gxx-1.))/(4.*gxx-1.)   ! MeV/fm**3
c
c
c
c
c        presx=epsx/3.-4.*bbagzz/3.   !MeV/fm**3
c
c  if (kdebug.eq.0) go to 270
c       calcen=gxx*epsx+(gxx-1.)*presx
c       amomc=gxx*(gxx-1.)*(epsx+presx)**2
c       write (156,*) 'calc E',calcen,erat*e(ijk),calcen-erat*e(ijk),amomc
c      1 ,amom*erat*erat
c   270 continue
c
c
c
c        if (presx .lt. 0.) presx=0.

c      if (presx .gt. 0.) epsx=epsx
c      if (presx .eq. 0.) epsx=3.*bbagde*((epsx+bbagde)/(4.*bbagde))**(4.
c     1 /3.)

      call eqsqcd3 (anq,epsx) ! epsx is already in MeV/fm3 and anq in fm-3

      presx=pq
c     ADD ON FEBRUARY 19, 2024
      if (presx .ge. epsx) then
      write (*,*) 'pijk>e (line5390):', 'pijk= ',presx,'eijk= ',epsx
c       pause
      endif
c       if (kdebug.gt.0) write (156,760) iphs3,i,j,k,tq,sq,amuq,gxx,epsx
c      1 ,presx
      eps=epsx
      pijk=presx/erat
      bbagzz=bbag
      usqrd=1.-1./gxx
      uijk=0.
      vijk=0.
      wijk=0.
c       if (gxx.gt.fast) write (156,*) ' faster',gxx,fast
      if (gxx.gt.fast) fast=gxx
      if (amom.eq.0.) go to 280
      uijk=dsqrt(usqrd)*xmom(ijk)/amom
      vijk=dsqrt(usqrd)*ymom(ijk)/amom
      wijk=dsqrt(usqrd)*zmom(ijk)/amom
      vel2=uijk**2+vijk**2+wijk**2
      if (isnan(vel2).or.vel2.ge.1.d0) then
      write (*,*) 'v^2>=1 corrected (5426) cyc/i/j/k=',ncyc,i,j,k,
     1 ' vel2=',vel2
      if (isnan(vel2).or.vel2.le.0.d0) then
      uijk=0.d0
      vijk=0.d0
      wijk=0.d0
      else
      vfac=dsqrt(0.9999d0/vel2)
      uijk=uijk*vfac
      vijk=vijk*vfac
      wijk=wijk*vfac
      endif
      endif
  280 if (iphs3.eq.1) go to 290
      if (amom.ne.0.) uij=dsqrt(usqrd)*xmom(ijk)/amom
      wij=dsqrt(usqrd)*zmom(ijk)/amom
      if (uijk.gt.0.) uijkx=amin1(uijk,uij)
      if (uijk.le.0.) uijkx=-amin1(dabs(uijk),abs(uij))
      if (iphs3.eq.1) go to 290
      uij=0.
      if (amom.ne.0.) uij=dsqrt(usqrd)*xmom(ijk)/amom
      wij=dsqrt(usqrd)*zmom(ijk)/amom
      if (uijk.gt.0.) uijkx=amin1(uijk,uij)
      if (uijk.le.0.) uijkx=-amin1(dabs(uijk),abs(uij))
      uijk=uijkx
c
  290 continue
c                                                                     |
c  |___________________________________________________________________|
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      koost(ijk0)=0
c
      go to 310
  300 continue
      iplc=ijkret(ijk)
      red1=(volnrm*retvol(iplc))
      if (red1.gt.1.) red1=1.
      red=1./red1
      e(ijk)=e(ijk)*red
      zmom(ijk)=zmom(ijk)*red
      ymom(ijk)=ymom(ijk)*red
      xmom(ijk)=xmom(ijk)*red
      capn(ijk)=capn(ijk)*red
      write (*,*) 'interx'
c
c  ____________________________________________________________________
c  |                                                                   |
      call interx (i,j,k,ijk,1,adens,ib4,numlvv,uijk,vijk,wijk,pijk)
c  |___________________________________________________________________|
c
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c       if (i.eq.50.and.k.eq.107) write (156,*) 'uijk',i,j,k,ijk,uijk,
c      1 uijkx
c       if (i.eq.50.and.k.eq.132) write (156,*) 'uijk',i,j,k,ijk,uijk,
c      1 uijkx
c       if (iphs3.eq.2) write (156,*) 'interx',i,j,k,ijk0
      e(ijk)=e(ijk)*red1
      zmom(ijk)=zmom(ijk)*red1
      ymom(ijk)=ymom(ijk)*red1
      xmom(ijk)=xmom(ijk)*red1
      capn(ijk)=capn(ijk)*red1
      velsq=c/b/b
      kextrn(ijk0)=0
c       if (dabs(wijk).lt.0.0000001) write (156,*) ' wijk '
  310 continue
c
      if (temp.le.0.) temp=0.001
c
      tempp=temp
      if (j.eq.1) eedns(i,k)=tq
c                Reaction plane  temp: here for _PLOT_ before 4019
      gamma=dsqrt(dabs(1.-uijk*uijk-vijk*vijk-wijk*wijk))
      gammaf=1./gamma
      gamsq=gammaf*gammaf
      tijk=temp
      if (keqost.le.0) tijk=tempzz
c      xni=xi*0.1453
c      zprint=0.1453*939.5731*eps*0.001
      xni=xi ! baryon density ?????
      zprint=eps ! MeV/fm**3 ?????
      if (j.eq.1) aplot(i,k)=tq
      if (j.eq.1) bplot(i,k)=sq
      if (j.eq.1) cplot(i,k)=pq
      if (pijk.lt.0.) pijk=0.
      if (kext.ne.0) go to 330
      ko=ko+1
      kstrmx=895195
      if (ko.gt.kstrmx) go to 320
      storev(ko)=uijk
      ko=ko+1
      storev(ko)=vijk
      ko=ko+1
      storev(ko)=wijk
      ko=ko+1
      storev(ko)=pijk
      go to 330
  320 continue
      if (kstrz.eq.0) rewind 889
      kstrz=kstrz+1
      write (889,*) ko,ijk,uijk,vijk,wijk,pijk,kstrz
  330 continue
c      preijk=0.1453*939.5731*pijk
      preijk=pijk ! MeV/fm3
c      if (iphs3.eq.2) xni=0.1453*a*gamma
      if (iphs3.eq.2) xni=a*gamma
c      zphys=0.1453*0.9395731*e(ijk)
      zphys=e(ijk) ! in MeV/fm3
c
      utrans=dsqrt(uijk*uijk+vijk*vijk)
      if (tq.eq.0.) tq=0.01
      if (iphs3.eq.1) captot=captot+capn(ijk)
      if (iphs3.eq.1) eptot=eptot+e(ijk)
      if (iphs3.eq.2) captoe=captoe+capn(ijk)
      if (iphs3.eq.2) eptote=eptote+e(ijk)
      if (xni.gt.1.e-8.and.modcyc.eq.0.and.kext.eq.0) go to 340
c      if (ncyc.eq.1) go to 340
      go to 350
  340 continue
c      if (ioptos4.ne.0) write (10,820) i,j,k,zphys,xni,pq,preijk,sq,uijk
c     1 ,vijk,wijk,tq

c      ------------------------------------------------------------------
c   ADDED ON FEBRUARY 14, 2024

c  Ensure velocities are non-zero if momenta are not
      if (wijk.ne.0..or.uijk.ne.0..or.vijk.ne.0.) go to 111
      if(e(ijk).gt.0.) then

      if(xmom(ijk).ne.0..or.ymom(ijk).ne.0..or.zmom(ijk).ne.0.) then
      anorm=sqrt(xmom(ijk)**2+ymom(ijk)**2+zmom(ijk)**2)
      uijk=xmom(ijk)/anorm
c      u(ijk)=uijk
      vijk=ymom(ijk)/anorm
c      v(ijk)=vijk
      wijk=zmom(ijk)/anorm
c      w(ijk)=wijk
      endif

      elseif(e(ijk).lt.0.) then
c       write(*,*)'probel t00 < 0',i,j,k,e(ijk)
      endif

  111 continue

      if (uijk**2+vijk**2+wijk**2 .ge. 1.) then
c       write (777,*) 'v^2>1 in line 5554',i,j,k,uijk,vijk,wijk,
c      1'v=',uijk**2+vijk**2+wijk**2
c       write (777,*) 'v^2>1 in line 5554',i,j,k,e(ijk),xmom(ijk),
c      1 ymom(ijk),zmom(ijk),capn(ijk)

      anorm=sqrt(xmom(ijk)**2+ymom(ijk)**2+zmom(ijk)**2)
      uijk=0.9999*xmom(ijk)/anorm
c      u(ijk)=uijk
      vijk=0.9999*ymom(ijk)/anorm
c      v(ijk)=vijk
      wijk=0.9999*zmom(ijk)/anorm
c      w(ijk)=wijk

c       write (777,*) 'New v^2 in line 5554',i,j,k,uijk,vijk,wijk,
c      1'v=',uijk**2+vijk**2+wijk**2
c       write (*,*) 'New v^2 in line 5554',i,j,k,uijk,vijk,wijk,
c      1'v=',uijk**2+vijk**2+wijk**2
       endif

c     added Angel 2024
c      u(ijk)=uijk
c      v(ijk)=vijk
c      w(ijk)=wijk
c     added Angel 2024


      if (ioptos4.ne.0) write (10,820)i,j,k,e(ijk),capn(ijk),
     1 pijk,wijk,uijk,vijk,zmom(ijk),ymom(ijk),xmom(ijk),tq,sq,
     2 amuq
      if (ioptos4.ne.0) then
       write (118,820)i,j,k,e(ijk),capn(ijk),
     1 pijk,wijk,uijk,vijk,zmom(ijk),ymom(ijk),xmom(ijk),tq,sq,
     2 amuq
       endif
c       if (sqrt(1.-wijk**2-uijk**2-vijk**2)*capn(ijk).lt.0.2) then
c       capn(ijk)=0.2/sqrt(1.-wijk**2-uijk**2-vijk**2)
c       endif
c       if ((e(ijk)+pijk)*(1.-wijk**2-uijk**2-vijk**2)-pijk .lt. 13.0)then
c       e(ijk)=(13.0+pijk)/(1.-wijk**2-uijk**2-vijk**2)-pijk
c       endif
      if (capn(ijk).gt.huge(1.0d0)) stop 58544
      if (isnan(capn(ijk))) stop 58555
      if (isnan(wijk)) stop 58566
      if (isnan(uijk)) stop 58577
      if (isnan(vijk)) stop 58588


      bardensrfpos(i,j,k)=0.
      bardensrfpos(i,j,k)=sqrt(1.-wijk**2-uijk**2-vijk**2)*
     1 capn(ijk)
     
      bardenslfpos(i,j,k)=0.
      bardenslfpos(i,j,k)=capn(ijk)

      if (bardensrfpos(i,j,k).gt.huge(1.0d0)) stop 5859

      if (isnan(bardensrfpos(i,j,k))) then
      write(*,*) '5858=',i,j,k,bardensrfpos(i,j,k),capn(ijk),
     1 wijk,uijk,vijk,wijk**2+uijk**2+vijk**2
      stop 585858
      endif

      if (e(ijk).gt.huge(1.0d0)) stop 5854
      if (isnan(e(ijk))) stop 5855



c      if (ncyc = 1) pause



      enerdenslfpos(i,j,k)=0.
      enerdenslfpos(i,j,k)=e(ijk)
      enerdensrfpos(i,j,k)=0.
      enerdensrfpos(i,j,k)=(1.-wijk**2-uijk**2-vijk**2)*(e(ijk)+
     1 pijk)-pijk


      if (e(ijk) .ne. 0. .or. capn(ijk) .ne. 0.)  then
      write (334,*) i,j,k,enerdenslfpos(i,j,k),bardenslfpos(i,j,k)
      endif

      if (enerdenslfpre(i,j,k) .ne. 0. .or. bardenslfpre(i,j,k) .ne. 0.)
     & then
      write (333,*) i,j,k,enerdenslfpre(i,j,k),bardenslfpre(i,j,k)
      endif


     
      if (isnan(enerdensrfpos(i,j,k))) stop 5866
      if (enerdensrfpos(i,j,k).gt.huge(1.0d0)) stop 5867

      if (enerdensrfpos(i,j,k).le.0) then
c       write(*,*)'e0=0 en POS',i,j,k, enerdensrfpos(i,j,k)
      enerdensrfpos(i,j,k)=0.
      endif
      xmomdenslfpos(i,j,k)=0.
      xmomdenslfpos(i,j,k)=xmom(ijk)
      ymomdenslfpos(i,j,k)=0.
      ymomdenslfpos(i,j,k)=ymom(ijk)
      zmomdenslfpos(i,j,k)=0.
      zmomdenslfpos(i,j,k)=zmom(ijk)

      if (isnan(xmomdenslfpos(i,j,k))) stop 5881
      if (xmomdenslfpos(i,j,k).gt.huge(1.0d0)) stop 5882
      if (isnan(ymomdenslfpos(i,j,k))) stop 5883
      if (ymomdenslfpos(i,j,k).gt.huge(1.0d0)) stop 5884
      if (isnan(zmomdenslfpos(i,j,k))) stop 5885
      if (zmomdenslfpos(i,j,k).gt.huge(1.0d0)) stop 5886

      vxpos(i,j,k)=0.
      vxpos(i,j,k)=uijk
      vypos(i,j,k)=0.
      vypos(i,j,k)=vijk
      vzpos(i,j,k)=0.
      vzpos(i,j,k)=wijk

      if (isnan(vxpos(i,j,k))) stop 5887
      if (vxpos(i,j,k).gt.huge(1.0d0)) stop 5888
      if (isnan(vypos(i,j,k))) stop 5890
      if (vypos(i,j,k).gt.huge(1.0d0)) stop 5891
      if (isnan(vzpos(i,j,k))) stop 5892
      if (vzpos(i,j,k).gt.huge(1.0d0)) stop 5893

      presspos(i,j,k) = 0.
      presspos(i,j,k) = pijk
      temppos(i,j,k)  = 0.
      temppos(i,j,k)  = temptab(enerdensrfpos(i,j,k),bardensrfpos(i
     & ,j,k))
      xmubpos(i,j,k)   = 0.
      xmubpos(i,j,k)  = chemtab(enerdensrfpos(i,j,k),bardensrfpos(i
     & ,j,k))
c       if (i==50 .and. j==50 .and. k==160) then
c       write(*,*) ncyc,i,j,k,vzpos(i,j,k)
c       pause
c       endif


      
     
      sumtemp=sumtemp+tq
      sumchem=sumchem+amuq
      sumentro=sumentro+sq
      nlpc=nlpc+1
     
c      if (ioptos4.ne.0) write (387,*) sq,tq
c ,zmom,xmom,ymom
      gamma=1./dsqrt(1.-uijk**2-vijk**2-wijk**2)
      if (j.eq.1.and.ioptos3.ne.0) write (26,710) ncyc,i,j,k,capn(ijk)
     1 ,capn(ijk)*dxdydz/pnu(1),xni,gamma*xni,capn(ijk)/xni/gamma,capn
     2 (ijk)/xni/gamma
c       if (kdebug.ne.0.and.abs(pq-preijk).gt.0.01*pq) write (156,720) i,
c      1 j,k,ijk,iphs3,zphys,xni,pq,preijk,sq,uijk,vijk,wijk,tq
c       if (preijk.ne.0..and.kdebug.ne.0) write (156,680) i,j,k,pq,preijk
      vel2=uijk**2+vijk**2+wijk**2
      gamma=1./dsqrt(1.-vel2)
      numcll=numcll+1
      dentot=dentot+xni
      entot=entot+eps
      tempav=tempav+tq
      if (kdebug.ne.0) then
        write (14,830) i,j,k,zprint,xni,pq,tq
        if (j.eq.1) write (13,820) i,j,k,zprint,xni,pq,preijk,sq,uijk
     1  ,vijk,wijk,tq
c
        if (iphs3.eq.1) write (11,820) i,j,k,zprint,xni,pq,preijk,sq
     1  ,uijk,vijk,wijk,tq
        if (iphs3.eq.2) write (12,820) i,j,k,zprint,xni,pq,preijk,sq
     1  ,uijk,vijk,wijk,tq
c
        if (j.eq.1) write (15,840) i,k,zprint,xni,preijk,utrans,wijk,tq
        if (j.eq.3) write (16,840) i,k,zprint,xni,preijk,utrans,wijk,tq
      endif
c
      totn=totn+2*xni*dxdydz*gammaf

  350 continue
  360 continue
  370 continue

c     shift from central values of cells to vertice values
c     for Cornelius



      do 2370 k=1,kbar
      do 2360 j=1,jbar
      do 2350 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      posNx(i,j,k)=0.
      posNy(i,j,k)=0.
      posNz(i,j,k)=0.
      posN0(i,j,k)=0.
      if (vxpos(i,j,k)**2+vypos(i,j,k)**2
     1 +vzpos(i,j,k)**2.lt.1.) then
      gammapos=1./sqrt(1.-vxpos(i,j,k)**2-vypos(i,j,k)**2
     1 -vzpos(i,j,k)**2)
c       posNx(i,j,k)=bardensrfpos(i,j,k)*gammapos*vxpos(i,j,k)
c       posNy(i,j,k)=bardensrfpos(i,j,k)*gammapos*vypos(i,j,k)
c       posNz(i,j,k)=bardensrfpos(i,j,k)*gammapos*vzpos(i,j,k)
c       posN0(i,j,k)=bardensrfpos(i,j,k)*gammapos
      posNx(i,j,k)=capn(ijk)*vxpos(i,j,k)
      posNy(i,j,k)=capn(ijk)*vypos(i,j,k)
      posNz(i,j,k)=capn(ijk)*vzpos(i,j,k)
      posN0(i,j,k)=capn(ijk)
      else
      posNx(i,j,k)=0.
      posNy(i,j,k)=0.
      posNz(i,j,k)=0.
      posN0(i,j,k)=0.
      vxpos(i,j,k)=0.
      vypos(i,j,k)=0.
      vzpos(i,j,k)=0.
      bardensrfpos(i,j,k)=0.
      enerdensrfpos(i,j,k)=0.
      presspos(i,j,k)=0.
      temppos(i,j,k)  = 0.
      xmubpos(i,j,k)   = 0.
      endif

      preNx(i,j,k)=0.
      preNy(i,j,k)=0.
      preNz(i,j,k)=0.
      preN0(i,j,k)=0.

      if (vxpre(i,j,k)**2+vypre(i,j,k)**2
     1 +vzpre(i,j,k)**2.lt.1.) then
      gammapre=1./sqrt(1.-vxpre(i,j,k)**2-vypre(i,j,k)**2
     1 -vzpre(i,j,k)**2)

      preNx(i,j,k)=bardensrfpre(i,j,k)*gammapre*vxpre(i,j,k)
      preNy(i,j,k)=bardensrfpre(i,j,k)*gammapre*vypre(i,j,k)
      preNz(i,j,k)=bardensrfpre(i,j,k)*gammapre*vzpre(i,j,k)
      preN0(i,j,k)=bardensrfpre(i,j,k)*gammapre

      else
      preNx(i,j,k)=0.
      preNy(i,j,k)=0.
      preNz(i,j,k)=0.
      preN0(i,j,k)=0.
      vxpre(i,j,k)=0.
      vypre(i,j,k)=0.
      vzpre(i,j,k)=0.
      bardensrfpre(i,j,k)=0.
      enerdensrfpre(i,j,k)=0.
      presspre(i,j,k)=0.
      temppre(i,j,k)=0.
      xmubpre(i,j,k)=0.
      endif

 2350 continue
 2360 continue
 2370 continue

      do 1370 k=2,kbar-1
      do 1360 j=2,jbar-1
      do 1350 i=2,ibar-1

      posNNx=(posNx(i,j,k)+posNx(i-1,j,k)
     1 +posNx(i,j-1,k)+posNx(i,j,k-1)
     2 +posNx(i-1,j-1,k)+posNx(i-1,j,k-1)
     3 +posNx(i,j-1,k-1)+posNx(i-1,j-1,k-1))/8.

      if (isnan(posNNx)) stop 5963
      if (posNNx.gt.huge(1.0d0)) stop 6011

      posNNy=(posNy(i,j,k)+posNy(i-1,j,k)
     1 +posNy(i,j-1,k)+posNy(i,j,k-1)
     2 +posNy(i-1,j-1,k)+posNy(i-1,j,k-1)
     3 +posNy(i,j-1,k-1)+posNy(i-1,j-1,k-1))/8.

      if (isnan(posNNy)) stop 5968
      if (posNNx.gt.huge(1.0d0)) stop 6019

      posNNz=(posNz(i,j,k)+posNz(i-1,j,k)
     1 +posNz(i,j-1,k)+posNz(i,j,k-1)
     2 +posNz(i-1,j-1,k)+posNz(i-1,j,k-1)
     3 +posNz(i,j-1,k-1)+posNz(i-1,j-1,k-1))/8.

      if (isnan(posNNz)) stop 5973
      if (posNNx.gt.huge(1.0d0)) stop 6027

      posNN0=(posN0(i,j,k)+posN0(i-1,j,k)
     1 +posN0(i,j-1,k)+posN0(i,j,k-1)
     2 +posN0(i-1,j-1,k)+posN0(i-1,j,k-1)
     3 +posN0(i,j-1,k-1)+posN0(i-1,j-1,k-1))/8.

      if (posNN0.lt.0.) posNN0=0.
      if (posNN0.gt.huge(1.0d0)) stop 6035
      if (isnan(posNN0)) stop 5979

      posT0x=(xmomdenslfpos(i,j,k)+xmomdenslfpos(i-1,j,k)
     1 +xmomdenslfpos(i,j-1,k)+xmomdenslfpos(i,j,k-1)
     2 +xmomdenslfpos(i-1,j-1,k)+xmomdenslfpos(i-1,j,k-1)
     3 +xmomdenslfpos(i,j-1,k-1)+xmomdenslfpos(i-1,j-1,k-1))/8.

      if (isnan(posT0x)) stop 6033

      posT0y=(ymomdenslfpos(i,j,k)+ymomdenslfpos(i-1,j,k)
     1 +ymomdenslfpos(i,j-1,k)+ymomdenslfpos(i,j,k-1)
     2 +ymomdenslfpos(i-1,j-1,k)+ymomdenslfpos(i-1,j,k-1)
     3 +ymomdenslfpos(i,j-1,k-1)+ymomdenslfpos(i-1,j-1,k-1))/8.

      if (isnan(posT0y)) stop 6040

      posT0z=(zmomdenslfpos(i,j,k)+zmomdenslfpos(i-1,j,k)
     1 +zmomdenslfpos(i,j-1,k)+zmomdenslfpos(i,j,k-1)
     2 +zmomdenslfpos(i-1,j-1,k)+zmomdenslfpos(i-1,j,k-1)
     3 +zmomdenslfpos(i,j-1,k-1)+zmomdenslfpos(i-1,j-1,k-1))/8.

      if (isnan(posT0z)) stop 6047

      posT00=(enerdenslfpos(i,j,k)+enerdenslfpos(i-1,j,k)
     1 +enerdenslfpos(i,j-1,k)+enerdenslfpos(i,j,k-1)
     2 +enerdenslfpos(i-1,j-1,k)+enerdenslfpos(i-1,j,k-1)
     3 +enerdenslfpos(i,j-1,k-1)+enerdenslfpos(i-1,j-1,k-1))/8.

      if (posT00.lt.0.) posT00=0.
      if (isnan(posT00)) stop 6054

      vxpost(i,j,k)=0.
      vypost(i,j,k)=0.
      vzpost(i,j,k)=0.

      if (isnan(vxpost(i,j,k))) stop 5990000
      if (isnan(vypost(i,j,k))) stop 5991000
      if (isnan(vzpost(i,j,k))) stop 5992000

      if(posNN0.eq.0.and.posT00.gt.0.)then ! Landau's definition
      vxpost(i,j,k)=posT0x/posT00
      vypost(i,j,k)=posT0y/posT00
      vzpost(i,j,k)=posT0z/posT00
      endif

      if (isnan(vxpost(i,j,k))) stop 6071
      if (isnan(vypost(i,j,k))) stop 6072
      if (isnan(vzpost(i,j,k))) stop 6073

      if(posNN0.ne.0.)then ! Eckart's definition
      vxpost(i,j,k)=posNNx/posNN0
      vypost(i,j,k)=posNNy/posNN0
      vzpost(i,j,k)=posNNz/posNN0
      endif

c       if(posNN0.gt.0.0)then
c       if (posNNx.ne.0.) vxpost(i,j,k)=posNNx/posNN0
c       if (posNNy.ne.0.) vypost(i,j,k)=posNNy/posNN0
c       if (posNNz.ne.0.) vzpost(i,j,k)=posNNz/posNN0
c       endif
c
c       if(posT00.gt.0.)then
c       if (posT0x.ne.0.) vxpost(i,j,k)=posT0x/posT00
c       if (posT0y.ne.0.) vypost(i,j,k)=posT0y/posT00
c       if (posT0z.ne.0.) vzpost(i,j,k)=posT0z/posT00
c       endif

      if (isnan(vxpost(i,j,k))) then
      write(*,*)'5990 stop',i,j,k,vxpost(i,j,k),posNNx,posNN0
      stop 5990
      endif

      if (isnan(vypost(i,j,k))) stop 5991
      if (isnan(vzpost(i,j,k))) stop 5992

      bardensrfpost(i,j,k)=(bardensrfpos(i,j,k)+bardensrfpos(i-1,j,k)
     1 +bardensrfpos(i,j-1,k)+bardensrfpos(i,j,k-1)
     2 +bardensrfpos(i-1,j-1,k)+ bardensrfpos(i-1,j,k-1)
     3 +bardensrfpos(i,j-1,k-1)+bardensrfpos(i-1,j-1,k-1))/8.

      if (bardensrfpost(i,j,k).lt.0) then
c       write(*,*)'bardensrfpost<=0 ',i,j,k,
c      1 bardensrfpost(i,j,k)
      bardensrfpost(i,j,k)=0.
      endif

      enerdensrfpost(i,j,k)=(enerdensrfpos(i,j,k)+enerdensrfpos(i-1,j,k)
     1 +enerdensrfpos(i,j-1,k)+enerdensrfpos(i,j,k-1)
     2 +enerdensrfpos(i-1,j-1,k)+ enerdensrfpos(i-1,j,k-1)
     3 +enerdensrfpos(i,j-1,k-1)+enerdensrfpos(i-1,j-1,k-1))/8.

      if (enerdensrfpost(i,j,k).lt.0) then
c       write(*,*)'enerdensrfpost<=0 in line 5898',i,j,k,
c      1 enerdensrfpost(i,j,k)
      enerdensrfpost(i,j,k)=0.
      endif

      if (vxpost(i,j,k)**2+vypost(i,j,k)**2
     1 +vzpost(i,j,k)**2.ge.1.) then
c       write(*,*) vxpost(i,j,k)**2+vypost(i,j,k)**2
c      1 +vzpost(i,j,k)**2
c       write(*,*)'RF',bardensrfpost(i,j,k),enerdensrfpost(i,j,k)
c       write(*,*)'N^mu',posNN0,posNNx,posNNy,posNNz
c       write(*,*)'T^mu nu',posT00,posT0x,posT0y,posT0z
       vxpost(i,j,k)=0.
       vypost(i,j,k)=0.
       vzpost(i,j,k)=0.
       bardensrfpost(i,j,k)=0.
       enerdensrfpost(i,j,k)=0.
c       stop 5798
      endif

      if (isnan(vxpost(i,j,k))) stop 6029
      if (isnan(vypost(i,j,k))) stop 6030
      if (isnan(vzpost(i,j,k))) stop 6031

c       enerdenslfpost(i,j,k)=(enerdenslfpos(i,j,k)+enerdenslfpos(i-1,j,k)
c      1 +enerdenslfpos(i,j-1,k)+enerdenslfpos(i,j,k-1)
c      2 +enerdenslfpos(i-1,j-1,k)+ enerdenslfpos(i-1,j,k-1)
c      3 +enerdenslfpos(i,j-1,k-1)+enerdenslfpos(i-1,j-1,k-1))/8.
c
c       xmomdenslfpost(i,j,k)=(xmomdenslfpos(i,j,k)+xmomdenslfpos(i-1,j,k)
c      1 +xmomdenslfpos(i,j-1,k)+xmomdenslfpos(i,j,k-1)
c      2 +xmomdenslfpos(i-1,j-1,k)+ xmomdenslfpos(i-1,j,k-1)
c      3 +xmomdenslfpos(i,j-1,k-1)+xmomdenslfpos(i-1,j-1,k-1))/8.
c
c       ymomdenslfpost(i,j,k)=(ymomdenslfpos(i,j,k)+ymomdenslfpos(i-1,j,k)
c      1 +ymomdenslfpos(i,j-1,k)+ymomdenslfpos(i,j,k-1)
c      2 +ymomdenslfpos(i-1,j-1,k)+ ymomdenslfpos(i-1,j,k-1)
c      3 +ymomdenslfpos(i,j-1,k-1)+ymomdenslfpos(i-1,j-1,k-1))/8.
c
c       zmomdenslfpost(i,j,k)=(xmomdenslfpos(i,j,k)+xmomdenslfpos(i-1,j,k)
c      1 +xmomdenslfpos(i,j-1,k)+xmomdenslfpos(i,j,k-1)
c      2 +xmomdenslfpos(i-1,j-1,k)+ xmomdenslfpos(i-1,j,k-1)
c      3 +xmomdenslfpos(i,j-1,k-1)+xmomdenslfpos(i-1,j-1,k-1))/8.

      z=enerdensrfpost(i,j,k)
      az=bardensrfpost(i,j,k)

      call eqsqcd3(az,z)

      presspos(i,j,k) = pq
      temppos(i,j,k)  = tq
      xmubpos(i,j,k)   = amuq

c      ---------------------------- PRE

      preNNx=(preNx(i,j,k)+preNx(i-1,j,k)
     1 +preNx(i,j-1,k)+preNx(i,j,k-1)
     2 +preNx(i-1,j-1,k)+preNx(i-1,j,k-1)
     3 +preNx(i,j-1,k-1)+preNx(i-1,j-1,k-1))/8.

      if (isnan(preNNx)) stop 6058

      preNNy=(preNy(i,j,k)+preNy(i-1,j,k)
     1 +preNy(i,j-1,k)+preNy(i,j,k-1)
     2 +preNy(i-1,j-1,k)+preNy(i-1,j,k-1)
     3 +preNy(i,j-1,k-1)+preNy(i-1,j-1,k-1))/8.

      if (isnan(preNNy)) stop 6063

      preNNz=(preNz(i,j,k)+preNz(i-1,j,k)
     1 +preNz(i,j-1,k)+preNz(i,j,k-1)
     2 +preNz(i-1,j-1,k)+preNz(i-1,j,k-1)
     3 +preNz(i,j-1,k-1)+preNz(i-1,j-1,k-1))/8.

      if (isnan(preNNz)) stop 6068

      preNN0=(preN0(i,j,k)+preN0(i-1,j,k)
     1 +preN0(i,j-1,k)+preN0(i,j,k-1)
     2 +preN0(i-1,j-1,k)+preN0(i-1,j,k-1)
     3 +preN0(i,j-1,k-1)+preN0(i-1,j-1,k-1))/8.

      if (preNN0.lt.0.) preNN0=0.
      if (isnan(preNN0)) stop 6080

      preT0x=(xmomdenslfpre(i,j,k)+xmomdenslfpre(i-1,j,k)
     1 +xmomdenslfpre(i,j-1,k)+xmomdenslfpre(i,j,k-1)
     2 +xmomdenslfpre(i-1,j-1,k)+xmomdenslfpre(i-1,j,k-1)
     3 +xmomdenslfpre(i,j-1,k-1)+xmomdenslfpre(i-1,j-1,k-1))/8.

      if (isnan(preT0x)) stop 6181

      preT0y=(ymomdenslfpre(i,j,k)+ymomdenslfpre(i-1,j,k)
     1 +ymomdenslfpre(i,j-1,k)+ymomdenslfpre(i,j,k-1)
     2 +ymomdenslfpre(i-1,j-1,k)+ymomdenslfpre(i-1,j,k-1)
     3 +ymomdenslfpre(i,j-1,k-1)+ymomdenslfpre(i-1,j-1,k-1))/8.

      if (isnan(preT0y)) stop 6188

      preT0z=(zmomdenslfpre(i,j,k)+zmomdenslfpre(i-1,j,k)
     1 +zmomdenslfpre(i,j-1,k)+zmomdenslfpre(i,j,k-1)
     2 +zmomdenslfpre(i-1,j-1,k)+zmomdenslfpre(i-1,j,k-1)
     3 +zmomdenslfpre(i,j-1,k-1)+zmomdenslfpre(i-1,j-1,k-1))/8.

      if (isnan(preT0z)) stop 6195

      preT00=(enerdenslfpre(i,j,k)+enerdenslfpre(i-1,j,k)
     1 +enerdenslfpre(i,j-1,k)+enerdenslfpre(i,j,k-1)
     2 +enerdenslfpre(i-1,j-1,k)+enerdenslfpre(i-1,j,k-1)
     3 +enerdenslfpre(i,j-1,k-1)+enerdenslfpre(i-1,j-1,k-1))/8.

      if (preT00.lt.0.) preT00=0.
      if (isnan(preT00)) stop 6203

      vxprev(i,j,k)=0.
      vyprev(i,j,k)=0.
      vzprev(i,j,k)=0.

      if(preNN0.eq.0.and.preT00.gt.0.)then ! Landau's definition
      vxprev(i,j,k)=preT0x/preT00
      vyprev(i,j,k)=preT0y/preT00
      vzprev(i,j,k)=preT0z/preT00
      endif

      if(preNN0.ne.0.)then ! Eckart's definition
      vxprev(i,j,k)=preNNx/preNN0
      vyprev(i,j,k)=preNNy/preNN0
      vzprev(i,j,k)=preNNz/preNN0
      endif

c       if(preNN0.gt.0.0)then
c       if (preNNx.ne.0.) vxprev(i,j,k)=preNNx/preNN0
c       if (preNNy.ne.0.) vyprev(i,j,k)=preNNy/preNN0
c       if (preNNz.ne.0.) vzprev(i,j,k)=preNNz/preNN0
c       endif

c       if(preT00.gt.0.)then
c       if (preT0x.ne.0.) vxprev(i,j,k)=preT0x/preT00
c       if (preT0y.ne.0.) vyprev(i,j,k)=preT0y/preT00
c       if (preT0z.ne.0.) vzprev(i,j,k)=preT0z/preT00
c       endif

      if (isnan(vxprev(i,j,k))) stop 6100
      if (isnan(vyprev(i,j,k))) stop 6101
      if (isnan(vzprev(i,j,k))) stop 6102

      bardensrfprev(i,j,k)=(bardensrfpre(i,j,k)+bardensrfpre(i-1,j,k)
     1 +bardensrfpre(i,j-1,k)+bardensrfpre(i,j,k-1)
     2 +bardensrfpre(i-1,j-1,k)+ bardensrfpre(i-1,j,k-1)
     3 +bardensrfpre(i,j-1,k-1)+bardensrfpre(i-1,j-1,k-1))/8.

      if (bardensrfprev(i,j,k).lt.0) then
c       write(*,*)'bardensrfprev<=0 in line 6015',i,j,k,
c      1 bardensrfprev(i,j,k)
      bardensrfprev(i,j,k)=0.
      endif

      enerdensrfprev(i,j,k)=(enerdensrfpre(i,j,k)+enerdensrfpre(i-1,j,k)
     1 +enerdensrfpre(i,j-1,k)+enerdensrfpre(i,j,k-1)
     2 +enerdensrfpre(i-1,j-1,k)+ enerdensrfpre(i-1,j,k-1)
     3 +enerdensrfpre(i,j-1,k-1)+enerdensrfpre(i-1,j-1,k-1))/8.

      if (enerdensrfprev(i,j,k).lt.0) then
c       write(*,*)'enerdensrfprev<=0 in line 6015',i,j,k,
c      1 enerdensrfprev(i,j,k)
      enerdensrfprev(i,j,k)=0.
      endif

      if (vxprev(i,j,k)**2+vyprev(i,j,k)**2
     1 +vzprev(i,j,k)**2.ge.1.) then
c       write(*,*) vxprev(i,j,k)**2+vyprev(i,j,k)**2
c      1 +vzprev(i,j,k)**2
       vxprev(i,j,k)=0.
       vyprev(i,j,k)=0.
       vzprev(i,j,k)=0.
       bardensrfprev(i,j,k)=0.
       enerdensrfprev(i,j,k)=0.
      endif

      if (isnan(vxprev(i,j,k))) stop 6137
      if (isnan(vyprev(i,j,k))) stop 6138
      if (isnan(vzprev(i,j,k))) stop 6139

c       enerdenslfprev(i,j,k)=(enerdenslfpre(i,j,k)+enerdenslfpre(i-1,j,k)
c      1 +enerdenslfpre(i,j-1,k)+enerdenslfpre(i,j,k-1)
c      2 +enerdenslfpre(i-1,j-1,k)+ enerdenslfpre(i-1,j,k-1)
c      3 +enerdenslfpre(i,j-1,k-1)+enerdenslfpre(i-1,j-1,k-1))/8.
c
c       xmomdenslfprev(i,j,k)=(xmomdenslfpre(i,j,k)+xmomdenslfpre(i-1,j,k)
c      1 +xmomdenslfpre(i,j-1,k)+xmomdenslfpre(i,j,k-1)
c      2 +xmomdenslfpre(i-1,j-1,k)+ xmomdenslfpre(i-1,j,k-1)
c      3 +xmomdenslfpre(i,j-1,k-1)+xmomdenslfpre(i-1,j-1,k-1))/8.
c
c       ymomdenslfprev(i,j,k)=(ymomdenslfpre(i,j,k)+ymomdenslfpre(i-1,j,k)
c      1 +ymomdenslfpre(i,j-1,k)+ymomdenslfpre(i,j,k-1)
c      2 +ymomdenslfpre(i-1,j-1,k)+ ymomdenslfpre(i-1,j,k-1)
c      3 +ymomdenslfpre(i,j-1,k-1)+ymomdenslfpre(i-1,j-1,k-1))/8.
c
c       zmomdenslfprev(i,j,k)=(xmomdenslfpre(i,j,k)+xmomdenslfpre(i-1,j,k)
c      1 +xmomdenslfpre(i,j-1,k)+xmomdenslfpre(i,j,k-1)
c      2 +xmomdenslfpre(i-1,j-1,k)+ xmomdenslfpre(i-1,j,k-1)
c      3 +xmomdenslfpre(i,j-1,k-1)+xmomdenslfpre(i-1,j-1,k-1))/8.

      z=enerdensrfprev(i,j,k)
      az=bardensrfprev(i,j,k)

      call eqsqcd3(az,z)
c     añadir common
      presspre(i,j,k) = pq
      temppre(i,j,k)  = tq
      xmubpre(i,j,k)   = amuq

 1350 continue
 1360 continue
 1370 continue

c
c  check external cells that have no internal neighbors and store the
c  stored velocities
      ko=0
c       if (kdebug.ne.0) write (156,*) ' Begin Phase 4, ncyc=',ncyc,
c      1 ' ___________________'
      if (kstrz.ne.0) rewind 889
c       if (kstrz.ne.0) write (156,*) ' kstrz',kstrz,' iphs3',iphs3
      do 410 k=1,kbar
      do 400 j=1,jbar
      do 390 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 390
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      call koccupy (i,j,k,iow,io6w,jdim1)
      if (iphs3.eq.1.and.iow.eq.0) go to 390
      if (iphs3.eq.2.and.iow.ne.0) go to 390
      kext=kextrn(ijk0)
      if (kext.ne.0) go to 390
      ko=ko+1
c  should make the dimension consistent with common
      if (ko.gt.kstrmx) go to 380
      u(ijk)=storev(ko)
      ko=ko+1
      v(ijk)=storev(ko)
      ko=ko+1
      w(ijk)=storev(ko)
      ko=ko+1
      p(ijk)=storev(ko)
      if (p(ijk).lt.0.) p(ijk)=0.
      if (e(ijk).gt.0..and.p(ijk).gt.e(ijk)) p(ijk)=e(ijk)
      if (u(ijk)**2+v(ijk)**2+w(ijk)**2 .gt. 1.) then
      write (*,*) 'v^2>1 in line 5554',i,j,k,u(ijk),v(ijk),w(ijk)
      wijk=0.
      uijk=0.
      vijk=0.
      w(ijk)=0.
      u(ijk)=0.
      v(ijk)=0.
c       pause 5692
      endif
      go to 390
  380 continue
      read (889,*) ko,ijkz,uijk,vijk,wijk,pijk
      if (ijk.ne.ijkz) write (*,*) ' problems with storev',ijk,ijkz
c       if (ijk.ne.ijkz) write (156,*) ' problems with storev',ijk,ijkz
      u(ijk)=uijk
      v(ijk)=vijk
      w(ijk)=wijk
      p(ijk)=pijk
      if (p(ijk).lt.0.) p(ijk)=0.
      if (e(ijk).gt.0..and.p(ijk).gt.e(ijk)) p(ijk)=e(ijk)
  390 continue
  400 continue
  410 continue
C
C     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c      sumN0=0
c      sumT00=0
c      sumT0Z=0
c      sumpress=0
c      sumwijk=0
c      do i=1,ibar
c      do j=1,jbar
c      do k=1,kbar
c      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c      sumN0=sumN0+capn(ijk)
c      sumT00=sumT00+e(ijk)
c      sumT0Z=sumT0Z+zmom(ijk)
c      sumpress=sumpress+p(ijk)
c      sumwijk=sumwijk+abs(w(ijk))
c      enddo
c      enddo
c      enddo
c      write(*,*) 'N0phas3=',sumN0*dxdydz,'T00=',sumT00*dxdydz,
c     1 'sumT0Z=',sumT0Z*dxdydz,sumpress*dxdydz,sumwijk
      if (iphs3.eq.1) koi=ko
c       if (ncyc.gt.690) write (156,*) '++++ ko=',ko
  420 continue
c       if (kdebug.ne.0) write (156,770) captot,captoe,(captot+captoe)
c      1 *dxdydz,eptot*dxdydz,eptote*dxdydz,(eptot+eptote)*dxdydz,koi,ko
c      2 ,numext
c     iphs3=3
c
c  _______________ now deal with the external cells _______________
c
      kphs=4
      if (interp.eq.2) go to 610 !interp=1
      if (numext.eq.0) go to 610 ! numext=0
      kmmm=0
      ko=0
      iexternl=3
      if (modcyc.eq.0) then
        if (kdebug.ne.0) write (12,850) totn,keqost,ibar,jbar,kbar
        if (kdebug.ne.0) write (11,810) ncyc,iphs3
      endif
c
c   ______________________________
      do 590 k=1,kbar
      do 580 j=1,jbar
      do 570 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 570
      kext=kextrn(ijk0)
      if (kext.eq.0) go to 570
      kmmm=kmmm+1
      ijkm=ijk-nqijb
      imjkm=ijkm-nq
      ipjkm=ijkm+nq
      ijmkm=ijkm-nqib
      imjmkm=ijmkm-nq
      ipjmkm=ijmkm+nq
      ijpkm=ijkm+nqib
      imjpkm=ijpkm-nq
      ipjpkm=ijpkm+nq
      a=capn(ijk)*ren0 ! N0 in fm-3, ren0=1.
      b=e(ijk)*rem0n0 ! T00 in MeV/fm3, rem0n0=1.
      cmom=xmom(ijk)**2+ymom(ijk)**2+zmom(ijk)**2
      c=cmom*remnsq
c      if (c.gt.bsq.and.ncyc.le.3) write (156,*) 'before aver',i,j,k,c,bsq
      ra=1./a
      bsq=b*b
      call aver (i,j,k,ijk,ijmkm,ijpkm,imjkm,ipjkm,xi,bsq,x,y,jdim1)
      if (j.eq.1.and.i.eq.imid) write (*,*) 'after aver  k=',k,'  eps='
     1 ,eps,'  xi=',xi
      koost(ijk0)=0
      if (xi.eq.0.) go to 440
      if (compmy.ge.xi) go to 430
      icmax=i
      jcmax=j
      kcmax=k
      compmy=xi
  430 if (y.lt.compe) go to 450
      compe=y
      go to 450
  440 xi=0.
      x=0.
c       write (156,870) i,j,k,ijk,a,b,c,capn(ijk),e(ijk),cmom
      y=0.
      capn(ijk)=0.
      e(ijk)=0.
      xmom(ijk)=0.
      ymom(ijk)=0.
      zmom(ijk)=0.
      u(ijk)=0.
      v(ijk)=0.
      w(ijk)=0.
      p(ijk)=0.
      go to 570
  450 continue
      eps=em0n0*y ! y is obtained from subroutine aver (energy density in the rest frame
      if (j.ne.1) go to 460
      adens(i,k)=xi ! xi obtained from subroutine aver (baryon desnsity in the rest frame)
  460 continue
c  note that hxi is not defined here
      cpntot=cpntot+a
      if (keqost.le.0) go to 470
      hxi=eqost(xi)
      pijk=eosgm1*eps+em0n0*hxi
      tijk=10.
      go to 500
  470 continue


      z=eps  ! eps energy density in the rest frame in MeV/fm3
      az=xi ! baryon density in the rest frame

c Angel, 2024, not needed with new EoS

c       presx=eps/3.-4.*bbagde/3. !pressure in MeV/fm3
c       if (presx .lt. 0.) presx=0.
c       if (presx .gt. 0.) z=eps
c       if (presx .eq. 0.) z=3.*bbagde*((eps+bbagde)/(4.*bbagde))**(4./
c      1 3.)
      if (keqost.eq.-2) go to 480
      if (keqost.eq.-3) go to 480
      if (keqost.eq.-9) go to 490
c                     CALLING SUBROUTINE eqsqcd3


 480  continue
c     ADD ON FEBRUARY 19, 2024
      if (z.le.0.) then
c       write (*,*) 'e<0 (line5810):','eijk= ',z
c       write (*,*)e(ijk),xmom(ijk),ymom(ijk),zmom(ijk),capn(ijk),p(ijk)
      z=0.
      eps=0.
      endif


      call eqsqcd3 (az,z)
c      
      pijk=pq !/(939.5731*0.1453)
      tijk=tq
      pres=pq
      temp=tq
      chem=amuq
      entro=sq

      if (pijk .gt. z) then
      write (*,*) 'pijk>e (line5810):', 'pijk= ',pijk,'eijk= ',z
      endif

c c     ADD ON FEBRUARY 19, 2024
c       if (pijk .ge. z) then
c       write (*,*) 'pijk>e (line5810):', 'pijk= ',pijk,'eijk= ',z
c       z=0.
c       az=0.
c
c       pijk=0.
c       uijk=0.
c       vijk=0.
c       wijk=0.
c
c       e(ijk)=0.
c       xmom(ijk)=0.
c       ymom(ijk)=0.
c       zmom(ijk)=0.
c       capn(ijk)=0.
c
c        tijk=0
c        pres=0.
c        temp=0.
c        chem=0.
c        entro=0.
c c       pause
c       endif
c       if (az.eq.0.) write (156,*) ' az=0',i,j,k,az,z
      go to 500
c_                                      __________________prxx
  490 call prxx (az,z)
      pijk=pres !/(939.5731*0.1453)
      tijk=tempzz
  500 continue
      p(ijk)=pijk
      if (p(ijk).lt.0.) p(ijk)=0.
      if (e(ijk).gt.0..and.p(ijk).gt.e(ijk)) p(ijk)=e(ijk)



c            FREEZE OUT Condition
c      write(*,*)freez, kfreeze, kfrzcy
c      write(*,*)'press=',press,'p=',p(ijk)
      if (xi.lt.freez.and.kfreez.ne.0.and.ncyc.gt.kfrzcy) go to 510



c      epsc=epn(xi)-epn0
      epsc=0.
      if (xi.ne.0.) xth=(eps/xi-eusual)*econ-epsc
c xth is the thermal energy density
      if (xth.lt.0.) xth=0.
      if (c.lt.bsq) x=c/bsq
      if (p(ijk).gt.press) go to 530
      p(ijk)=0.
      go to 530




c
c    remove particles from system  FREEZE OUT
c
  510 CONTINUE
      write(*,*) 'enter into hydro freeze out'
      if (c.lt.bsq) x=c/bsq
      call counms (i,j,k,numout,numou2,1)
      write (*,*) 'stop in freeze out'
      stop 8522
      numouu=numouu+numout
      numouu=numouu+numout+numou2
c      epsc=epn(xi)-epn0
      epsc=0.
      if (xi.ne.0.) xth=(eps/xi-eusual)*econ-epsc
c xth is the thermal energy density
      if (xth.lt.0.) xth=0.
      xarg=0.086673*xth/xi**.66666666
c      epsc=epn(xi)-epn0
      epsc=0.
      xth=(eps/xi-eusual)*econ-epsc
      if (j.ne.1) go to 520
      adens(i,k)=en
c       write (156,*) ' freeze',i,k,en
      numfz1=numfz1+1
      adens(i,k)=-99999999.
      edensz(numfz1)=xi
      idensz(numfz1)=100*i+k
  520 continue
      numcst=numcst+1
      capn(ijk)=0.
      e(ijk)=0.
      xmom(ijk)=0.
      ymom(ijk)=0.
      zmom(ijk)=0.
      u(ijk)=0.
      v(ijk)=0.
      w(ijk)=0.
      p(ijk)=0.
      go to 560





  530 continue
      uijk=0.
      vijk=0.
      wijk=0.
      if (xmom(ijk).eq.0..and.ymom(ijk).eq.0..and.zmom(ijk).eq.0.) go to
     1 560
      if (cmom.lt.1.e-14) go to 560

      if (x.lt.-1.e-5) then
      x=c/b**2
      if (x.lt.-1.e-5)  go to 620
      endif

      x=amax1(x,0.d+00)
      repx=dsqrt(x/cmom)
      uijk=repx*xmom(ijk)
      vijk=repx*ymom(ijk)
      wijk=repx*zmom(ijk)
      gammaf=1./dsqrt(dabs(1.-uijk*uijk-vijk*vijk-wijk*wijk))
      gamma=dsqrt(1.-x)  !inverso de gamma!
      xi=a*gamma
      if (uijk**2+vijk**2+wijk**2 .gt. 1.) then
      write (*,*) 'v^2>1'
      stop 5943
      endif
c      z=0.1453*939.5731*eps*0.001
c       z=0.1453*939.5731*eps

c     Volodya, but no effect
      z=(e(ijk)+p(ijk))*gamma**2 - p(ijk) !inverso de gamma!
      if(z.ne.eps) then
c       write (*,*) 'e0 correct=',z,'eps=',eps
c       write (*,*) i,j,k,e(ijk),xi,p(ijk),uijk,vijk,wijk,xmom(ijk),
c     1 ymom(ijk),zmom(ijk)
      endif
c      preijk=0.1453*939.5731*pijk
      preijk=pijk
c      if (iphs3.eq.2) xni=0.1453*a*gamma
      if (iphs3.eq.2) xni=a*gamma
      if (tq.gt.0.01) go to 540
      if (ncyc.gt.40) go to 540
      kr=k+1
      ijkr=nq*(ijbar*(kr-1)+ibar*(j-1)+i-1)+1
c       if (capn(ijkr).ne.0..and.capn(ijk).eq.0.) write (156,*) 'z<',i,j,k
  540 continue
      if (xni.gt.1.e-8.and.modcyc.eq.0) then
c Type 3 cells are not printed !!! *******************************
c     write (10,1330) i,j,k,z,xi,pq,preijk,sq   ,uijk,vijk,wijk,tq
c Type 3 cells are not printed !!! *******************************
        if (kdebug.ne.0) then
          write (14,830) i,j,k,z,xi,pq,tq
          write (12,820) i,j,k,z,xi,pq,preijk,sq,uijk,vijk,wijk,tq
          if (k.eq.kmid) write (113,*) i,j,xi,tq
c           write (156,780) i,j,k,ijk,z,xi,pq,preijk,sq,uijk,vijk,wijk,tq
          if (j.eq.1) write (13,820) i,j,k,z,xi,pq,preijk,sq,uijk,vijk
     1    ,wijk,tq
        endif
        totn=totn+2*xi*dxdydz*gammaf
      endif
      if (j.eq.1) aplot(i,k)=tq
      if (j.eq.1) bplot(i,k)=sq
      if (j.eq.1) cplot(i,k)=pq
c
      ko=ko+1
      if (ko.gt.895195) go to 550
      storev(ko)=float(ijk)+.1
      ko=ko+1
      storev(ko)=uijk
      ko=ko+1
      storev(ko)=vijk
      ko=ko+1
      storev(ko)=wijk
      numko=ko/4
      go to 560
  550 continue
  560 continue
  570 continue
  580 continue
  590 continue
      if (modcyc.eq.0) then
        if (kdebug.ne.0) write (14,850) totn,keqost,ibar,jbar,kbar
        if (kdebug.ne.0) write (12,850) totn,keqost,ibar,jbar,kbar
      endif
c
      kphs=5
      ko=0
      if (numko.ne.numext) write (*,*) ' phase3 NUMEXT error ',numext
     1 ,numko
c       if (numko.ne.numext) write (156,*) ' phase3 NUMEXT error ',numext
c      1 ,numko
c      if (numko.ne.numext) stop
      do 600 i=1,numext
      ko=ko+1
      if (ko.gt.895194) go to 600
      ijk=int(storev(ko))
      if (ijk.eq.0) go to 600
      ko=ko+1
      u(ijk)=storev(ko)
      ko=ko+1
      v(ijk)=storev(ko)
      ko=ko+1
      w(ijk)=storev(ko)
      if (u(ijk)**2+v(ijk)**2+w(ijk)**2 .gt. 1.) then
      write (777,*) 'v^2>1 in line 5554',i,j,k,u(ijk),v(ijk),w(ijk)
      wijk=0.
      uijk=0.
      vijk=0.
      w(ijk)=0.
      u(ijk)=0.
      v(ijk)=0.
c        pause 6040
      endif
  600 continue
c      write (156,*) 'fastest cell=',fast
  610 continue
c










      xmaxpress=0.
      xmaxenerf=0.
      xmaxenerlf=0.
      maxii=0
      maxjj=0
      maxkk=0
      maxpii=0
      maxpjj=0
      maxpkk=0

      do i=1,ibar
      do j=1,jbar
      do k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      xgamma=sqrt(1.-w(ijk)**2-u(ijk)**2-v(ijk)**2)
      enerf=(e(ijk)+p(ijk))*xgamma**2-p(ijk)
      if (xmaxenerf .lt. enerf) then
      xmaxenerf=enerf
      maxii=i
      maxjj=j
      maxkk=k
      endif
      if (xmaxpress .lt. p(ijk)) then
      xmaxpress=p(ijk)
      maxpii=i
      maxpjj=j
      maxpkk=k
      endif
      enddo
      enddo
      enddo

      ijk_e0max=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1)+maxii-1)+1

      ijk_pmax=nq*(ijbar*(maxpkk-1)+ibar*(maxpjj-1)+maxpii-1)+1
      xgamma=sqrt(1.-w(ijk_pmax)**2-u(ijk_pmax)**2-v(ijk_pmax)**2)
      enerf=(e(ijk_pmax)+p(ijk_pmax))*xgamma**2-p(ijk_pmax)

      write (733,*) 'REST FRAME',ncyc,maxii,maxjj,maxkk,xmaxenerf,
     1 p(ijk_e0max),'   '
     2 ,maxpii,maxpjj,maxpkk,enerf,xmaxpress,xmaxpress-p(ijk_pmax)

      write (732,*) 'REST FRAME',ncyc,maxii,maxjj,maxkk,xmaxenerf

c      ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1)+maxii-1)+1
c       write (732,*) ncyc,maxii,maxjj,maxkk,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       write (732,*) '------------------------------------------------'
c       write (732,*) '          NEIGHBOURS TRANSVERSE PLANE           '
c c     ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c       ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1)+maxii+1-1)+1
c       write (732,*) ncyc,maxii+1,maxjj,maxkk,ncyc,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj+1-1)+maxii+1-1)+1
c       write (732,*) ncyc,maxii+1,maxjj+1,maxkk,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1-1)+maxii+1-1)+1
c       write (732,*) ncyc,maxii+1,maxjj-1,maxkk,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1-1)+maxii-1)+1
c       write (732,*) ncyc,maxii,maxjj-1,maxkk,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj+1-1)+maxii-1)+1
c       write (732,*) ncyc,maxii,maxjj+1,maxkk,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1-1)+maxii-1-1)+1
c       write (732,*) ncyc,maxii-1,maxjj-1,maxkk,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1)+maxii-1-1)+1
c       write (732,*) ncyc,maxii-1,maxjj,maxkk,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj+1-1)+maxii-1-1)+1
c       write (732,*) ncyc,maxii-1,maxjj+1,maxkk,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       write (732,*) '                                                '
c       write (732,*) '                 REACTION PLANE                 '
c       write (732,*) '                                                '
c c      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c       ijk=nq*(ijbar*(maxkk+1-1)+ibar*(maxjj-1)+maxii-1-1)+1
c       write (732,*) ncyc,maxii-1,maxjj,maxkk+1,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk+1-1)+ibar*(maxjj-1)+maxii-1)+1
c       write (732,*) ncyc,maxii,maxjj,maxkk+1,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk+1-1)+ibar*(maxjj-1)+maxii+1-1)+1
c       write (732,*) ncyc,maxii+1,maxjj,maxkk+1,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk-1-1)+ibar*(maxjj-1)+maxii-1-1)+1
c       write (732,*) ncyc,maxii-1,maxjj,maxkk-1,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk-1-1)+ibar*(maxjj-1)+maxii-1)+1
c       write (732,*) ncyc,maxii,maxjj,maxkk-1,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       ijk=nq*(ijbar*(maxkk-1-1)+ibar*(maxjj-1)+maxii+1-1)+1
c       write (732,*) ncyc,maxii+1,maxjj,maxkk-1,e(ijk),capn(ijk),
c      1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
c       write (732,*) '#################################################'

      maxii=maxpii
      maxjj=maxpjj
      maxkk=maxpkk

      ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1)+maxii-1)+1
      write (732,*) ncyc,maxii,maxjj,maxkk,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      write (732,*) '------------------------------------------------'
      write (732,*) '          NEIGHBOURS TRANSVERSE PLANE           '
c     ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1)+maxii+1-1)+1
      write (732,*) ncyc,maxii+1,maxjj,maxkk,ncyc,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj+1-1)+maxii+1-1)+1
      write (732,*) ncyc,maxii+1,maxjj+1,maxkk,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1-1)+maxii+1-1)+1
      write (732,*) ncyc,maxii+1,maxjj-1,maxkk,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1-1)+maxii-1)+1
      write (732,*) ncyc,maxii,maxjj-1,maxkk,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj+1-1)+maxii-1)+1
      write (732,*) ncyc,maxii,maxjj+1,maxkk,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1-1)+maxii-1-1)+1
      write (732,*) ncyc,maxii-1,maxjj-1,maxkk,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj-1)+maxii-1-1)+1
      write (732,*) ncyc,maxii-1,maxjj,maxkk,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk-1)+ibar*(maxjj+1-1)+maxii-1-1)+1
      write (732,*) ncyc,maxii-1,maxjj+1,maxkk,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      write (732,*) '                                                '
      write (732,*) '                 REACTION PLANE                 '
      write (732,*) '                                                '
c      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      ijk=nq*(ijbar*(maxkk+1-1)+ibar*(maxjj-1)+maxii-1-1)+1
      write (732,*) ncyc,maxii-1,maxjj,maxkk+1,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk+1-1)+ibar*(maxjj-1)+maxii-1)+1
      write (732,*) ncyc,maxii,maxjj,maxkk+1,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk+1-1)+ibar*(maxjj-1)+maxii+1-1)+1
      write (732,*) ncyc,maxii+1,maxjj,maxkk+1,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk-1-1)+ibar*(maxjj-1)+maxii-1-1)+1
      write (732,*) ncyc,maxii-1,maxjj,maxkk-1,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk-1-1)+ibar*(maxjj-1)+maxii-1)+1
      write (732,*) ncyc,maxii,maxjj,maxkk-1,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      ijk=nq*(ijbar*(maxkk-1-1)+ibar*(maxjj-1)+maxii+1-1)+1
      write (732,*) ncyc,maxii+1,maxjj,maxkk-1,e(ijk),capn(ijk),
     1 p(ijk),w(ijk),u(ijk),v(ijk),zmom(ijk),ymom(ijk),xmom(ijk)
      write (732,*) '#################################################'



c      sumN0=0
c      sumT00=0
c      sumT0Z=0
c      sumpress=0
c      sumwijk=0
c      do i=1,ibar
c      do j=1,jbar
c      do k=1,kbar
c      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c      sumN0=sumN0+capn(ijk)
c      sumT00=sumT00+e(ijk)
c      sumT0Z=sumT0Z+zmom(ijk)
c      sumpress=sumpress+p(ijk)
c      sumwijk=sumwijk+abs(w(ijk))
c      enddo
c      enddo
c      enddo
c      write(*,*) 'N0phas3=',sumN0*dxdydz,'T00=',sumT00*dxdydz,'sumT0Z=',
c     1 sumT0Z*dxdydz,sumpress*dxdydz,sumwijk
c      pause
c     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      return
  620 continue
      write (*,880) kstop,t,ncyc,i,j,k,numit,a,b,c,xi,x
c       write (156,880) kstop,t,ncyc,i,j,k,numit,a,b,c,xi,x
      stop 6
c
  630 format (/,'phase3: capn<0   i=',i6,i8,3i4,i9,'  capn=',f10.7,' npt
     1=',f9.6,6f10.6,/)
  640 format ('phase3: capn=0, Missing e',i8,2f12.7,i7,2f12.7)
  650 format ('Phase3: error-koost1',3i4,i8,i9)
  660 format ('Phase3: error in b',3f12.8,' a,c=',2f12.8)
  670 format ('Phase3: error in c',3f12.8)
  680 format (3i8,' Phase3: pq=',f12.7,'  preijk=',f12.7)
  690 format (3i4,' Phase3: pijk',5f12.6)
  700 format ('KOAST-Phase3:',i4,2i10,6f20.12)
  710 format (4i4,6f12.5,f7.1,7f12.5)
  720 format ('P DIFF-phase3',3i4,i8,i2,2f10.4,1x,2f12.3,8f10.4)
  730 format (' rtmi',f13.9,' xv=',3f13.9,f15.11)
  740 format (' gxx',3f15.10,' eps',2f14.9,' P',2f14.8,2f15.12)
  750 format ('    usq',f12.9,' gg',f12.8,' eps=',f12.8,' P=',f11.8,5f13
     1 .9)
  760 format (/i2,3i4,' T=',f7.2,' S=',f6.2,f7.2,'  g2',f9.6,'  eps',f9.
     1 4,' P',f9.3,/)
  770 format ('captot=',2f9.3,f10.5,'  eptot',3f11.3,'   ko,numext',4i6)
  780 format (' ext',3i4,i9,10f11.5)
  790 format (' new cell ',i6,3f10.6,e15.8)
  800 format (' SUM of changes ',5f15.10)
  810 format ('# Time(ncyc)= ',i4,' Cell-type: ',i2/
     1 '# i-x j-y k-z     T00        N0            p','           vz
     2         vx          vy        T0Z         T0Y             T0X        
     3       T           s            mub            '/
     4'#',10x,'   [MeV/fm3]','    [1/fm3]    [MeV/fm3]        [c]      
     5     [c]         [c]     [MeV/fm3]    [MeV/fm3]      [MeV/fm3]
     6  [MeV]        [fm-3]      [MeV]')
c  820 format (3i4,f10.4,f10.5,3f9.1,' ',3f8.4,f7.1,4f9.3)
  820 format (3i4,' ',12(E12.5,' '))
  821 format (3i4,' ',6(E12.5,' '))
  830 format (3i4,2f9.4,2f9.1,3f8.1,' ',3f9.4,f7.1,f9.3,i4)
  840 format (2i4,5f12.5,3f8.1)
  850 format ('# Tot.N = ',f10.3,'  keqost= ',i4,3x,' xmax= ',i4,' ymax=
     1 ',i4,' zmax= ',i4)
  860 format (' expansion has begun')
  870 format (' empty cell  ??',3i3,i6,6f12.8)
  880 format (' stop no.',i2,' at t=',e12.5,' cyc',i4,' cell (',i3,',
     1',i3,',',i3,')...',/,' numit=',i4,' a=',e16.9,' b=',e16.9,' c=,'
     2 ,e16.9,' xi=',e16.9,' x=',e16.9)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine phase2 (mleave,zpmin,zpmax,enmaxc)
c
c   Taken from HYDRO8E.6.F   20 June 2007   Frankfurt
c   Tidied 21 June 2007
c   NOTE THAT FROM SUBROUTINE phase1, ENERGY AND MOMENTA ARE GIVEN IN 
C   MeV AND capn IS MULTIPLIED BY dxdydz. b IS GIVEN IN MeV*fm-3 AND 
C   cmom AND c ARE GIVEN IN MeV2*fm-3
c
      implicit real*8 (a-h,o-z)
      real*4 start
      integer*2 koost
c
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfz
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
c
      common /gss/ sigx, sigy, sigz, radx, rady, radz, kspher
c
      common /transp/ xte, yte, zte, hpx, hmx, hpy, hmy, hpz, hmz, i, j,
     1 k, l, ueff, veff, weff, cneff, xmeff, ymeff, zmeff, eneff, koko
c
      common /midle/ imid, kmid, kend, imidpt, kmidpt, klmin, krmax,
     1 kdisp, ipmin, ipmax, imid1
c
      common /numcmm/ numbk, numcel, numcst, numext, numfz1, numlv,
     1 numout, numouu, numou2, jdim1, kfreez, kfrzcy, press, freez,
     2 compmy, koast, kext, cpntot, compmx, compe,iundee(58)
c
      common /other/ eusual, epn0, t0, epion, ppion, edelt, delta, pdelt
     1 , temp, tempkn, eps, xth, epsc, dnsmax, ncycpt, ncycrp, ncycaf
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233), px(85877055), py(85877055), volm(8), volmn(8), pz(85877055)
c
      common /bucket/ lablbk(1527120), xmomcl(1527120), ymomcl(1527120),
     1 zmomcl(1527120), encell(1527120), capncl(1527120), cnvvl(1527120)
c
      common /spcl/ delx, delz
c
      common /bccket/ lablbc(1527120), celvol(1527120)
      common /zavr/ zavmid, zavall
c
      common /backet/ lbblbk(22), cpnncl(22), convol(22)
c
      common /volbuk/ lblvol(1527120), retvol(1527120)
      common /leaves/ xmomlv(300), ymomlv(300), zmomlv(300)
c
      common /xposn/ itrck(30000000)
c
      common /plcm7/ kextrn(30542411), koost(30542411)
c
      common /timing/ diffp1, diffp2, diffp3, diffp4, diffp5, start
c
      dimension irx(8), irxbk(8), irx2(8), irxbk2(8)
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
c
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c
c       if (kdebug.ne.0) write (156,*) ' Begin Phase 2, ncyc=',ncyc,' ______
c      1_____________'
      diff2a=secnds(start)
      kphs=2
      iwr=0
      fimid=imid
      fkmid=kmid
      aimid=0
      akmid=0.
      zpart=0.
      cpmax=0.
      kleave=0
      zpmin=1000.
      zpmax=-1000.
      xpmin=1000.
      xpmax=-1000.
      ypmax=-1000.
      xmax=0.
      xmo=0.
      ymo=0.
      zmo=0.
      do 10 i=1,1527120
      lablbk(i)=-1
      lablbc(i)=-1
      xmomcl(i)=0.
      ymomcl(i)=0.
      zmomcl(i)=0.
      encell(i)=0.
      capncl(i)=0.
   10 continue
      do 40 i=1,ibar
      do 30 j=1,jbar
      do 20 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
   20 continue
   30 continue
   40 continue
      ihalf=imid-ktnc/2+2
   50 continue
c   ___________________________ coding for turbulence figures
c
      do 80 ichy=1,2
      ko=-2
   60 ko=ko+3
      if (ko.gt.kp) go to 70
      kopx=ko/3+1
      if (itrck(kopx).ne.ichy) go to 60
      xte=px(ko)
      yte=py(ko)
      zte=pz(ko)
      if (yte.gt.1.0) go to 60
      go to 60
   70 continue
   80 continue
      diff2b=secnds(start)
      diffp1=diffp1+diff2b-diff2a
      kpn=1
      move=0
      sumdel=0.
      smx=0.
      smz=0.
      kol=0
      komx=0
      rtmx=0.
      zaver=0.
      zaver1=0.
      izaver=0
      izave1=0
      zave=0.
      izcl=0
c   _______________________________________________
      ko=-2
   90 ko=ko+3
c
      if (ko.gt.kp) go to 150
      xte=px(ko)
      yte=py(ko)
      zte=pz(ko)
c   ___________________________ coding for turbulence figures
      i=xte+1.
      j=yte+1.
      k=zte+1.
      if (j.eq.1) zaver=zaver+zte
      if (i.eq.imid.and.j.eq.1) zaver1=zaver1+zte
      if (j.eq.1) izaver=izaver+1
      if (i.eq.imid.and.j.eq.1) izave1=izave1+1
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (ijk.gt.lcmb) go to 290
      l=1
      if (ko.gt.kpbas(1).and.ko.le.kpbas(2)) l=2
      if (ko.gt.kpbas(3).and.ko.le.kpbas(4)) l=2
      if (xte.gt.900) go to 140
      ax=xte-aint(xte)
      ay=yte-aint(yte)
      az=zte-aint(zte)
      azo=zte-aint(zte)
      xn=xte
      yn=yte
      zn=zte
      koko=ko
      if (kslab.ne.5) go to 100
      rrx=((xte-float(imid))/sigx)**2+((zte-float(kmid))/sigz)**2
c
      if (rrx.le.rtmx) go to 100
      komx=ko
      rtmx=rrx
  100 continue
c
c       if (capn(ijk).eq.0.) write (156,*) 'capn zero-first',ijk,i,j,k,ko
c      1 ,kp
c
      diff3a=secnds(start)
      call wtnwt (jdim1,interv)
      diff3b=secnds(start)
      diffp3=diffp3+diff3b-diff3a
      ratold=pnu(l)/capn(ijk) ! dimensionless since capn and pnu=dpt are
      num=1./ratold+.15
      dele=ratold*e(ijk) ! [MeV]
      delmx=ratold*xmom(ijk) ! [MeV]
      delmy=ratold*ymom(ijk) ! [MeV]
      delmz=ratold*zmom(ijk) ! [MeV]
      xtm=xte-fimid
      ytm=yte
      ztm=zte-fkmid
      aimid=aimid+xtm
      akmid=akmid+ztm
c
c   Compute angular momentum
c
      if (xtm.lt.0.0.or.ztm.lt.0.0) go to 110
      xmo=xmo+(ymom(ijk)*ztm-zmom(ijk)*ytm)*ratold
      ymo=ymo+(zmom(ijk)*xtm-xmom(ijk)*ztm)*ratold
      zmo=zmo+(delmy*xtm-delmx*ytm)
      zpart=delmy*xtm-delmx*ytm
      zport=zport+zpart
  110 continue
      if (dabs(weff).lt.0.0000001) weff=0.
      xtecl=(-xte+fimid)*dx/convpu
      ytecl=(yte)*dy/convpu
      ztecl=(-zte+fkmid+1.)*dz/convpu
      if (yte.le.0.999) zave=zave+ztecl
      if (yte.le.0.999) izcl=izcl+1
      if (ioptos1.ne.0) then
        if (xtecl.gt.xmax) xmax=xtecl
      endif
      xn=xte+ueff*dtodx
      yn=yte+veff*dtody
      zn=zte+weff*dtodz
      veloc=ueff*ueff+veff*veff+weff*weff
      in=xn+1.
      jn=yn+1.
      kn=zn+1.
      ijkn=nq*(ijbar*(kn-1)+ibar*(jn-1)+in-1)+1
      ktrak=0
  120 continue
      if (dabs(zn-zte).gt.zmovemx) zmovemx=abs(zn-zte)
      if (dabs(xn-xte).gt.xmovemx) xmovemx=abs(xn-xte)
      if (capn(ijk).gt.cpmax) cpmax=capn(ijk)
      if (zpmin.gt.zn) zpmin=zn
      if (zpmax.lt.zn) zpmax=zn
      if (xpmin.gt.xn) xpmin=xn
      if (xpmax.lt.xn) xpmax=xn
      if (ypmax.lt.yn) ypmax=yn
      ax=abs(xn-aint(xn))
      ay=abs(yn-aint(yn))
      az=abs(zn-aint(zn))
      if (yn.lt.0.) yn=abs(yn)+1.e-4
      if (xn.lt.0..or.xn.gt.fibar) go to 160
      if (yn.gt.fjbar) go to 160
      if (zn.lt.0..or.zn.gt.fkbar) go to 160
      dist=yte*yte+(xte-imid)**2
      sdist=dsqrt(dist)
      anumb1=capn(ijk)/pnu(1)
      xdif=xn-x
      if (interp.eq.2) go to 190
      if (ijkn.ne.ijk) go to 170
  130 continue
      px(kpn)=xn
      py(kpn)=yn
      pz(kpn)=zn
  140 continue
      kpn=kpn+3
      go to 90
  150 continue
      diff2c=secnds(start)
      diffp2=diffp2+diff2c-diff2b
      kp=kpn-1
      go to 280
  160 leave=1
      mleave=mleave+1
c       if (mleave.le.100) write (156,*) ' LEAVE SET',ko,i,j,k,xte,yte,zte
c      1 ,xn,yn,zn,fibar,fjbar,fkbar
      if (mleave.le.100) write (*,*) ' LEAVE SET',ko,i,j,k,xte,yte,zte
     1 ,xn,yn,zn
c       if (mleave.le.199.and.xn.lt.0..or.xn.gt.fibar) write (156,*) '
c      1 xn,fibar ',xn,fibar
c       if (mleave.le.199.and.yn.gt.fjbar) write (156,*) ' yn,fjbar ',yn
c      1 ,fjbar
c       if (mleave.le.199.and.zn.lt.0..or.zn.gt.fkbar) write (156,*) '
c      1 zn,fkbar ',zn,fkbar,mleave
  170 continue
      if (capn(ijk).ne.0.) go to 180
c       if (e(ijk).eq.0.) write (156,380) i,j,k,ko,capn(ijk),e(ijk)
      if (e(ijk).eq.0.) write (*,380) i,j,k,ko,capn(ijk),e(ijk)
      go to 130
  180 ratold=pnu(l)/capn(ijk)
      num=1./ratold+.15
      dele=ratold*e(ijk)
      delmx=ratold*xmom(ijk)
      delmy=ratold*ymom(ijk)
      delmz=ratold*zmom(ijk)
      call strmom (ijk,idzm,-delmx,-delmy,-delmz,-dele,-pnu(l))
      if (mleave.eq.1) go to 260
      call strmom (ijkn,idzm,delmx,delmy,delmz,dele,pnu(l))
      go to 130
c
c  interp=2 - linear assignment
  190 continue
      do 200 ml=1,19
      lbblbk(ml)=-1
      cpnncl(ml)=0.
  200 continue
      iplc=ijkret(ijk)
      done=dfloat(1)
      write(*,*) 'into area1'
      call area (-done,volm,irx,irxbk,i,j,k,l,xte,yte,zte)
      write(*,*) 'ino area2'
c      pause 4764
c      write (156,1470) (volm(ijkl),ijkl=1,8)
c      write (156,1450) (irx(ijkl),ijkl=1,8)
c      write (156,*) xn,yn,zn,weff
      call area (done,volmn,irx2,irxbk2,in,jn,kn,l,xn,yn,zn)
      su=volmn(1)+volmn(2)+volmn(3)+volmn(4)+volmn(5)+volmn(6)+volmn(7)
     1 +volmn(8)
      sumad=0.
      sumlv=0.
      sump=0.
      sume=0.
      sumx=0.
      sumy=0.
      sumz=0.
      sumpv=0.
      do 230 ml=1,19
      if (cpnncl(ml).eq.0.) go to 230
      if (cpnncl(ml).lt.0.) go to 210
      iplcml=1
      sumad=sumad+cpnncl(ml)
      go to 220
c
c  Identify the amount of matter leaving each cell.
c  Distribute according to the ratio of volumes
c
  210 continue
      iplcml=ijkret(lbblbk(ml))
      sumlv=sumlv+cpnncl(ml)
      ratio=cpnncl(ml)/retvol(iplcml)
      sump=sump+ratio*capn(lbblbk(ml))
      sume=sume+ratio*e(lbblbk(ml))
      sumz=sumz+ratio*zmom(lbblbk(ml))
      sumy=sumy+ratio*ymom(lbblbk(ml))
      sumx=sumx+ratio*xmom(lbblbk(ml))
  220 continue
  230 continue
c
c  Move the quantities associated with the motion
c
      summov=0.
      do 250 ml=1,19
      if (cpnncl(ml).eq.0.) go to 250
      if (cpnncl(ml).lt.0.) go to 240
c  Add
      ratio=-cpnncl(ml)/sumad
c       write (156,330) ijkm,lbblbk(ml),ratio,cpnncl(ml),ratio*sump
      ijkm=lbblbk(ml)
      call strmom (ijkm,irz,ratio*sumx,ratio*sumy,ratio*sumz,ratio*sume
     1 ,ratio*sump)
      summov=summov+ratio*sump
      go to 250
c
c  Subtract
  240 continue
      iplcml=ijkret(lbblbk(ml))
      ratio=cpnncl(ml)/retvol(iplcml)
      qn=ratio*capn(lbblbk(ml))
      qe=ratio*e(lbblbk(ml))
      qz=ratio*zmom(lbblbk(ml))
      qy=ratio*ymom(lbblbk(ml))
      qx=ratio*xmom(lbblbk(ml))
      ijkm=lbblbk(ml)
      call strmom (ijkm,irz,qx,qy,qz,qe,qn)
      summov=summov+ratio*sump
  250 continue
      go to 130
c
c  --- particle leaves system; calculate edge loss ---
c
  260 game=1./dsqrt(1.-(ueff**2+veff**2+weff**2))
      xmomot=xmomot+delmx
      zmomot=zmomot+delmz
      ek=amin1(emid3,econ*(game-1.))
      eps=ek/game+eusual
      xi=pnu(1)
c      epsc=epn(xi)-epn0
      epsc=0.
      xth=(eps/xi-eusual)*econ-epsc
c  xth is the thermal energy density
      if (xth.lt.0.) xth=0.
      xarg=0.086673*xth/xi**.6666666666666
      capnlf=capnlf+pnu(l)
      cnth=pnu(l)
      kleave=kleave+1
      if (kleave.eq.41) write (*,*) ' kleave exceeded 40',kleave,ncyc
c
c  calculate the contribution of the particle to the rapidity
c
      uijk=delmx/dele
      vijk=delmy/dele
      wijk=delmz/dele
      gamma=1./dsqrt(dabs(1.-uijk*uijk-vijk*vijk-wijk*wijk))
      xi=pnu(l)/gamma*.145*rdxdydz
      y=(1.+wijk)/(1.-wijk)
      y=0.5*dlog(y)
      if (uijk**2+vijk**2+wijk**2 .gt. 1.) then
      write (*,*) 'v^2>1'
      stop 6670
      endif
c       if (kleave.lt.4) write (156,370) pnu(l),xn,yn,zn,wijk,y,xi,gamma,
c      1 wijk
      if (kleave.gt.300) go to 270
      xmomlv(kleave)=delmx*rdxdydz
      ymomlv(kleave)=delmy*rdxdydz
      zmomlv(kleave)=delmz*rdxdydz
  270 continue
      x=ueff**2+veff**2+weff**2
      cmom=delmx**2+delmy**2+delmz**2
c       if (kleave.lt.10) write (156,390) kleave,i,j,k,pnu(l),xmomot,
c      1 zmomot,dele,delmx,delmy,delmz,ratold,wijk,gamma
      leave=0
      go to 90
c
c  ___________________________________________________________
c
  280 continue
c      if (kslab.eq.5) write(*,*) 'komx=',komx
c       if (kslab.eq.5) write (156,350) komx,dsqrt(rtmx),px(komx),pz(komx)
c      1 ,px(komx)-imid,pz(komx)-kmid
c
      if(izaver.ne.0.)zavmid=zaver/float(izaver)
      if(izcl.ne.0)zzzav=zave/float(izcl)
      zlave=0.
      if(izave1.ne.0)zlave=zaver1/float(izave1)
c       if (mod(ncyc,50).eq.1) write (156,340) zaver,izaver,zaver1,izave1
c      1 ,zavmid,dx/convpu/2.,zlave,zave,izcl,zzzav
c       if (kdebug.ne.0.and.mod(ncyc,5).eq.1) write (156,*) 'xmax ',xmax
c       if (kdebug.ne.0) write (156,300) float(ncyc)*dt,xmo,ymo,zmo,dt
      return
  290 continue
      write (156,*) 'Phase2 error: ijk too large',i,j,k,ko,ijk,ijkmax
      write (*,*) 'Phase2 error: ijk too large',i,j,k,ko,ijk,ijkma
      stop 9206
c
  300 format (10x,'Ang. Momentum (particle)',6f13.5,/)
  310 format (2f9.3)
  320 format (3i4,2f9.2,f9.2)
  330 format ('strmom2',2i9,5f10.5)
  340 format (3x,' zaver',f13.2,i8,f12.2,i7,3f8.3,' zave',f8.2,i8,f8.3)
  350 format ('greatest extent',i8,' rtmx',5f9.3)
  360 format (2f10.4,i9)
  370 format (' dN for particle leaving mesh',10f10.4)
  380 format (' E zero in phase2 ',3i4,' ',i8,2e20.9)
  390 format (' leave =',i4,3i3,10e11.4)
  400 format (2f6.2,f7.1,f7.3,4f8.3)
      end
   
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine phase1 (mleave,kpbar)
c
      implicit real*8 (a-h,o-z)
      integer*2 koost
c
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfz
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
c
      common /midle/ imid, kmid, kend, imidpt, kmidpt, klmin, krmax,
     1 kdisp, ipmin, ipmax, imid1
c
      common /numcmm/ numbk, numcel, numcst, numext, numfz1, numlv,
     1 numout, numouu, numou2, jdim1, kfreez, kfrzcy, press, freez,
     2 compmy, koast, kext, cpntot, compmx, compe,iundee(58)
c
      common /other/ eusual, epn0, t0, epion, ppion, edelt, delta, pdelt
     1 , temp, tempkn, eps, xth, epsc, dnsmax, ncycpt, ncycrp, ncycaf
c
      common /thestf/ rap00, Tfin, dxv, dyv, dzv, epsilon00, b0, rho00,
     1 e00
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c
      common /newcrp/ entot, dentot, tempav, numcll
c
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233), px(85877055), py(85877055), pz(85877055)
c
      common /plcm7/ kextrn(30542411), koost(30542411)
c
      common /eqqcd/ tq, eq1, pq, sq, anq, amuq
c
      common /miscnw/ cnfst, numfst
c
      common /mincl/ imincl, imaxcl, jmincl, jmaxcl, kmincl, kmaxcl
      
      
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
c
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
      external feqcd3
c
      xmomtt=0.
      ymomtt=0.
      xmomta=0.
      numfst=0
      cnfst=0.
      cntot=0.
c
      do 30 i=1,ibar
      do 20 j=1,jbar
      do 10 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
   10 continue
   20 continue
   30 continue
c
      q0=0.
      imincl=1000
      imaxcl=-1000
      jmincl=1000
      jmaxcl=0
      kmincl=1000
      kmaxcl=0
      if (kslab.lt.8) call chksym
      el=0.
      ymomta=0.
      eaddx=0.
      sume=0.
      numcel=0
      if (ioptos3.ne.0) write (50,*) ' ncyc=',ncyc,'  N, E, p, w, u, v,
     1zmom, ...'
      do 60 i=1,ibar
      do 50 j=1,jbar
      do 40 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 40
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      numcel=numcel+1
      sume=sume+e(ijk)
      if (i.lt.imincl) imincl=i
      if (i.gt.imaxcl) imaxcl=i
      if (k.lt.kmincl) kmincl=k
      if (k.gt.kmaxcl) kmaxcl=k
      if (j.lt.jmincl) jmincl=j
      if (jmaxcl.lt.j) jmaxcl=j
      if (j.gt.1) go to 40
      gamma=1./dsqrt(1.-(u(ijk)**2+v(ijk)**2+w(ijk)**2))
      xnii=capn(ijk)/gamma ! maybe baryon density in the rest frame?
c     I am not sure what capn(ijk) is. It should be N0 but I have not
c     found any place in the code where capn(ijk) is defined. (december,
c     2024)
c  ioptos3 __
   40 continue
   50 continue
   60 continue
      sume=sume*dxdydz
c       if (kdebug.gt.0) write (156,670) imincl,imaxcl,ibar-imaxcl,jmaxcl
c      1 ,jbar-jmaxcl,kmincl,kmaxcl,kbar-kmaxcl
      orige=sume
c       if (kdebug.gt.0) write (156,*) ' Begin Phase 1, ncyc=',ncyc,
c      1 ' _____'
c
      do 120 i=1,ibar
      do 110 j=1,jbar
      do 100 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 100
      b=e(ijk)*rem0n0 ! [MeV*fm-3] rem0n0=1
      cmom=xmom(ijk)**2+ymom(ijk)**2+zmom(ijk)**2 ! [Mev2*fm-6]
      c=cmom*remnsq ! remnsq=1
c       if (ncyc.le.0.and.c.gt.b*b) write (156,660) ncyc,i,j,k,ijk,b,c,b*b
  100 continue
  110 continue
  120 continue
      if (kdisp.eq.1) go to 170
      if (mleave.gt.0) go to 170
      if (kslab.eq.8) go to 170 ! kslab=8
c   Perform checks
      do 160 i=1,ibar
      do 150 j=1,jbar
      do 140 k=1,kbar/2
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      ib=ibar+1-i
      kb=kbar+3-k
      if (kslab.eq.4.or.kslab.eq.5.or.kslab.eq.7) kb=kb-2
c  check symmetry
      ijkb=nq*(ijbar*(kb-1)+ibar*(j-1)+ib-1)+1
      if (capn(ijk).eq.0..and.capn(ijkb).eq.0.) go to 140
c
c       if (dabs(capn(ijk)-capn(ijkb)).gt.1.e-8) write (156,510) i,j,k,ib,
c      1 kb,ijk,ijkb,capn(ijk),capn(ijkb),ibar,jbar,kbar
      if (dabs(capn(ijk)-capn(ijkb)).gt.1.e-8) stop 13
c
c       if (dabs(p(ijk)-p(ijkb)).gt.1.e-9*p(ijk)) write (156,680) i,ib,j,k
c      1 ,kb,p(ijk),p(ijkb),p(ijk)-p(ijkb),ijk,ijkb
c
c       if (dabs(e(ijk)-e(ijkb)).gt.e(ijk)*.000000000001) write (156,500)
c      1 i,j,k,kb,e(ijk),e(ijkb),ijk,ijkb,abs(e(ijk)-e(ijkb)),e(ijk)*
c      2 .00000001,abs(e(ijk)-e(ijkb))-e(ijk)*.00000001
c
      if (kslab.eq.8) go to 130
      if (dabs(xmom(ijk)+xmom(ijkb)).gt..000001*abs(xmom(ijk))) write (
     1 156,530) i,j,k,kb,xmom(ijk),xmom(ijkb),xmom(ijk)+xmom(ijkb)
c
      if (dabs(ymom(ijk)-ymom(ijkb)).gt..00000001*abs(ymom(ijk))) write
     1 (6,520) i,j,k,ib,ymom(ijk),ymom(ijkb),ymom(ijk)-ymom(ijkb)
c
c       if (dabs(zmom(ijk)+zmom(ijkb)).gt.1.e-9) write (156,540) i,j,k,kb
c      1 ,zmom(ijk),zmom(ijkb),zmom(ijk)+zmom(ijkb),ijk,ijkb,e(ijk),e(ijkb
c      2 )
c
      if (e(ijk).lt.0.) write (*,*) 'phs1 XXXXX1  E<0',ijk,i,j,k,e(ijk)
c       if (e(ijk).lt.0.) write (156,*) 'phs1 XXXXX1  E<0',ijk,i,j,k,
c      1 e(ijk)
c
c       if (dabs(w(ijk)+w(ijkb)).gt.1.e-11) write (156,690) i,ib,j,k,kb,w
c      1 (ijk),w(ijkb),abs(w(ijk)+w(ijkb)),ijk,ijkb,zmom(ijk),zmom(ijkb),u
c      2 (ijk),u(ijkb)
c
c       if (dabs(u(ijk)+u(ijkb)).gt.1.e-9) write (156,700) i,ib,j,k,kb,
c      1 u(ijk),u(ijkb),abs(u(ijk)+u(ijkb)),ijk,ijkb,xmom(ijk),xmom(ijkb)
c
  130 continue
c
      xdiff=0.5*(xmom(ijk)+xmom(ijkb))
      ydiff=0.5*(ymom(ijk)-ymom(ijkb))
      zdiff=0.5*(zmom(ijk)+zmom(ijkb))
      zmom(ijk)=zmom(ijk)-zdiff
      zmom(ijkb)=zmom(ijkb)+zdiff
      ediff=0.5*(e(ijk)-e(ijkb))
      e(ijk)=e(ijk)-ediff
      e(ijkb)=e(ijkb)+ediff
      pdiff=0.5*(p(ijk)-p(ijkb))
      p(ijk)=p(ijk)-pdiff
      p(ijkb)=p(ijkb)+pdiff
c       if (dabs(pdiff).gt.1.e-10) write (156,710) i,ib,j,k,kb,pdiff,
c      1 p(ijk),p(ijkb),ijk,ijkb
      if (b0.ne.0.) go to 140
c      ib=ibar+1-i
      ijkb=nq*(ijbar*(k-1)+ibar*(j-1)+ib-1)+1
c       if (dabs(capn(ijk)-capn(ijkb)).gt.1.e-9) write (156,*) 'phase1:
c      1 err capn2',i,j,k,' kb',kb,capn(ijk),capn(ijkb)
      if (dabs(capn(ijk)-capn(ijkb)).gt.1.e-9) write (*,*) 'phase1: err
     1capn2',i,j,k,' kb',kb,capn(ijk),capn(ijkb)
c
c       if (dabs(e(ijk)-e(ijkb)).gt.e(ijk)*.000001) write (156,600) i,ib,j
c      1 ,k,e(ijk),e(ijkb),capn(ijk),capn(ijkb),ijk,ijkb
c
      if (e(ijk).lt.0.) write (*,*) 'PHASE1:  E<0',ijk,i,j,k,e(ijk)
c       if (e(ijk).lt.0.) write (156,*) 'PHASE1:  E<0',ijk,i,j,k,e(ijk)
c
c       if (dabs(p(ijk)-p(ijkb)).gt.1.e-10) write (156,640) i,ib,j,k,
c      1 p(ijk),p(ijkb),capn(ijk),capn(ijkb),ijk,ijkb
      if (kslab.eq.8) go to 140
      if (dabs(xmom(ijk)+xmom(ijkb)).gt..000001*abs(xmom(ijk))) write (
     1 156,570) i,ib,j,k,xmom(ijk),xmom(ijkb),capn(ijk),capn(ijkb),ijk,
     2 ijkb
      if (dabs(ymom(ijk)-ymom(ijkb)).gt..000001*abs(ymom(ijk))) write (
     1 156,580) i,ib,j,k,ymom(ijk),ymom(ijkb),capn(ijk),capn(ijkb),ijk,
     2 ijkb
c
      if (dabs(zmom(ijk)-zmom(ijkb)).gt..000001*abs(zmom(ijk))) write (
     1 156,590) i,ib,j,k,zmom(ijk),zmom(ijkb),capn(ijk),capn(ijkb),ijk,
     2 ijkb
     2 ,abs(zmom(ijk)-zmom(ijkb)),zmom(ijk)-zmom(ijkb)
      if (dabs(zmom(ijk)-zmom(ijkb)).gt..000001*abs(zmom(ijk))) write (*
     1 ,590) i,ib,j,k,zmom(ijk),zmom(ijkb),capn(ijk),capn(ijkb),ijk,ijkb
     2 ,abs(zmom(ijk)-zmom(ijkb)),zmom(ijk)-zmom(ijkb)
      if (dabs(zmom(ijk)-zmom(ijkb)).gt..000001*abs(zmom(ijk))) stop 444
c
c       if (dabs(u(ijk)+u(ijkb)).gt.1.e-9) write (156,610) i,ib,j,k,u(ijk)
c      1 ,u(ijkb),capn(ijk),capn(ijkb),ijk,ijkb,xmom(ijk),xmom(ijkb)
c       if (dabs(v(ijk)-v(ijkb)).gt.1.e-9) write (156,620) i,ib,j,k,v(ijk)
c      1 ,v(ijkb),capn(ijk),capn(ijkb),ijk,ijkb
c       if (dabs(w(ijk)-w(ijkb)).gt.1.e-9) write (156,630) i,ib,j,k,w(ijk)
c      1 ,w(ijkb),capn(ijk),capn(ijkb),ijk,ijkb,abs(w(ijk)-w(ijkb))
  140 continue
  150 continue
  160 continue
  170 continue
      kphs=1
      if (kdebug.ne.0) write (156,*) ' kphs',kphs,ioptos2
      cpntot=0.
      smz=0.
      if (mod(ncyc,ncycrp).ne.0) go to 210
      if (ioptos3.eq.0) go to 210
      write (*,*) ' i=imid profiles  capn e p  ncyc=',ncyc
      j=1
      i=imid
      do 190 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 190
  190 continue
      j=1
      k=kmid
      do 200 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 200
  200 continue
  210 continue
c
c  Begin __________________________________________________
c
      numneg=0
      sumneg=0.
      edifft=0.
      do 420 k=1,kbar
      do 410 j=1,jbar
      do 400 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (e(ijk).le.EFLOOR.or.capn(ijk).le.BFLOOR) then
      e(ijk)=0.d0
      capn(ijk)=0.d0
      p(ijk)=0.d0
      xmom(ijk)=0.d0
      ymom(ijk)=0.d0
      zmom(ijk)=0.d0
      u(ijk)=0.d0
      v(ijk)=0.d0
      w(ijk)=0.d0
      go to 400
      endif
      if (capn(ijk).eq.0.) go to 400
      if (capn(ijk).lt.0.00001) then
c         if (e(ijk).gt.0.01) write (156,*) ' ENERGY
c      1  ERR, n=0 ',i,j,k,ijk,capn(ijk),e(ijk)
        if (e(ijk).gt.0.01) write (*,*) ' ENERGY
     1  ERR, n=0 ',i,j,k,ijk,capn(ijk),e(ijk)
      endif
      cntot=cntot+capn(ijk) ! Unknown what capn(ijk) really is
      b=e(ijk)*rem0n0 ! [MeV*fm-3] rem0n0=1
      cmom=xmom(ijk)**2+ymom(ijk)**2+zmom(ijk)**2 ! [MeV2*fm-6]
      cmomol=cmom ! [MeV2*fm-6]
      c=cmom*remnsq ! [MeV2*fm-6] remnsq=1
      numfil=numfil+1 ! numfil is set to zero in subroutine picrel, just
c     before calling this subroutine
      smz=smz+zmom(ijk) ! [MeV*fm-3]
      ipjk=ijk+nq ! nq is set to 3 in main program
      imjk=ijk-nq
      ijpk=ijk+nqib
      ijmk=ijk-nqib
      ijkp=ijk+nqijb
      ijkm=ijk-nqijb
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      ipjk0=ijk0+1
      imjk0=ijk0-1
      ijpk0=ijk0+ibar
      ijmk0=ijk0-ibar
      ijkp0=ijk0+ijbar
      ijkm0=ijk0-ijbar
      adjflx=1.
c  koost determines if a cell is coasting
      icb=koost(ijk0)
      if (icb.eq.1) adjflx=1./dsqrt(u(ijk)**2+v(ijk)**2+w(ijk)**2+.00000
     1 1)
      uijk=u(ijk)*adjflx
      vijk=v(ijk)*adjflx
      wijk=w(ijk)*adjflx
      if (e(ijk).le.EFLOOR.or.capn(ijk).le.BFLOOR) then
      uijk=0.d0
      vijk=0.d0
      wijk=0.d0
      vel2=0.d0
      else
      vel2=uijk**2+vijk**2+wijk**2
      if (isnan(vel2).or.vel2.ge.1.d0) then
      write (*,*) 'v^2>=1 corrected (7034) cyc/i/j/k=',ncyc,i,j,k,
     1 ' vel2=',vel2,' u,v,w=',uijk,vijk,wijk
      if (isnan(vel2).or.vel2.le.0.d0) then
      uijk=0.d0
      vijk=0.d0
      wijk=0.d0
      else
      vfac=dsqrt(0.9999d0/vel2)
      uijk=uijk*vfac
      vijk=vijk*vfac
      wijk=wijk*vfac
      endif
      endif
      endif
      pip=0.
      pim=0.
      pjp=0.
      pjm=0.
      pkp=0.
      pkm=0.
      qkp=0.
      qkm=0.
      qip=0.
      qim=0.
      gamma=dsqrt(dabs(1.-uijk*uijk-vijk*vijk-wijk*wijk)) 
      sound=.0
      if (i.eq.ibar) go to 220
      adjflx=1.
      ncb=koost(ipjk0)
      if (ncb.eq.1) adjflx=1./dsqrt(u(ipjk)**2+v(ipjk)**2+w(ipjk)**2+.00
     1 001)
      uip=.5*(uijk+u(ipjk)*adjflx)
      if (capn(ipjk).eq.0.) go to 230
      pip=.5*(p(ijk)+p(ipjk))
      go to 230
  220 uip=uijk
      pip=p(ijk)
  230 if (i.eq.1) go to 240
      adjflx=1.
      ncb=koost(imjk0)
      if (ncb.eq.1) adjflx=1./dsqrt(u(imjk)**2+v(imjk)**2+w(imjk)**2+.00
     1 01)
      uim=.5*(uijk+u(imjk)*adjflx)
      if (capn(imjk).eq.0.) go to 250
c???
      pim=.5*(p(ijk)+p(imjk))
      go to 250
  240 uim=uijk
      pim=p(ijk)
  250 if (j.eq.jbar) go to 260
      adjflx=1.
      ncb=koost(ijpk0)
      if (ncb.eq.1) adjflx=1./dsqrt(u(ijpk)**2+v(ijpk)**2+w(ijpk)**2+.00
     1 01)
      vjp=.5*(vijk+v(ijpk)*adjflx)
      if (capn(ijpk).eq.0.) go to 270
      pjp=.5*(p(ijk)+p(ijpk))
      go to 270
  260 vjp=vijk
      pjp=p(ijk)
  270 if (j.eq.1) go to 280
      adjflx=1.
      ncb=koost(ijmk0)
      if (ncb.eq.1) adjflx=1./dsqrt(u(ijmk)**2+v(ijmk)**2+w(ijmk)**2+.00
     1 01)
      vjm=.5*(vijk+v(ijmk)*adjflx)
      if (capn(ijmk).eq.0.) go to 290
      pjm=.5*(p(ijk)+p(ijmk))
      go to 290
  280 vjm=0.
      if (jdim1.ne.1) pjm=p(ijk)
  290 if (k.eq.kbar) go to 300
      adjflx=1.
      ncb=koost(ijkp0)
      if (ncb.eq.1) adjflx=1./dsqrt(u(ijkp)**2+v(ijkp)**2+w(ijkp)**2+.00
     1 01)
      wkp=.5*(wijk+w(ijkp)*adjflx)
      if (capn(ijkp).eq.0.) go to 310
      pkp=.5*(p(ijk)+p(ijkp))
c      qkp=q0*dz*sound*capn(ijk)*(wijk-w(ijkp))
c      if (qkp.lt.0.) qkp=0.
      go to 310
  300 wkp=wijk
      pkp=p(ijk)
  310 if (k.eq.1) go to 320
      adjflx=1.
      ncb=koost(ijkm0)
      if (ncb.eq.1) adjflx=1./dsqrt(u(ijkm)**2+v(ijkm)**2+w(ijkm)**2+.00
     1 01)
      wkm=.5*(wijk+w(ijkm)*adjflx)
c
      if (capn(ijkm).eq.0.) go to 330
      pkm=.5*(p(ijk)+p(ijkm))
c      qkm=q0*dz*sound*capn(ijk)*(w(ijkm)-wijk)
      go to 330
  320 wkm=wijk
      pkm=p(ijk)
  330 continue
      xmold=xmom(ijk)
      ymold=ymom(ijk)
      zmold=zmom(ijk)
      eold=e(ijk)
c       if (jdim1.ne.3.and.ymom(ijk).ne.0.) write (156,*) 'ymom',i,j,k,pjm
c      1 ,pjp,pjm-pjp
      zorig=zmom(ijk)
      zdiff=pkm-pkp
      if (dabs(zdiff).lt.1.e-11) zdiff=0.
      eold=e(ijk)
      if (i.eq.imid.and.k.eq.kmid) ymomtt=ymomtt+dtody*(pjm-pjp)
      if (i.le.imid.and.k.eq.kmid) ymomta=ymomta+dtody*(pjm-pjp)
      if (j.eq.1.and.i.le.imid) xmomta=xmomta+dtodx*(pim-pip)
      if (j.eq.1.and.k.eq.kmid.and.i.le.imid) xmomtt=xmomtt+dtodx*(pim
     1 -pip)
      ijkmn=ijk
      ijkmir=ijk
      imn=i
      imx=i
      ediff=dtodx*(pim*uim-pip*uip)+dtody*(pjm*vjm-pjp*vjp)+dtodz*(pkm
     1 *wkm-pkp*wkp)
c      nkx=n2*ijbar
      delex=dtodx*(pim*uim-pip*uip)
      deley=dtody*(pjm*vjm-pjp*vjp)
      delez=dtodz*(pkm*wkm-pkp*wkp)
      deleto=delex+deley+delez
      delx=dtodx*(pim-pip)
      dely=dtody*(pjm-pjp)
      delz=dtodz*(pkm-pkp)
      xmomnw=xmom(ijk)+delx
      ymomnw=ymom(ijk)+dely
      zmomnw=zmom(ijk)+delz
c--
      edifft=edifft+ediff
c--
      enew=e(ijk)+ediff
      b=enew
      cmomol=xmom(ijk)**2+ymom(ijk)**2+zmom(ijk)**2
      cmom=xmomnw**2+ymomnw**2+zmomnw**2
      if (cmom.lt.b*b) go to 380
      numfst=numfst+1
      cnfst=cnfst+capn(ijk)
      go to 380
      cmhalf=dsqrt(cmom)
      cvv=dsqrt(u(ipjk)**2+v(ipjk)**2+w(ipjk)**2)
      xnw=cmom/enew**2
      xol=cmomol/eold**2
      write (*,470) i,j,k,cmom,enew,xnw,cmomol,eold,xol,capn(ijk)
      write (*,490) dtodx*pim,dtodx*pip,xmom(ijk),dtodx*pjm,dtodx*pjp
     1 ,ymom(ijk),dtodx*pkm,dtodx*pkp,zmom(ijk)
      write (*,480) delx,dely,delz,ediff
      adjflx=1.
      if (capn(ipjk).eq.0.) go to 340
      chkp=(xmom(ipjk)**2+ymom(ipjk)**2+zmom(ipjk)**2)
      write (*,*) 'chkp',chkp,e(ipjk)**2,chkp/e(ipjk)**2
      if (chkp.lt.e(ipjk)**2) go to 340
      uhkp=(u(ipjk)**2+v(ipjk)**2+w(ipjk)**2)
      adjflx=1./dsqrt(uhkp+.0001)
  340 uip=.5*(uijk+u(ipjk)*adjflx)
      write (*,*) 'uip',uip,adjflx
      adjflm=1.
      if (capn(imjk).eq.0.) go to 350
      chkm=(xmom(imjk)**2+ymom(imjk)**2+zmom(imjk)**2)
      write (*,*) 'chkm',chkm,e(imjk)**2,chkm/e(imjk)**2
      if (chkm.lt.e(imjk)**2) go to 350
      uhkm=(u(imjk)**2+v(imjk)**2+w(imjk)**2)
      adjflm=1./dsqrt(uhkm+.0001)
  350 uim=.5*(uijk+u(imjk)*adjflm)
      write (*,*) 'uim',uim,adjflm
c   z__
      adjflx=1.
      if (capn(ijkp).eq.0.) go to 360
      chkp=(xmom(ijkp)**2+ymom(ijkp)**2+zmom(ijkp)**2)
      write (*,*) 'chkp',chkp,e(ijkp),chkp/e(ijkp)**2
      if (chkp.lt.e(ijkp)**2) go to 360
      uhkp=(u(ijkp)**2+v(ijkp)**2+w(ijkp)**2)
      adjflx=1./dsqrt(uhkp+.0001)
  360 ukp=.5*(uijk+u(ijkp)*adjflx)
      write (*,*) 'ukp',ukp,adjflx
      adjflm=1.
      if (capn(ijkm).eq.0.) go to 370
      chkm=(xmom(ijkm)**2+ymom(ijkm)**2+zmom(ijkm)**2)
      write (*,*) 'chkm',chkm,e(ijkm)**2,chkm/e(ijkm)**2
      if (chkm.lt.e(ijkm)**2) go to 370
      uhkm=(u(ijkm)**2+v(ijkm)**2+w(ijkm)**2)
      adjflm=1./dsqrt(uhkm+.0001)
  370 ukm=.5*(uijk+u(ijkm)*adjflm)
      write (*,*) 'ukm',ukm,adjflm
      stop 335
c
c--
c--
  380 continue
      xmom(ijk)=xmomnw
      ymom(ijk)=ymomnw
      zmom(ijk)=zmomnw
      e(ijk)=enew
      if (e(ijk).lt.0.) e(ijk)=0.  !Angel
      if (e(ijk).ge.0.) go to 390
c
c       if (kdebug.ne.0) write (156,650) numneg,ijk,i,j,k,e(ijk),eold,
c      1 ediff,capn(ijk)*dxdydz/pnu(1)
      if (kdebug.ne.0) write (*,650) numneg,ijk,i,j,k,e(ijk),eold,ediff
     1 ,capn(ijk)*dxdydz/pnu(1)
      numneg=numneg+1
      sumneg=sumneg+e(ijk)
c
      if (numneg.lt.115) go to 390
c       write (156,*) ' Too many E<0 in PHASE1.F',numneg
      write (*,*) ' Too many E<0 in PHASE1.F',numneg
      stop 55
c
  390 continue
  400 continue
  410 continue
  420 continue
c    
      sume=0.
      if (kdebug.ne.0.and.numneg.gt.0) write (*,*) ' numneg',numneg,' su
     1mneg',sumneg
c       if (kdebug.ne.0.and.numneg.gt.0) write (156,*) ' numneg',numneg,'
c      1 sumneg',sumneg
      sume=0.
c
      do 450 i=1,ibar
      do 440 j=1,jbar
      do 430 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 430
      b=e(ijk)*rem0n0 ! [MeV*fm-3]
      cmom=xmom(ijk)**2+ymom(ijk)**2+zmom(ijk)**2 ! [MeV2*fm-6]
      c=cmom*remnsq ! [MeV2*fm-6]
      capn(ijk)=dxdydz*capn(ijk)
      e(ijk)=dxdydz*e(ijk) ! [MeV]
      sume=sume+e(ijk)
      icb=0
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      koost(ijk0)=icb
      xmom(ijk)=dxdydz*xmom(ijk) ! [MeV]
      ymom(ijk)=dxdydz*ymom(ijk) ! [MeV]
      zmom(ijk)=dxdydz*zmom(ijk) ! [MeV]
      sume=sume+e(ijk)
  430 continue
  440 continue
  450 continue
      return
c
  460 format (2i4,i9,3f10.5)
  470 format (3i4,1x,3f12.6,3x,3f12.6,f12.5)
  480 format (13x,4f12.6,/)
  490 format ((33f10.5,2x))
  500 format ('phase1 err e:k',3i4,' kb=',i4,2f13.8,2i8,3f13.8)
  510 format ('phase1:err capn1 ',3i4,' b:',2i4,3f13.8,3i8)
  520 format ('phase1:err ymom1',3i4,' ib',i4,3f14.9)
  530 format ('phase1:err xmom1',3i4,' kb',i4,3f14.9)
  540 format ('phase1:err zmomk',3i4,' kb',i4,3f14.9,2i8,2f13.8)
  550 format (2f5.0,f10.3,3f10.6)
  560 format (4i4,i10,f12.8,f12.5,f12.6,3f16.12,f7.1,f14.7,5f17.10)
  570 format ('error:phase1-xmom2 i=',2i4,' j,k=',i3,i4,2f12.8,' N',2f12
     1 .8,2i11)
  580 format ('error:phase1-ymom2 i=',2i4,' j,k=',i3,i4,2f12.8,' N',2f12
     1 .8,2i11)
  590 format ('error:phase1-zmom2 i=',2i4,' j,k=',i3,i4,2f15.8,' N',2f12
     1 .8,2i11,2e15.6)
  600 format ('error:phase1-  e   i=',2i4,' j,k=',i3,i4,2f12.8,' N',2f12
     1 .8,2i11,2f12.8)
  610 format ('error:phase1-  u i i=',2i4,' j,k=',i3,i4,2f12.8,' N',2f12
     1 .8,2i11,2f12.8)
  620 format ('error:phase1-  v i i=',2i4,' j,k=',i3,i4,2f12.8,' N',2f12
     1 .8,2i11)
  630 format ('error:phase1-  w i i=',2i4,' j,k=',i3,i4,2f12.8,' N',2f12
     1 .8,2i11,e15.7)
  640 format ('error:phase1-  p i i=',2i4,' j,k=',i3,i4,2f12.8,' N',2f12
     1 .8,2i11)
  650 format ('| phs1 ERROR E<0',i3,i9,' (',i3,i2,i4,')  E=',f10.5,' Eol
     1d=',f9.5,'  Edif=',f9.5,'  #',f4.1,' |')
  660 format ('koast-phs1a',4i4,i10,2f20.13,' >',f20.13)
  670 format (' imincl=',3i4,' jmaxcl=',2i3,' kmincl=',3i4)
  680 format ('phase1:err p:k  i=',2i3,' j=',i2,' k=',2i4,3f20.13,2i10)
  690 format ('phase1:err w:k  i=',2i3,' j=',i2,' k=',2i4,3e14.6,2i10,2e
     1 14.6,' zmom',2e14.6)
  700 format ('phase1:err u:k  i=',2i3,' j=',i2,' k=',2i4,3e14.6,2i10,'
     1xmom',4e14.6)
  710 format ('phase1:err P:k  i=',2i3,' j=',i2,' k=',2i4,3e14.6,2i10)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine outfin (timfmc,centrz,ktemp)
c
      implicit real*8 (a-h,o-z)
c
      parameter (BFLOOR_OUT=1.d-4)

      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
c
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
c
      common /mssss/ outn, porton, npxz, npt, kout
c
      common /numcmm/ numbk, numcel, numcst, numext, numfz1, numlv,
     1 numout, numouu, numou2, jdim1, kfreez, kfrzcy, press, freez,
     2 compmy, koast, kext, cpntot, compmx, compe,iundee(58)
c
      common /other/ eusual, epn0, t0, epion, ppion, edelt, delta, pdelt
     1 , temp, tempkn, eps, xth, epsc, dnsmax, ncycpt, ncycrp, ncycaf
c
      common /matric/ storev(929801), idensz(500), edensz(500)
c
      common /pscm2/ cycen, istop, isume
c
      common /pscmz/ zz
c
      common /foclose/ EPSFOE, EPSFOB, IAdaptiveFO, NInactiveFO,
     1 NInactiveReq, NMinFOCheck, ICycleExtendStep, ICycleMaxHard,
     2 IClosedFO
      
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      dimension entrpy(100), spherc(6), flux(6), evectr(9), flx(6),
     1 spher(6), fluxkp(100), fluzkp(100)
c
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233)
c
c      dimension  vv1(820521), vv2(820521), ptM(820521)
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
c
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c
      equivalence (nx(1),nxt), (nx(2),nxp), (ny(1),nyt), (ny(2),nyp),
     1 (nz(1),nzt), (nz(2),nzp), (nplot(1),nplott), (nplot(2),nplotp),
     2 (gam(1),gamt), (gam(2),gamp), (smaj(1),smajt), (smaj(2),smajp),
     3 (smin(1),smint), (smin(2),sminp), (fx(1),fxt), (fx(2),fxp), (fz(1
     4 ),fzt), (fz(2),fzp), (w0(1),wt), (w0(2),wp), (sumnu(1),sumnut),
     5 (sumnu(2),sumnup)
c
      data x, ueff /0.,0./
c
      if (t+1.e-7.gt.tout) go to 50
      if (ncyc.lt.2) go to 60
   10 if (IAdaptiveFO.ne.0.and.ncyc.ge.NMinFOCheck.and.IClosedFO.eq.1)
     1 isume=1
      if (ncyc.ge.icysum) then
      if (IAdaptiveFO.ne.0.and.IClosedFO.eq.0.and.ncyc.ge.NMinFOCheck
     1 .and.icysum.lt.ICycleMaxHard) then
      icysum=min0(icysum+ICycleExtendStep,ICycleMaxHard)
      write (*,*) 'Extending icysum to',icysum,' at cyc=',ncyc,
     1 ' (FO not closed yet)'
      else
      isume=1
      if (IAdaptiveFO.ne.0.and.IClosedFO.eq.0.and.
     1 icysum.ge.ICycleMaxHard) then
      write (*,*) 'Reached ICycleMaxHard=',ICycleMaxHard,
     1 ' without robust FO closure'
      endif
      endif
      endif
      if (isume.eq.1) go to 100
c      if (mod(ncyc,5).eq.0.and.ncyc.ge.50) go to 100
c
   20 if (istop.eq.1) go to 40
   30 t=t+dt
      ncyc=ncyc+1
      k10=3*(npt/10)
      k20=3*(npt/20)
      return
   40 continue
      numlv=-1
      numtd=numtd+1
      go to 160
c
   50 tout=tout+dto
   60 ix1=mod(ncyc,ncycpt)
      if (ncyc.lt.5) go to 70
      if (ix1.ne.0) go to 90
   70 continue
      if (ncyc.eq.1) go to 80
      i=-1
   80 continue
c
   90 continue
      iprt=2
      go to 10
c               ____________________________700: from after 510
  100 cycen=cycen+1.
      if (cycen.ge.cycenx) istop=1
      
c       write(*,*) 'ncyc=',ncyc,'cycen=',cycen,'cycenx=',cycenx
      
      if (ncyc.lt.icysum) istop=1
      do 110 i=1,6
      spherc(i)=0.
      flux(i)=0.
  110 continue
      cz1=0.
      cz2=0.
c       magas v1, v2
      ncount=0
      ncount1=0
      v1=0.
      v2=0.
c       /vladimir
      do 150 k=1,kbar
      do 140 j=1,jbar
      do 130 i=1,ibar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).le.BFLOOR_OUT) go to 130
      ijkz=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      cnth=capn(ijk)*dxdydz
      cz1=cz1+cnth
      cz2=cz2+capn(ijk)
      spherc(1)=spherc(1)+2*xmom(ijk)*xmom(ijk)
      spherc(2)=spherc(2)+2*xmom(ijk)*ymom(ijk)
      spherc(4)=spherc(4)+2*xmom(ijk)*zmom(ijk)
      spherc(3)=spherc(3)+2*ymom(ijk)*ymom(ijk)
      spherc(5)=spherc(5)+2*ymom(ijk)*zmom(ijk)
      spherc(6)=spherc(6)+2*zmom(ijk)*zmom(ijk)
      flux(1)=flux(1)+2*xmom(ijk)*xmom(ijk)/cnth
      flux(2)=flux(2)+2*xmom(ijk)*ymom(ijk)/cnth
      flux(4)=flux(4)+2*xmom(ijk)*zmom(ijk)/cnth
      flux(3)=flux(3)+2*ymom(ijk)*ymom(ijk)/cnth
      flux(5)=flux(5)+2*ymom(ijk)*zmom(ijk)/cnth
      flux(6)=flux(6)+2*zmom(ijk)*zmom(ijk)/cnth
      uvwijk=u(ijk)*u(ijk)+v(ijk)*v(ijk)+w(ijk)*w(ijk)
c       if (uvwijk.ge.1.) write (156,390) u(ijk),v(ijk),w(ijk),i,j,k
      if (isnan(uvwijk).or.uvwijk.ge.1.d0) then
      write (*,*) 'v^2>=1 corrected (2156) cyc/i/j/k=',ncyc,i,j,k,
     1 ' vel2=',uvwijk,' u,v,w=',u(ijk),v(ijk),w(ijk)
      if (isnan(uvwijk).or.uvwijk.le.0.d0) then
      u(ijk)=0.d0
      v(ijk)=0.d0
      w(ijk)=0.d0
      else
      vfac=dsqrt(0.9999d0/uvwijk)
      u(ijk)=u(ijk)*vfac
      v(ijk)=v(ijk)*vfac
      w(ijk)=w(ijk)*vfac
      endif
      uvwijk=u(ijk)*u(ijk)+v(ijk)*v(ijk)+w(ijk)*w(ijk)
      endif
c       magas v1, v2
      ptMijk=dsqrt(xmom(ijk)**2+ymom(ijk)**2)
      if (ptMijk.eq.0.) go to 130
      vv1ijk=xmom(ijk)/ptMijk
      vv2ijk=(xmom(ijk)/ptMijk)**2-(ymom(ijk)/ptMijk)**2
      if (j.ne.1) koffff=2
      if (j.eq.1) koffff=1
c       v1 is averaging for only one side of OX
      if (xmom(ijk).lt.0) go to 120
      v1=v1+koffff*vv1ijk
      ncount1=ncount1+koffff
  120 v2=v2+koffff*vv2ijk
      ncount=ncount+koffff
  130 continue
  140 continue
  150 continue
c       Magas v1, v2
      if (ncount1.ne.0) v1=-v1/ncount1
      if (ncount.ne.0) v2=v2/ncount
      if (dabs(1.-cycenx).lt.0.1) cycen=1.
      if (ncyc.ge.9999) go to 160
      go to 20
c   _______________________________________________________________
c  !________________________________________ final output _________!
c
  160 continue
c       write (156,380)
      write (*,380)
      compe=compe*.145*939.5
      kpzz=2
      if (kpzz.eq.2) return
c       write (156,360) compe
c  stuff removed
      k=0
c       write (156,400)
      do 180 i=1,3
      do 170 j=1,i
      k=k+1
      spherc(k)=spherc(k)+spher(k)
      flux(k)=flux(k)+flx(k)
c       write (156,450) i,j,spherc(k),flux(k),spher(k),flx(k)
  170 continue
  180 continue
c  set y comps. of y to zero from symmetry
c
      flux(2)=0.
      flux(5)=0.
      flx(2)=0.
      flx(5)=0.
      spherc(2)=0.
      spherc(5)=0.
      spher(2)=0.
      if (jdim1.eq.1) return
      spher(5)=0.
      call eigen (flux,evectr,3,0)
      xrat=flux(1)/flux(6)
      yrat=flux(3)/flux(6)
      anglr=dasin(evectr(9))
      anglr=180.*anglr/pi
c       write (156,410)
c       write (156,460) flux(6),flux(3),flux(1),xrat,yrat,anglr
c       write (156,430)
      do 190 i=1,3
      n1=3*(i-1)+1
      n2=3*i
c       write (156,480) i,(evectr(j),j=n1,n2)
  190 continue
      call eigen (spherc,evectr,3,0)
c       write (156,470) spherc(1),spherc(3),spherc(6)
      do 200 i=1,3
      n1=3*(i-1)+1
      n2=3*i
c       write (156,480) i,(evectr(j),j=n1,n2)
  200 continue
c       write (156,430)
c
c       write (156,420)
c
      if (ncyc.lt.icysum) go to 30
      ymax=0.
      do 210 i=1,ncyc
      if (entrpy(i).gt.ymax) ymax=entrpy(i)
  210 continue
      k=ymax
      ymax=k+1
      do 220 i=1,100
      edensz(i)=float(i)
  220 continue
      k=ncyc+1
      if (ncyc.ge.100) go to 240
      do 230 i=k,100
      entrpy(i)=0.
  230 continue
  240 continue
      ncyc=min0(ncyc,100)
      numlv=-1
      rewind 1
c       write (156,440)
      sumde2=0.
      if (ncyc.gt.100) go to 340
      do 270 i=1,ibar
      do 260 j=1,jbar
      do 250 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 250
      cnth=capn(ijk)*rdxdydz
  250 continue
  260 continue
  270 continue
      npt=3*npt
      ymax=0.
      do 280 i=1,ncyc
      if (ymax.lt.fluxkp(i)) ymax=fluxkp(i)
  280 continue
      k=ncyc+1
      do 290 i=k,100
      fluxkp(i)=0.
  290 continue
      if (ymax.eq.0.) go to 300
  300 ymax=0.
      do 310 i=1,ncyc
      if (ymax.lt.fluzkp(i)) ymax=fluzkp(i)
  310 continue
      if (ymax.eq.0.) go to 330
      k=ncyc+1
      do 320 i=k,100
      fluzkp(i)=0.
  320 continue
  330 continue
c         write (156,410)
c       write (156,460) flux(6),flux(3),flux(1),xrat,yrat,anglr
c       write (156,380)
c       write (156,370)
c       write (156,380)
      numlv=-1
  340 continue
c   ________________________________________________________________
c  !___________________________________ End final output ___________!
c
      return
c
  350 format (2x,E12.6,5x,E12.6,5x,i3)
  360 format (//,' maximum temperature ',f8.2,'  energy ',f8.2/)
  370 format ('    End of Output on Unit ***  6  *** ')
  380 format (/,45(2h__),/)
  390 format (' u v w ',3f12.5,3i4)
  400 format (1h1,/,'    spherc flux',/)
  410 format (//,13x,' flux components',20x,'1/3    2/3    theta',/)
  420 format (1h1)
  430 format (/)
  440 format (/,'1(cyc   F-xx    F-yy    F-zz    F-1+    F-2+[mN
     1]  F-3+    ',/)
  450 format (' sphericity flux  ',2i2,4f13.5)
  460 format (' flux',3f13.5,3x,2f7.2,f10.2)
  470 format (' spherc',3f13.5)
  480 format (i3,1x,3f10.5)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine opnfil
      implicit real*8 (a-h,o-z)
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      open (unit=156,file='fort.6',status='unknown')
      open (unit=118,file='fort_118.dat',status='unknown')
      open (unit=333,file='fort.333',status='unknown')
      open (unit=221,file='oscfil.dat',status='unknown')
      open (unit=334,file='fort.334',status='unknown')
      open (unit=336,file='fort.336',status='unknown')
      open (unit=337,file='fort.337',status='unknown')
      open (unit=338,file='fort.338',status='unknown')
      open (unit=339,file='fort.339',status='unknown')
      open (unit=889,file='fort.889',status='unknown')
      return
   30 continue
      stop 1030
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function mylenc (string)
      character temp*1
      character string*(*)
      character kstrng
      common /myl/ kstrng(100)
      j=len(string)
      do 10 i=j,1,-1
      temp=char(0)
      if ((string(i:i).ne.' ').and.(string(i:i).ne.temp)) go to 20
   10 continue
      i=0
   20 mylenc=i
      do 30 k=1,i
      kstrng(k)=string(k:k)
   30 continue
      return
c
      end
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine markr1
c
      implicit real*8 (a-h,o-z)
c
c   Generate the marker particles using npartp, npatz algorithm
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
c
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233)
c
      common /xposn/ itrck(30000000)
c
      common /matric/ zmin(300,300), zmax(300,300), energ(300,300), velz
     1 (300,300), dens(300,300), zsmz(200), adens(300,700), kmove(700),
     2 Rpl(300,300), Rmi(300,300), chyz(300,300), dummy(1)
c
      common /thestf/ rap00, Tfin, dxv, dyv, dzv, epsilon00, b0, rho00,
     1 e00
c
      common /mssss/ outn, porton, npxz, npt, kout
c
      common /midle/ imid, kmid, kend, imidpt, kmidpt, klmin, krmax,
     1 kdisp, ipmin, ipmax, imid1
c
      dimension px(85877055), py(85877055), pz(85877055), zminp(300,300)
     1 , zmaxp(300,300)
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c
c   Now distrbute the marker particles   _______________________________
c   Define how much each marker particle will carry
c
      dpt=dxdydz/dfloat(nparp*nparp*nparpz)
      npxz=0
      esum=0.
      asum=0.
      asumot=0.
      esumot=0.
      outn=0.
      ko=-2
      mzmax=0
      jom=0
      iom=0
      kout=0
      ymin=1000.
      ymax=-1.
      parpi=1./dfloat(nparp)
      parpzi=1./dfloat(nparpz)
      ypmin=1./dfloat(nparp)
      ypmin=0.35
      if (nparp.eq.1) parpi=0.
      if (nparpz.eq.1) parpzi=0.
      if (nparp.gt.8) write (*,*) ' nparp>8'
      if (nparpz.gt.8) write (*,*) ' nparpz>8'
c      if (nparp.gt.8) stop 1236
c      if (nparpz.gt.8) stop 1238
c      if (nparp.eq.6) stop 1236
      ess=0.
      us=0.
      ps=0.
      do 30 i=1,ibar
      do 20 j=1,jbar
      zminp(i,j)=10000.
      zmaxp(i,j)=0.
      do 10 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (u(ijk).eq.0.) go to 10
      ess=ess+e(ijk)
      us=us+u(ijk)
      ps=ps+p(ijk)
c
   10 continue
   20 continue
   30 continue
c       if (kdebug.ne.0) write (156,*) 'ess',ess,'  us',us,'  ps',ps,ess/
c      1 ps
c
c    Assign marker particles to cells __________________________________
c
      ibar2=(ibar)/2
      kbar2=(kbar+3)/2
      ibarx=ibar
      if (b0.eq.0.) ibarx=ibar2
      bsum=0.
      besum=0.
      rrmax=0.
      numpl=0
      nummn=0
c
c  Note that we take advantage of symmetry for b=0 reactions
c
      do 380 i=1,ibarx
      do 370 j=1,jbar
      do 360 k=1,kbar2
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (ijk.gt.ijkmax) go to 460
      kb=kbarmx+1-k
      ijkb=nq*(ijbar*(kb-1)+ibar*(j-1)+i-1)+1
      if (u(ijk).eq.0.) go to 350
      bsum=bsum+p(ijk)
      besum=besum+e(ijk)
      if (krandm.eq.1.or.krandm.eq.0) go to 40
      call random_seed
   40 continue
c
c   mz = # of marker particles in this cell
c
      temp=p(ijk)/(dpt)+.5
c      temp=capn(ijk)/(dpt)+.5
      mz=temp
      if (mz.ne.0) go to 50
      asumot=asumot+p(ijk)
      esumot=esumot+e(ijk)
      e(ijk)=0.
      p(ijk)=0.
      xmom(ijk)=0.
      ymom(ijk)=0.
      zmom(ijk)=0.
      capn(ijk)=0.
      go to 340
   50 continue
      if (mz.gt.mzmax) mzmax=mz
      lji=0
      laccpt=0
      laccp2=0
   60 lji=lji+1
   
c	NOW we have many marker particles 29-09-2023
c       if (lji.le.1000) go to 70
c      outn=outn+p(ijk)
c      go to 350
c   70 continue
c   Distribute marker particles


	if (krandm.eq.0) go to 80
c
      call random_number (x)
      call random_number (y)
      call random_number (z)
      if (krandm.eq.1) go to 220
      if (krandm.eq.2) go to 220
      xrd=x
      yrd=y
      zrd=z
      if (krandm.eq.-1) go to 80
      f=1./dfloat(krandm-1)
      x=x*f
      y=y*f
      z=z*f
      go to 220
   80 continue
c
c   no randomization
      x=.5
      y=.5
      z=.5
      if (nparp.eq.1) go to 90
      if (nparp.eq.2) go to 100
      if (nparp.eq.3) go to 110
      if (nparp.eq.4) go to 120
      if (nparp.eq.5) go to 130
      if (nparp.eq.7) go to 140
      write (*,*) ' nparp problem',nparp
      stop 123
   90 xdf=0
      ydf=0.
      go to 150
  100 xdf=parpi*mod((lji-1)/2,2)-.25
      ydf=parpi*mod((lji-1)/(nparp*nparpz),2)-.25
      go to 150
  110 xdf=parpi*(mod((lji-1)/nparpz,3)-1)
      ydf=parpi*(mod((lji-1)/(nparp*nparpz),3)-1)
      go to 150
  120 ydf=.25*mod((lji-1)/(nparp*nparpz),4)-.375
      xdf=.25*mod((lji-1)/nparpz,4)-.375
      go to 150
  130 ydf=.2*mod((lji-1)/(nparp*nparpz),5)-.4
      xdf=.2*mod((lji-1)/nparpz,5)-.4
      go to 150
  140 ydf=(mod((lji-1)/(nparp*nparpz),7)-3)/7.
      xdf=(mod((lji-1)/nparpz,7)-3)/7.
c
  150 if (i.lt.imid) xdf=-xdf
      if (nparpz.eq.1) go to 160
      if (nparpz.eq.2) go to 170
      if (nparpz.eq.3) go to 180
      if (nparpz.eq.4) go to 190
      if (nparpz.eq.5) go to 200
      write (*,*) ' Need coding here'
      stop 9871
  160 zdf=0.
      go to 210
  170 zdf=0.5*mod(lji,2)-.25
      go to 210
  180 zdf=.3333333333*(mod(lji,3)-1)
      go to 210
  190 zdf=.25*mod((lji-1),4)-.375
      go to 210
  200 zdf=.2*mod((lji-1),5)-.4
  210 x=xdf+x
      y=ydf+y
      z=zdf+z
  220 continue
      xc=i-1.+x
      yc=j-1.+y
      zc=k-1.+z
      thalf=0.
c???
      if (krandm.ge.0) go to 270
      xcmidl=.5
      xcmidr=.5
c???
      if (xc.gt.dfloat(imid)) go to 230
      plrad=yc*yc+(xc-xcmidr)**2
      xq=yc*yc+(xc-xcmidl)**2
      if (plrad.ge.thalf) go to 270
      go to 240
  230 plrad=yc*yc+(xc-xcmidl)**2
      if (plrad.ge.thalf) go to 270
  240 xc=xc-x+xrd
      yc=yc-y+yrd
      zc=zc-z+zrd
      plrad=yc*yc+(xc-xcmidr)**2
      if (plrad.lt.thalf) go to 250
      xc=xc-xrd+x
      yc=yc-yrd+y
      zc=zc-zrd+z
      go to 260
  250 x=xrd
      y=yrd
      z=zrd
  260 continue
c
c   To which of the original cells does this correspond?
c
c   CHANGE  8 DEC 2010
  270 continue
      jo=yc*dfloat(nreduy)+1
      if (jo.gt.jom) jom=jo
      io=dfloat(nredux)*(xc-imid-.1)+dfloat(imid1)+.5
      kentr=0
      if (i.gt.imid1.and.kentr.eq.0) io=io+1
      if (io.gt.iom) iom=io
      zko=(zc-dfloat(kmid)-1.)*dzvfm
      if (io.le.0) stop 'io < 0'
      if (io.le.0) go to 60
c
      if (ymin.gt.yc) ymin=yc
      if (ymax.lt.yc) ymax=yc
c
      if (zko.gt.zmax(io,jo)) kout=kout+1
      if (zko.lt.zmin(io,jo)) kout=kout+1
      if (zko.gt.zmax(io,jo)) go to 60
      if (zko.lt.zmin(io,jo)) go to 60
c   CHANGE  8 DEC 2010
      if (j.eq.1) npxz=npxz+1
      kkk=zc+1.
      dist=yc*yc+(xc-imid)**2
      rdist=dsqrt(dist)
c       if (i.eq.imid.and.j.eq.0) write (156,520) i,j,k,xc,zc,rdist,skt
      if (kslab.eq.4.and.rdist.gt.tnc) go to 60
      sdist=dsqrt(dist)
      skt=ktnc/2
      if (b0.lt.0.01) go to 280
      if (b0.lt.0.01 .and. sdist .gt. skt) go to 60
      if (b0.lt.0.01.and.sdist.gt.13..and.ktnc.eq.16) go to 60
      if (b0.lt.0.01.and.sdist.gt.7..and.ktnc.eq.12) go to 60
      if (kslab.eq.4.and.u(ijk).eq.0) go to 60
  280 continue
c   Calculate polar radius
      if (krandm.ge.0) go to 300
c   __
      if (tnc*tnc-yc*yc.lt.0.) go to 290
      x1=xcmidr-dsqrt(tnc*tnc-yc*yc)
      plrad=dsqrt(yc*yc+(xc-xcmidr)**2)
c
      if (k.eq.kmid.and.plrad.gt.tnc) komit=komit+1
c       if (k.eq.kmid.and.plrad.gt.tnc) write (156,*) ' omitz',komit,xc,yc
c      1 ,zc,plrad,k
      if (plrad.gt.tnc) go to 330
  290 continue
  300 continue
      ko=ko+3
      if (e(ijk).eq.0.) write (*,*) 'E=0',i,j,k,ijk,' N=',capn(ijk)
      if (e(ijk).eq.0.) go to 330
      capn(ijk)=capn(ijk)+dpt
      asum=asum+dpt
      if (ko.gt.kpmax) go to 450
      laccpt=laccpt+1
      px(ko)=xc
      py(ko)=yc
      pz(ko)=zc
      kopxt=ko/3+1
      if (zc.lt.zminp(i,j)) zminp(i,j)=zc
      if (zc.gt.zmaxp(i,j)) zmaxp(i,j)=zc
      xcd=xc-dfloat(imid)
      if (xcd.lt.0.) itrck(kopxt)=1
      if (xcd.ge.0.) itrck(kopxt)=2
      if (xcd.ge.0.) numpl=numpl+1
      if (xcd.lt.0.) nummn=nummn+1
c       if (zko.ge.zmax(io,jo)) write (156,*) 'error ',xc,yc,zc
c       if (zko.le.zmin(io,jo)) write (156,*) 'error ',xc,yc,zc
c       if (py(ko).le.0.) write (156,*) ' error in y',y
      if (zko.gt.zmax(io,jo)) go to 60
      if (zko.lt.zmin(io,jo)) go to 60
      ib=ibar+1-i
      kb=kbar+3-k
      ijkb=nq*(ijbar*(kb-1)+ibar*(j-1)+ib-1)+1
      ko=ko+3
      if (e(ijkb).eq.0.) write (*,480) ib,j,kb,ijkb,capn(ijkb),p(ijkb),i
     1 ,j,k,ijk,capn(ijk),e(ijk),p(ijk)
      capn(ijkb)=capn(ijkb)+dpt
      asum=asum+dpt
      if (ko.gt.kpmax) go to 450
c
      xc=ibar+1-i-x
      yc=j-1.+y
      zc=kbar+3-k-z
      px(ko)=xc
      py(ko)=yc
      pz(ko)=zc
      if (zc.lt.zminp(ib,j)) zminp(ib,j)=zc
      if (zc.gt.zmaxp(ib,j)) zmaxp(ib,j)=zc
      kopxt=ko/3+1
      xcd=xc-dfloat(imid)
      if (xcd.lt.0.) itrck(kopxt)=1
      if (xcd.ge.0.) itrck(kopxt)=2
      if (xcd.ge.0.) numpl=numpl+1
      if (xcd.lt.0.) nummn=nummn+1
      kkk=zc+1.
      dist=yc*yc+(xc-imid)**2
      sdist=dsqrt(dist)
      if (b0.ne.0.) go to 330
c
c   Use symmetry  - reflect about imid line
      ijkb=nq*(ijbar*(k-1)+ibar*(j-1)+ib-1)+1
      if (e(ijkb).ne.0.) go to 310
      write (*,500) i,ib,j,k,ijk,ijkb,capn(ijkb),capn(ijk),e(ijk)
      e(ijkb)=0.5*e(ijk)
      e(ijk)=0.5*e(ijk)
      write (*,*) i,j,k,ijk,' ref'
  310 ko=ko+3
      capn(ijkb)=capn(ijkb)+dpt
      if (i.ge.42) write (*,*) 'm8',i,j,k,ijkb,capn(2121703),p(2121703)
     1 ,capn(2160115),p(2160115),capn(2160355)
      asum=asum+dpt
      if (i.eq.35.and.k.eq.1) write (*,490) i,ib,j,k,kb,ijkb,capn(ijk),e
     1 (ijk),p(ijk)
      if (ko.gt.kpmax) go to 450
c
      xc=ibar+1-i-x
      yc=j-1.+y
      zc=k-1.+z
      px(ko)=xc
      py(ko)=yc
      pz(ko)=zc
      if (zc.lt.zminp(ib,j)) zminp(ib,j)=zc
      if (zc.gt.zmaxp(ib,j)) zmaxp(ib,j)=zc
      kopxt=ko/3+1
      xcd=xc-dfloat(imid)
      if (xcd.lt.0.) itrck(kopxt)=1
      if (xcd.ge.0.) itrck(kopxt)=2
      if (xcd.ge.0.) numpl=numpl+1
      if (xcd.lt.0.) nummn=nummn+1
c
c   Use symmetry  - reflect about kmid line
c
      ijkb=nq*(ijbar*(kb-1)+ibar*(j-1)+i-1)+1
      if (e(ijkb).ne.0.) go to 320
      write (*,510) i,j,k,kb,ijk,ijkb,capn(ijkb),capn(ijk),e(ijk)
      e(ijkb)=0.5*e(ijk)
      e(ijk)=0.5*e(ijk)
  320 ko=ko+3
      capn(ijkb)=capn(ijkb)+dpt
      asum=asum+dpt
      if (ko.gt.kpmax) go to 450
c
      xc=i-1.+x
      yc=j-1.+y
      zc=kbar+3-k-z
      px(ko)=xc
      py(ko)=yc
      pz(ko)=zc
      if (zc.lt.zminp(i,j)) zminp(i,j)=zc
      if (zc.gt.zmaxp(i,j)) zmaxp(i,j)=zc
      kopxt=ko/3+1
      xcd=xc-dfloat(imid)
      if (xcd.lt.0.) itrck(kopxt)=1
      if (xcd.ge.0.) itrck(kopxt)=2
      if (xcd.ge.0.) numpl=numpl+1
      if (xcd.lt.0.) nummn=nummn+1
c
      kkk=zc+1.
      dist=yc*yc+(xc-imid)**2
      sdist=dsqrt(dist)
      if (rrmax.lt.sqrt(dist)) rrmax=dist
  330 if (laccpt.lt.mz) go to 60
  340 continue
  350 continue
      kb=kbarmx+1-k
      ijkb=nq*(ijbar*(kb-1)+ibar*(j-1)+i-1)+1
  360 continue
  370 continue
  380 continue
      npt=ko
      do 420 i=1,ibar
      do 410 j=1,jbar
      do 400 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).ne.0.) go to 390
      u(ijk)=0.
      v(ijk)=0.
      e(ijk)=0.
      zmom(ijk)=0.
      go to 400
  390 continue
c      if(kreat1.ne.0)write(458,103) i,j,k,ijk,capn(ijk),e(ijk),
c     1 zmom(ijk),w(ijk),p(ijk),xmom(ijk),ymom(ijk),u(ijk),v(ijk)
  400 continue
  410 continue
  420 continue
      i=-1
      return
c
  450 continue
c       write (156,*) ' too many particles',i,j,k,ko,'>',kpmax
      write (*,*) ' too many particles',i,j,k,ko,'>',kpmax
      stop 41
  460 continue
c       write (156,*) ' number of cells too large ',ijk
      write (*,*) ' number of cells too large ',ijk
      stop 42
c
  470 format (2i4,2f11.4)
  480 format ('markr: E1=0:',3i4,i8,2f12.8,3i4,i8,'  N=',f12.8,' E=',3f1
     1 2.8)
  490 format (5i4,i8,3f12.8)
  500 format ('markr: E2=0',i4,i3,' j=',i2,i4,2i11,' N=',2f11.6,' E=',2f
     1 11.6)
  510 format ('markr: E3=0',i4,' j=',i2,2i4,2i11,' N=',2f11.6,' E=',f11.
     1 6)
  520 format (3i3,2f7.2,' R=',f5.2,f5.2,f5.0)
  530 format (/,'In the setup, the following cells were omitted because'
     1 ,/,' the amount of matter was smaller than one-half of',f8.5,/)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine koccupy (iq,jq,kq,io,io6,jdim1)
c
      implicit real*8 (a-h,o-z)
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfz
c
      common /plcm2/ aa2(91627237)
c
      common /occp/ numcel, numce2
c
      dimension capn(91627233)
c
      equivalence (capn,aa2(1))
c
      io=0
      io6=0
      koun=0
      isp=nq*(ijbar*(kq-1)+ibar*(jq-1)+iq-1)+1
      i=iq
      j=jq
      k=kq
      nqk=nq*ijbar
      numcel=0
      numce2=0
      if (j.gt.1) go to 90
c
c      first do j=1
c      we don't have to do the ijm check since that cell is the
c      the reflected cell: isp.
c
      isp=nq*(i-1)+1
      j=1
      isp=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      imj=isp-nq
      ipj=isp+nq
      ijp=isp+nqib
      ijku=isp-nqk
      ijkd=isp+nqk
      if (capn(isp).eq.0.) go to 80
      if (iq.le.1) go to 10
      if (capn(imj).ne.0.) koun=koun+1
   10 if (iq.ge.ibar) go to 20
      if (capn(ipj).ne.0.) koun=koun+1
   20 if (jdim1.ne.1.and.capn(ijp).ne.0.) koun=koun+1
      if (kq.le.1) go to 30
      if (capn(ijku).ne.0.) koun=koun+1
   30 continue
      if (kq.ge.kbar) go to 40
      if (capn(ijkd).ne.0.) koun=koun+1
   40 continue
      numcel=numcel+1
      if (koun.gt.1) io6=1
      if (koun.gt.2) io=1
      if (kq.le.1) go to 50
      if (kq.ge.kbar) go to 50
      if (iq.le.1) go to 50
      if (iq.ge.ibar) go to 50
      return
   50 continue
      if (koun.ge.0) io6=1
      if (koun.ge.2) io=1
      return
c
c    check lower level
c
c      i+1,j,k-1
      i1=ijku+nq
      if (capn(i1).eq.0.) go to 70
c      i-1,j,k-1
      i1=ijku-nq
      if (capn(i1).eq.0.) go to 70
c      i,j+1,k-1
      i1=ijku+nqib
      if (jdim1.gt.1.and.capn(i1).eq.0.) go to 70
c
c    check upper level
c
c      i+1,j,k+1
      i1=ijkd+nq
      if (capn(i1).eq.0.) go to 70
c      i-1,j,k+1
      i1=ijkd-nq
      if (capn(i1).eq.0.) go to 70
c      i,j+1,k+1
      if (jdim1.eq.1) go to 60
      i1=ijkd+nqib
      if (capn(i1).eq.0.) go to 70
c
c    check centre level
c
c      i-1,j+1,k
      i1=ijp-nq
      if (capn(i1).eq.0.) go to 70
c      i+1,j+1,k
      i1=ijp+nq
      if (capn(i1).eq.0.) go to 70
c      i-1,j+1,k-1
      i1=ijkd-nq+nqib
c    if (capn(i1).eq.0.) go to 40
c      i-1,j+1,k+1
      i1=ijku-nq+nqib
c    if (capn(i1).eq.0.) go to 40
c
c    i+1,j+1,k-1
      i1=ijkd+nq+nqib
c    if (capn(i1).eq.0.) go to 40
c    i+1,j+1,k+1
      i1=ijku+nq+nqib
c    if (capn(i1).eq.0.) go to 40
c
c      an interior cell
c
   60 continue
      io=1
   70 continue
   80 continue
      return
   90 continue
      koun=0
c
c  _______________________________________________________________
c
      if (jdim1.eq.1) return
      isp=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      imj=isp-nq
      ipj=isp+nq
      ijm=isp-nqib
      ijp=isp+nqib
      ijku=isp-nqk
      ijkd=isp+nqk
      if (capn(isp).eq.0.) return
      if (iq.le.1) go to 100
      if (capn(imj).ne.0.) koun=koun+1
  100 if (iq.ge.ibar) go to 110
      if (capn(ipj).ne.0.) koun=koun+1
  110 if (capn(ijm).ne.0.) koun=koun+1
      if (capn(ijp).ne.0.) koun=koun+1
      if (kq.le.1) go to 120
      if (capn(ijku).ne.0.) koun=koun+1
  120 continue
      if (kq.ge.kbar) go to 130
      if (capn(ijkd).ne.0.) koun=koun+1
  130 continue
      if (kq.le.1) go to 140
      if (kq.ge.kbar) go to 140
      if (iq.le.1) go to 140
      if (iq.ge.ibar) go to 140
      if (koun.gt.1) io6=1
      if (koun.gt.3) io=1
      return
  140 continue
      if (koun.ge.1) io6=1
      if (koun.ge.3) io=1
      return
c
c    check lower level
c
c      i+1,j,k-1
      i1=ijku+nq
      if (capn(i1).eq.0.) go to 160
c      i-1,j,k-1
      i1=ijku-nq
      if (capn(i1).eq.0.) go to 160
c      i,j+1,k-1
      i1=ijku+nqib
      if (capn(i1).eq.0.) go to 160
c      i,j-1,k-1
      i1=ijku-nqib
      if (capn(i1).eq.0.) go to 160
c
c    check upper level
c
c      i+1,j,k+1
      i1=ijkd+nq
      if (capn(i1).eq.0.) go to 160
c      i-1,j,k+1
      i1=ijkd-nq
      if (capn(i1).eq.0.) go to 160
c      i,j+1,k+1
      i1=ijkd+nqib
      if (capn(i1).eq.0.) go to 160
c      i,j-1,k+1
      i1=ijkd-nqib
      if (capn(i1).eq.0.) go to 160
c
c    check centre level
c
c    i-1,j-1,k
      i1=ijm-nq
      if (capn(i1).eq.0.) go to 160
c      i+1,j-1
      i1=ijm+nq
      if (capn(i1).eq.0.) go to 160
c      i-1,j+1,k
      i1=ijp-nq
      if (capn(i1).eq.0.) go to 160
c      i+1,j+1,k
      i1=ijp+nq
      if (capn(i1).eq.0.) go to 160
      kchose=0
c    ???
      if (kchose.eq.0) go to 150
c
c    next few lines examine the corner cells
c
c    upper level
c
c    i+1,j+1,k+1
      i1=ijku+nq+nqib
      if (capn(i1).eq.0.) go to 160
c      i+1,j-1,k+1
      i1=ijku+nq-nqib
      if (capn(i1).eq.0.) go to 160
c      i-1,j-1,k+1
      i1=ijku-nq-nqib
      if (capn(i1).eq.0.) go to 160
c      i-1,j+1,k+1
      i1=ijku-nq+nqib
      if (capn(i1).eq.0.) go to 160
c
c    lower level
c
c    i+1,j+1,k-1
      i1=ijkd+nq+nqib
      if (capn(i1).eq.0.) go to 160
c      i+1,j-1,k-1
      i1=ijkd+nq-nqib
      if (capn(i1).eq.0.) go to 160
c      i-1,j-1,k-1
      i1=ijkd-nq-nqib
      if (capn(i1).eq.0.) go to 160
c      i-1,j+1,k-1
      i1=ijkd-nq+nqib
      if (capn(i1).eq.0.) go to 160
  150 continue
      io=1
c
c      an interior cell
c
  160 continue
      return
c       write (156,*) 'i,j, or k at boundaries',iq,jq,kq,' increase CELL
c      1 number limit in the code or',' increase kbarex in the input
c      2    *** !!!!!!!!!! ***'
      write (*,*) 'i,j, or k at boundaries',iq,jq,kq
      stop 1234
c
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine interx (i,j,k,ijk,iphs3,adens,ib4,numlvv,uijk,vijk,wijk
     1 ,pijk)
c
      implicit real*8 (a-h,o-z)
      integer*2 koost
c      real*4 aa1, u, v, w
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
c
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
c
      common /pscmz/ zz
c
      common /yestemp/ jtemp
c
      common /juet/ dxyz1, pwq, amnW, beta2, gammaf, betat, cosfi, gjut,
     1 y0, hj, sumths
c
      common /prstbl/ prs(99,102), trs(99,102), pres, tempzz, chem,
     1 amurs(99,102), entrs(99,102), entro, ipres
c
      common /other/ eusual, epn0, t0, epion, ppion, edelt, delta, pdelt
     1 , temp, tempkn, eps, xth, epsc, dnsmax, ncycpt, ncycrp, ncycaf
c
      common /numcmm/ numbk, numcel, numcst, numext, numfz1, numlv,
     1 numout, numouu, numou2, jdim1, kfreez, kfrzcy, press, freez,
     2 compmy, koast, kext, cpntot, compmx, compe,iundee(58)
c
      common /newcrp/ entot, dentot, tempav, numcll
c
      common /rot/ a, b, c
c
      common /densav/ densav
c
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c
      common /hxi33/ hxi3a, xi, hxi
c
      common /plasma/ sumplast, iplasma, kplasma
c
      common /eqqcd/ tq, eq1, pq, sq, anq1, amuq
c
      common /plcm7/ kextrn(30542411), koost(30542411)
c
      common /cnuclr/ hvc, hv2, pi0, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
c     
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233), px(85877055), py(85877055), pz(85877055), numlvv(50),
     3 adens(300,700)
c
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
c
      external fct1, fun, fct2
c
c      entrop(x)=0.5213+1.5*dlog(x+dexp(-2.*0.5213/3.))+1.809*dsqrt(x)/(
c     1. +1.139*dsqrt(x)+1.417*x+1.014*x*dsqrt(x))
      entrop(x)=0.5213+1.5*dlog(dble(x)+dexp(dble(-2.*0.5213/3.)))+1.809
     1 *dsqrt(dble(x))/(1.+1.139*dsqrt(dble(x))+1.417*x+1.014*x*dsqrt
     2 (dble(x)))
c
c      tempfn(x)=1.5*x/(x+dexp(-2.*.5213/3.))+0.9045*dsqrt(x)*(1.-1.417*
c     1-2.028*x*dsqrt(x))/(1.+1.139*dsqrt(x)+1.417*x+1.014*x*dsqrt(x))**
c
      tempfn(x)=1.5*x/(x+dexp(dble(-2.*.5213/3.)))+0.9045*dsqrt(dble(x))
     1 *(1.-1.417*x-2.028*x*dsqrt(dble(x)))/(1.+1.139*dsqrt(dble(x))+1.4
     2 17*x+1.014*x*dsqrt(dble(x)))**2
c
      kext=0
      pijk=0.
      a=capn(ijk)*ren0
      b=e(ijk)*rem0n0
      cmom=xmom(ijk)**2+ymom(ijk)**2+zmom(ijk)**2
      c=cmom*remnsq
      ra=1./a
      bsq=b*b
      if (c.le.bsq) go to 10
c
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (c.gt.bsq) koost(ijk0)=1
      call intern (i,j,k,ijk,2,adens,ib4,numlvv,uijk,vijk,wijk,pijk
     1 ,dfloat(1))
      return
   10 continue
c
c      if here, then i,j,k is an interior cell
c
      xi=0.
      iconv=0
      numit=0
   20 continue
      xmaxa=0.99999999*a
      deps=0.000001
      call rtmi (xx,froot,fun,dfloat(0),xmaxa,deps,60,ier)
c       write (156,*) ' interior ',i,j,k,xx
c       if (ier.ne.0) write (156,*) ' ier not zero in INTERX ',a,b,c,froot
      xv=1.-xx*xx/(a*a)
      yy=b-dsqrt(c*xv)
      xi=xx
      y=yy
      x=xv
      if (keqost.le.0) go to 80
      hxi=eqost(xi)
      xioasq=(xi*ra)**2
      x=1.-xioasq
      yden=1.+bcoef-bcoef*xioasq
      if (iconv.eq.1) go to 70
      ww=c*yden**2
      z=((1.+bcoef)*b+hxi)**2
      woz=ww/z
      if (woz.gt.1.) go to 30
      xinew=a*dsqrt(1.-woz)
      delxi=abs(xinew-xi)
      xi=xinew
      numit=numit+1
      if (delxi.lt.1.e-5) go to 60
      if (numit.lt.500) go to 20
   30 xi=0.
      numit=0
   40 xioasq=(xi*ra)**2
      x=1.-xioasq
      if (x.le.0.) go to 50
      xiterm=eqost(xi)
      yden=5.-2.*xioasq
      ww=dsqrt(c)*yden
      xinew=0.5*((1.+bcoef)*b-ww/dsqrt(x)+bcoef*xi+xiterm)
      delxi=abs(xinew-xi)
      xi=xinew
      numit=numit+1
      if (delxi.lt.1.d-5) go to 60
      if (xi.lt.0.) go to 50
      if (numit.lt.500) go to 40
   50 call rtmi (xi,froot,fct1,dfloat(0),a,1.d-6,23,ier)
      if (ier.eq.0) go to 60
      call rtmi (xi,froot,fct1,dfloat(0),a,1.d-3,23,ier)
      if (ier.eq.0) go to 60
      call rtmi (xi,froot,fct1,dfloat(0),a,5.d-3,23,ier)
      if (ier.eq.0) go to 60
      call rtmi (xi,froot,fct2,dfloat(0),a,5.d-3,35,ier)
      if (ier.eq.0) go to 60
      call rtmi (xi,froot,fct2,0.5*a,a,5.d-3,35,ier)
      if (ier.eq.0) go to 60
      call rtmi (xi,froot,fct2,dfloat(0),0.6*a,5.d-3,35,ier)
      if (ier.eq.0) go to 60
      stop 21
      stop 157
   60 iconv=1
      go to 20
   70 y=((b*xioasq-x*hxi))/yden
   80 en=en0*xi
      if (xi.le.compmy) go to 90
      compmy=xi
      icmax=i
      jcmax=j
      kcmax=k
      if (xi.le.compmx) go to 90
      compmx=xi
c       write (156,*) ' greater compression ',compmx,i,j,k
   90 continue
      eps=em0n0*y
      if (compe.gt.y) go to 100
      compe=y
  100 continue
      if (j.ne.1) go to 110
      adens(i,k)=en
  110 continue
      if (keqost.le.0) go to 120
      hxi=eqost(xi)
      pijk=eosgm1*eps+em0n0*hxi
      go to 150
  120 continue
c      z=0.1453*939.5731*eps
c      az=0.1453*xi
c      z=1000.*eps
      z=eps ! eps is already in MeV/fm3 in the rest frame
      az=xi

c Angel, 2024, not needed with new EoS
c       presx=eps/3.-4.*bbag/3. !pressure in MeV/fm3
c       if (presx .lt. 0.) presx=0.
c       if (presx .gt. 0.) z=eps
c       if (presx .eq. 0.) z=3.*bbag*((eps+bbag)/(4.*bbag))**(4./3.)
      
      if (keqost.eq.-2) go to 130
      if (keqost.eq.-3) go to 130
      if (keqost.eq.-9) go to 140
c_                                      __________________eqsqcd3
  130 call eqsqcd3 (az,z)
      pijk=pq !/(939.5731*0.1453)
      pres=pq
c       write (156,*) 'line 1123',az,z,pq
      temp=tq
      chem=amuq
      entro=sq

c     ADD ON FEBRUARY 19, 2024
      if (pijk .ge. z) then
      write (*,*) 'pijk>e (line8553):', 'pijk= ',pijk,'eijk= ',z
c       pause
      endif


      open (unit=333,file='check12281',status='unknown')
c      write(333,*) i,j,k,pijk,press,temp,chem,entro
c       if (az.eq.0.) write (156,*) ' az=0-x',i,j,k,az,z
c       if (pijk.lt.0.) write (156,*) ' interx',pijk,pq,az,z
      go to 150
c_                                      __________________prxx
  140 call prxx (az,z)
      pijk=pres 
  150 continue
      prsijk=pijk
      epsc=0.
      if (xi.ne.0.) xth=(eps/xi-eusual)*econ-epsc
c xth is the thermal energy density
      if (xth.lt.0.) xth=0.
      xarg=0.
      if (xi.ne.0.) xarg=0.086673*xth/xi**.66666666
      temp=0.
      kplasma=0
      stop 'temp1'
      if (temp1.ne.0.) temp=xth/temp1
      if (ipres.lt.0.and.ipres.gt.-50.and.jtemp.eq.1) then
        iplasma=iplasma+1
        sumplast=sumplast+tempzz
        kplasma=1
c  kplasma is one only if QGP, otherwise zero
      endif
      if (keqost.le.0) temp=tempzz
      cpntot=cpntot+a
      xid=xi
      if (xi.gt.dnsmax) xid=dnsmax
      if (temp.eq.0.) temp=0.001
      if (x.lt.-1.e-5) go to 160
      x=amax1(x,0.d+00)
      repx=dsqrt(x/cmom)
      uijk=repx*xmom(ijk)
      vijk=repx*ymom(ijk)
      wijk=repx*zmom(ijk)
      if (uijk**2+vijk**2+wijk**2 .gt. 1.) then
      write (*,*) 'v^2>1'
      stop 8642
      endif
  160 continue
      return
c
  170 format (' ier=',i2,'  a=',5f15.8)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine intern (i,j,k,ijk,iphs3,adens,ib4,numlvv,uijk,vijk,wijk
     1 ,pijk,reqwe)
c
      implicit real*8 (a-h,o-z)
      integer*2 koost
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
c
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfz
c
      common /pscmz/ zz
c
      common /yestemp/ jtemp
c
      common /maxinf/ icmax, jcmax, kcmax
c     
      common /prstbl/ prs(99,102), trs(99,102), pres, tempzz, chem,
     1 amurs(99,102), entrs(99,102), entro, ipres
c
      common /other/ eusual, epn0, t0, epion, ppion, edelt, delta, pdelt
     1 , temp, tempkn, eps, xth, epsc, dnsmax, ncycpt, ncycrp, ncycaf
c
      common /matric/ storev(929801), idensz(500), edensz(500)
c
      common /numcmm/ numbk, numcel, numcst, numext, numfz1, numlv,
     1 numout, numouu, numou2, jdim1, kfreez, kfrzcy, press, freez,
     2 compmy, koast, kext, cpntot, compmx, compe,iundee(58)
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      common /thestf/ rap00, Tfin, dxv, dyv, dzv, epsilon00, b0, rho00,
     1 e00
c
      common /newcrp/ entot, dentot, tempav, numcll
c
      common /rot/ a, b, c
c
      common /densav/ densav
c
      common /plcm7/ kextrn(30542411), koost(30542411)
c
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c
      common /mssss/ outn, porton, npxz, npt, kout
c
      common /hxi33/ hxi3a, xi, hxi
c
      common /plasma/ sumplast, iplasma, kplasma
c
      common /eqqcd/ tq, eq1, pq, sq, anq1, amuq
c      
      common /cnuclr/ hvc, hv2, pi0, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
c
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233), px(85877055), py(85877055), pz(85877055), numlvv(50),
     3 adens(300,700), koun(27), akoun(27)
c
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
c
      external fct1, funx4, fun, funwr
c
c      entrop(x)=0.5213+1.5*dlog(x+dexp(-2.*0.5213/3.))+1.809*dsqrt(x)/(
c     1 1.139*dsqrt(x)+1.417*x+1.014*x*dsqrt(x))
c
      tempfn(x)=1.5*x/(x+dexp(dble(-2.*.5213/3.)))+0.9045*dsqrt(dble(x))
     1 *(1.-1.417*x-2.028*x*dsqrt(dble(x)))/(1.+1.139*dsqrt(dble(x))+1.4
     2 17*x+1.014*x*dsqrt(dble(x)))**2
c
c     After calling phase1 e(ijk),zmom(ijk),ymom(ijk) and xmom(ijk) are
c     given in MeV. capn(ijk) is multiplied by dxdydz and pressure is
c     given in MeV/fm3 b is given in MeV*fm-3, cmom and c in MeV2*fm-6
      kext=0
      if (capn(ijk).eq.0.) go to 690
      pijk=0.
      a=capn(ijk)   ! N0*V
      b=e(ijk) ! T00 in MeV
      cmom=(xmom(ijk))**2+(ymom(ijk))**2+(zmom(ijk))**2 ! MeV2
      c=cmom*remnsq ! MeV2*fm-6
      a=reqwe*a ! reqwe=1
      b=reqwe*b ! reqwe=1
      c=reqwe*reqwe*c ! reqwe=1
      ra=1./a
      bsq=b*b ! MeV2
      num=0
      koun(1)=0
      do 10 iq=1,27
      koun(iq)=0
      akoun(iq)=0.
   10 continue
      call koccupy (i,j,k,iow,io6w,jdim1)
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (c.gt.bsq) koost(ijk0)=1
      if (iow.eq.1.and.iphs3.eq.2) return
      if (iow.eq.0.and.iphs3.eq.1) return
      if (iow.eq.0) go to 160
c
c      if here, then i,j,k is an interior cell
c
      xi=0.
      iconv=0
      numit=0
      if (c.le.bsq) go to 90
c
c  this section computes the average of the neighbouring cells
c  velocities and pressures for the case of a coasting cell (which
c  occurs when c>b*b or p/e=u>1)
c
   20 ipjk=ijk+nq
      imjk=ijk-nq
      ijpk=ijk+nqib
      ijmk=ijk-nqib
      ijkp=ijk+nqijb
      ijkm=ijk-nqijb
c
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      koost(ijk0)=1
      cells=0.
      velsum=0.
c     koast=1 - changed to eliminate pressure difference!!!!
      koast=0
      prssum=0.
      x=u(ijk)**2+v(ijk)**2+w(ijk)**2
c     write(*,*) 'a',i,j,k,u(ijk),v(ijk),w(ijk),x
      if (i.eq.ibar.or.capn(ipjk).eq.0.) go to 30
      cells=cells+1.
      velsum=velsum+u(ipjk)**2+v(ipjk)**2+w(ipjk)**2
c     write(*,*) 'b',i,j,k,u(ipjk),v(ipjk),w(ipjk),velsum
      prssum=p(ipjk) ! in MeV*fm-3
   30 if (i.eq.1.or.capn(imjk).eq.0.) go to 40
      cells=cells+1.
      velsum=velsum+u(imjk)**2+v(imjk)**2+w(imjk)**2
      prssum=prssum+p(imjk)
   40 if (j.eq.jbar.or.capn(ijpk).eq.0.) go to 50
      cells=cells+1.
      velsum=velsum+u(ijpk)**2+v(ijpk)**2+w(ijpk)**2
      prssum=prssum+p(ijpk)
   50 if (j.eq.1.or.capn(ijmk).eq.0.) go to 60
      cells=cells+1.
      velsum=velsum+u(ijmk)**2+v(ijmk)**2+w(ijmk)**2
      prssum=prssum+p(ijmk)
   60 if (k.eq.kbar.or.capn(ijkp).eq.0.) go to 70
      cells=cells+1.
      velsum=velsum+u(ijkp)**2+v(ijkp)**2+w(ijkp)**2
      prssum=prssum+p(ijkp)
   70 if (k.eq.1.or.capn(ijkm).eq.0.) go to 80
      cells=cells+1.
      velsum=velsum+u(ijkm)**2+v(ijkm)**2+w(ijkm)**2
      prssum=prssum+p(ijkm)
   80 if (cells.gt.0.) x=velsum/cells
      if (cells.gt.0.) prssum=prssum/(cells+1.)
      if (x.gt.1.) x=.9
      go to 570
   90 continue
c
      xmaxa=a
      deps=0.000000001
      xlower=0.00000000001
      call rtmi (xx,froot,fun,xlower,xmaxa,deps,18,ier)
      if (ier.eq.0) go to 100
      call rtmi (xx,froot,fun,xlower,xmaxa,100.*deps,20,ier)
      if (ier.eq.0) go to 100
      if (ier.ne.0) write (*,*) 'ier=',ier,ncyc,i,j,k
      deps=0.00001
      xlower=0.001
c
      call rtmi (xx,froot,fun,xlower,xmaxa,deps,24,ier)
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.95*xmaxa,deps,25,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.90*xmaxa,deps,25,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.85*xmaxa,deps,30,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.80*xmaxa,deps,30,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.75*xmaxa,deps,30,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.70*xmaxa,deps,30,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.65*xmaxa,deps,30,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.60*xmaxa,deps,30,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.55*xmaxa,deps,30,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.50*xmaxa,deps,30,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.45*xmaxa,deps,30,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.40*xmaxa,deps,30,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,xlower,.35*xmaxa,deps,30,ier
     1 )
      if (ier.ne.0) call rtmi (xx,froot,fun,0.5*xmaxa,xmaxa,deps,30,ier)
      if (ier.ne.0) call rtmi (xx,froot,fun,0.85*xmaxa,xmaxa,10.*deps,30
     1 ,ier)
      if (ier.ne.0) call rtmi (xx,froot,fun,0.75*xmaxa,xmaxa,10.*deps,30
     1 ,ier)
c       if (ier.ne.0) write (156,*) ' rtmi error ',i,j,k,ier,xx,a,b,c,bsq
c      1 ,xmaxa
      if (ier.ne.0) call tablt
  100 xv=1.-xx*xx/(a*a)
      if (ier.ne.0) write (*,*) 'ier=',ier,xx,a
      yy=b-dsqrt(c*xv)
      xi=xx
      y=yy
      x=xv
      if (keqost.le.0) go to 590
      hxi=eqost(xi)
      xioasq=(xi*ra)**2
      x=1.-xioasq
      yden=1.+bcoef-bcoef*xioasq
c
      if (iconv.eq.1) go to 150
      ww=c*yden**2
      z=((1.+bcoef)*b+hxi)**2
      if (ww.eq.0.and.z.eq.0) then
        woz=1
      else
        woz=ww/z
      endif
      if (woz.gt.1.) go to 110
      xinew=a*dsqrt(1.-woz)
      delxi=abs(xinew-xi)
      xi=xinew
      numit=numit+1
      if (delxi.lt.1.e-5) go to 140
      if (numit.lt.500) go to 90
  110 xi=0.
      numit=0
  120 xioasq=(xi*ra)**2
      x=1.-xioasq
      if (x.le.0.) go to 130
      xiterm=eqost(xi)
      yden=5.-2.*xioasq
      ww=dsqrt(c)*yden
      xinew=0.5*((1.+bcoef)*b-ww/dsqrt(x)+bcoef*xi+xiterm)
      delxi=abs(xinew-xi)
      xi=xinew
      numit=numit+1
      if (delxi.lt.1.e-5) go to 140
      if (xi.lt.0.) go to 130
      if (numit.lt.500) go to 120
  130 call rtmi (xi,froot,fct1,dfloat(0),a,deps,23,ier)
      if (ier.eq.0) go to 140
      call rtmi (xi,froot,fct1,dfloat(0),a,1000.*deps,23,ier)
      if (ier.eq.0) go to 140
      call rtmi (xi,froot,fct1,dfloat(0),a,5000.*deps,23,ier)
      if (ier.eq.0) go to 140
      if (ier.ne.0) go to 670
      xi=a*dsqrt(1.-x)
  140 iconv=1
      go to 90
  150 y=((b*xioasq-x*hxi))/yden
c       if (dabs(xi-xx).gt.0.00001) write (156,730) a,b,c,xi,xx
      go to 590
c
c      i,j,k is an exterior cell  ______________
c
c      The original version of the code set the density to 1.0 for all
c      exterior cells.  However, after many iterations essentially all
c      (>95%) of the occupied cells were exterior.  We thus replace
c      the original algorithm for finding the new density by:
c
c       1. If one or more of the neighbouring cells are interior cells
c         replace the density by the average of the densities of the
c         neighbouring interior cells.  This has the property that the
c         density is correct for the initial iterations (and thus
c         reproduces that aspect of the original algorithm) and,
c         moreover, is nearer the correct physical situation.
c
c       2. If none of the neighbouring cells is an interior cell,
c         then follow the procedure below:
c         store the coordinates of the cell i,j,k in array kmove,
c         (this is no longer necessary, as bit 1 of ymom determines
c            whether the cell is external)
c         and procede to calculate the density in the
c         remaining cells.  After one has attempted to
c         calculate the density in all the cells, go back to the
c         external cells having no neighbours.  For each of these
c         cells, the routine AVER is called.  First, try to calculate
c         the density as the average of the densities of external
c         cells having at least one internal neighbour.
c         If this fails, treat the cell as an internal cell.
c
  160 xi=1.0
c  change
      ijkm=ijk-nqijb
      imjkm=ijkm-nq
      ipjkm=ijkm+nq
      ijmkm=ijkm-nqib
      imjmkm=ijmkm-nq
      ipjmkm=ijmkm+nq
      ijpkm=ijkm+nqib
      imjpkm=ijpkm-nq
      ipjpkm=ijpkm+nq
      inter=1
      x=1.0
      num=0
      sum=0.
c  change
      if (c.gt.bsq) go to 20
c
c      (i-1,j,k)
      if (i.eq.1) go to 170
      call koccupy (i-1,j,k,iow,io6w,jdim1)
      ics=1
      if (iow.eq.0) go to 170
      num=num+1
      ijkz=imjkm+nqijb
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c      (i+1,j,k)
  170 if (i.eq.ibar) go to 180
      call koccupy (i+1,j,k,iow,io6w,jdim1)
      ics=2
      if (iow.eq.0) go to 180
      num=num+1
      ijkz=ipjkm+nqijb
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i,j-1,k)
  180 if (jdim1.eq.1) go to 220
      if (j.eq.1) go to 190
      call koccupy (i,j-1,k,iow,io6w,jdim1)
      ics=3
      ijkz=ijmkm+nqijb
      if (iow.eq.0) go to 210
      num=num+1
      ics=3
      ijkz=ijmkm+nqijb
      go to 200
c
c      following lines take into account the plane of
c      symmetry about y = 0.
c
  190 continue
      call koccupy (i,j,k,iow,io6w,jdim1)
      if (iow.eq.0) go to 210
      num=num+1
      ics=3
      ijkz=ijk
  200 continue
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i,j+1,k)
  210 if (j.eq.jbar) go to 220
      call koccupy (i,j+1,k,iow,io6w,jdim1)
      ics=4
      if (iow.eq.0) go to 220
      num=num+1
      ics=4
      ijkz=ijpkm+nqijb
      if (ijkz.gt.nlc) go to 220
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i,j,k-1)
  220 if (k.eq.1) go to 230
      call koccupy (i,j,k-1,iow,io6w,jdim1)
      ics=5
      if (iow.eq.0) go to 230
      num=num+1
      ics=5
      ijkz=ijk-nqijb
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i,j,k+1)
  230 if (i.eq.1) go to 240
      call koccupy (i,j,k+1,iow,io6w,jdim1)
      ics=6
      if (iow.eq.0) go to 240
      num=num+1
      ics=6
      ijkz=ijk+nqijb
      iqw=nq*(ijbar*(k)+ibar*(j-1)+i-1)+1
      if (ijkz.gt.nlc) go to 240
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
  240 continue
c
c    now do the next nearest neighbours
c
c       could multiply num and sum here to weight the above cells
c
c      (i-1,j,k-1)
      if (i.eq.1) go to 250
      if (k.eq.1) go to 250
      call koccupy (i-1,j,k-1,iow,io6w,jdim1)
      if (iow.eq.0) go to 250
      num=num+1
      ics=7
      ijkz=imjkm
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i+1,j,k-1)
c
  250 if (k.eq.1) go to 260
      if (i.eq.ibar) go to 260
      call koccupy (i+1,j,k-1,iow,io6w,jdim1)
      if (iow.eq.0) go to 260
      num=num+1
      ics=8
      ijkz=ipjkm
      if (ijkz.gt.nlc) go to 260
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i-1,j,k+1)
c
  260 if (i.eq.1) go to 270
      if (k.eq.kbar) go to 270
      call koccupy (i-1,j,k+1,iow,io6w,jdim1)
      if (iow.eq.0) go to 270
      num=num+1
      ics=9
      ijkz=imjkm+nqijb+nqijb
      if (ijkz.gt.nlc) go to 270
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i+1,j,k+1)
c
  270 if (i.eq.ibar) go to 280
      if (k.eq.kbar) go to 280
      call koccupy (i+1,j,k+1,iow,io6w,jdim1)
      if (iow.eq.0) go to 280
      num=num+1
      ics=10
      ijkz=ipjkm+nqijb+nqijb
      if (ijkz.gt.nlc) go to 280
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      upper level  j+1
c
  280 if (j.eq.jbar) go to 320
c
c      (i-1,j+1,k)
c
      if (i.eq.1) go to 290
      call koccupy (i-1,j+1,k,iow,io6w,jdim1)
      if (iow.eq.0) go to 290
      num=num+1
      ics=11
      ijkz=imjpkm+nqijb
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i+1,j+1,k)
c
  290 if (i.eq.ibar) go to 300
      call koccupy (i+1,j+1,k,iow,io6w,jdim1)
      if (iow.eq.0) go to 300
      num=num+1
      ics=12
      ijkz=ipjpkm+nqijb
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i,j+1,k-1)
c
  300 if (k.eq.1) go to 310
      call koccupy (i,j+1,k-1,iow,io6w,jdim1)
      if (iow.eq.0) go to 310
      num=num+1
      ics=13
      ijkz=ijpkm
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i,j+1,k+1)
c
  310 if (k.eq.kbar) go to 320
      call koccupy (i,j+1,k+1,iow,io6w,jdim1)
      ics=14
      if (iow.eq.0) go to 320
      num=num+1
      ics=14
      ijkz=ijpkm+nqijb+nqijb
      if (ijkz.gt.nlc) go to 320
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c      lower level   must always check  for reflections
c
  320 jb=j-1
      if (jdim1.eq.1) go to 520
c
c   CHANGE 10/8/2001
      if (j.eq.1) jb=1
c
c      (i-1,j-1,k)
c
      if (i.eq.1) go to 350
      call koccupy (i-1,j-1,k,iow,io6w,jdim1)
      if (iow.eq.0) go to 350
      num=num+1
      ics=15
      if (j.eq.1) go to 330
      ijkz=imjmkm+nqijb
      go to 340
  330 ijkz=imjpkm+nqijb
  340 continue
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i+1,j-1,k)
c
  350 if (i.eq.ibar) go to 380
      call koccupy (i+1,j-1,k,iow,io6w,jdim1)
      if (iow.eq.0) go to 380
      num=num+1
      ics=16
      if (j.eq.1) go to 360
      ijkz=ipjmkm+nqijb
      go to 370
  360 ijkz=ipjpkm+nqijb
  370 continue
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i,j-1,k-1)
c
  380 if (k.eq.1) go to 410
      call koccupy (i,j-1,k-1,iow,io6w,jdim1)
      if (iow.eq.0) go to 410
      num=num+1
      ics=17
      if (j.eq.1) go to 390
      ijkz=ijmkm
      go to 400
  390 ijkz=ijpkm
  400 continue
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      (i,j-1,k+1)
c
  410 if (k.eq.kbar) go to 440
      call koccupy (i,j-1,k+1,iow,io6w,jdim1)
      if (iow.eq.0) go to 440
      num=num+1
      ics=18
      if (j.eq.1) go to 420
      ijkz=ijmkm+nqijb+nqijb
      if (ijkz.gt.nlc) go to 440
      go to 430
  420 ijkz=ijpkm+nqijb+nqijb
  430 continue
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
  440 continue
      if (sum.gt.0.01) go to 530
c
c    now check the corner cells
c
c      lower level
c
      k2pl=nqijb+nqijb
      if (j.eq.1) go to 480
c    need not check if j = 1 since these have already
c       been checked
c
c      i-1,j-1,k-1
c
      if (i.eq.1) go to 450
      if (k.eq.1) go to 450
      call koccupy (i-1,j-1,k-1,iow,io6w,jdim1)
      if (iow.eq.0) go to 450
      num=num+1
      ics=19
      ijkz=imjmkm
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      i-1,j-1,k+1
c
  450 if (i.eq.1) go to 460
      if (k.eq.kbar) go to 460
      call koccupy (i-1,j-1,k+1,iow,io6w,jdim1)
      if (iow.eq.0) go to 460
      num=num+1
      ics=20
      ijkz=imjmkm+k2pl
      if (ijkz.gt.nlc) go to 460
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      i+1,j-1,k-1
c
  460 if (k.eq.1) go to 470
      if (i.eq.ibar) go to 470
      call koccupy (i+1,j-1,k-1,iow,io6w,jdim1)
      if (iow.eq.0) go to 470
      num=num+1
      ics=21
      ijkz=ipjmkm
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      i+1,j-1,k+1
c
  470 if (i.eq.ibar) go to 480
      if (k.eq.kbar) go to 480
      call koccupy (i+1,j-1,k+1,iow,io6w,jdim1)
      if (iow.eq.0) go to 480
      num=num+1
      ics=22
      ijkz=ipjmkm+k2pl
      if (ijkz.gt.nlc) go to 480
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c    check upper level
c
  480 if (jbar.eq.j) go to 520
c
c      i-1,j+1,k-1
c
      if (i.eq.1) go to 490
      if (k.eq.1) go to 490
      call koccupy (i-1,j+1,k-1,iow,io6w,jdim1)
      if (iow.eq.0) go to 490
      num=num+1
      ics=23
      ijkz=imjpkm
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      i-1,j+1,k+1
c
  490 if (i.eq.1) go to 500
      if (k.eq.kbar) go to 500
      call koccupy (i-1,j+1,k+1,iow,io6w,jdim1)
      if (iow.eq.0) go to 500
      num=num+1
      ics=24
      ijkz=imjpkm+k2pl
      if (ijkz.gt.nlc) go to 500
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c      i+1,j+1,k-1
c
  500 if (k.eq.1) go to 510
      if (i.eq.ibar) go to 510
      call koccupy (i+1,j+1,k-1,iow,io6w,jdim1)
      if (iow.eq.0) go to 510
      num=num+1
      ics=25
      ijkz=ipjpkm
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
c      i+1,j+1,k+1
c
  510 if (i.eq.ibar) go to 520
      if (k.eq.kbar) go to 520
      call koccupy (i+1,j+1,k+1,iow,io6w,jdim1)
      if (iow.eq.0) go to 520
      num=num+1
      ics=26
      ijkz=ipjpkm+k2pl
      if (ijkz.gt.nlc) go to 520
      velsum=u(ijkz)*u(ijkz)+v(ijkz)*v(ijkz)+w(ijkz)*w(ijkz)
      gamm=dsqrt(1.-velsum)
      sum=sum+capn(ijkz)*gamm
      koun(num)=ics
      akoun(num)=gamm*capn(ijkz)
c
  520 continue
  530 continue
      if (iphs3.eq.2) pijk=eps/3.-4.*bbag/3.
      if (num.eq.0) go to 560
      xi=sum/float(num)
      if (xi.eq.0.) go to 550
      if (xi.gt.1.0) xi=(2.+sum)/float(num+2)
      if (xi.le.compmy) go to 540
      compmy=xi
      icmax=i
      jcmax=j
      kcmax=k
      if (xi.lt.compmx) go to 540
      compmx=xi
      icmax=i
      jcmax=j
      kcmax=k
  540 call enbkst (ijk,xi)
c  store xi in array ensort
      numbk=numbk+1
      go to 580
c
c  One gets here if num is non-zero, but sum is zero - should not be
c  possible
  550 continue
c       write (156,740) i,j,k,ijk,capn(ijk),capn(ijkz),num
c       write (156,710) num,(koun(iy),iy=1,num)
c       write (156,720) num,(akoun(iy),iy=1,num)
      write (*,*) ' '
      write (*,*) 'Stopping in inter - error, sum=0   num=',num,i,j,k
c       write (156,*) 'Stopping in inter - error, sum=0   num=',num,i,j,k
      stop 2222
c
c  if here, no neighbouring internal cells
c
  560 numext=numext+1
      kext=1
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      kextrn(ijk0)=1
      go to 660
  570 continue
      xi=a*dsqrt(1.-x) ! baryon density in the rest frame, x=vx2+vy2+vz2
      y=b-dsqrt(c*x) ! b=energy density in the rest frame in MeV*fm-3, c=totmom in MeV2
      go to 590
  580 continue
c
c  if here, i,j,k is an exterior cell with density xi
c
      if (c.gt.bsq) go to 20
c
      densav=xi
      call xiter (i,j,k,x,y)
c
  590 en=en0*xi ! baryon density in the rest frame
      if (xi.lt.compmy) go to 600
      compmy=xi
      icmax=i
      jcmax=j
      kcmax=k
      if (xi.lt.compmx) go to 600
      compmx=xi
      nummm=-1
  600 continue
      eps=em0n0*y ! energy density in the rest frame in MeV*fm-3
      if (compe.lt.y) compe=y
      if (j.ne.1) go to 610
      adens(i,k)=en
  610 continue
      if (keqost.le.0) go to 620
      hxi=eqost(xi)
      pijk=eosgm1*eps+em0n0*hxi
      go to 650
  620 continue
c      z=0.1453*939.5731*eps
c      az=0.1453*xi
      xi=xi
      z=eps ! eps is the energy density in the local rest frame in MeV
c      z=eps 
      az=xi ! N0*V in the rest frame.lt.1

c Angel, 2024, not needed with new EoS
c
c       presx=eps/3.-4.*bbag/3. !pressure in MeV/fm3
c       if (presx .lt. 0.) presx=0.
c       if (presx .gt. 0.) z=eps
c       if (presx .eq. 0.) z=3.*bbag*((eps+bbag)/(4.*bbag))**(4./3.)
      
      if (keqost.eq.-2) go to 630
      if (keqost.eq.-3) go to 630
      if (keqost.eq.-9) go to 640





c_                                      __________________eqsqcd3
  630 continue

c     ADD ON FEBRUARY 19, 2024
      if (z.le.0.) then
c       write (*,*) 'e<0 (line9473):','eijk= ',z
c       write (*,*)e(ijk),xmom(ijk),ymom(ijk),zmom(ijk),capn(ijk),p(ijk)
      z=0.
      eps=0.
      endif


      call eqsqcd3 (az,z)


      pijk=pq 
      pres=pq
      temp=tq
      chem=amuq
      entro=sq


      if (pijk .gt. z) then
      write (*,*) 'pijk>e (line9473):', 'pijk= ',pijk,'eijk= ',z
      endif





c       if (az.eq.0.) write (156,*) ' az=0-9',i,j,k,az,z,tq
      go to 650
c_                                      __________________prxx
  640 call prxx (az,z)
      pijk=pres 
      write (*,*) 'pres',pres,pijk
  650 continue
      prsijk=pijk
      if (koast.eq.1) pijk=prssum
      if (pijk.lt.0.) pijk=0.
      if (b.gt.0..and.pijk.gt.b) then
      if (pijk.gt.1.d12) write(*,*) 'WARNING p capped in intern cyc=',
     1 ncyc,' i,j,k=',i,j,k,' p=',pijk,' e_lab=',b,' iphs3=',iphs3
      pijk=b
      endif
      epsc=0.
      if (xi.ne.0.) xth=(eps/xi-eusual)*econ-epsc
c xth is the thermal energy density
      if (xth.lt.0.) xth=0.
      xarg=0.
      if (xi.ne.0.) xarg=0.086673*xth/xi**.66666666
      temp1=tempfn(xarg)
      temp=0.
      kplasma=0
      if (temp1.ne.0.) temp=xth/temp1
      if (ipres.lt.0.and.ipres.gt.-50.and.jtemp.eq.1) then
        iplasma=iplasma+1
        sumplast=sumplast+tempzz
        kplasma=1
c  kplasma is one only if QGP, otherwise zero
      endif
      if (keqost.le.0) temp=tempzz
      cpntot=cpntot+a
      if (temp.eq.0.) temp=0.001
      amattr=a*dxdydz
c
      uijk=0.
      vijk=0.
      wijk=0.
      if (xmom(ijk).eq.0..and.ymom(ijk).eq.0..and.zmom(ijk).eq.0.) go to
     1 660
      x=amax1(x,0.d+00) ! x is set to 1 in line 11035 of this subroutine
      repx=dsqrt(x/cmom) ! [MeV-1]
      uijk=repx*xmom(ijk)
      vijk=repx*ymom(ijk)
      wijk=repx*zmom(ijk)
      vel2=uijk**2+vijk**2+wijk**2
      if (isnan(vel2).or.vel2.ge.1.d0) then
      write (*,*) 'v^2>=1 corrected (9579) cyc/i/j/k=',ncyc,i,j,k,
     1 ' vel2=',vel2
      if (isnan(vel2).or.vel2.le.0.d0) then
      uijk=0.d0
      vijk=0.d0
      wijk=0.d0
      else
      vfac=dsqrt(0.9999d0/vel2)
      uijk=uijk*vfac
      vijk=vijk*vfac
      wijk=wijk*vfac
      endif
      endif
c       if (jdim1.ne.3.and.vijk.ne.0.) write (156,*) 'vijk',i,j,k,ijk,
c      1 iphs3,ymom(ijk),vijk
  660 continue
      return
  670 kstop=1
      go to 680
      kstop=4
  680 continue
c       write (156,750) kstop,t,ncyc,i,j,k,numit,a,b,c,xi,x
      write (*,750) kstop,t,ncyc,i,j,k,numit,a,b,c,xi,x
      write (*,*) 'Pausing-intern'
      stop 670
c      return
  690 continue
      write (*,*) 'Problem in INTERN, CAPN=0',i,j,k,ijk,capn(ijk),iphs3
c       write (156,*) 'Problem in INTERN, CAPN=0',i,j,k,ijk,capn(ijk),
c      1 iphs3
      stop 1239
c
  700 format ('KOAST-intern',4i4,2i10,6f20.12)
  710 format (i3,'  koun:',10i4)
  720 format (i3,'  akoun:',10f10.7)
  730 format (' error in calculation of xi',6f13.7)
  740 format (/,' External cell-error, sum=0.',3i4,i8,2f10.5,i5)
  750 format (' stop no.',i2,' in intern at t=',e12.5,' cyc',i4,' cell (
     1',i3,',',i3,',',i3,')'/' numit=',i4,' a=',f13.9,' b=',e16.9,' c='
     2 ,e16.9,' xi=',e16.9,' x=',e16.9)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function ijkret (ijk)
      implicit real*8 (a-h,o-z)
c
      common /volbuk/ lablbk(1527120), cn(1527120)
c
      data modnum /304139/
c
      k=mod(ijk,modnum)+1
      do 10 i=k,1527120
      if (lablbk(i).eq.ijk) go to 40
      if (lablbk(i).eq.-1) go to 30
   10 continue
c
      do 20 i=1,k
      if (lablbk(i).eq.ijk) go to 40
      if (lablbk(i).eq.-1) go to 30
   20 continue
   30 continue
c       write (156,50) ijk
      ijkret=1
      return
   40 continue
      ijkret=i
      return
c
   50 format (' Trouble in ijkret ',i7)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE gser (gamser,a,x,gln)
      implicit real*8 (a-h,o-z)
      INTEGER ITMAX
c      REAL a, gamser, gln, x, EPS
      PARAMETER (ITMAX=100, EPS=3.e-7)
CU    USES gammln
      INTEGER n
c      REAL ap, del, sum, gammln
      gln=gammln(a)
      if (x.le.0.) then
        if (x.lt.0.) stop 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 10 n=1,ITMAX
      ap=ap+1.
      del=del*x/ap
      sum=sum+del
      if (dabs(del).lt.dabs(sum)*EPS) go to 20
   10 continue
      stop 'a too large, ITMAX too small in gser'
   20 gamser=sum*dexp(-x+a*log(x)-gln)
      return
      END

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE gcf (gammcf,a,x,gln)
c
      implicit real*8 (a-h,o-z)
c
      INTEGER ITMAX
c      REAL a, gammcf, gln, x, EPS, FPMIN
      PARAMETER (ITMAX=100, EPS=3.e-7, FPMIN=1.e-30)
CU    USES gammln
      INTEGER i
c      REAL an, b, c, d, del, h, gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      done=dfloat(1)
      do 10 i=1,ITMAX
      an=-i*(i-a)
      b=b+dfloat(2)
      d=an*d+b
      if (dabs(d).lt.FPMIN) d=FPMIN
      c=b+an/c
      if (dabs(c).lt.FPMIN) c=FPMIN
      d=1./d
      del=d*c
      h=h*del
      if (dabs(del-done).lt.EPS) go to 20
   10 continue
      stop 'a too large, ITMAX too small in gcf'
   20 gammcf=dexp(-x+a*log(x)-gln)*h
      return
      END

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function gamnix1 (x)
      implicit real*8 (a-h,o-z)
      common /gamsol/ en, bbagy, rrr, amom, bbagde
      common /redund/ a, b, c
      fac=(4.-x*x)/(4.*b-4.*bbagy/139.14)
      gamnix1=x-a*dsqrt(dfloat(1)-c*fac*fac)
c       write (156,10) x,a,b,c,fac,gamnix1,a*dsqrt(dfloat(1)-c*fac*fac),
c      1 1.-c*fac*fac
      return
c
   10 format ('nix1 ',f10.7,3f10.4,6f10.5)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function gamnix (xi)
      implicit real*8 (a-h,o-z)
      common /gamsol/ en, bbagy, rrr, amom, bbagde
      common /redund/ a, b, c
      fac=(4.-xi*xi)/(4.*b-4.*bbagde)
      gamnix=xi-a*dsqrt(dfloat(1)-c*fac*fac)
c       if (dabs(gamnix).lt.0.1) write (156,10) xi,a,b,c,fac,gamnix,c*(
c      1 4.-xi*xi)**2,1./xi,bbagde
      return
c
   10 format ('nix ',f10.7,3f10.4,2f10.5,2e20.12,f10.5)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function gammp (a,x)
      implicit real*8 (a-h,o-z)
c      REAL a, gammp, x
CU    USES gcf,gser
c      REAL gammcf, gamser, gln
      if (x.lt.0..or.a.le.0.) stop 'bad arguments in gammp'
      if (x.lt.a+1.) then
        call gser (gamser,a,x,gln)
        gammp=gamser
      else
        call gcf (gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      END

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function gammln (xx)
      implicit real*8 (a-h,o-z)
c      REAL gammln, xx
      INTEGER j
      DOUBLE PRECISION ser, stp, tmp, x, y, cof(6)
      SAVE cof, stp
      DATA cof, stp /76.18009172947146d0,-86.50532032941677d0,24.0140982
     1 4083091d0,-1.231739572450155d0,.1208650973866179d-2,-.53952393849
     2 53d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 10 j=1,6
      y=y+1.d0
      ser=ser+cof(j)/y
   10 continue
      gammln=tmp+log(stp*ser/x)
      return
      END

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c       function gammh (x)
c       implicit real*8 (a-h,o-z)
c c      implicit double precision (a-h,o-z)
c       common /gamsol/ en, bbag, rrr, amom, bbagde
c c      bbagz=bbag/erat
c       fac=16.*(x-1.)*x
c       erat=136.14
c       en31=en*erat
c       eps=(3.*en31+4.*(x-1.)*bbag)/(4.*x-1.)
c       eps=0.75*erat*dsqrt(amom/(x*(x-1.)))+bbag
c       write (156,*) dsqrt(amom/(x*(x-1.))),0.75*erat*dsqrt(amom/(x*(x-
c      1 1.)))
c       bb=bbag/eps
c       am2=fac*(1.-bb)**2/((4.*x-1.-4.*bb*(x-1.))**2)
c       gammh=am2-rrr
c       pres=eps/3.-4.*bbag/3.
c       return
c c
c       end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c       function gammfp (x)
c       implicit real*8 (a-h,o-z)
c       common /gamsol/ en, bbagy, rrr, amomz, bbagde
c       bbag=bbagy
c       fac=16.*(x-1.)*x
c       en31=en*136.14
c       eps=(3.*en31+4.*(x-1.)*bbag)/(4.*x-1.)
c       bb=bbag/eps
c       am2=fac*(1.-bb)**2/((4.*x-1.-4.*bb*(x-1.))**2)
c       gammfp=am2-rrr
c       pres=eps/3.-4.*bbag/3.
c       return
c c
c       end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c       function gammf (x)
c       implicit real*8 (a-h,o-z)
c c
c       common /gamsol/ en, bbag, rrr, amom, bbagde
c       erat=136.14
c       fac=16.*(x-1.)*x
c       en31=en*erat
c       eps=(3.*en31+4.*(x-1.)*bbag)/(4.*x-1.)
c       bb=bbag/eps
c       am2=fac*(1.-bb)**2/((4.*x-1.-4.*bb*(x-1.))**2)
c       gammf=am2-rrr
c       pres=eps/3.-4.*bbag/3.
c       return
c c
c       end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function funx4 (x)
      implicit real*8 (a-h,o-z)
c
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      common /prstbl/ prs(99,102), trs(99,102), pres, temp, chem, amurs(
     1 99,102), entrs(99,102), entro, ipres
c
      common /rot/ a, b, c
c
      common /densav/ densav
c
      common /eqqcd/ tq, eq1, pq, sq, anq1, amuq
c
      common /cnuclr/ hvc, hv2, pi0, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
c
      y=b-dsqrt(c*x) ! energy density in the local rest frame in MeV/fm3
      if (x.lt.0.999) y1=y*densav/(a*dsqrt(1.-x))


c      z=0.1453*939.5731*y1
c      az=0.1453*densav
c      z=1000.*y1
      z=y1 ! energy density in the local rest frame in MeV/fm3
      az=densav

      if (keqost.le.0) go to 10
      hxi=eqost(densav)
c    y=(densav*b/a*dsqrt(1.-x)-x*hxi)/(1.+bcoef*x)
c  both y and y1 give the same results
      pres=bcoef*y1+hxi
      go to 40
   10 continue

c Angel, 2024, not needed with new EoS
c       presx=y1/3.-4.*bbag/3. !pressure in MeV/fm3
c       if (presx .lt. 0.) presx=0.
c       if (presx .gt. 0.) z=y1
c       if (presx .eq. 0.) z=3.*bbag*((y1+bbag)/(4.*bbag))**(4./3.)

      if (keqost.eq.-2) go to 20
      if (keqost.eq.-3) go to 20
      if (keqost.eq.-9) go to 30
c_                                      __________________eqsqcd3
   20 call eqsqcd3 (az,z)
      pijk=pq !/(939.5731*0.1453)
      pres=pq !/(939.5731*0.1453)
      temp=tq
      chem=amuq
      entro=sq

c     ADD ON FEBRUARY 19, 2024
      if (pijk .ge. z) then
      write (*,*) 'pijk>e (line9849):', 'pijk= ',pijk,'eijk= ',z
c       pause
      endif
c       if (az.eq.0.) write (156,*) ' az=0-7',az,z
      go to 40
c_                                      __________________prxx
   30 call prxx (az,z)
      pres=pres !/(0.1453*939.5731)
   40 continue
      funx4=x-c/(b+a*pres*dsqrt(1.-x)/densav)**2
      return
c
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c       function funwr (an)
c c
c       implicit real*8 (a-h,o-z)
c c
c       common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
c      1 icysum, kresrt
c c
c       common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
c      1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
c      2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
c      3 kbarmx, modcyc, ioptos4
c c
c       common /prstbl/ prs(99,102), trs(99,102), pres, temp, chem, amurs(
c      1 99,102), entrs(99,102), entro, ipres
c c
c       common /rot/ a, b, c
c c
c       common /eqqcd/ tq, eq1, pq, sq, anq1, amuq
c       common /rtrt/ iqwe
c c
c       common /transp/ xte, yte, zte, hpx, hmx, hpy, hmy, hpz, hmz, i, j,
c      1 k, l, ueff, veff, weff, cneff, xmeff, ymeff, zmeff, eneff, koko
c c
c       common /cnuclr/ hvc, hv2, pi0, amp, amn, amnuc, ak, an0, bin, a2j,
c      1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
c      2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
c      3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
c c
c       x=1.-an*an/(a*a)
c       iqwe=1
c       y=b-dsqrt(c*x) ! energy density in the local rest frame in MeV/fm3
c       if (keqost.le.0) go to 10
c       hsi=eqost(an)
c       pres=bcoef*y+hsi
c       z=0.1453*939.5731*y
c       az=0.1453*an
c       go to 40
c    10 continue
c c      z=0.1453*939.5731*y
c c      az=0.1453*an
c
c c      z=1000.*y
c       z=y ! energy density in the local rest frame in MeV/fm3
c       az=an
c c Angel, 2024, not needed with new EoS
c c       presx=y/3.-4.*bbag/3. !pressure in MeV/fm3
c c       if (presx .lt. 0.) presx=0.
c c       if (presx .gt. 0.) z=y
c c       if (presx .eq. 0.) z=3.*bbag*((y+bbag)/(4.*bbag))**(4./3.)
c
c       if (keqost.eq.-2) go to 20
c       if (keqost.eq.-3) go to 20
c       if (keqost.eq.-9) go to 30
c c_                                      __________________eqsqcd3
c    20 call eqsqcd3 (az,z)
c       pijk=pq !/(939.5731*0.1453)
c       iqwe=0
c       pres=pq !/(939.5731*0.1453)
c       temp=tq
c       chem=amuq
c       entro=sq
c       open (unit=333,file='check11560',status='unknown')
c       write(333,*) i,j,k,pijk,press,temp,chem,entro
c c       if (az.eq.0.) write (156,*) ' az=0-6',i,j,k,az,z
c       go to 40
c c_                                      __________________prxx
c    30 if (az.gt.0.) call prxx (az,z)
c       if (az.eq.0.) pres=z/3.
c c_   Inserted on 5/10/00
c       pr=pres
c       pres=pres !/(0.1453*939.5731)
c    40 funwr=c-(b+pres)*(b+pres)*x
c       gsq=1./(1.-x)
c       write (156,50) an,funwr,a,b,c,pres,x,y,gsq,dsqrt(c/(gsq*gsq*x))
c      1 ,dsqrt(c/(gsq*gsq*x))-pres-y
c       return
c c
c    50 format ('anw',f12.7,f14.7,'  a,b,c',3f10.4,' p',3f11.6,3x,f11.7,2f
c      1 16.7)
c       end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c       function funalt (x)
c       implicit real*8 (a-h,o-z)
c c      implicit double precision (a-h,o-z)
c       common /gamsol/ en, bbag, rrr, amom, bbagde
c       common /redund/ a, b, c
c c   Here, x = xi^2
c       right=1.-c*(4.-x)**2/(16.*(b-bbagde)**2)
c       funalt=x-a*a*right
c       xv=1.-x/a**2
c       anum=(b*(1.-xv)-4.*xv*bbagde/3.)
c       eps=anum/(1.+xv*b)
c       pres=eps/3.-4.*bbag/3./(0.1453*939.5731)
c       if (dabs(funalt).lt.0.1) write (156,10) x,xv,eps,pres,right,funalt
c      1 ,bbagde
c       return
c c
c    10 format (' funalt',f13.8,' xv=',f12.9,' eps=',2f13.8,f10.5,2f12.8)
c       end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function fun (an)
c
      implicit real*8 (a-h,o-z)
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c
      common /prstbl/ prs(99,102), trs(99,102), pres, temp, chem, amurs(
     1 99,102), entrs(99,102), entro, ipres
c
      common /rot/ a, b, c
c
      common /eqqcd/ tq, eq1, pq, sq, anq1, amuq
c
      common /transp/ xte, yte, zte, hpx, hmx, hpy, hmy, hpz, hmz, i, j,
     1 k, l, ueff, veff, weff, cneff, xmeff, ymeff, zmeff, eneff, koko
c
      common /cnuclr/ hvc, hv2, pi0, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
c
      x=1.-an*an/(a*a)
      if (x.eq.1.) x=0.999999999999
      y=b-dsqrt(c*x) ! energy density in the local rest frame in MeV/fm3 
      if (keqost.le.0) go to 10
      hsi=eqost(an)
      pres=bcoef*y+hsi
c      z=0.1453*939.5731*y
c      az=0.1453*an

      z=1000.*y
      az=an
      go to 40
   10 continue
c      z=0.1453*939.5731*y
c      az=0.1453*an
c      z=1000.*y
      z=y ! energy density in the rest frame in MeV/fm3
      az=an

c Angel 2024 not need with new EoS
c       presx=y/3.-4.*bbag/3. !pressure in MeV/fm3
c       if (presx .lt. 0.) presx=0.
c       if (presx .gt. 0.) z=y
c       if (presx .eq. 0.) z=3.*bbag*((y+bbag)/(4.*bbag))**(4./3.)
      
      if (keqost.eq.-2) go to 20
      if (keqost.eq.-3) go to 20
      if (keqost.eq.-9) go to 30
c_                                      __________________eqsqcd3
      
   20 continue
c       ADD ON FEBRUARY 19, 2024
      if (z.le.0.) then
c       write (*,*) 'e<0 (line 10042):','eijk= ',z
c       write (*,*) a,b,pijk
      z=0.
      y=0.
      endif

      call eqsqcd3 (az,z)
   

      pijk=pq
      pres=pq
      temp=tq
      chem=amuq
      entro=sq

      if (pijk .gt. z) then
      write (*,*) 'pijk>e (line10042):', 'pijk= ',pijk,'eijk= ',z
      endif


c c     ADD ON FEBRUARY 19, 2024
c       if (pijk .ge. z) then
c       write (*,*) 'pijk>e (line10042):', 'pijk= ',pijk,'eijk= ',z
c       pijk=0.
c       z=0.
c       az=0.
c        b=0. !b --> e
c        a=0. !a --> capn
c        c=0. !c --> mom**2
c        pres=0.
c        temp=0.
c        chem=0.
c        entro=0.
c c       pause
c       endif

      go to 40
c_                                      __________________prxx
   30 if (az.gt.0.) call prxx (az,z)
      if (az.eq.0.) pres=z/3.
c_   Inserted on 5/10/00
      pr=pres
      pres=pres !/(0.1453*939.5731)
   40 fun=c-(b+pres)*(b+pres)*x
      gsq=1./(1.-x)
      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




      function feqcd3 (t)
      implicit real*8 (a-h,o-z)
c
      common /cnuclr/ hvc, hv2, pi, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      common /eqqcd/ tq, eq1, pq, sq, anq, amuq
c
      common /fq3/ op9anq, anq13
c
      common /uniq/ mqcdwt
c
      t2=t*t ! in MeV2
      amuq=0.
      if (anq13.eq.0) go to 10
      aa=t2*op9anq/3. ! op9anq in MeV-2, aa adimensional
      rr=(0.5+dsqrt(dble(aa**3+0.25)))**(1./3.)
      amuq=(rr-aa/rr)*anq13*qq9 ! qq9 in MeV*fm, amuq in MeV
   10 amuq2=amuq*amuq
      eq=(aq3*amuq2*amuq2+alam4+aq1*t2*t2+aq2*amuq2*t2)*ph3
      feqcd3=eq1-eq
      tq=t
      return
c
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function fct2 (x)
      implicit real*8 (a-h,o-z)
c
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c
      common /rot/ a, b, c
c
      s=eqost(x)
      t=1.+bcoef-bcoef*x*x/(a*a)
      u=s+b*(1.+bcoef)
      v=c*t*t/(u*u)
      if (v.gt.1.) go to 10
      fct2=x-a*dsqrt(1.-v)
c       write (156,20) x,fct2,a,b,c,s,t,u,v
      return
   10 continue
      fct2=1.
c       write (156,20) x,fct2,a,b,c,s,t,u,v
      return
   20 format (' fct2 ',4f10.5,f11.5,4f10.5)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      function fct1 (x)
      implicit real*8 (a-h,o-z)
c
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c
      common /rot/ a, b, c
c
      s=eqost(x)
      t=1.+bcoef-bcoef*x*x/(a*a)
      u=s+b*(1.+bcoef)
      v=c*t*t/(u*u)
      if (v.gt.1.) go to 10
      fct1=x-a*dsqrt(1.-v)
      return
   10 continue
      fct1=1.
      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      function erf (x)
      implicit real*8 (a-h,o-z)
      half=dfloat(1)/dfloat(2)
      if (x.lt.0.) then
        erf=-gammp(half,x**2)
      else
        erf=gammp(half,x**2)
      endif
      return
      END

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine eqsqcd3 (anq,eq)
      implicit real*8 (a-h,o-z)
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      common /rtrt/ iqwe
c
      common /uniq/ mqcdwt
c
      common /cnuclr/ hvc, hv2, pi, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
c
      common /eqqcd/ tq, eq1, pq, sq, anq1, amuq
c
      common /transp/ vdummy(9), i, j, k, l, udummy(9)
c
      common /fq3/ op9anq, anq13
c
      external feqcd3
!      open (unit=5,file='check',status='unknown')

      xn0 = 0.15d0 ! fm-3
      ene0 = 146d0 ! MeV/fm3
      if (anq.le.0.) go to 20 !anq is baryon density in the rest frame
      pq=presstab(eq,anq)
      if (pq .lt. 0.) pq=0.
      tq=temptab(eq,anq)
      amuq=chemtab(eq,anq)
      if (tq > 0.0) then
      	  sq = (eq + pq - amuq * anq) / tq;
      else
          sq = 0.
      endif
      return
   20 continue
       
      tq=0.
      sq=0.
      pq=0.
      amuq=0.
      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function eqost (x) ! x= baryon density ib the rest frame
      implicit real*8 (a-h,o-z)
c
      common /newpar/ a, b, c, an0, bn0, cn0, a1, a2, a0, b0, c0
c  function to return e. o. s.
c  x is the particle number density  n/n0
c
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c     Now keqost=-2
      if (iabs(keqost).eq.2) go to 10
      if (iabs(keqost).eq.4) go to 40
      if (iabs(keqost).eq.3) go to 30
      if (iabs(keqost).eq.5) go to 80
      if (iabs(keqost).eq.6) go to 110
      if (iabs(keqost).eq.7) go to 120
      if (iabs(keqost).eq.8) go to 130
      if (iabs(keqost).eq.9) go to 140
      if (iabs(keqost).ne.1) stop 41
      y=x**.66666666666666666
      eqost=0.02134*(0.666666666666667-bcoef)*x*y-0.07409*x*x*(1.-bcoef)
     1 +0.03592*x*x*y*(1.6666666666666667-bcoef)-bcoef*x
      return
c
c  arne's e.o.s
c
c       5/3=1.66666666
   10 if (x.gt.epsnew) go to 20
      eqost=0.02226*x**(1.66666666666666667)*(0.6666666666666-bcoef)
     1 -barn*x*x*(1.-bcoef)-bcoef*x ! bcoef=1./3.
      return
   20 eqost=0.000238563436*comprs*(x*(1.-bcoef)+dsqrt(x)*(2.*bcoef-1)
     1 -bcoef)+bcoef*0.008588283692-bcoef ! bcoef=1./3.
      eqost=x*eqost
      return
   30 y=x**.6666666666666
      eqost=-0.02470*x*x+0.03592*x*x*y-bcoef*x
      return
   40 if (x.gt.epsnew) go to 50
      eqost=0.02226*x**(1.66666666666666667)*(0.6666666666666-bcoef)
     1 -barn*x*x*(1.-bcoef)-bcoef*x
      return
   50 if (x.gt.a1) go to 60
      eqost=0.000238563436*a*(x*(1.-bcoef)+dsqrt(x)*(2.*bcoef-1)-bcoef)
     1 +bcoef*0.008588283692-bcoef
      eqost=x*eqost
c
      return
   60 if (x.gt.a2) go to 70
      eqost=0.0002385*b*(x*bn0*(1.-bcoef)+dsqrt(x)*(2.*bcoef-1.)*bn0**1.
     1 5-bcoef*bn0*bn0)-bcoef-bcoef*b0/931.57311
      eqost=x*eqost
      return
   70 eqost=0.0002385*c*(x*cn0*(1.-bcoef)+dsqrt(x)*(2.*bcoef-1.)*cn0**1.
     1 5-bcoef*cn0*cn0)-bcoef-bcoef*c0/931.57311
      eqost=x*eqost
      return
c
c  kapusta's e.o.s
c
   80 if (x.gt.epsnew) go to 90
      eqost=0.02226*x**(1.66666666666666667)*(0.6666666666666-bcoef)
     1 -barn*x*x*(1.-bcoef)-bcoef*x
      return
   90 if (x.gt.1.) go to 100
      eqost=0.000238563436*comprs*(x*(1.-bcoef)+dsqrt(x)*(2.*bcoef-1)
     1 -bcoef)+bcoef*0.008588283692-bcoef
      eqost=x*eqost
c
c  0.000238563436 = 2./9./931.57311
c
      return
  100 eqost=comprs*0.000039761*(-5.+5.*x-2.*x*dlog(x))-0.660941144*x
      return
  110 eqost=comprs*(x-1.)*(2.*x-bcoef*(x-1.))/18./939.5731-bcoef*(1.-.00
     1 85882836)
      eqost=eqost*x
      return
  120 eqost=comprs*(x*(1.-bcoef)+2.*bcoef-(1.+bcoef)/x)/18./939.5731
     1 -bcoef*(1.-.0085882836)
      eqost=x*eqost
      return
  130 if (x.lt.an0) go to 110
      eqost=a*x*(1.-bcoef)/939.5731-c*(1.+bcoef)/x/939.5731-bcoef*(1.+b/
     1 939.5731)
      eqost=x*eqost
      return
  140 eqost=0.
      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine enbkst (ijk,v)
c
      implicit real*8 (a-h,o-z)
c
      common /bucket/ lablbk(1527120), ensort(1527120), eza(1527120), ezb(1527120), ezc(1527120), ezd(1527120), cnzzz(1527120)
      data modnum /304139/
c
      k=mod(ijk,modnum)+1
      do 10 i=k,1527120
      if (lablbk(i).ne.-1) go to 10
      lablbk(i)=ijk
      ensort(i)=v
      go to 30
   10 continue
c
      do 20 i=1,k
      if (lablbk(i).ne.-1) go to 20
      lablbk(i)=ijk
      ensort(i)=v
      go to 30
   20 continue
c       write (156,40) ijk,v
      write (*,40) ijk,v
      stop 43
   30 continue
      return
c
   40 format (' out of space in enbkst ',i9,f10.6)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c    subroutine eigen
c
c       purpose
c         compute eigenvalues and eigenvectors of a real symmetric
c         matrix
c
c       usage
c         call eigen(a,r,n,mv)
c
c       description of parameters
c         a - original matrix (symmetric), destroyed in computation.
c            resultant eigenvalues are developed in diagonal of
c            matrix a in descending order.
c         r - resultant matrix of eigenvectors (stored columnwise,
c            in same sequence as eigenvalues)
c         n - order of matrices a and r
c         mv- input code
c               0  compute eigenvalues and eigenvectors
c               1  compute eigenvalues only (r need not be
c                   dimensioned but must still appear in calling
c                   sequence)
c
c       remarks
c         original matrix a must be real symmetric (storage mode=1)
c         matrix a cannot be in the same location as matrix r
c
c       subroutines and function subprograms required
c         none
c
c       method
c         diagonalization method originated by jacobi and adapted
c         by von neumann for large computers as found in 'mathematical
c         methods for digital computers', edited by a. ralston and
c         h.s. wilf, john wiley and sons, new york, 1962, chapter 7
c
      subroutine eigen (a,r,n,mv)
      implicit real*8 (a-h,o-z)
c
      dimension a(1), r(1)
c
c       ...............................................................
c
c       if a double precision version of this routine is desired, the
c       c in column 1 should be removed from the double precision
c       statement which follows.
c
c    double precision a,r,anorm,anrmx,thr,x,y,sinx,sinx2,cosx,
c   1              cosx2,sincs,range
c
c       the c must also be removed from double precision statements
c       appearing in other routines used in conjunction with this
c       routine.
c
c       the double precision version of this subroutine must also
c       contain double precision fortran functions.  sqrt in statements
c       40, 68, 75, and 78 must be changed to dsqrt.  abs in statemen
c       62 must be changed to d abs. the constant in statement 5 should
c       be changed to 1.0d-12.
c
c       ...............................................................
c
c       generate identity matrix
c
      range=1.0d-15
c      if (mv-1) 10,40,10
      if (mv.eq.0) go to 30
      iq=-n
      do 20 j=1,n
      iq=iq+n
      do 10 i=1,n
      ij=iq+i
      r(ij)=0.0
      if (i.ne.j) go to 10
c      if (i-j) 30,20,30
      r(ij)=1.0
   10 continue
   20 continue
c
c       compute initial and final norms (anorm and anormx)
c
   30 anorm=0.0
      do 50 i=1,n
      do 40 j=i,n
c      if (i-j) 50,60,50
      if (i.eq.j) go to 40
      ia=i+(j*j-j)/2
      anorm=anorm+a(ia)*a(ia)
   40 continue
   50 continue
c      if (anorm) 320,320,70
      if (anorm.le.0.) go to 210
      anorm=1.414*dsqrt(anorm)
      anrmx=anorm*range/float(n)
c
c       initialize indicators and compute threshold, thr
c
      ind=0
      thr=anorm
   60 thr=thr/float(n)
   70 l=1
   80 m=l+1
c
c       compute sin and cos
c
   90 mq=(m*m-m)/2
      lq=(l*l-l)/2
      lm=l+mq
      if (dabs(a(lm))-thr.lt.0.) go to 170
c      if(dabs(a(lm))-thr) 250,120,120
      ind=1
      ll=l+lq
      mm=m+mq
      x=0.5*(a(ll)-a(mm))
      y=-a(lm)/dsqrt(a(lm)*a(lm)+x*x)
      if (x.gt.0.) go to 100
c      if (x) 130,140,140
      y=-y
  100 sinx=y/dsqrt(2.0*(1.0+(dsqrt(1.0-y*y))))
      sinx2=sinx*sinx
      cosx=dsqrt(1.0-sinx2)
      cosx2=cosx*cosx
      sincs=sinx*cosx
c
c       rotate l and m columns
c
      ilq=n*(l-1)
      imq=n*(m-1)
      do 160 i=1,n
      iq=(i*i-i)/2
      if (i-l.eq.0) go to 150
      if (i-m.eq.0) go to 150
      if (i-m.gt.0) go to 110
c      if (i-l) 150,220, 150
c  150 if (i-m) 160,220,170
      im=i+mq
      go to 120
  110 im=m+iq
  120 if (i-l.ge.0) go to 130
c  180 if (i-l) 190,200,200
      il=i+lq
      go to 140
  130 il=l+iq
  140 x=a(il)*cosx-a(im)*sinx
      a(im)=a(il)*sinx+a(im)*cosx
      a(il)=x
  150 if (mv-1.eq.0) go to 160
c  220 if (mv-1) 230,240,230
      ilr=ilq+i
      imr=imq+i
      x=r(ilr)*cosx-r(imr)*sinx
      r(imr)=r(ilr)*sinx+r(imr)*cosx
      r(ilr)=x
  160 continue
      x=2.0*a(lm)*sincs
      y=a(ll)*cosx2+a(mm)*sinx2-x
      x=a(ll)*sinx2+a(mm)*cosx2+x
      a(lm)=(a(ll)-a(mm))*sincs+a(lm)*(cosx2-sinx2)
      a(ll)=y
      a(mm)=x
c
c       tests for completion
c
c       test for m = last column
c
  170 if (m-n.eq.0) go to 180
c  250 if (m-n) 260,270,260
      m=m+1
      go to 90
c
c       test for l = second from last column
c
  180 if (l-(n-1).eq.0) go to 190
c  270 if (l-(n-1)) 280,290,280
      l=l+1
      go to 80
c  290 if (ind-1) 310,300,310
  190 if (ind-1.ne.0) go to 200
      ind=0
      go to 70
c
c       compare threshold with final norm
c
c  310 if (thr-anrmx) 320,320,80
  200 if (thr-anrmx.gt.0.) go to 60
c
c       sort eigenvalues and eigenvectors
c
  210 iq=-n
      do 240 i=1,n
      iq=iq+n
      ll=i+(i*i-i)/2
      jq=n*(i-2)
      do 230 j=i,n
      jq=jq+n
      mm=j+(j*j-j)/2
      if (a(ll)-a(mm).ge.0.) go to 230
c      if (a(ll)-a(mm)) 330,360,360
      x=a(ll)
      a(ll)=a(mm)
      a(mm)=x
      if (mv-1.eq.0) go to 230
c      if (mv-1) 340,360,340
      do 220 k=1,n
      ilr=iq+k
      imr=jq+k
      x=r(ilr)
      r(ilr)=r(imr)
      r(imr)=x
  220 continue
  230 continue
  240 continue
      return
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine dmill (num,i100,i0,i)
c   routine to separate a number into millions, thousands, and units
      i100=num/1000000
      k=num-1000000*i100
      i0=k/1000
      i=k-1000*i0
      return
      end
      

c      ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine counms (iq,jq,kq,num,num2,kfrz)
c
      implicit real*8 (a-h,o-z)
c      real*4 aa1
c
c  kfrz=1 -> freeze out
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
c
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
c
      common /other/ eusual, epn0, t0, epion, ppion, edelt, delta, pdelt
     1 , temp, tempkn, eps, xth, epsc, dnsmax, ncycpt, ncycrp, ncycaf
c
      common /extrmes/ xextrn, xextrx, yextrx, zextrn, zextrx, kox
c
      dimension px(85877055), py(85877055), pz(85877055)
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c
      num=0
      num2=0
c
c   Following are the extreme positions of particles in this cell
      xextrx=0.
      yextrx=0.
      zextrx=0.
      xextrn=1110.
      zextrn=1110.
      ko=-2
   10 ko=ko+3
      if (ko.gt.kp) go to 80
      xte=px(ko)
      yte=py(ko)
      zte=pz(ko)
c
      l=1
      if (ko.gt.kpbas(1).and.ko.le.kpbas(2)) l=2
      if (ko.gt.kpbas(3).and.ko.le.kpbas(4)) l=2
      i=xte+1.
      j=yte+1.
      k=zte+1.
      if (i.ne.iq) go to 70
      if (j.ne.jq) go to 70
      if (k.ne.kq) go to 70
c if here, particle is in cell
      if (l.eq.2) num2=num2+1
      if (l.eq.1) num=num+1
      if (kfrz.lt.0) go to 30
      if (kfrz.gt.0) go to 60
c
c    Search for smallest, largest x and corresponding y,z
      if (xte.gt.xextrn) go to 20
c   new minimum
      xextrn=xte
      yextrn=yte
      zextrn=zte
      zextrx=zte
      kox=ko
c      go to 20
   20 continue
      if (xte.lt.xextrx) go to 70
c   new maximum
      xextrx=xte
      yextrx=yte
      zextrn=zte
      zextrx=zte
      kox=ko
      go to 70
c
c    Search for smallest, largest z and corresponding x,y
   30 continue
      if (kfrz.eq.-2) go to 50
      if (zte.gt.zextrn) go to 40
c   new minimum
      zextrn=zte
      xextrn=xte
      xextrx=xte
      yextrn=yte
      kox=ko
c      go to 20
   40 continue
      if (zte.lt.zextrx) go to 70
c   new maximum
      xextrn=xte
      xextrx=xte
      yextrx=yte
      zextrx=zte
      kox=ko
      go to 70
   50 continue
      if (yte.lt.yextrx) go to 70
c   new maximum
      xextrx=xte
      xextrn=xte
      zextrn=zte
      yextrx=yte
      zextrx=zte
      kox=ko
      go to 70
   60 continue
c
c  give the particle a large x value so as to remove it from the mesh
      xn=999.
      px(ko)=xn
   70 continue
      go to 10
   80 continue
      return
c
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine cnucl
c
      implicit real*8 (a-h,o-z)
c
c    calculation of  nuclear constants  and  eos
c
      common /cnuclr/ hvc, hv2, pi, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
c
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c
c  units in mev, fm      c=k=1
c
      common /eqtype/ kcomp, kther, kdelt, kpion, kqcd
c
c---------------------------------------------------------------
c
c  type of the equation of state is regulated by the following
c  parameters :
c
c  kcomp  see subroutine compr1
c  kther  see subroutine therm1
c  kdelt  see subroutine delta1
c  kpion  not implemented  massless pions are standard
c  kqcd  not implemented  bag eqoation of state is standard
c
c       constants
c
      hvc=197.3285751 ! MeV*fm
      hv2=hvc*hvc ! MeV2*fm2
      hv3=hvc*hv2
      pi=3.1415926535897932384626
      pi2=pi*pi
      amp=938.279627
      amn=939.5731
      amnuc=0.5*(amp+amn)
c
c    compressibility
c
      ak=210.
c
c    binding energy
c
      bin=-8.0
c
c    normal nuclear density (fm-3)
c
      an0=0.1453
c
c    low t fermi gas expansion parameter
c
      bfer=1.809
c
c    compression energy parameters from j,kapusta prc29(84)1735
c
      a2j=21.1
      a3j=-38.3
      a4j=-26.7
      a5j=35.9
c
c for parameters  of QGP
c
      deg=4
      ent=(hv2*2.*pi/amn)**(3./2.)
      ph3=1./hv3
      ampi=139.5688
      ampi0=134.9645
      zq01=81.*pi2*hv3
      aaq01=1./(hvc*pi**0.666666666)
c
c                        bag constant and critical temperature
c
      alam=235.
      alam4=alam**4.
      bbag=alam4*ph3
      tcr00=90./pi2
      if (keqost.eq.-2) go to 20
      if (keqost.eq.-3) go to 10
      if (keqost.eq.-9) go to 20
      go to 30
   10 continue
c  For Nf=3  See L.Csernai: Intro... Wiley 94, p. 88, 128, 138 & 173
      anflav=3.
      aq1=19.*pi2/12.
      aq2=1./2.
      aq3=1./(36.*pi2)
      aq4=1./3.
      aq5=19.*pi2/9.
      zq1=zq01/3.
      aaq1=9.*aaq01
      alim=27.
      amin=-19.
      atq1=3.*hv3
      atq2=1./pi2
      qq9=(dfloat(27)*pi2*hv3)**(1./3.) ! MeV*fm
      op9=pi**(2./3.)/hv2
      tcr=alam*(tcr00/(47.5-3.))**0.25
      go to 30
   20 continue
c  For Nf=2  See L.Csernai: Intro... Wiley 94, p. 88, 128, 138 & 173
      anflav=2.
      aq1=37.*pi2/30.
      aq2=1./3.
      aq3=1./(54.*pi2)
      aq4=2./9.
      aq5=74.*pi2/45.
      zq1=zq01/2.
      aaq1=15.*0.666666666**(0.333333333)*aaq01
      alim=45.
      amin=-37.
      atq1=9.*hv3/2.
      atq2=1./pi2
      qq9=(81.*pi2*hv3/2.)**(1./3.)
      op9=(2.*pi/3.)**(2./3.)/hv2 ! in (MeV*fm)-2
      tcr=alam*(tcr00/(37.-3.))**0.2
   30 continue
c
c          parameters for delta1
c
      amdelta=1232.0
      degdelt=16.
      entdelt=(hv2*2.*pi/amdelta)**(3./2.)
      ddd=amdelta-amnuc
      dcoef=4.*(amdelta/amnuc)**(3./2.)
c
c          compr parameters for "grant-kapusta" prc32(85)663
c
      agk=400.
      ak=200.
c
c   type for start (parameters set in cnucl!!!)
c
      kcomp=2
      kther=1
      kdelt=1
      kpion=0
      kqcd=0
c
c                   "quadratic", "boltzmann",  "no delta"
c
c  **************************************************************
c
c    variable parameters may be defined or modified here by read
c
      an0=0.1453
      ak=240.
      bin=-0.
      amnuc=amn
      ddd=amdelta-amnuc
      return
c
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine chksym
c
c   Check file is symmetric
c
      implicit real*8 (a-h,o-z)
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
c
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfz
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      dimension px(85877055), py(85877055), pz(85877055)
c
      common /midle/ imid, kmid, kend, imidpt, kmidpt, klmin, krmax,
     1 kdisp, ipmin, ipmax, imid1
c
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233)
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
      if (kslab.eq.8) return
c
c   find mid line
      j=1
      imidhr=(ipmin+ipmax)/2
      numerr=0
c       if (kdebug.ne.0) write (156,*) ' chksym',imid,kmid,klmin,krmax,
c      1 kdisp,ipmin,ipmax
      k=kmid
      kb=kmid+1
      j=1
      do 10 i=ipmin,ipmax
      ij=nq*(ibar*(j-1)+i-1)+1
   10 continue
      do 50 j=1,30
c
      do 40 k=klmin,kmid
      kb=krmax-k+klmin
      do 30 i=ipmin,imid
c
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      ib=ipmax-i+ipmin
      ijkb=nq*(ijbar*(kb-1)+ibar*(j-1)+i-1)+1
      if1=nq*ijbar*(kb-1)
      if2=nq*ijbar*(k-1)
      ig1=nq*(ijbar*(kb-1)+ibar*(j-1))
      ig2=nq*(ijbar*(kb-1)+ibar*(j-1))
      ih1=nq*ijbar*(kb-1)+ibar*(j-1)+nq*(i-1)
      ih2=nq*ijbar*(k-1)+ibar*(j-1)+nq*(i-1)
      ij1=ih1+nq*(i-1)
      ij2=ih2+nq*(i-1)
      ibjk=nq*(ijbar*(k-1)+ibar*(j-1)+ib-1)+1
      ibjkb=nq*(ijbar*(kb-1)+ibar*(j-1)+ib-1)+1
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
c
      if (capn(ijk).eq.0..and.capn(ibjkb).eq.0.) go to 30
c
      diffn=capn(ijk)-capn(ibjkb)
      diffn1=capn(ibjk)-capn(ijkb)
      diffe=e(ijk)-e(ibjkb)
      diffe1=e(ibjk)-e(ijkb)
      diffz=zmom(ijk)+zmom(ibjkb)
      diffz1=zmom(ibjk)+zmom(ijkb)
      diffw=w(ijk)+w(ibjkb)
      diffw1=w(ibjk)+w(ijkb)
      diffp=p(ijk)-p(ibjkb)
      diffp1=p(ibjk)-p(ijkb)
      diff=abs(capn(ijk)-capn(ijkb))
      if (diff.gt.0.00000001) kdisp=1
c
c       if (kdisp.eq.0.and.diff.gt.0.000000001) write (156,70) i,ib,j,k,kb
c      1 ,ijk,ibjkb,diff,capn(ijk),capn(ibjkb)
      if (kdisp.eq.0.and.diff.gt.0.000000001) write (*,70) i,ib,j,k,kb
     1 ,ijk,ibjkb,diff,capn(ijk),capn(ibjkb)
c
      if (dabs(diffn).gt.0.00000001.or.abs(diffn1).gt.0.00000001) numerr
     1 =numerr+1
      if (dabs(diffn).gt.0.00000001.or.abs(diffn1).gt.0.00000001) write
     1 (*,80) i,ib,j,k,kb,diffn,capn(ijk),capn(ibjkb),diffn1,capn(ijkb)
     2 ,capn(ibjk),ijk,ibjkb,ibjk,ijkb
      if (dabs(diffe).gt.0.00000001.or.abs(diffe1).gt.0.00000001) write
     1 (*,90) i,ib,j,k,kb,diffe,e(ijk),e(ibjkb),e(ijkb),e(ibjk),diffe1
     2 ,ijk,ibjkb,ibjk,ijkb
c N
      if (dabs(diffn).gt.0.00000001.or.abs(diffn1).gt.0.00000001) write
     1 (6,80) i,j,k,ib,kb,diffn,capn(ijk),capn(ibjkb),diffn1,capn(ijkb)
     2 ,capn(ibjk),ijk,ibjkb,ibjk,ijkb
c E
      if (dabs(diffe).gt.0.00000001.or.abs(diffe1).gt.0.00000001) write
     1 (6,90) i,j,k,ib,kb,diffe,e(ijk),e(ibjkb),e(ijkb),e(ibjk),diffe1
     2 ,ijk,ibjkb,ibjk,ijkb
c Z
      if (dabs(diffz).gt.0.00000001.or.abs(diffz1).gt.0.00000001) then
c         write (156,100) i,ib,j,k,kb,diffz,zmom(ijk),zmom(ibjkb),diffz1
c      1  ,zmom(ibjk),w(ijkb),ijk,ibjkb,ibjk,ijkb,capn(ijk),capn(ibjkb)
c      2  ,capn(ibjk),capn(ijkb)
        kmin=999
        do 20 kq=1,kbar
        kqb=krmax-kq+klmin
        ijkq=nq*(ijbar*(kq-1)+ibar*(j-1)+i-1)+1
        ibjkbq=nq*(ijbar*(kqb-1)+ibar*(j-1)+ib-1)+1
        if (capn(ijkq).eq.0..and.capn(ibjkbq).eq.0.) go to 20
        zd=zmom(ijkq)+zmom(ibjkbq)
c         write (156,60) kq,kqb,ijkq,capn(ijkq),zmom(ijkq),ibjkbq,capn
c      1  (ibjkbq),zmom(ibjkbq),zd
   20   continue
      endif
c W
      if (dabs(diffw).gt.0.00000001.or.abs(diffw1).gt.0.00000001) write
     1 (6,110) i,j,k,ib,kb,diffw,w(ijk),w(ibjkb),diffw1,w(ibjk),w(ijkb)
     2 ,ijk,ibjkb,ibjk,ijkb
c P
      if (dabs(diffp).gt.0.00000001.or.abs(diffp1).gt.0.00000001) write
     1 (6,120) i,j,k,ib,kb,diffp,p(ijk),p(ibjkb),diffp1,p(ibjk),p(ijkb)
     2 ,ijk,ibjkb,ibjk,ijkb
      if (dabs(diffp).gt.0.00000001.or.abs(diffp1).gt.0.00000001) write
     1 (6,*) i,j,k,ib,kb,' E=',e(ijk),e(ibjkb),e(ibjk),e(ijkb)
      if (dabs(diffp).gt.0.00000001.or.abs(diffp1).gt.0.00000001) write
     1 (6,*) i,j,k,ib,kb,' N=',capn(ijk),capn(ibjkb),capn(ibjk),capn
     2 (ijkb)
      if (numerr.ge.10) stop 'in CHKSYM NUMERR >10'
c
   30 continue
   40 continue
   50 continue
      return
c
   60 format (2i4,2(i10,2f13.7),f14.8)
   70 format (83x,'|',/,'For i=',i3,i4,' j=',i2,' k=',i3,i4,2i8,' capn d
     1iff=',3f12.8,' NOT SYMMETRIC about imid',3x,'|',/,83x,'|')
   80 format (5i4,' N Error ',2(' err=',f11.8,'  N=',2f12.8),2x,2i8,2x,2
     1 i8/)
   90 format (5i4,' E Error ',2(' err=',f13.6,'  E=',2f12.6),2x,2i8,2x,2
     1 i8/)
  100 format (5i4,' Z Error ',2(' err=',f11.5,'  Z=',2f12.6),2x,2i8,2x,2
     1 i8,/,30x,'CAPN:',4f12.8,/)
  110 format (5i4,' W Error ',2(' err=',f11.8,'  W=',2f12.8),2x,2i8,2x,2
     1 i8/)
  120 format (5i4,' P Error ',2(' err=',f11.8,'  P=',2f12.8),2x,2i8,2x,2
     1 i8/)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine checks
c
      implicit real*8 (a-h,o-z)
c      real*4 aa1, u, v, w
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
c
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
c
      common /numcmm/ numbk, numcel, numcst, numext, numfz1, numlv,
     1 numout, numouu, numou2, jdim1, kfreez, kfrzcy, press, freez,
     2 compmy, koast, kext, cpntot, compmx, compe,iundee(58)
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      common /thestf/ rap00, Tfin, dxv, dyv, dzv, epsilon00, b0, rho00,
     1 e00
c
      common /misca/ kmove(700), numlvv(52)
c
      common /dmisca/ adens(300,700), eedns(300,700), aplot(300,700),
     1 bplot(300,700), cplot(300,700)
c
      common /midle/ imid, kmid, kend, imidpt, kmidpt, klmin, krmax,
     1 kdisp, ipmin, ipmax, imid1
c
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233), px(85877055), py(85877055), pz(85877055)
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
c
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c
      kmin=1000
c       if (kdebug.ne.0) write (156,*) 'checks',kdisp
      kmax=0
      j=1
      do 20 i=1,ibar
      do 10 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 10
      if (k.gt.kmax) kmax=k
      if (k.lt.kmin) kmin=k
   10 continue
   20 continue
      klimit=(kmax+kmin)/2-1
c      write(*,*) 'kmin',kmin,kmax,klimit,kbar/2
c
      sume=0.
      do 50 i=1,ibar
      do 40 j=1,jbar
      do 30 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 30
      sume=sume+e(ijk)
   30 continue
   40 continue
   50 continue
      do 80 i=1,ibar
      do 70 j=1,jbar
      do 60 k=1,klimit
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      ib=ibar+1-i
      kb=kbar+3-k
      if (kslab.eq.4.or.kslab.eq.5.or.kslab.eq.7) kb=kb-2
      ijkb=nq*(ijbar*(kb-1)+ibar*(j-1)+ib-1)+1
      if (k.gt.kb) write (*,*) ' k,kb',k,kb,ijk,ijkb,i,j
      if (ijkb.le.0) go to 60
c
      if (capn(ijk).eq.0..and.capn(ijkb).eq.0.) go to 60
c       if (dabs(capn(ijk)-capn(ijkb)).gt.1.e-7) write (156,170) i,ib,j,k
c      1 ,kb,capn(ijk),capn(ijkb),ijk,ijkb
      if (dabs(xmom(ijk)+xmom(ijkb)).gt..00001*abs(xmom(ijk))) write (
     1 156,190) i,ib,j,k,kb,xmom(ijk),xmom(ijkb)
      zdiffk=zmom(ijk)+zmom(ijkb)
      if (dabs(zmom(ijk)+zmom(ijkb)).gt..0001*abs(zmom(ijk))) write (
     1 156,290) i,ib,j,k,kb,zmom(ijk),zmom(ijkb),zdiffk,ibar,kbar,ijk,
     2 ijkb
c       if (e(ijk).lt.0.) write (156,160) ijk,i,j,k,e(ijk)
c       if (dabs(p(ijk)-p(ijkb)).gt.1.e-9*p(ijk)) write (156,240) i,ib,j,k
c      1 ,kb,p(ijk),p(ijkb),dabs(p(ijk)-p(ijkb)),ijk,ijkb
c       if (dabs(w(ijk)+w(ijkb)).gt.1.e-9) write (156,270) i,ib,j,k,kb,
c      1 w(ijk),w(ijkb),abs(w(ijk)+w(ijkb)),ijk,ijkb
c       if (dabs(u(ijk)+u(ijkb)).gt.1.e-7) write (156,260) i,ib,j,k,kb,
c      1 u(ijk),u(ijkb),abs(u(ijk)+u(ijkb)),ijk,ijkb
c       if (dabs(e(ijk)-e(ijkb)).gt.abs(e(ijk))*.000001) write (156,280) i
c      1 ,ib,j,k,kb,e(ijk),e(ijkb),abs(e(ijk)-e(ijkb)),ijk,ijkb
      if (b0.ne.0.) go to 60
      if (kdisp.ne.0) go to 60
c --- new ijkb  Block 2
c
      kb=k
      ijkb=nq*(ijbar*(k-1)+ibar*(j-1)+ib-1)+1
c       if (dabs(capn(ijk)-capn(ijkb)).gt.1.e-7) write (156,180) i,ib,j,k,
c      1 kb,capn(ijk),capn(ijkb),ijk,ijkb
c       if (dabs(xmom(ijk)+xmom(ijkb)).gt..0001*abs(xmom(ijk))) write (
c      1 156,200) i,ib,j,k,kb,xmom(ijk),xmom(ijkb),ijk,ijkb
c       if (dabs(zmom(ijk)-zmom(ijkb)).gt..001*abs(zmom(ijk))) write (156
c      1 ,210) i,ib,j,k,kb,zmom(ijk),zmom(ijkb),zmom(ijk)-zmom(ijkb)
      if (e(ijk).lt.0.) write (*,*) 'XXXXX4  E<0',ijk,i,j,k,e(ijk)
c       if (e(ijk).lt.0.) write (156,*) 'XXXXX4  E<0',ijk,i,j,k,e(ijk)
c       if (dabs(e(ijk)-e(ijkb)).gt.e(ijk)*.001) write (156,220) i,ib,j,k,
c      1 kb,e(ijk),e(ijkb),abs(e(ijk)-e(ijkb))
c       if (dabs(w(ijk)-w(ijkb)).gt.1.e-9) write (156,230) i,ib,j,k,kb,
c      1 w(ijk),w(ijkb),ijk,ijkb
c       if (dabs(u(ijk)+u(ijkb)).gt.1.e-7) write (156,*) 'err u',i,k,ib,u
c      1 (ijk),u(ijkb)
c       if (dabs(p(ijk)-p(ijkb)).gt.1.e-10) write (156,250) i,ib,j,k,kb,p
c      1 (ijk),p(ijkb),dabs(p(ijk)-p(ijkb)),ijk,ijkb
   60 continue
   70 continue
   80 continue
      sumv=0.
      if (kdebug.eq.0) go to 150
      do 110 j=1,6
      if (jdim1.eq.1.and.j.gt.1) go to 110
      do 100 i=1,ibar
      do 90 k=1,kbar
      eedns(i,k)=0.
      adens(i,k)=0.
      aplot(i,k)=0.
      cplot(i,k)=0.
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 90
      eedns(i,k)=zmom(ijk)
      adens(i,k)=capn(ijk)
      u(ijk)=0.
   90 continue
  100 continue
  110 continue
      sume=0.
      do 140 i=1,ibar
      do 130 j=1,jbar
      do 120 k=1,kbar
      ijk=nq*(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      if (capn(ijk).eq.0.) go to 120
      sume=sume+e(ijk)
  120 continue
  130 continue
  140 continue
c       write (156,*) 'out of checks'
  150 continue
      return
c
  160 format (' Checks  E<0',i10,' (',i3,i2,i4,')',f10.5)
  170 format ('checks:err capn1 i=',i3,i4,' j=',i2,' k=',i3,i4,2f13.8,2i
     1 10)
  180 format ('checks:err capn2 i=',i3,i4,' j=',i2,' k=',i3,i4,2f13.8,2i
     1 10)
  190 format ('checks:err xmom1 i=',i3,i4,' j=',i2,' k=',i3,i4,2f13.8)
  200 format ('checks:err xmom2 i=',i3,i4,' j=',i2,' k=',i3,i4,2f13.8,2i
     1 10)
  210 format ('checks:err zmom2 i=',i3,i4,' j=',i2,' k=',i3,i4,3f13.8)
  220 format ('checks:err   E2  i=',i3,i4,' j=',i2,' k=',i3,i4,3f13.8)
  230 format ('checks:err   w2: i=',i3,i4,' j=',i2,' k=',i3,i4,2f13.8,2i
     1 10)
  240 format ('checks:err  P1:  i=',i3,i4,' j=',i2,' k=',i3,i4,3f13.8,2i
     1 10)
  250 format ('checks:err  P2:  i=',i3,i4,' j=',i2,' k=',i3,i4,3f13.8,2i
     1 10,/)
  260 format ('checks:err   u1  i=',i3,i4,' j=',i2,' k=',i3,i4,3f13.8,2i
     1 10)
  270 format ('checks:err   w1  i=',i3,i4,' j=',i2,' k=',i3,i4,3f13.8,2i
     1 10)
  280 format ('checks:err  E1   i=',i3,i4,' j=',i2,' k=',i3,i4,3f13.8,2i
     1 10)
  290 format ('checks:err zmom1 i=',i3,i4,' j=',i2,' k=',2i4,3f16.12,4i1
     1 0)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     xi=capn/gamma baryon density density in the rest frame
c     x=u(ijk)**2+v(ijk)**2+w(ijk)**2
c     y=(3.*(b*xioasq-x*hxi))/yden
c     y=b-dsqrt(c*x) enrgy density in the rest frame
      subroutine aver (i,j,k,ijk,ijmkm,ijpkm,imjkm,ipjkm,xi,bsq,x,y
     1 ,jdim1)
c
      implicit real*8 (a-h,o-z)
      integer*2 koost
c      real*4 aa1, u, v, w
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
c
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
c
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233), koun(30), akoun(30)
c
      common /rot/ a, b, c
c
      common /densav/ densav
c
      common /plcm7/ kextrn(30542411), koost(30542411)
c
      external fct1, fun, funx4
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
c
      numm=-1
      if (c.le.bsq) go to 90
      ijk0=(ijbar*(k-1)+ibar*(j-1)+i-1)+1
      koost(ijk0)=1
   10 ipjk=ijk+nq
      imjk=ijk-nq
      ijpk=ijk+nqib
      ijmk=ijk-nqib
      ijkp=ijk+nqijb
      ijkm=ijk-nqijb
      call counms (i,j,k,numm,numn,0)
      if (numm+numn.ne.0) go to 20
      x=0.
      y=0.
      a=0.
      b=0.
      c=0.
      bsq=0.
      xi=0.
      return
   20 continue
      cells=0.
      velsum=0.
      x=u(ijk)**2+v(ijk)**2+w(ijk)**2
      if (i.eq.ibar.or.capn(ipjk).eq.0.) go to 30
      cells=cells+1.
      velsum=velsum+u(ipjk)**2+v(ipjk)**2+w(ipjk)**2
   30 if (i.eq.1.or.capn(imjk).eq.0.) go to 40
      cells=cells+1.
      velsum=velsum+u(imjk)**2+v(imjk)**2+w(imjk)**2
   40 if (j.eq.jbar.or.capn(ijpk).eq.0.) go to 50
      cells=cells+1.
      velsum=velsum+u(ijpk)**2+v(ijpk)**2+w(ijpk)**2
   50 if (j.eq.1.or.capn(ijmk).eq.0.) go to 60
      cells=cells+1.
      velsum=velsum+u(ijmk)**2+v(ijmk)**2+w(ijmk)**2
   60 if (k.eq.kbar.or.capn(ijkp).eq.0.) go to 70
      cells=cells+1.
      velsum=velsum+u(ijkp)**2+v(ijkp)**2+w(ijkp)**2
   70 if (k.eq.1.or.capn(ijkm).eq.0.) go to 80
      cells=cells+1.
      velsum=velsum+u(ijkm)**2+v(ijkm)**2+w(ijkm)**2
   80 if (cells.gt.0.) x=velsum/cells
      if (cells.eq.0.) go to 500
      if (x.ge.1.0) x=0.9999
      xi=a*dsqrt(1.-x) ! a=capn, xi energy density in the rest frame
      go to 480
   90 continue
      num=0
      imjmkm=ijmkm-nq
      ipjmkm=ijmkm+nq
      imjpkm=ijpkm-nq
      ipjpkm=ijpkm+nq
      sum=0.
c     for Mac shift
      ib1=ishft(1,i-2)
      ib2=ishft(1,i-1)
      ib3=ishft(1,i)
c      (i-1,j,k)
      if (i.eq.1) go to 100
      call koccupy (i-1,j,k,iow,io6w,jdim1)
      if (io6w.eq.0) go to 100
      ics=1
      ijkz=imjkm+nqijb
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 100
      sum=sum+energ
      num=num+1
      koun(num)=ics
c      (i+1,j,k)
  100 if (i.eq.ibar) go to 110
      call koccupy (i+1,j,k,iow,io6w,jdim1)
      if (io6w.eq.0) go to 110
      ics=2
      ijkz=ipjkm+nqijb
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 110
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c      (i,j-1,k)
  110 if (jdim1.eq.1) go to 150
      if (j.eq.1) go to 120
      call koccupy (i,j-1,k,iow,io6w,jdim1)
      if (io6w.eq.0) go to 140
      ics=3
      ijkz=ijmkm+nqijb
      go to 130
c
c      following lines take into account the plane of
c      symmetry about y = 0.
c
  120 continue
      call koccupy (i,j-1,k,iow,io6w,jdim1)
      ics=3
      if (io6w.eq.0) go to 140
      ics=3
      ijkz=ijk
  130 continue
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 140
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i,j+1,k)
  140 if (j.eq.jbar) go to 150
      call koccupy (i,j+1,k,iow,io6w,jdim1)
      if (io6w.eq.0) go to 150
      ics=4
      ijkz=ijpkm+nqijb
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 150
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i,j,k-1)
  150 if (k.eq.1) go to 160
      call koccupy (i,j,k-1,iow,io6w,jdim1)
      if (io6w.eq.0) go to 160
      ics=5
      ijkz=ijk-nqijb
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 160
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i,j,k+1)
  160 if (i.eq.1) go to 170
      call koccupy (i,j,k+1,iow,io6w,jdim1)
      ics=6
      if (io6w.eq.0) go to 170
      ics=6
      ijkz=ijk+nqijb
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 170
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
  170 continue
c
c    now do the next nearest neighbours
c
c       could multiply num and sum here to weight the above cells
c
c      (i-1,j,k-1)
c
      if (i.eq.1) go to 180
      if (k.eq.1) go to 180
      call koccupy (i-1,j,k-1,iow,io6w,jdim1)
      if (io6w.eq.0) go to 180
      ics=7
      ijkz=imjkm
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 180
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i+1,j,k-1)
c
  180 if (k.eq.1) go to 190
      if (i.eq.ibar) go to 190
      call koccupy (i+1,j,k-1,iow,io6w,jdim1)
      if (io6w.eq.0) go to 190
      ics=8
      ijkz=ipjkm
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 190
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i-1,j,k+1)
c
  190 if (i.eq.1) go to 200
      if (k.eq.kbar) go to 200
      call koccupy (i-1,j,k+1,iow,io6w,jdim1)
      if (io6w.eq.0) go to 200
      ics=9
      ijkz=imjkm+nqijb+nqijb
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 200
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i+1,j,k+1)
c
  200 if (i.eq.ibar) go to 210
      if (k.eq.kbar) go to 210
      call koccupy (i+1,j,k+1,iow,io6w,jdim1)
      if (io6w.eq.0) go to 210
      ics=10
      ijkz=ipjkm+nqijb+nqijb
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 210
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      upper level  j+1
c
  210 if (jdim1.eq.1) go to 370
      if (j.eq.jbar) go to 250
c
c      (i-1,j+1,k)
c
      if (i.eq.1) go to 220
      call koccupy (i-1,j+1,k,iow,io6w,jdim1)
      if (io6w.eq.0) go to 220
      ics=11
      ijkz=imjpkm+nqijb
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 220
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i+1,j+1,k)
c
  220 if (i.eq.ibar) go to 230
      call koccupy (i+1,j+1,k,iow,io6w,jdim1)
      if (io6w.eq.0) go to 230
      ics=12
      ijkz=ipjpkm+nqijb
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 230
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i,j+1,k-1)
c
  230 if (k.eq.1) go to 240
      call koccupy (i,j+1,k-1,iow,io6w,jdim1)
      if (io6w.eq.0) go to 240
      ics=13
      ijkz=ijpkm
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 240
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i,j+1,k+1)
c
  240 if (k.eq.kbar) go to 250
      call koccupy (i,j+1,k+1,iow,io6w,jdim1)
      if (io6w.eq.0) go to 250
      ics=14
      ijkz=ijpkm+nqijb+nqijb
      energ=retrev(ijkz)
      if (energ.eq.0.) go to 250
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      lower level   must always check  for reflections
c
  250 jb=j-1
      if (j.eq.1) jb=1
c
c      (i-1,j-1,k)
c
      if (i.eq.1) go to 280
      call koccupy (i-1,j-1,k,iow,io6w,jdim1)
      if (io6w.eq.0) go to 280
      ics=15
      if (j.eq.1) go to 260
      ijkz=imjmkm+nqijb
      energ=retrev(ijkz)
      go to 270
  260 ijkz=imjpkm+nqijb
      energ=retrev(ijkz)
  270 continue
      if (energ.eq.0.) go to 280
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i+1,j-1,k)
c
  280 if (i.eq.ibar) go to 310
      call koccupy (i+1,j-1,k,iow,io6w,jdim1)
      if (io6w.eq.0) go to 310
      ics=16
      if (j.eq.1) go to 290
      ijkz=ipjmkm+nqijb
      energ=retrev(ijkz)
      go to 300
  290 ijkz=ipjpkm+nqijb
      energ=retrev(ijkz)
  300 continue
      if (energ.eq.0.) go to 310
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i,j-1,k-1)
c
  310 if (k.eq.1) go to 340
      call koccupy (i,j-1,k-1,iow,io6w,jdim1)
      if (io6w.eq.0) go to 340
      ics=17
      if (j.eq.1) go to 320
      ijkz=ijmkm
      energ=retrev(ijkz)
      go to 330
  320 ijkz=ijpkm
      energ=retrev(ijkz)
  330 continue
      if (energ.eq.0.) go to 340
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
c
c      (i,j-1,k+1)
c
  340 if (k.eq.kbar) go to 370
      call koccupy (i,j-1,k+1,iow,io6w,jdim1)
      if (io6w.eq.0) go to 370
      ics=18
      if (j.eq.1) go to 350
      ijkz=ijmkm+nqijb+nqijb
      energ=retrev(ijkz)
      go to 360
  350 ijkz=ijpkm+nqijb+nqijb
      energ=retrev(ijkz)
  360 continue
      if (energ.eq.0.) go to 370
      sum=sum+energ
      num=num+1
      koun(num)=ics
      akoun(num)=energ
  370 continue
      if (num.eq.0) go to 430
      if (sum.lt.0.0001) go to 430
      xi=sum/float(num)
      densav=xi
      deps=0.00001
      call rtmi (x,froot,funx4,dfloat(0),dfloat(1),deps,40,ier)
      if (ier.ne.0) call rtmi (x,froot,funx4,dfloat(0),dfloat(1),100.
     1 *deps,40,ier)
      if (ier.ne.0) go to 380
      y=b-dsqrt(c*x) ! energy density inn the rest frame
      return
c
c  no soln - try again
  380 continue
c       write (156,510) ier,a,b,c,xi,sum,num
      if (xi.lt..9.and.a.gt.1.15) go to 430
      hxi=eqost(xi)
      hxi3a=3.*a*hxi
      hxi3a=hxi3a/xi
      numit=0
      x=0.
  390 xnew=c*((3.+2.*x)/(5.*b+hxi3a*dsqrt(1.-x)))**2
      delx=abs(xnew-x)
      x=xnew
      numit=numit+1
      if (delx.lt.1.e-5) go to 420
      if (numit.eq.200) go to 400
      if (x.lt.1.0) go to 390
  400 x=1.0
      numit=0
  410 xnew=1.-(((3.+2.*x)*dsqrt(c/x)-5.*b)/hxi3a)**2
      delx=abs(xnew-x)
      x=xnew
      numit=numit+1
      if (delx.lt.1.e-5) go to 420
      if (numit.eq.200) go to 490
      if (x.gt.0.) go to 410
      go to 490
  420 y=(xi*b/a*dsqrt(1.-x)-x*hxi)/(1.+2./3.*x)
      return
  430 continue
c  changes  4/27/87 & 5/26/01 & 10/2/22
      xlower=0.0000001
      deps=0.000001
      call rtmi (xi,froot,fun,xlower,0.9999999*a,deps,30,ier)
      if (ier.eq.0) return
      if (ier.ne.0) call rtmi (xi,froot,fun,xlower,0.9999999*a,10.*deps,
     1 30,ier)
      if (ier.eq.0) return
      zxero=dfloat(0)
      deps=0.0001
      if (ier.ne.0) call rtmi (xi,froot,fun,zxero,a,deps,30,ier)
      if (ier.ne.0) call rtmi (xi,froot,fun,zxero,0.9999999*a,deps,20
     1 ,ier)
      if (ier.ne.0) call rtmi (xi,froot,fun,zxero,0.9999999*a,10.*deps,2
     1 0,ier)
      if (ier.ne.0) call rtmi (xi,froot,fun,zxero,0.9999999*a,100.*deps,
     1 20,ier)
      if (ier.ne.0) go to 440
      x=1.-xi*xi/(a*a)
      y=b-dsqrt(c*x)
      return
  440 continue
c       write (156,520) ier,a,b,c,i,j,k
      xi=0.
      ra=1./a
      numit=0
  450 hxi=eqost(xi)
      xioasq=(xi*ra)**2
      x=1.-xioasq
      if (x.lt.0.) go to 460
      if (x.gt.1.0) go to 460
      yden=5.-2.*xioasq
      ww=c*yden**2
      z=(5.*b+3.*hxi)**2
      woz=ww/z
      if (woz.gt.1.) go to 460
      xinew=a*dsqrt(1.-woz)
      delxi=abs(xinew-xi)
      xi=xinew
      numit=numit+1
      if (delxi.lt.1.e-5) go to 480
      if (numit.lt.500) go to 450
  460 xi=0.
      numit=0
  470 xioasq=(xi*ra)**2
      x=1.-xioasq
      if (x.le.0.) go to 10
      xiterm=eqost(xi)
      yden=5.-2.*xioasq
      ww=dsqrt(c)*yden
      xinew=0.5*(5.*b+3.*xi*xiterm-ww/dsqrt(x))
      delxi=abs(xinew-xi)
      xi=xinew
      if (xi.lt.0.) go to 10
      numit=numit+1
      if (delxi.lt.1.e-5) go to 480
      if (numit.lt.500) go to 470
      dep=0.002
      call rtmi (xi,froot,fct1,dfloat(0),a,dep,500,ier)
  480 continue
      hxi=eqost(xi)
      ra=1./a
      xioasq=(xi*ra)**2 ! xioasq=1-v2
      x=1.-xioasq ! uijk2+vijk2+wijk2
      if (x.gt.1.) go to 10
      yden=5.-2.*xioasq
      y=(3.*(b*xioasq-x*hxi))/yden ! b energy density T00
      
      return
  490 continue
c       write (156,530) i,j,k
      xi=1.0
      go to 420
c
c    calculate en assuming velocities are known (from previous round)
  500 x=u(ijk)*u(ijk)+v(ijk)*v(ijk)+w(ijk)*w(ijk)
      gamm=dsqrt(1.-x)
      xi=capn(ijk)*gamm
      hxi=eqost(xi)
      go to 420
c
  510 format (' no soln in aver-external',i4,5f12.6,i6)
  520 format (' no soln. in AVER (interior)',i4,3f16.7,3i3)
  530 format (/,' no solution for x in aver ',4i4,/)
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine arne1
c
      implicit real*8 (a-h,o-z)
c
      common /option/ cycenx, bcoef, eta, comprs, epsnew, barn, keqost,
     1 icysum, kresrt
c
      common /newstf/ vlxyz, ioptosq, ismin1, ismax, nxv, nyv, nredux,
     1 nreduy, nreduz, ioptos1, ioptos2, ioptos3, kw, kdebug, ktnc,
     2 kslab, interp, interv, npart, npartz, nparp, nparpz, krandm,
     3 kbarmx, modcyc, ioptos4
c
      data a /20.73242/
      data h0 /-8./
      c=3.*a/(2.*comprs+9.*h0)
      beta=2.*comprs/(2.*comprs+9.*h0)
      iter=0
      epsnew=1./(1.+c)**1.92
   10 iter=iter+1
      epsold=epsnew
      e3=epsold**(.666666666666667)
      e2=dsqrt(epsold)
      del=6.*epsold*(1.-c*e3-beta*e2)/(4.*c*e3+3.*beta*e2)
      epsnew=epsold+del
      if (dabs(del).lt.0.00000001) go to 20
      if (iter.gt.99) go to 30
      go to 10
   20 barn=a*epsnew**(-.333333333333333)-2.*comprs*(1.-2./dsqrt(epsnew)+
     1 1./epsnew)/9.-h0/epsnew
c
c  normalize by dividing by mo
      barn=barn/931.5731
c       if (keqost.eq.9) write (156,40) epsnew,barn,iter
      return
   30 continue
c       write (156,50)
      stop 42
c
   40 format (/,' matching xi,  b',2f12.6,i5)
   50 format (' iter')
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine area (addr,volm,irx,irxbk,i,j,k,l,xte,yte,zte)
c
      implicit real*8 (a-h,o-z)
c      real*4 aa1, u, v, w
c
      common /plcm1/ aa1(91627237)
      common /plcm2/ aa2(91627237)
      common /plcm3/ aa3(91627237)
      common /plcm4/ aa4(85877065)
c
      common /pscm1/ aa(1), ap, at, day, dt, dto, dtodx, dtody, dtodz,
     1 dx, dxdydz, dy, dz, econ, emevp, emevt, emid3, em0, em0n0, en0,
     2 eosgm1, eps0, e0(2), fibar, fjbar, fkbar, fx(2), fz(2), gam(2),
     3 gevp, gevt, hxieq1, ibar, id(20), ijbar, ijkbar, ijkmax, iprt,
     4 jbar, kbar, kp, kpbas(4), kpmax, lcmb, ncyc, nlc, nplot(2), nq,
     5 nqib, nqijb, nsc, numit, numtd, nx(2), ny(2), nz(2), pi, pmove,
     6 pnu(2), rdx, rdxdydz, rdy, rdz, remnsq, rem0n0, ren0, smaj(2),
     7 smin(2), sumnu(2), t, time, tnc, tout, twopi, w0(2), xdisp, zmom0
     8 (2), convpu, dtmod, dzvfm
c
      dimension u(91627233), v(91627233), w(91627233), capn(91627233), p
     1 (91627233), e(91627233), xmom(91627233), ymom(91627233), zmom(91627233)
c
      common /backet/ lbblbk(22), cpnncl(22), convol(22)
c
      dimension volm(8), irx(8), irxbk(8)
c
      equivalence (aa1(1),u), (aa1(2),v), (aa1(3),w), (aa2(1),capn),
     1 (aa2(2),p), (aa2(3),e), (aa3(1),xmom), (aa3(2),ymom), (aa3(3)
     2 ,zmom)
      equivalence (aa4(1),px), (aa4(2),py), (aa4(3),pz)
c
      sumvol=0.
      sumu=0.
      sumv=0.
      sumw=0.
      leave=0
      hpx=aint(xte+1.5)-.5-xte
      hmx=1.-hpx
      hpy=aint(yte+1.5)-.5-yte
      hmy=1.-hpy
      hpz=aint(zte+1.5)-.5-zte
      hmz=1.-hpz
      ir=hpx+.5+float(i)
      jr=hpy+.5+float(j)
      kr=hpz+.5+float(k)
      ir1=nq*(ijbar*(kr-1)+ibar*(jr-1)+ir-1)+1
      ir2=ir1-nq
      ir3=ir2-nqijb
      ir4=ir1-nqijb
      ir5=ir1-nqib
      ir6=ir5-nq
      ir7=ir6-nqijb
      ir8=ir5-nqijb
      do 10 ml=1,8
      volm(ml)=0.
      irx(ml)=1
   10 continue
c  Assumption that index 1 will not correspond to a physical cell
      if (i.eq.1.or.i.eq.ibar) go to 200
      if (j.eq.1.or.j.eq.jbar) go to 200
      if (k.eq.1.or.k.eq.kbar) go to 200
   20 continue
      kr1=50
      kr2=70
      kr3=90
      kr4=110
      kr5=130
      kr6=150
      kr7=170
   30 continue
      vol=hmx*hmy*hmz
      sumvol=sumvol+vol
      kcc=1
      volm(1)=vol
      cq=vol*addr
      if (dabs(cq).lt.0.0000001) go to 40
      call strmoma (ir1,irz,cq)
      if (addr.gt.0.) call strmomb (ir1,irzz,cq)
      convol(irz)=convol(irz)+vol
   40 irx(1)=ir1
      irxbk(1)=irz
      if (kr1.eq.50) go to 50
      if (kr1.eq.90) go to 90
      if (kr1.eq.110) go to 110
      if (kr1.eq.1900) go to 190
      stop 'kr1 error'
c      go to kr1, (50,90,110,190)
   50 continue
      vol=hpx*hmy*hmz
      sumvol=sumvol+vol
      kcc=2
      volm(2)=vol
      cq=vol*addr
      if (dabs(cq).lt.0.0000001) go to 60
      call strmoma (ir2,irz,cq)
      if (addr.gt.0.) call strmomb (ir2,irzz,cq)
      convol(irz)=convol(irz)+vol
c c      if (vol.gt.0.0000001) write (156,610) ir2,irz,kcc,vol,convol(irz)
   60 irx(2)=ir2
      irxbk(2)=irz
      if (kr2.eq.70) go to 70
      if (kr2.eq.110) go to 110
      if (kr2.eq.130) go to 130
      if (kr2.eq.190) go to 190
      stop 'kr2 problem'
   70 continue
      vol=hpx*hmy*hpz
      sumvol=sumvol+vol
      kcc=3
      volm(3)=vol
      cq=vol*addr
      if (dabs(cq).lt.0.0000001) go to 80
      call strmoma (ir3,irz,cq)
      if (addr.gt.0.) call strmomb (ir3,irzz,cq)
      convol(irz)=convol(irz)+vol
   80 irx(3)=ir3
      irxbk(3)=irz
      if (kr3.eq.90) go to 90
      if (kr3.eq.130) go to 130
      if (kr3.eq.150) go to 150
      if (kr3.eq.190) go to 190
      stop 'kr3 problems'
c      go to kr3, (90,130,150,190)
   90 continue
      vol=hmx*hmy*hpz
      sumvol=sumvol+vol
      kcc=4
      volm(4)=vol
      cq=vol*addr
      if (dabs(cq).lt.0.0000001) go to 100
      call strmoma (ir4,irz,cq)
      if (addr.gt.0.) call strmomb (ir4,irzz,cq)
      convol(irz)=convol(irz)+vol
  100 irx(4)=ir4
      irxbk(4)=irz
      if (kr4.eq.110) go to 110
      if (kr4.eq.150) go to 150
      if (kr4.eq.170) go to 170
      if (kr4.eq.190) go to 190
      stop 'kr4 problems'
  110 continue
      vol=hmx*hpy*hmz
      sumvol=sumvol+vol
      kcc=5
      volm(5)=vol
      cq=vol*addr
      if (dabs(cq).lt.0.0000001) go to 120
      call strmoma (ir5,irz,cq)
      if (addr.gt.0.) call strmomb (ir5,irzz,cq)
      convol(irz)=convol(irz)+vol
  120 irx(5)=ir5
      irxbk(5)=irz
      if (kr5.eq.130) go to 130
      if (kr5.eq.150) go to 150
      if (kr5.eq.170) go to 170
      if (kr5.eq.190) go to 190
      stop 'kr5 problems'
  130 continue
      vol=hpx*hpy*hmz
      sumvol=sumvol+vol
      kcc=6
      volm(6)=vol
      cq=vol*addr
      if (dabs(cq).lt.0.0000001) go to 140
      call strmoma (ir6,irz,cq)
      if (addr.gt.0.0000001) call strmomb (ir6,irzz,cq)
      convol(irz)=convol(irz)+vol
  140 irx(6)=ir6
      irxbk(6)=irz
      if (kr6.eq.150) go to 150
      if (kr6.eq.190) go to 190
      stop 'kr6 problem'
c      go to kr6, (150,190)
  150 continue
      vol=hpx*hpy*hpz
      sumvol=sumvol+vol
      kcc=7
      volm(7)=vol
      cq=vol*addr
      if (dabs(cq).lt.0.0000001) go to 160
      call strmoma (ir7,irz,cq)
      if (addr.gt.0.) call strmomb (ir7,irzz,cq)
      convol(irz)=convol(irz)+vol
  160 irx(7)=ir7
      irxbk(7)=irz
      if (kr7.eq.170) go to 170
      if (kr7.eq.190) go to 190
  170 continue
      vol=hmx*hpy*hpz
      sumvol=sumvol+vol
      kcc=8
      volm(8)=vol
      cq=vol*addr
      if (dabs(cq).lt.0.0000001) go to 180
      call strmoma (ir8,irz,cq)
      convol(irz)=convol(irz)+vol
      if (addr.gt.0.) call strmomb (ir8,irzz,cq)
  180 irx(8)=ir8
      irxbk(8)=irz
  190 if (j.eq.1.and.hpy.lt.0.5) sumv=2.*sumv*yte
      return
  200 if (hpy.lt.0.5) go to 210
      jtest=jbar
      koct=1
      go to 220
  210 jtest=1
      koct=2
  220 if (hpx.gt.0.5) go to 230
      if (hpz.ge.0.5) go to 280
      go to 310
  230 if (hpz.lt.0.5) go to 340
      if (i.ne.ibar.and.k.eq.kbar) go to 250
      if (i.eq.ibar.and.k.eq.kbar) go to 260
      if (i.eq.ibar.and.k.ne.kbar) go to 270
  240 if (j.ne.jtest) go to 20
      go to (570,580), koct
  250 if (j.ne.jtest) go to 560
      go to (450,510), koct
  260 if (j.ne.jtest) go to 520
      go to (420,380), koct
  270 if (j.ne.jtest) go to 600
      go to (540,490), koct
  280 if (i.ne.1.and.k.eq.kbar) go to 250
      if (i.eq.1.and.k.eq.kbar) go to 290
      if (i.eq.1.and.k.ne.kbar) go to 300
      go to 240
  290 if (j.ne.jtest) go to 430
      go to (170,390), koct
  300 if (j.ne.jtest) go to 550
      go to (440,470), koct
  310 if (i.eq.1.and.k.ne.1) go to 300
      if (i.eq.1.and.k.eq.1) go to 320
      if (i.ne.1.and.k.eq.1) go to 330
      go to 240
  320 if (j.ne.jtest) go to 480
      go to (400,360), koct
  330 if (j.ne.jtest) go to 590
      go to (530,460), koct
  340 if (i.eq.ibar.and.k.ne.1) go to 270
      if (i.eq.ibar.and.k.eq.1) go to 350
      if (i.ne.ibar.and.k.eq.1) go to 330
      go to 240
  350 if (j.ne.jtest) go to 500
      go to (410,370), koct
  360 kr1=190
c  360 assign 190 to kr1
      go to 30
  370 kr2=190
c     assign 190 to kr2
      go to 50
  380 kr3=190
c     assign 190 to kr3
      go to 70
  390 kr4=190
c     assign 190 to kr4
      go to 90
  400 kr5=190
c     assign 190 to kr5
      go to 110
  410 kr6=190
c assign 190 to kr6
      go to 130
  420 kr7=190
c  assign 190 to kr7
      go to 150
  430 kr4=170
c  assign 170 to kr4
      go to 90
  440 kr5=170
c  assign 170 to kr5
      go to 110
  450 kr7=170
c  assign 170 to kr7
      go to 150
  460 kr1=50
c  assign 50 to kr1
      kr2=190
c  assign 190 to kr2
      go to 30
  470 kr1=90
c  assign 90 to kr1
      kr4=190
c  assign 190 to kr4
      go to 30
  480 kr1=110
c  assign 110 to kr1
      kr5=190
c  assign 190 to kr5
      go to 30
  490 kr2=70
c  assign 70 to kr2
      kr3=190
c  assign 190 to kr3
      go to 50
  500 kr2=130
c  assign 130 to kr2
      kr6=190
c  assign 190 to kr6
      go to 50
  510 kr3=90
c  assign 90 to kr3
      kr4=190
c  assign 190 to kr4
      go to 70
  520 kr3=150
c  assign 150 to kr3
      kr7=190
c  assign 190 to kr7
      go to 70
  530 kr5=130
c  assign 130 to kr5
      kr6=190
c  assign 190 to kr6
      go to 110
  540 kr6=150
c  assign 150 to kr6
      kr7=190
c  assign 190 to kr7
      go to 130
  550 kr1=90
c  assign 90 to kr1
      kr4=110
c  assign 110 to kr4
      kr5=170
c  assign 170 to kr5
      go to 30
  560 kr3=90
c  assign 90 to kr3
      kr4=150
c  assign 150 to kr4
      kr7=170
c  assign 170 to kr7
      go to 70
  570 kr5=130
c  assign 130 to kr5
      kr6=150
c  assign 150 to kr6
      kr7=170
c  assign 170 to kr7
      go to 110
  580 kr1=50
c  assign 50 to kr1
      kr2=70
c  assign 70 to kr2
      kr3=90
c  assign 90 to kr3
      kr4=190
c  assign 190 to kr4
      go to 30
  590 kr1=50
c  assign 50 to kr1
      kr2=110
c  assign 110 to kr2
      kr5=130
c  assign 130 to kr5
      kr6=190
c  assign 190 to kr6
      go to 30
  600 kr2=70
c  assign 70 to kr2
      kr3=130
c  assign 130 to kr3
      kr6=150
c  assign 150 to kr6
      kr7=190
c  assign 190 to kr7
      go to 50
c
      end
