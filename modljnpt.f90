MODULE LJNPTMOD
     IMPLICIT NONE
     INTEGER, SAVE :: N_ATOM,NEQ,MCNT,MSIZE,N_MOL,LRGMSZE,NREPS,EVFREQ
     INTEGER, SAVE :: GRIDSIZE,FILEIND,VFILEIND,QFILEIND,VOLFIND,FRCFIND,PWFRCFIND,EVALFIND,ARRYSIZE,mynode
     INTEGER, PARAMETER :: NUMOBSERV=3,MAXGAUSS=20,MAXTYPE=20,ISIZE=5
     DOUBLE PRECISION, SAVE :: ROECURR,ROENEW,ENOT,BL,LS12,LS6,EAVE,PCURR,PAVE,BLDF,VFAC6,VFAC12
     DOUBLE PRECISION, SAVE :: VNEW,VCURR,VMCONV,LJ12,LJ6,LJNEW12,LJNEW6,PVCURR,VDIFF,EAVEFD,COE12,COE6
     DOUBLE PRECISION, PARAMETER :: ATOMICMASS=0.020614788876D0,PI=3.1415926D0,EVTOL=0.0001D0
     DOUBLE PRECISION, PARAMETER :: ANGTOCM=1.0D-24,NAVOG=6.02214179D23,PCONV=0.00733707D0
     DOUBLE PRECISION, ALLOCATABLE, SAVE :: ZTOT(:),ZCURR(:),BETA(:),BETALEFT(:),BETARGHT(:),VDATA(:)
     DOUBLE PRECISION, ALLOCATABLE, SAVE :: PMAT12(:,:),PMAT6(:,:),PARRY12(:),PARRY6(:),EAVE2(:),ZTV2(:)
     DOUBLE PRECISION, ALLOCATABLE, SAVE :: PRESSURE(:),PRESRGHT(:),PRESLEFT(:)
     CHARACTER*100 :: CNST,FNAME,VFNAME,QFNAME,VOLFNAME,FRCFNAME,PWFRCFNAME,EVALFNAME
     LOGICAL :: PBC,UPDATED,EXTENDED,PGRID,NRECSNT,PTINIT,PCALC,NVT,NPT,FDP
END MODULE LJNPTMOD
  
SUBROUTINE CLEANUP_LJNPT
  USE LJNPTMOD
  DEALLOCATE(ZTOT,ZCURR,BETA,BETALEFT,BETARGHT,VDATA,PMAT12,PMAT6,PARRY12,PARRY6,PRESSURE,PRESRGHT,PRESLEFT,ZTV2,EAVE2)
  CLOSE(FILEIND)
  CLOSE(VFILEIND)
  CLOSE(QFILEIND)
  CLOSE(VOLFIND)
  CLOSE(FRCFIND)
  CLOSE(PWFRCFIND)
  CLOSE(EVALFIND)
END SUBROUTINE CLEANUP_LJNPT

SUBROUTINE INITIALIZE_LJNPTMOD(N,TMAX,TMIN,PMIN,PMAX,BLX,NP,GSIZE,NODE,ENT,EG,BLDIF,EFREQ)
  USE LJNPTMOD
  IMPLICIT NONE
  INTEGER :: I,J,K,N,NSETS,ITYPE,MSIZE1,MSIZE2,NG,CP,CNT,ACNT,NAM,IND,NM,LSZE
  INTEGER :: NP,GSIZE,PCNT,NODE,EG, EFREQ
  DOUBLE PRECISION :: TMAX,TMIN,BLX,BLY,BLZ,COE,AT,ALPHA,TI,TARRAY(NP),ENT,PMIN,PMAX,RCT,BLDIF
  CHARACTER*100 :: DUMMY

  N_ATOM=N
  NREPS=NP                                                  ! Number of parallel replicas
  GRIDSIZE=GSIZE
  NRECSNT=.TRUE.
  PTINIT=.TRUE.                                             ! Parallel tempering initialization
  mynode=NODE                                               ! node (0 = lowest)
  VMCONV=(ANGTOCM*NAVOG)/N_ATOM                             ! Molar volume conversion - NAVOG=Avogadros, ANGTOCM=Angstroms^3 to cm^3
  ENOT=ENT                                                  ! E_not for scaling energy
  EAVE=0.0D0
  PCALC=.FALSE.
  PCURR=0.0D0
  PAVE=0.0D0
  BLDF=BLDIF                                                ! Box length difference for finite diff pressure           
  EVFREQ=EFREQ                                              ! Frequency for computing eigenvalues

  COE12 = -4.0D0*48.0D0
  COE6 = 4.0D0*24.0D0

  IF(BLX.GT.0) THEN                                         ! Periodic boundary conditions ?
     PBC=.TRUE.
     BL=BLX
  ELSE
     PBC=.FALSE.
     BL=-1.0D0
     PMAX=0.0D0
  ENDIF
  
  ALLOCATE(PMAT12(N_ATOM,N_ATOM),PMAT6(N_ATOM,N_ATOM),PARRY12(N_ATOM),PARRY6(N_ATOM))       ! Data matrix for observables, see VGWFUPDATE
  
  CALL GRIDINIT(TMAX,TMIN,PMIN,PMAX,EG)
 
END SUBROUTINE INITIALIZE_LJNPTMOD

SUBROUTINE GRIDINIT(TMAX,TMIN,PMIN,PMAX,EG2)
  USE LJNPTMOD
  IMPLICIT NONE
  INTEGER :: I,NODE,EG2
  DOUBLE PRECISION :: TMAX,TMIN,PMIN,PMAX,CTE,T1,T2,TSTEP,RTEMP(NREPS+2)
  
  ARRYSIZE=GRIDSIZE+1
  NODE=mynode+1
  
  ALLOCATE(BETA(ARRYSIZE),BETALEFT(ARRYSIZE),BETARGHT(ARRYSIZE),ZTOT(ARRYSIZE),ZCURR(ARRYSIZE))
  ALLOCATE(VDATA(ARRYSIZE),PRESSURE(ARRYSIZE),PRESRGHT(ARRYSIZE),PRESLEFT(ARRYSIZE),ZTV2(ARRYSIZE),EAVE2(ARRYSIZE))                   ! Set-up tau grid (beta/2)
  
  IF(PMIN.GT.0) THEN
     PGRID=.TRUE.
     BETARGHT=1.0D0/TMIN
     BETA=BETARGHT
     BETALEFT=BETARGHT

     TSTEP=(PMAX-PMIN)/NREPS
     DO I=0,NREPS+1
        RTEMP(I+1)=PMIN+I*TSTEP
     ENDDO
     
     T1=RTEMP(NODE+2)
     T2=RTEMP(NODE+1)
     TSTEP=(T2-T1)/GRIDSIZE

     DO I=0,GRIDSIZE
        PRESRGHT(I+1)=T1+I*TSTEP
     ENDDO

     T1=T2
     T2=RTEMP(NODE)
     TSTEP=(T2-T1)/GRIDSIZE

     DO I=0,GRIDSIZE
        PRESSURE(I+1)=T1+I*TSTEP
     ENDDO

     IF(NODE.NE.1) THEN
        T1=RTEMP(NODE)
        T2=RTEMP(NODE-1)
        TSTEP=(T2-T1)/GRIDSIZE
        DO I=0,GRIDSIZE
           PRESLEFT(I+1)=T1+I*TSTEP               ! Beta for left node
        ENDDO
     ENDIF
  ELSE
     PGRID=.FALSE.
     IF (TMIN < 0.000001D0) TMIN=0.000001D0       ! to avoid devision by zero
     
     IF(EG2.EQ.1) THEN      
        TSTEP=(TMAX-TMIN)/NREPS
        DO I=0,NREPS+1
           RTEMP(I+1)=TMIN+I*TSTEP
        ENDDO
     ELSE
        CTE=(LOG(TMAX/TMIN))/(NREPS-1)
        CTE=EXP(CTE)
        DO I=0,NREPS+1
           RTEMP(I+1)=TMIN*CTE**I
        ENDDO
     ENDIF
     
     ZTOT=0.0D0
     ZCURR=0.0D0
     VDATA=0.0D0                                  ! Volume data
         
     T1=1.0D0/RTEMP(NODE+2)
     T2=1.0D0/RTEMP(NODE+1)
     TSTEP=(T2-T1)/GRIDSIZE
     
     DO I=0,GRIDSIZE
        BETARGHT(I+1)=T1+I*TSTEP                  ! Beta (inverse temperature) for right node
     ENDDO
     
     T1=T2
     T2=1.0D0/RTEMP(NODE)
     TSTEP=(T2-T1)/GRIDSIZE
     
     DO I=0,GRIDSIZE
        BETA(I+1)=T1+I*TSTEP                      ! Beta for current (MC) node 
     ENDDO
     
     IF(NODE.NE.1) THEN
        T1=1.0D0/RTEMP(NODE)
        T2=1.0D0/RTEMP(NODE-1)
        TSTEP=(T2-T1)/GRIDSIZE
        DO I=0,GRIDSIZE
           BETALEFT(I+1)=T1+I*TSTEP               ! Beta for left node
        ENDDO
     ENDIF
     
     PRESSURE(:)=PMAX              !*PCONV        ! Pressure array
     PRESRGHT(:)=PMAX
     PRESLEFT(:)=PMAX

  ENDIF
  
END SUBROUTINE GRIDINIT

SUBROUTINE LJNPTSETUP(Q,CONT)
  USE LJNPTMOD
  USE LJNPTMC
  IMPLICIT NONE
  INTEGER :: I,NODE,CONT,NLINES,LBLK,PBLK,LLINE,PLINE,EOF
  DOUBLE PRECISION :: Q(3,N_ATOM)
  CHARACTER*7 :: BLK
  CHARACTER*10 :: TMP
  CHARACTER*13 :: BUFF
 
  NODE=mynode+1
  FILEIND=100+NODE
  VFILEIND=800+NODE
  VOLFIND=900+NODE
  FRCFIND=1100+NODE
  PWFRCFIND=1200+NODE
  EVALFIND=1400+NODE

  NLINES=0
  EOF=0
  PBLK=0
  LLINE=0

  CALL SYSTEM("mkdir vmdcoords")
  CALL SYSTEM("mkdir latestcoords")
  CALL SYSTEM("mkdir vcurr")
  CALL SYSTEM("mkdir forces")
  CALL SYSTEM("mkdir pwforces")
  CALL SYSTEM("mkdir evals")

  IF(NODE.LT.10) THEN
     WRITE(FNAME,"(A6,I1)") "state.",NODE
     WRITE(VFNAME,"(A21,I1,A4)")"./vmdcoords/vmdcoords",NODE,".xyz"
     WRITE(QFNAME,"(A23,I1)") "./latestcoords/lcoords.", NODE
     WRITE(VOLFNAME,"(A16,I1)") "./vcurr/currvol.", NODE
     WRITE(FRCFNAME,"(A16,I1)") "./forces/forces.", NODE
     WRITE(PWFRCFNAME,"(A20,I1)") "./pwforces/pwforces.", NODE
     WRITE(EVALFNAME,"(A14,I1)") "./evals/evals.", NODE
  ELSE
     WRITE(FNAME,"(A6,I2)") "state.",NODE
     WRITE(VFNAME,"(A21,I2,A4)")"./vmdcoords/vmdcoords.",NODE,".xyz"
     WRITE(QFNAME,"(A23,I2)") "./latestcoords/lcoords.", NODE
     WRITE(VOLFNAME,"(A16,I2)") "./vcurr/currvol.", NODE
     WRITE(FRCFNAME,"(A16,I2)") "./forces/forces.", NODE
     WRITE(PWFRCFNAME,"(A20,I2)") "./pwforces/pwforces.", NODE
     WRITE(EVALFNAME,"(A14,I2)") "./evals/evals.", NODE
  ENDIF

  IF(CONT.EQ.1) THEN
     OPEN(UNIT=FILEIND,FILE=FNAME,STATUS='OLD')
     OPEN(UNIT=QFILEIND,FILE=QFNAME,STATUS='OLD')
     OPEN(UNIT=VOLFIND,FILE=VOLFNAME,STATUS='OLD')
     OPEN(UNIT=FRCFIND,FILE=FRCFNAME,STATUS='OLD')
     OPEN(UNIT=PWFRCFIND,FILE=PWFRCFNAME,STATUS='OLD')
     OPEN(UNIT=EVALFIND,FILE=EVALFNAME,STATUS='OLD')
     DO WHILE(.NOT.EOF)
        READ(FILEIND,*,IOSTAT=EOF) BLK,TMP
        IF(.NOT.EOF) THEN
           NLINES=NLINES+1
        ENDIF
        IF(BLK.EQ.'<BLOCK>') THEN
           PBLK=LBLK
           PLINE=LLINE
           BLCK=PBLK
           READ(TMP,*) LBLK
           LLINE=NLINES
        ENDIF
     ENDDO

     IF(LBLK.GT.1) THEN
        REWIND(FILEIND)
        DO I=1,PLINE
           READ(FILEIND,*)
        ENDDO
        DO I=1,N_ATOM
           READ(QFILEIND,*) Q(1,I),Q(2,I),Q(3,I)
        ENDDO
        READ(FILEIND,*) BUFF, NTOTAL
        READ(FILEIND,*) BUFF, NTACCPT
        READ(FILEIND,*) BUFF, NTRANS  
        READ(FILEIND,*)
        READ(FILEIND,*) BUFF, STPLEN
        READ(FILEIND,*) BUFF, BL
        READ(FILEIND,*) BUFF, VSLEN, NVACCT, NVSTEP
        READ(FILEIND,*) 
        READ(FILEIND,*) BUFF, ENOT
        READ(FILEIND,*)
        READ(FILEIND,*) 
        READ(FILEIND,*) BUFF, NSWAPS(1), NSWAPS(2)
        READ(FILEIND,*) BUFF, NSWPATT(1), NSWPATT(2)       
     ELSE
        WRITE(*,*) "Only one block. Did not read previous data."
     ENDIF

     NSTEP=NSTEP+NTOTAL

     REWIND(FILEIND)
     
     DO I=1,LLINE-1
        READ(FILEIND,*)
     ENDDO

  ELSE
     OPEN(UNIT=FILEIND,FILE=FNAME)
     OPEN(UNIT=QFILEIND,FILE=QFNAME)
     OPEN(UNIT=VOLFIND,FILE=VOLFNAME)
     OPEN(UNIT=FRCFIND,FILE=FRCFNAME)
     OPEN(UNIT=PWFRCFIND,FILE=PWFRCFNAME)
     OPEN(UNIT=EVALFIND,FILE=EVALFNAME)
  ENDIF

  NVT=.FALSE.
  NPT=.FALSE.
  CHNGV=.FALSE.
  PCALC=.FALSE.
  FDP=.FALSE.

  VCURR=BL**3

  IF((VSLEN.GT.0).AND.(PBC).AND.(VSFREQ.GT.0)) THEN
     CHNGV=.TRUE.
     NPT=.TRUE.
  ELSE
     PRESSURE=0.0D0
     PRESLEFT=0.0D0
     PRESRGHT=0.0D0
     IF(PBC.AND.(VSFREQ.GT.0)) THEN
        PCALC=.TRUE.
        NVT=.TRUE.
     ELSE
        BL=-1.0D0
        VCURR=1.0D0
     ENDIF
  ENDIF
  
  PVCURR=PRESSURE(ARRYSIZE)*VCURR
  
  Q=Q/ABS(BL)
  
  LS12=4.0D0/(BL**12)
  LS6=4.0D0/(BL**6)

  IF((BLDF.GT.0).AND.NVT) THEN
     FDP=.TRUE.
     BLDF=BLDF+BL
     VDIFF=(BLDF**3)-(BL**3)
     VFAC6=(BL/BLDF)**6
     VFAC12=(BL/BLDF)**12
     IF(mynode.eq.0) then 
          write(3455,*) VDIFF,BLDF
          call flush(3455)
     ENDIF
  ENDIF

  IF(ENOT.GT.0) CALL ALLENERGY(Q,ENOT)

  CALL PAIRSINIT(Q)

 ! IF(mynode.eq.0) THEN
 !    OPEN(UNIT=97,FILE='initialdata')
 !    WRITE(97,*) "Initial Energy: ", LJ12*LS12 + LJ6*LS6
 !    WRITE(97,*) "Inititial Box Length: ", BL
 !    WRITE(97,*) " "
 !    DO I=1,N_ATOM        
 !       WRITE(97,*) Q(1,I),Q(2,I),Q(3,I)
 !    ENDDO
 !    CLOSE(97)
 ! ENDIF

  OPEN(UNIT=VFILEIND,FILE=VFNAME)

END SUBROUTINE LJNPTSETUP
    
SUBROUTINE PENERGY(Q,IND)
  USE LJNPTMOD
  USE LJNPTMC
  IMPLICIT NONE
  DOUBLE PRECISION :: Q(3,N_ATOM),QIJ,RSQ,TENRG
  INTEGER :: I,J,IND
  
  LJNEW12=LJ12
  LJNEW6=LJ6

  DO I=1,N_ATOM
     LJNEW12=LJNEW12-PMAT12(IND,I)
     LJNEW6=LJNEW6-PMAT6(IND,I)
  ENDDO
  
  PARRY12(IND)=0.0D0
  PARRY6(IND)=0.0D0
  
  DO I=1,N_ATOM
     IF(I.NE.IND) THEN
        RSQ=0.0D0
        DO J=1,3
           QIJ=Q(J,IND)-Q(J,I)
           IF(PBC) THEN                                            ! Periodic boundary conditions
              IF(QIJ.GT.0.5D0) QIJ=QIJ-1.0D0
              IF(QIJ.LT.-0.5D0) QIJ=QIJ+1.0D0  
           ENDIF
           RSQ=RSQ+QIJ**2
        ENDDO
        RSQ=1/RSQ
        PARRY12(I)=LS12*(RSQ**6)      
        PARRY6(I)=-LS6*(RSQ**3)           
        LJNEW12=LJNEW12+PARRY12(I)         
        LJNEW6=LJNEW6+PARRY6(I)
     ENDIF
  ENDDO

  TENRG=LJNEW12 + LJNEW6 - ENOT
  ROENEW = -BETA(ARRYSIZE)*(TENRG + PRESSURE(ARRYSIZE)*VCURR) !PVCURR)
  
END SUBROUTINE PENERGY

SUBROUTINE CALCFORCES(Q,FRCS,FRCSR,FRCSA)
  USE LJNPTMOD
  USE LJNPTMC
  IMPLICIT NONE
  DOUBLE PRECISION :: Q(3,N_ATOM),FRCS(3,N_ATOM),RIJ(3),RSQ,GFACT
  DOUBLE PRECISION :: FRCSR(3,N_ATOM),FRCSA(3,N_ATOM),FREP,FATT
  INTEGER :: I,J,K

  FRCS=0.0D0
  FRCSR=0.0D0
  FRCSA=0.0D0

  DO I=1,N_ATOM-1
     DO J=I+1,N_ATOM
        RSQ=0.0D0
        DO K=1,3
           RIJ(K)=Q(K,I)-Q(K,J)
           IF(PBC) THEN
              IF(RIJ(K).GT.0.5D0) RIJ(K)=RIJ(K)-1.0D0
              IF(RIJ(K).LT.-0.5D0) RIJ(K)=RIJ(K)+1.0D0
           ENDIF
           RSQ=RSQ+(BL*RIJ(K))**2
        ENDDO
        FREP = COE12/RSQ**7
        FATT = COE6/RSQ**4
        GFACT= FATT + FREP                     
        DO K=1,3
           RSQ=GFACT*(BL*RIJ(K))
           FRCS(K,I)=FRCS(K,I)-RSQ
           FRCS(K,J)=FRCS(K,J)+RSQ

           RSQ=FREP*(BL*RIJ(K))
           FRCSR(K,I)=FRCSR(K,I)-RSQ
           FRCSR(K,J)=FRCSR(K,J)+RSQ

           RSQ=FATT*(BL*RIJ(K))
           FRCSA(K,I)=FRCSA(K,I)-RSQ
           FRCSA(K,J)=FRCSA(K,J)+RSQ
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE CALCFORCES

SUBROUTINE CALCFORCESPW(Q,FRCS,NPAIRS,RCUT)
  USE LJNPTMOD
  USE LJNPTMC
  IMPLICIT NONE
  DOUBLE PRECISION :: Q(3,N_ATOM),FRCS(3,(N_ATOM*N_ATOM)/2),RIJ(3)
  DOUBLE PRECISION :: RSQ,GFACT,RCUT,RC2,FREP,FATT
  INTEGER :: I,J,K, NPAIRS

  FRCS=0.0D0
 
  NPAIRS = 0

  RC2 = RCUT**2

  DO I=1,N_ATOM-1
     DO J=I+1,N_ATOM
        RSQ=0.0D0
        DO K=1,3
           RIJ(K)=Q(K,I)-Q(K,J)
           IF(PBC) THEN
              IF(RIJ(K).GT.0.5D0) RIJ(K)=RIJ(K)-1.0D0
              IF(RIJ(K).LT.-0.5D0) RIJ(K)=RIJ(K)+1.0D0
           ENDIF
           RSQ=RSQ+(BL*RIJ(K))**2
        ENDDO
        IF(RSQ.LE.RC2) THEN
           NPAIRS = NPAIRS + 1
           FREP = COE12/RSQ**7
           FATT = COE6/RSQ**4
           GFACT = FATT + FREP     

           DO K=1,3
              RSQ=GFACT*(BL*RIJ(K))
              FRCS(1,NPAIRS) = FRCS(1,NPAIRS) + RSQ*RSQ
              
              RSQ=FREP*(BL*RIJ(K))
              FRCS(2,NPAIRS) = FRCS(2,NPAIRS) + RSQ*RSQ
              
              RSQ=FATT*(BL*RIJ(K))
              FRCS(3,NPAIRS) = FRCS(3,NPAIRS) + RSQ*RSQ
           ENDDO

           DO K=1,3
              FRCS(K,NPAIRS) = SQRT(FRCS(K,NPAIRS))
           ENDDO
              
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE CALCFORCESPW

SUBROUTINE PRCALC(Q)
  USE LJNPTMOD
  USE LJNPTMC
  IMPLICIT NONE
  DOUBLE PRECISION :: Q(3,N_ATOM),GRAD(3,N_ATOM),RIJ(3),RSQ,GFACT
  DOUBLE PRECISION :: EC,EP
  INTEGER :: I,J,K

  !IF(FDP) THEN
   !  EP=LJ12*VFAC12+LJ6*VFAC6
    ! EC=LJ12+LJ6
    ! PCURR=(EP-EC)/VDIFF    
  !ELSE    
 
     GRAD=0.0D0
     
     DO I=1,N_ATOM-1
        DO J=I+1,N_ATOM
           RSQ=0.0D0
           DO K=1,3
              RIJ(K)=Q(K,I)-Q(K,J)
              IF(PBC) THEN
                 IF(RIJ(K).GT.0.5D0) RIJ(K)=RIJ(K)-1.0D0
                 IF(RIJ(K).LT.-0.5D0) RIJ(K)=RIJ(K)+1.0D0
              ENDIF
              RSQ=RSQ+(BL*RIJ(K))**2
           ENDDO
           GFACT=24.0D0/RSQ**4 - 48.0D0/RSQ**7
           DO K=1,3
              RSQ=GFACT*(BL*RIJ(K))
              GRAD(K,I)=GRAD(K,I)-RSQ
              GRAD(K,J)=GRAD(K,J)+RSQ
           ENDDO
        ENDDO
     ENDDO
     
     GFACT=0.0D0
     
     DO I=1,N_ATOM
        DO J=1,3
           GFACT=GFACT+GRAD(J,I)*(BL*Q(J,I))
        ENDDO
     ENDDO
     
     PCURR=(1/BL**3)*(N_ATOM/BETA(ARRYSIZE)+GFACT/3.0D0)
  
END SUBROUTINE PRCALC

SUBROUTINE HESSIAN(Q)
  USE LJNPTMOD
  USE QUENCHMOD
  IMPLICIT NONE
  INTEGER, PARAMETER :: MIT=1500
  DOUBLE PRECISION :: Q(3,N_ATOM),QTMP(3,N_ATOM),DDU(3*N_ATOM,3*N_ATOM),RSQ,RSQ4,RSQ5,RSQ7,RSQ8
  DOUBLE PRECISION :: QIJ(3),EVALS(3*N_ATOM),WORK(9*N_ATOM-1),COE1,COE2,COE3,COE4,DFACT,GRD
  INTEGER I,J,K,L,LWORK,INFO,NIT(1)
  CHARACTER*1 :: JZ, UL

  DO I=1,N_ATOM
     DO J=1,3
        QTMP(J,I) = Q(J,I)
     ENDDO
  ENDDO

  GRD = QUENCHF(QTMP,N_ATOM,MIT,NIT)

  write(mynode+601,*) GRD, NIT(1)

  JZ = 'V'
  UL = 'U'
  LWORK = 9*N_ATOM-1

  DDU = 0.0D0

  COE1 = 24.0D0
  COE2 = 48.0D0
  COE3 = 192.0D0
  COE4 = 672.0D0

  DO I=1,N_ATOM-1
    DO J=I+1,N_ATOM
      RSQ=0.0D0
      DO K=1,3
        QIJ(K) = QTMP(K,I)-QTMP(K,J)
        RSQ = RSQ + QIJ(K)**2
      ENDDO
      
      RSQ4 = RSQ**4
      RSQ8 = RSQ4*RSQ4
      RSQ5 = RSQ4*RSQ
      RSQ7 = RSQ8/RSQ
     
      DO K=1,3
        DFACT =  COE1/RSQ4 - COE2/RSQ7 - (QIJ(K)**2)*(COE3/RSQ5 - COE4/RSQ8)             ! Diagonal  elements
        DDU(3*(I-1)+K,3*(I-1)+K) = DDU(3*(I-1)+K,3*(I-1)+K) + DFACT
        DDU(3*(J-1)+K,3*(J-1)+K) = DDU(3*(J-1)+K,3*(J-1)+K) + DFACT
        DDU(3*(J-1)+K,3*(I-1)+K) = -1.0D0*DFACT
        DDU(3*(I-1)+K,3*(J-1)+K) = DDU(3*(J-1)+K,3*(I-1)+K)                            ! Transpose
        DO L=K+1,3
           DFACT = (COE4/RSQ8 - COE3/RSQ5)*(QIJ(K)*QIJ(L))                             ! Off diagonal elements
           DDU(3*(I-1)+L,3*(I-1)+K) = DDU(3*(I-1)+L,3*(I-1)+K) + DFACT                 
           DDU(3*(J-1)+L,3*(J-1)+K) = DDU(3*(J-1)+L,3*(J-1)+K) + DFACT    
           DDU(3*(I-1)+K,3*(I-1)+L) = DDU(3*(I-1)+K,3*(I-1)+L) + DFACT                 ! Self transpose
           DDU(3*(J-1)+K,3*(J-1)+L) = DDU(3*(J-1)+K,3*(J-1)+L) + DFACT
           DDU(3*(J-1)+L,3*(I-1)+K) = -1.0* DFACT                                      ! Jth element
           DDU(3*(J-1)+K,3*(I-1)+L) = -1.0* DFACT  
           DDU(3*(I-1)+L,3*(J-1)+K) = DDU(3*(J-1)+L,3*(I-1)+K)                         ! Transpose
           DDU(3*(I-1)+K,3*(J-1)+L) = DDU(3*(J-1)+K,3*(I-1)+L)                         ! Transpose
        ENDDO
     ENDDO
  ENDDO
ENDDO

!WRITE(444+mynode,*) " "

!DO I=1,N_ATOM
!   DO J=1,N_ATOM
!      WRITE(444+mynode,"(F12.5)") DDU(3*J,3*I)
!   ENDDO
!   WRITE(444+mynode,*) " "
!ENDDO

!WRITE(444+mynode,*) "***********"
!WRITE(444+mynode,*) " "

!CALL FLUSH(11+mynode)

CALL DSYEV(JZ,UL,3*N_ATOM,DDU,3*N_ATOM,EVALS,WORK,LWORK,INFO)

DO I=7,3*N_ATOM
   IF((ABS(EVALS(I)).GE.EVTOL).AND.(EVALS(I).NE.'NaN')) THEN
      WRITE(EVALFIND,"(F12.5)") SQRT(EVALS(I))
   ENDIF
ENDDO

!WRITE(EVALFIND, *) ""

CALL FLUSH(EVALFIND)

END SUBROUTINE HESSIAN

SUBROUTINE PAIRSINIT(QC)
  USE LJNPTMOD
  IMPLICIT NONE
  DOUBLE PRECISION :: QC(3,N_ATOM),RSQ,QIJ,ECURR
  INTEGER :: I,J,K
  
  LJ12=0.0D0
  LJ6=0.0D0
  PMAT12=0.0D0
  PMAT6=0.0D0

  DO I=1,N_ATOM-1
     DO J=I+1,N_ATOM
        RSQ=0.0D0
        DO K=1,3
           QIJ=QC(K,I)-QC(K,J)
           IF(PBC) THEN                                            ! Periodic boundary conditions
              IF(QIJ.GT.0.5D0) QIJ=QIJ-1.0D0
              IF(QIJ.LT.-0.5D0) QIJ=QIJ+1.0D0         
           ENDIF
           RSQ=RSQ+QIJ**2
        ENDDO
        RSQ=1/RSQ
        PMAT12(I,J)=LS12*(RSQ**6)
        PMAT6(I,J)=-LS6*(RSQ**3)
        PMAT12(J,I)=PMAT12(I,J)
        PMAT6(J,I)=PMAT6(I,J)
        LJ12=LJ12+PMAT12(I,J)
        LJ6=LJ6+PMAT6(I,J)
     ENDDO
  ENDDO

  ECURR=LJ12 + LJ6 - ENOT

  DO I=1,ARRYSIZE
     ZCURR(I)=-BETA(I)*(ECURR + PRESSURE(I)*VCURR)
  ENDDO
 
  ROECURR=ZCURR(ARRYSIZE)
   
 END SUBROUTINE PAIRSINIT

SUBROUTINE ALLENERGY(Q,CENRG)
  USE LJNPTMOD
  IMPLICIT NONE
  DOUBLE PRECISION :: Q(3,N_ATOM),RSQ,U12,U6,CENRG,QIJ       
  INTEGER :: I,J,K

  U12=0.0D0
  U6=0.0D0

  DO I=1,N_ATOM-1
     DO J=I+1,N_ATOM
        RSQ=0.0D0
        DO K=1,3
           QIJ=Q(K,I)-Q(K,J)
           IF(PBC) THEN                                          ! Periodic boundary conditions
               IF(QIJ.GT.0.5D0) QIJ=QIJ-1.0D0
               IF(QIJ.LT.-0.5D0) QIJ=QIJ+1.0D0          
            ENDIF
           RSQ=RSQ+QIJ**2
        ENDDO
        RSQ=1/RSQ
        U12=U12+RSQ**6
        U6=U6-RSQ**3
     ENDDO
  ENDDO   

  CENRG=U12*LS12 + U6*LS6

END SUBROUTINE ALLENERGY


