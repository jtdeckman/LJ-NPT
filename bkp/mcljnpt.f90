MODULE LJNPTMC
  IMPLICIT NONE
  INTEGER*8, SAVE :: NUPDATES,NEQUIL,NSTEP,NTOTAL,NACCPT,NTRANS
  INTEGER*8, SAVE :: NTACCPT,NACC,NVSTEP,NVACC,NVACCT,NSWPATT(2),NSWAPS(2)
  INTEGER, SAVE :: UFREQ,BLCK,VSFREQ,REFRESH,MCUF,WFREQ,NVOUT,SWPINT(2),CSTEP(2)
  DOUBLE PRECISION, SAVE :: STPLEN,VSLEN,RCON,MAXBL
  LOGICAL, SAVE :: CHNGV

  CONTAINS
     FUNCTION RANF (DUMMY)
       INTEGER  ::   L, C, M
       PARAMETER (L = 1029,C = 221591,M = 1048576)
       DOUBLE PRECISION :: RANF 
       INTEGER ::  SEED
       DOUBLE PRECISION :: DUMMY
       SAVE        SEED
       DATA        SEED / 0 /
       
       SEED = MOD (SEED*L + C, M)
       RANF = REAL (SEED)/M
       
       RETURN
     END FUNCTION RANF

END MODULE LJNPTMC

SUBROUTINE INITIALIZE_LJNPTMC(NSTP,UFR,WRF,MCF,SI,SLEN,VLEN,RCN,Q,CONT,N,VSF,RFRSH,MBL,VOF)
  USE LJNPTMC
  IMPLICIT NONE
  INTEGER*8 :: NSTP,NE
  INTEGER :: UFR,WRF,MCF,SI,CONT,N,VSF,VOF
  DOUBLE PRECISION :: SLEN,SFAC,RCN,VLEN,RFRSH,Q(3,N),MBL

  NSTEP=NSTP                    ! Number of MC steps to do
  WFREQ=WRF                     ! Observable write frequency (every n updates)
  UFREQ=UFR                     ! Observable update frequency
  MCUF=MCF                      ! Scale step lengths every MCUF MC steps
  SWPINT=SI                     ! PT swap interval
  STPLEN=SLEN                   ! Step length for translational steps
  VSLEN=VLEN                    ! Step length for volume changes
  VSFREQ=VSF                    ! Attempt volume step frequency
  RCON=RCN                      ! Constraining radius for non-PBC

  NACCPT=0                      ! Number of accepted MC steps
  NTOTAL=0                      ! Number of total steps
  NTRANS=0                      ! Number of translational steps
  NVSTEP=0                      ! Number of volume change steps
  NSWPATT=0                     ! Number of PT swap attempts
  NSWAPS=0                      ! Number of swaps made (1 = left node, 2=right node)
  CSTEP=0                       ! Current step (used for PT swapping; swap every CSTEP=SWPINT) 
  BLCK=0                        ! Current block of observable data
  NACC=0                        ! Number of accepted translational steps
  NVACC=0                       ! Number of accepted volume change steps
  NVACCT=0                      ! Total accepted volume steps
  NTACCPT=0                     ! Accepted translational steps
  CHNGV=.FALSE.                 ! Change volume (do volume steps)
  REFRESH=RFRSH                 ! Energy refresh frequency
  MAXBL=MBL                     ! Maximum box length
  NVOUT=VOF                     ! Volume output frequency

  CALL LJNPTSETUP(Q,CONT)

END SUBROUTINE INITIALIZE_LJNPTMC

SUBROUTINE MCSTEP(Q,VWRF,TAG)
  USE LJNPTMOD
  USE LJNPTMC
  IMPLICIT NONE
  INTEGER :: I,J,STYPE,VWRF
  DOUBLE PRECISION :: Q(3,N_ATOM),QVEC(3)
  LOGICAL :: INSIDE,ATEST,TFLAG,UD
  CHARACTER*2 :: TAG
  
  TFLAG=.FALSE.

  CALL RANDOM_SEED    
  
  DO WHILE(.NOT.TFLAG) 
     
     CALL TSTEP(Q)
     
     IF((MOD(NTRANS,VSFREQ).EQ.0).AND.CHNGV) CALL VSTEP(Q)
     
    ! IF((MOD(NTRANS,VSFREQ).EQ.0).AND.PCALC) CALL PRCALC(Q)

     IF(MOD(NTOTAL,REFRESH).EQ.0) CALL PAIRSINIT(Q)
     
     CALL ATTSWAP(Q,TFLAG)
     
     IF((MOD(NTOTAL,VWRF).EQ.0).AND.(VWRF.GT.0)) CALL VMDCOORDS(Q,TAG)
     
  ENDDO

END SUBROUTINE MCSTEP

SUBROUTINE TSTEP(Q)
  USE LJNPTMOD
  USE LJNPTMC
  IMPLICIT NONE
  INTEGER :: I,J,IND
  DOUBLE PRECISION :: Q(3,N_ATOM),QT(3,N_ATOM),RNUM,ACCR,ECURR,DMY
  LOGICAL :: INSIDE

  CALL RANDOM_NUMBER(RNUM)
 ! RNUM=RANF(DMY)
 
  IND=INT(RNUM*N_ATOM)+1 
  IF(IND.GE.N_ATOM) IND=N_ATOM
  
  QT=Q

  NTRANS=NTRANS+1
  IF(mynode.GT.0) CSTEP(1)=CSTEP(1)+1
  IF(mynode.LT.(NREPS-1)) CSTEP(2)=CSTEP(2)+1

  INSIDE=.TRUE.
 
  DO J=1,3
     CALL RANDOM_NUMBER(RNUM)
     !RNUM=RANF(DMY)
     QT(J,IND)=Q(J,IND)+STPLEN*(RNUM-0.5D0)
     IF(PBC) THEN
        IF(QT(J,IND).GE.(0.5D0)) QT(J,IND)=QT(J,IND)-1.0D0
        IF(QT(J,IND).LE.(-0.5D0)) QT(J,IND)=QT(J,IND)+1.0D0
        IF(ABS(QT(J,IND)).GT.(0.5D0)) INSIDE=.FALSE.
     ENDIF
  ENDDO

  IF(.NOT.PBC) CALL RCENTER(QT,N_ATOM,RCON,INSIDE)

  IF(INSIDE) THEN
     CALL PENERGY(QT,IND)
     CALL RANDOM_NUMBER(RNUM)
     IF((RNUM < EXP(ROENEW-ROECURR))) THEN
        NACC=NACC+1
        NTACCPT=NTACCPT+1
        ROECURR=ROENEW                             ! Set current density to new density (roe)
        LJ12=LJNEW12
        LJ6=LJNEW6
        Q=QT
        DO I=1,N_ATOM
           PMAT12(IND,I)=PARRY12(I)
           PMAT12(I,IND)=PARRY12(I)
           PMAT6(IND,I)=PARRY6(I)
           PMAT6(I,IND)=PARRY6(I)
        ENDDO
        ECURR=LJ12+LJ6-ENOT
        DO I=1,ARRYSIZE
           ZCURR(I)=-BETA(I)*(ECURR + PRESSURE(I)*VCURR)
        ENDDO 
     ENDIF             
  ENDIF

  CALL UPDATE(Q)

  IF(MOD(NTRANS,MCUF)==0) THEN
     ACCR=FLOAT(NACC)/FLOAT(MCUF)
     NACC=0
     IF(STPLEN.GE.(0.30D0)) THEN
      !  STPLEN=STPLEN/1.10779652734191D0
         STPLEN=0.30D0
     ELSE 
        IF(ACCR.GT.0.4D0) STPLEN=STPLEN*1.12799165273419D0
        IF(ACCR.LT.0.3D0) STPLEN=STPLEN/1.10779652734191D0
     ENDIF
  ENDIF
  
END SUBROUTINE TSTEP

SUBROUTINE VSTEP(Q)
  USE LJNPTMOD
  USE LJNPTMC
  IMPLICIT NONE
  INTEGER :: I,J,VEC
  DOUBLE PRECISION :: BLTMP,Q(3,N_ATOM),ENEW,RNUM,ACCR,LSTMP12,LSTMP6,PVNEW,FAC12,FAC6,LJN12,LJN6,DMY

  NVSTEP=NVSTEP+1

  CALL RANDOM_NUMBER(RNUM)
!  RNUM=RANF(DMY)
  BLTMP=BL+VSLEN*(RNUM-0.5D0)
 ! VNEW=VCURR+VSLEN*(RNUM-0.5D0)
 ! BLTMP=VNEW**(1.0D0/3.0D0)
  IF(BLTMP.GT.MAXBL) THEN 
     BLTMP=MAXBL
 !    VNEW=BLTMP**3
  ENDIF

  VNEW=BLTMP**3
  PVNEW=PRESSURE(ARRYSIZE)*VNEW

  FAC12=(BL/BLTMP)**12
  FAC6=(BL/BLTMP)**6

  LJN12=LJ12*FAC12
  LJN6=LJ6*FAC6

  ENEW=LJN12+LJN6-ENOT
  ROENEW=-BETA(ARRYSIZE)*(ENEW + PRESSURE(ARRYSIZE)*VNEW)

  CALL RANDOM_NUMBER(RNUM)
 ! RNUM=RANF(DMY)

  IF(RNUM < EXP(ROENEW - ROECURR + N_ATOM*LOG(VNEW/VCURR) + 2.0D0*LOG(BLTMP/BL)).AND.(BLTMP.GT.0)) THEN 
     NVACC=NVACC+1
     NVACCT=NVACCT+1
     ROECURR=ROENEW
     LJ12=LJN12
     LJ6=LJN6
     PMAT12=PMAT12*FAC12
     PMAT6=PMAT6*FAC6
     BL=BLTMP
     LS12=4.0D0/(BL**12)
     LS6=4.0D0/(BL**6)
     VCURR=VNEW
     PVCURR=PVNEW
     DO I=1,ARRYSIZE
        ZCURR(I)=-BETA(I)*(ENEW + PRESSURE(I)*VNEW)
     ENDDO
  ENDIF

  CALL UPDATE(Q)

  IF(MOD(NVSTEP,MCUF).EQ.0) THEN
     ACCR=FLOAT(NVACC)/FLOAT(MCUF)
     NVACC=0
  !   IF(VSLEN.GE.(0.5D0)) THEN
     IF(VSLEN.GE.(BL/2.0D0)) THEN
        VSLEN=VSLEN/1.10779652734191D0
     ELSE
        IF(ACCR.GT.0.4D0) VSLEN=VSLEN*1.12799165273419D0
        IF(ACCR.LT.0.3D0) VSLEN=VSLEN/1.10779652734191D0
     ENDIF
  ENDIF

END SUBROUTINE VSTEP

SUBROUTINE RCENTER(Q,N,RC,INSIDE)
  IMPLICIT NONE
  INTEGER :: I,J,N,RN
  DOUBLE PRECISION :: RC,Q(3,N),QC(3),RMAX,RCURR
  LOGICAL :: INSIDE

   DO I=1,N                                  ! Determine center of mass vector
      DO J=1,3
         QC(J)=QC(J)+Q(J,I)                    
      ENDDO
   ENDDO
    
   QC=QC/N
  
   DO I=1,N
      DO J=1,3
         Q(J,I)=Q(J,I)-QC(J)                 ! Zero center of mass
      ENDDO
   ENDDO

   RMAX=-1.0D0

   DO I=1,N
      RCURR=0.0D0
      DO J=1,3
         RCURR=RCURR+Q(J,I)**2
      ENDDO
      IF(RCURR.GT.RMAX) RMAX=RCURR
   ENDDO
   
   IF(RMAX.GT.RC**2) THEN
      INSIDE=.FALSE.
   ELSE
      INSIDE=.TRUE.
   ENDIF

 END SUBROUTINE RCENTER

 SUBROUTINE PBCHECK(Q,N,BL,TVAL)
   IMPLICIT NONE
   INTEGER :: I,J,N
   DOUBLE PRECISION :: Q(3,N),BL,BL2
   LOGICAL :: TVAL
   
   TVAL=.TRUE.

   BL2=BL/2.0D0

   DO I=1,N
      DO J=1,3
         IF(Q(J,I).GE.BL2) Q(J,I) = Q(J,I) - BL
         IF(Q(J,I).LE.-BL2) Q(J,I) = Q(J,I) + BL
         IF(ABS(Q(J,I).GT.BL2)) TVAL=.FALSE.                           ! If its still outside the box then reject config
      ENDDO
   ENDDO
   
 END SUBROUTINE PBCHECK


