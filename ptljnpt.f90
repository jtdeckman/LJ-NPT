SUBROUTINE PTSWAP(Q,SENDER)
  USE LJNPTMOD
  USE LJNPTMC
  IMPLICIT NONE
  INCLUDE "mpif.h"
  INTEGER :: I,EXCH,ISTAT(MPI_STATUS_SIZE),ierr
  DOUBLE PRECISION :: Q(3*N_ATOM),QTMP(3*N_ATOM),LNP(ARRYSIZE+4),SDATA(ARRYSIZE+4)
  DOUBLE PRECISION :: PMTEMP12(N_ATOM,N_ATOM),PMTEMP6(N_ATOM,N_ATOM),RNUM,W,ECURR,VTMP
  LOGICAL :: SENDER,UD,NR,TFLAG

  ECURR=LJ12+LJ6-ENOT

  IF(SENDER) THEN
     DO I=1,ARRYSIZE
         LNP(I)=-BETARGHT(I)*(ECURR + PRESRGHT(I)*VCURR)
     ENDDO
     CALL MPI_RECV(SDATA(1:ARRYSIZE+4),ARRYSIZE+4,MPI_DOUBLE_PRECISION,mynode+1,7,MPI_COMM_WORLD,ISTAT,ierr)   
     CALL RANDOM_NUMBER(RNUM)
     NSWPATT(2)=NSWPATT(2)+1
     W=EXP(LNP(ARRYSIZE)+SDATA(ARRYSIZE)-ROECURR-SDATA(ARRYSIZE+1))                  ! Weight
     IF(W.GT.RNUM) THEN    
        EXCH=1
        CALL MPI_SEND(EXCH,1,MPI_INTEGER,mynode+1,8,MPI_COMM_WORLD,ierr)
        IF(CHNGV) THEN                                                        ! For volume changes (npt ensemble)
           LNP(ARRYSIZE+2)=BL
           BL=SDATA(ARRYSIZE+2)
           LS12=4.0D0/(BL**12)
           LS6=4.0D0/(BL**6)
           VCURR=BL**3
           PVCURR=PRESSURE(ARRYSIZE)*VCURR
        ENDIF
        LNP(ARRYSIZE+3)=LJ12
        LNP(ARRYSIZE+4)=LJ6     
        CALL MPI_SEND(LNP(1:ARRYSIZE+4),ARRYSIZE+4,MPI_DOUBLE_PRECISION,mynode+1,9,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(Q(1:3*N_ATOM),3*N_ATOM,MPI_DOUBLE_PRECISION,mynode+1,10,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(PMAT12(1:N_ATOM,1:N_ATOM),N_ATOM**2,MPI_DOUBLE_PRECISION,mynode+1,11,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(PMAT6(1:N_ATOM,1:N_ATOM),N_ATOM**2,MPI_DOUBLE_PRECISION,mynode+1,12,MPI_COMM_WORLD,ierr)
        CALL MPI_RECV(Q(1:3*N_ATOM),3*N_ATOM,MPI_DOUBLE_PRECISION,mynode+1,13,MPI_COMM_WORLD,ISTAT,ierr)
        CALL MPI_RECV(PMAT12(1:N_ATOM,1:N_ATOM),N_ATOM**2,MPI_DOUBLE_PRECISION,mynode+1,14,MPI_COMM_WORLD,ISTAT,ierr)
        CALL MPI_RECV(PMAT6(1:N_ATOM,1:N_ATOM),N_ATOM**2,MPI_DOUBLE_PRECISION,mynode+1,15,MPI_COMM_WORLD,ISTAT,ierr)
        LJ12=SDATA(ARRYSIZE+3)
        LJ6=SDATA(ARRYSIZE+4)
        ZCURR(1:ARRYSIZE)=SDATA(1:ARRYSIZE)
        ROECURR=ZCURR(ARRYSIZE)  
        NSWAPS(2)=NSWAPS(2)+1    
     ELSE
        EXCH=0
        CALL MPI_SEND(EXCH,1,MPI_INTEGER,mynode+1,8,MPI_COMM_WORLD,ierr)
     ENDIF
  ELSE
     DO I=1,ARRYSIZE
        LNP(I)=-BETALEFT(I)*(ECURR + PRESLEFT(I)*VCURR)
     ENDDO
     LNP(ARRYSIZE+1)=ROECURR
     LNP(ARRYSIZE+2)=BL
     LNP(ARRYSIZE+3)=LJ12
     LNP(ARRYSIZE+4)=LJ6
     CALL MPI_SEND(LNP(1:ARRYSIZE+4),ARRYSIZE+4,MPI_DOUBLE_PRECISION,mynode-1,7,MPI_COMM_WORLD,ierr)
     CALL MPI_RECV(EXCH,1,MPI_INTEGER,mynode-1,8,MPI_COMM_WORLD,ISTAT,ierr)
     NSWPATT(1)=NSWPATT(1)+1
     IF(EXCH.EQ.1) THEN    
        CALL MPI_RECV(LNP(1:ARRYSIZE+4),ARRYSIZE+4,MPI_DOUBLE_PRECISION,mynode-1,9,MPI_COMM_WORLD,ISTAT,ierr)
        CALL MPI_RECV(QTMP(1:3*N_ATOM),3*N_ATOM,MPI_DOUBLE_PRECISION,mynode-1,10,MPI_COMM_WORLD,ISTAT,ierr)
        CALL MPI_RECV(PMTEMP12(1:N_ATOM,1:N_ATOM),N_ATOM**2,MPI_DOUBLE_PRECISION,mynode-1,11,MPI_COMM_WORLD,ISTAT,ierr)
        CALL MPI_RECV(PMTEMP6(1:N_ATOM,1:N_ATOM),N_ATOM**2,MPI_DOUBLE_PRECISION,mynode-1,12,MPI_COMM_WORLD,ISTAT,ierr)
        CALL MPI_SEND(Q(1:3*N_ATOM),3*N_ATOM,MPI_DOUBLE_PRECISION,mynode-1,13,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(PMAT12(1:N_ATOM,1:N_ATOM),N_ATOM**2,MPI_DOUBLE_PRECISION,mynode-1,14,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(PMAT6(1:N_ATOM,1:N_ATOM),N_ATOM**2,MPI_DOUBLE_PRECISION,mynode-1,15,MPI_COMM_WORLD,ierr)
        IF(CHNGV) THEN
           BL=LNP(ARRYSIZE+2)
           LS12=4.0D0/(BL**12)
           LS6=4.0D0/(BL**6)
           VCURR=BL**3
           PVCURR=PRESSURE(ARRYSIZE)*VCURR
        ENDIF        
        Q=QTMP
        PMAT12=PMTEMP12
        PMAT6=PMTEMP6
        LJ12=LNP(ARRYSIZE+3)
        LJ6=LNP(ARRYSIZE+4)
        ZCURR(1:ARRYSIZE)=LNP(1:ARRYSIZE)
        ROECURR=ZCURR(ARRYSIZE)
        NSWAPS(1)=NSWAPS(1)+1
     ENDIF
  ENDIF

  CALL UPDATE(Q)

END SUBROUTINE PTSWAP

SUBROUTINE ATTSWAP(Q,TFLAG)
  USE LJNPTMOD
  USE LJNPTMC
  IMPLICIT NONE
  INCLUDE "mpif.h"
  INTEGER :: I
  DOUBLE PRECISION :: Q(3*N_ATOM)
  LOGICAL :: TFLAG
  LOGICAL, SAVE :: SENDNR,VSENDNR
  INTEGER, SAVE :: STATUS(MPI_STATUS_SIZE,100),TSTAT(MPI_STATUS_SIZE)
  INTEGER, SAVE :: REQ(100),REQT(100),REQV(100),TBUF,SBUF,RBUF,VBUF,ierr

  IF(PTINIT) THEN
     TBUF=0
     SBUF=0
     VBUF=0
     TFLAG=.FALSE.
     SENDNR=.TRUE.
     VSENDNR=.TRUE.
     PTINIT=.FALSE.
     IF(mynode.GT.0) CALL MPI_IRECV(TBUF,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,REQ(1),ierr)
  ENDIF
  
  IF(NVOUT.GT.0) THEN
     IF(mynode.EQ.0) THEN
        IF(MOD(NTOTAL,NVOUT).EQ.0) THEN
           DO I=1,NREPS-1
              CALL MPI_ISEND(VBUF,1,MPI_INTEGER,I,2,MPI_COMM_WORLD,REQV(I),ierr)   ! Send terminate signal
           ENDDO
           CALL MPI_WAITALL(NREPS-1,REQV,STATUS,ierr)
           WRITE(VOLFIND,"(I15,F15.5)") NTOTAL,BL
           CALL FLUSH(VOLFIND)
        ENDIF
     ELSE
        IF(VSENDNR) THEN
           CALL MPI_IRECV(VBUF,1,MPI_INTEGER,0,2,MPI_COMM_WORLD,REQ(2),ierr)      
           VSENDNR=.FALSE.
        ENDIF
        
        CALL MPI_TEST(REQ(2),VSENDNR,TSTAT,ierr)
     
        IF(VSENDNR) THEN
           WRITE(VOLFIND,"(I15,F15.5)") NTOTAL,BL   ! Output current volume of node to file
           CALL FLUSH(VOLFIND)
        ENDIF
     ENDIF
  ENDIF

  IF(mynode.EQ.0) THEN
    IF(NTOTAL.GE.NSTEP) THEN
        DO I=1,NREPS-1
           CALL MPI_ISEND(TBUF,1,MPI_INTEGER,I,1,MPI_COMM_WORLD,REQT(I),ierr)             ! Send terminate signal
        ENDDO
        CALL MPI_WAITALL(NREPS-1,REQT,STATUS,ierr)
        TFLAG=.TRUE.
     ENDIF
  ELSE             
     CALL MPI_TEST(REQ(1),TFLAG,STATUS,ierr)
      
     IF((CSTEP(1).GE.SWPINT(1)).AND.(.NOT.TFLAG)) THEN        
  
       ! IF(SENDNR) THEN                                                                ! Send new receive
       !    CALL MPI_IRECV(RBUF,1,MPI_INTEGER,mynode-1,3,MPI_COMM_WORLD,REQ(3),ierr)
       !    SENDNR=.FALSE.
       ! ENDIF

       ! CALL MPI_TEST(REQ(3),SENDNR,TSTAT,ierr)

       ! IF(SENDNR) THEN 
           CALL PTSWAP(Q,.FALSE.)
           CSTEP(1)=0
       ! ENDIF
     ENDIF
  ENDIF
           
  IF((CSTEP(2).GE.SWPINT(2)).AND.(.NOT.TFLAG).AND.(mynode.LT.(NREPS-1))) THEN
  !   CALL MPI_ISEND(SBUF,1,MPI_INTEGER,mynode+1,3,MPI_COMM_WORLD,REQ(4),ierr)  
  !   CALL MPI_WAIT(REQ(4),TSTAT,ierr)
     CALL PTSWAP(Q,.TRUE.)
     CSTEP(2)=0              
  ENDIF

END SUBROUTINE ATTSWAP

!SUBROUTINE ATTSWAPOLD(Q,ENDSIGNAL)
!  USE LJNPTMOD
!  USE LJNPTMC
!  IMPLICIT NONE
!  INCLUDE "mpif.h"
!  INTEGER :: I,mynode
!  DOUBLE PRECISION :: Q(3*N_ATOM)
!  LOGICAL :: ENDSIGNAL
!  LOGICAL, SAVE :: SENDNR,TFLAG
!  INTEGER, SAVE :: CCODE,RCDE,ISTAT(MPI_STATUS_SIZE),TSREQ(5),ierr
  
!  mynode=mynode

!  IF(PTINIT) THEN
!     TFLAG=.FALSE.
!     SENDNR=.TRUE.
!     PTINIT=.FALSE.
!  ENDIF
  
!  IF(NTOTAL.GE.NSTEP) TFLAG=.TRUE.
!     
!  IF(mynode.GT.0) THEN                     
!     IF(SENDNR.AND.(.NOT.TFLAG)) THEN                                                                           ! Send new receive
!        CALL MPI_IRECV(RCDE,1,MPI_INTEGER,mynode-1,MPI_ANY_TAG,MPI_COMM_WORLD,TSREQ(3),ierr)
!        SENDNR=.FALSE.
!     ENDIF
!     CALL MPI_TEST(TSREQ(3),SENDNR,ISTAT,ierr)
!     IF(SENDNR.AND.(.NOT.TFLAG)) CALL PTSWAP(Q,.FALSE.)
!  ENDIF
           
!  IF((CSTEP.GE.SWPINT).AND.(.NOT.TFLAG).AND.(mynode.LT.(NREPS-1))) THEN
!     CALL MPI_ISEND(1,1,MPI_INTEGER,mynode+1,0,MPI_COMM_WORLD,TSREQ(2),ierr)  
!     CALL MPI_WAIT(TSREQ(2),ISTAT,ierr)
!     CALL PTSWAP(Q,.TRUE.)
!     CSTEP=0              
!  ENDIF
    
!  ENDSIGNAL = TFLAG

! END SUBROUTINE ATTSWAPOLD
