     PROGRAM PBMCVGW
     IMPLICIT NONE
     INCLUDE "mpif.h"
     DOUBLE PRECISION, ALLOCATABLE :: Q(:,:)
     DOUBLE PRECISION :: TSTART,SL,NUM,TMAX,BLX,BLY,BLZ,TMIN,RCUTOFF,RCON,ENOT
     DOUBLE PRECISION :: VLEN,PMIN,PMAX,RFRSH,MAXBL,BLDIFF
     INTEGER :: I,EOF,NP,mynode,WF,UF,CCODE,VWRF,CONT,EG,VOF,ierr
     INTEGER*8 :: NSTEP
     INTEGER :: NDE,WRF,GRIDSIZE,ISTAT(MPI_STATUS_SIZE),VSFREQ,SWPINT,NATOMS=0
     INTEGER :: NEQL,FRCFRQ,PWFRQ
     CHARACTER*2 :: TAG
     LOGICAL :: INSIDE

     CALL MPI_INIT(ierr)
     CALL MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
     CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NP,ierr)

     CALL CPU_TIME(TSTART)

     OPEN(UNIT=17,FILE='coords',STATUS='OLD')
     DO 
        READ(17,*,IOSTAT=EOF)
        IF (EOF==0) THEN
           NATOMS=NATOMS+1
        ELSE
           EXIT
        ENDIF
     ENDDO
     
     REWIND(17)
   
     ALLOCATE(Q(3,NATOMS))
     
     DO I=1,NATOMS
        READ(17,*) Q(1,I),Q(2,I),Q(3,I)
     ENDDO

     CLOSE(17)
      
     OPEN(UNIT=8,FILE='inputf.dat',STATUS='OLD')        
        READ(8,*) TMIN                               ! Minimum temperature for PTMC
        READ(8,*) TMAX                               ! Maximum temperature  "   "
        READ(8,*) PMIN                               ! Minimum Pressure
        READ(8,*) PMAX                               ! Maximum Pressure
        READ(8,*) NSTEP                              ! Number of MC steps
        READ(8,*) NEQL                               ! Number of equilibration steps
        READ(8,*) BLX                                ! Box length x (for PBC)
        READ(8,*) SL                                 ! Step length
        READ(8,*) WF                                 ! Reset # accepted moves and scale step length after WF steps
        READ(8,*) UF                                 ! Observables update frequency
        READ(8,*) WRF                                ! Observable write freq.
        READ(8,*) RCON                               ! Constraining radius
        READ(8,*) GRIDSIZE                           ! Tau grid size
        READ(8,*) EG                                 ! Equidistant grid
        READ(8,*) SWPINT                             ! Swap interval
        READ(8,*) TAG                                ! Atom name label (for VMD)
        READ(8,*) VWRF                               ! VMD coords write freq
        READ(8,*) CONT                               ! Continue previous MC run (1 = yes)
        READ(8,*) ENOT                               ! E zero
        READ(8,*) VLEN                               ! Volume step length
        READ(8,*) VSFREQ                             ! Attempt volume change frequency (after every n MC steps)
        READ(8,*) RFRSH                              ! Refresh frequency
        READ(8,*) MAXBL                              ! Maximum box length
        READ(8,*) VOF                                ! Volume output (or Node dump) frequency
        READ(8,*) BLDIFF                             ! Compute pressure via finite differences, box length difference
        READ(8,*) FRCFRQ                             ! Frequency for writing net forces
        READ(8,*) PWFRQ                              ! Frequency for writing pair-wise forces
     CLOSE(8)   
   
     INSIDE=.TRUE.

     IF(.NOT.CONT) THEN
        IF(BLX.GT.0) THEN
           CALL RCENTER(Q,NATOMS,RCON,INSIDE)
           CALL PBCHECK(Q,NATOMS,BLX,INSIDE)
        ELSE
           CALL RCENTER(Q,NATOMS,RCON,INSIDE)
        ENDIF
     ENDIF

     IF(.NOT.INSIDE) THEN
        WRITE(*,*) "CLUSTER RADIUS EXCEEDS CONSRAINING RADIUS OR THAT OF BOX. EXITING."
        GOTO 70
     ELSE
        CALL INITIALIZE_LJNPTMOD(NATOMS,TMAX,TMIN,PMIN,PMAX,BLX,NP,GRIDSIZE,mynode,ENOT,EG,BLDIFF)
        CALL INITIALIZE_LJNPTMC(NSTEP,UF,WRF,WF,SWPINT,SL,VLEN,RCON,Q,CONT,NATOMS,VSFREQ,RFRSH,MAXBL,VOF,NEQL,FRCFRQ,PWFRQ)        
        CALL MCSTEP(Q,VWRF,TAG)
        CALL CLEANUP_LJNPT
     ENDIF

70   CALL MPI_BARRIER(MPI_COMM_WORLD)
        
     DEALLOCATE(Q)
          
     CALL MPI_FINALIZE(ierr)

  END PROGRAM PBMCVGW     
           
