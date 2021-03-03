!=======================================================================
! Generated by  : PSCAD v4.6.3.0
!
! Warning:  The content of this file is automatically generated.
!           Do not modify, as any changes made here will be lost!
!-----------------------------------------------------------------------
! Component     : Main
! Description   : Main Page
!-----------------------------------------------------------------------


!=======================================================================

      SUBROUTINE MainDyn(f)

!---------------------------------------
! Standard includes
!---------------------------------------

      INCLUDE 'nd.h'
      INCLUDE 'emtconst.h'
      INCLUDE 'emtstor.h'
      INCLUDE 's0.h'
      INCLUDE 's1.h'
      INCLUDE 's2.h'
      INCLUDE 's4.h'
      INCLUDE 'branches.h'
      INCLUDE 'pscadv3.h'
      INCLUDE 'fnames.h'
      INCLUDE 'radiolinks.h'
      INCLUDE 'matlab.h'
      INCLUDE 'rtconfig.h'

!---------------------------------------
! Function/Subroutine Declarations
!---------------------------------------


!---------------------------------------
! Variable Declarations
!---------------------------------------


! Subroutine Arguments
      REAL,    INTENT(IN)  :: f

! Electrical Node Indices
      INTEGER  NT_1(3)

! Control Signals
      INTEGER  Brk1, Flt1
      REAL     IFlt, Is(3), RT_1(3), IBrkA, IBrkB
      REAL     IBrkC, Ia, Ib, Ic, Vload(3)

! Internal Variables
      LOGICAL  LVD1_1
      INTEGER  IVD1_1, IVD1_2, IVD1_3, IVD1_4
      REAL     RVD1_1, RVD1_2, RVD1_3, RVD1_4
      REAL     RVD1_5, RVD1_6, RVD1_7

! Indexing variables
      INTEGER ICALL_NO                            ! Module call num
      INTEGER ISTOI, ISTOF, IT_0                  ! Storage Indices
      INTEGER IPGB                                ! Control/Monitoring
      INTEGER SS, INODE, IBRCH, IXFMR             ! SS/Node/Branch/Xfmr


!---------------------------------------
! Local Indices
!---------------------------------------

! Dsdyn <-> Dsout transfer index storage

      NTXFR = NTXFR + 1

      TXFR(NTXFR,1) = NSTOL
      TXFR(NTXFR,2) = NSTOI
      TXFR(NTXFR,3) = NSTOF
      TXFR(NTXFR,4) = NSTOC

! Define electric network subsystem number

      SS        = NODE(NNODE+1)

! Increment and assign runtime configuration call indices

      ICALL_NO  = NCALL_NO
      NCALL_NO  = NCALL_NO + 1

! Increment global storage indices

      ISTOI     = NSTOI
      NSTOI     = NSTOI + 2
      ISTOF     = NSTOF
      NSTOF     = NSTOF + 17
      IPGB      = NPGB
      NPGB      = NPGB + 7
      INODE     = NNODE + 2
      NNODE     = NNODE + 17
      IBRCH     = NBRCH(SS)
      NBRCH(SS) = NBRCH(SS) + 27
      IXFMR     = NXFMR
      NXFMR     = NXFMR + 3
      NCSCS     = NCSCS + 0
      NCSCR     = NCSCR + 0

!---------------------------------------
! Transfers from storage arrays
!---------------------------------------

      IFlt     = STOF(ISTOF + 2)
      Brk1     = STOI(ISTOI + 1)
      Flt1     = STOI(ISTOI + 2)
      IBrkA    = STOF(ISTOF + 9)
      IBrkB    = STOF(ISTOF + 10)
      IBrkC    = STOF(ISTOF + 11)
      Ia       = STOF(ISTOF + 12)
      Ib       = STOF(ISTOF + 13)
      Ic       = STOF(ISTOF + 14)

! Array (1:3) quantities...
      DO IT_0 = 1,3
         Is(IT_0) = STOF(ISTOF + 2 + IT_0)
         RT_1(IT_0) = STOF(ISTOF + 5 + IT_0)
         Vload(IT_0) = STOF(ISTOF + 14 + IT_0)
      END DO


!---------------------------------------
! Electrical Node Lookup
!---------------------------------------


! Array (1:3) quantities...
      DO IT_0 = 1,3
         NT_1(IT_0) = NODE(INODE + 0 + IT_0)
      END DO

!---------------------------------------
! Configuration of Models
!---------------------------------------

      IF ( TIMEZERO ) THEN
         FILENAME = 'Main.dta'
         CALL EMTDC_OPENFILE
         SECTION = 'DATADSD:'
         CALL EMTDC_GOTOSECTION
      ENDIF
!---------------------------------------
! Generated code from module definition
!---------------------------------------


! 10:[xfmr-3p2w] 3 Phase 2 Winding Transformer 'T1'
!  TRANSFORMER SATURATION SUBROUTINE
      IVD1_1 = NEXC
      CALL TSAT1_EXE( (IBRCH+10), (IBRCH+11), (IBRCH+12),SS,1.0,0)

! 40:[tbreak] Timed Breaker Logic (Obsolete) 
! Timed breaker logic
      CALL EMTDC_TBREAK(0.26, 0.31, 0, Brk1)

! 50:[tfault] Timed Fault Logic 
! Timed fault logic
      CALL EMTDC_TFAULT((0.25+0.05),0.25, 0, Flt1)

! 70:[breaker3] 3 Phase Breaker 'Brk1'
      IVD1_4 = NSTORI
      NSTORI = NSTORI + 3
! Three Phase Breaker
      CALL EMTDC_BREAKER1(SS, (IBRCH+22),0.01,1000000.0,RTCF(NRTCF),0,NI&
     &NT(1.0-REAL(Brk1)))
      CALL EMTDC_BREAKER1(SS, (IBRCH+23),0.01,1000000.0,RTCF(NRTCF),0,NI&
     &NT(1.0-REAL(Brk1)))
      CALL EMTDC_BREAKER1(SS, (IBRCH+24),0.01,1000000.0,RTCF(NRTCF),0,NI&
     &NT(1.0-REAL(Brk1)))
!
      IVD1_1 = 2*E_BtoI(OPENBR( (IBRCH+22),SS))
      IVD1_2 = 2*E_BtoI(OPENBR( (IBRCH+23),SS))
      IVD1_3 = 2*E_BtoI(OPENBR( (IBRCH+24),SS))
      NRTCF = NRTCF + 1
      IF (FIRSTSTEP .OR. (STORI(IVD1_4+0) .NE. IVD1_1)) THEN
         CALL PSCAD_AGI2(ICALL_NO,761528435,IVD1_1,"BOpen1")
      ENDIF
      IF (FIRSTSTEP .OR. (STORI(IVD1_4+1) .NE. IVD1_2)) THEN
         CALL PSCAD_AGI2(ICALL_NO,761528435,IVD1_2,"BOpen2")
      ENDIF
      IF (FIRSTSTEP .OR. (STORI(IVD1_4+2) .NE. IVD1_3)) THEN
         CALL PSCAD_AGI2(ICALL_NO,761528435,IVD1_3,"BOpen3")
      ENDIF
      STORI(IVD1_4+0) = 2*E_BtoI(OPENBR( (IBRCH+22),SS))
      STORI(IVD1_4+1) = 2*E_BtoI(OPENBR( (IBRCH+23),SS))
      STORI(IVD1_4+2) = 2*E_BtoI(OPENBR( (IBRCH+24),SS))

! 80:[tpflt] Three Phase Fault 
      CALL E3PHFLT1_EXE(SS, (IBRCH+16), (IBRCH+17), (IBRCH+18), (IBRCH+1&
     &9), (IBRCH+20), (IBRCH+21),0,Flt1,10,0.01)
      LVD1_1 = (OPENBR( (IBRCH+16),SS).AND.OPENBR( (IBRCH+17),SS).AND.OP&
     &ENBR( (IBRCH+18),SS).AND.OPENBR( (IBRCH+19),SS).AND.OPENBR( (IBRCH&
     &+20),SS).AND.OPENBR( (IBRCH+21),SS))
      IVD1_1 = E_BtoI(LVD1_1)
      IF(FIRSTSTEP .OR. (IVD1_1 .NE. STORI(NSTORI))) THEN
         CALL PSCAD_AGI2(ICALL_NO,2087521150,1-IVD1_1,"AG1")
         STORI(NSTORI) = IVD1_1
      ENDIF
      NSTORI = NSTORI + 1

! 90:[datamerge] Merges data signals into an array 
      RT_1(1) = IBrkA
      RT_1(2) = IBrkB
      RT_1(3) = IBrkC

! 100:[pgb] Output Channel 'Two Phase to Ground Fault'

      PGB(IPGB+4) = IFlt

! 110:[pgb] Output Channel '3 Phase Breaker Current'

      DO IVD1_1 = 1, 3
         PGB(IPGB+5+IVD1_1-1) = RT_1(IVD1_1)
      ENDDO

! 1:[source3] Three Phase Voltage Source Model 1 'Src'
!  3-Phase source: Src
      RVD1_1 = RTCF(NRTCF+12)
      RVD1_2 = RTCF(NRTCF+14)
      RVD1_3 = RTCF(NRTCF+13)
      CALL ESYS651_EXE(SS, (IBRCH+1), (IBRCH+2), (IBRCH+3),0,0,0, SS, NT&
     &_1(1),NT_1(2),NT_1(3), 0, RVD1_2, RVD1_1, 0.05, 1.0, 1.0, 1.0,RVD1&
     &_3, 1.0, 0.02, 0.05, 1.0, 0.02, 0.05, RVD1_4, RVD1_5, RVD1_6, RVD1&
     &_7)

!---------------------------------------
! Feedbacks and transfers to storage
!---------------------------------------

      STOF(ISTOF + 2) = IFlt
      STOI(ISTOI + 1) = Brk1
      STOI(ISTOI + 2) = Flt1
      STOF(ISTOF + 9) = IBrkA
      STOF(ISTOF + 10) = IBrkB
      STOF(ISTOF + 11) = IBrkC
      STOF(ISTOF + 12) = Ia
      STOF(ISTOF + 13) = Ib
      STOF(ISTOF + 14) = Ic

! Array (1:3) quantities...
      DO IT_0 = 1,3
         STOF(ISTOF + 2 + IT_0) = Is(IT_0)
         STOF(ISTOF + 5 + IT_0) = RT_1(IT_0)
         STOF(ISTOF + 14 + IT_0) = Vload(IT_0)
      END DO


!---------------------------------------
! Transfer to Exports
!---------------------------------------

!---------------------------------------
! Close Model Data read
!---------------------------------------

      IF ( TIMEZERO ) CALL EMTDC_CLOSEFILE
      RETURN
      END

!=======================================================================

      SUBROUTINE MainOut()

!---------------------------------------
! Standard includes
!---------------------------------------

      INCLUDE 'nd.h'
      INCLUDE 'emtconst.h'
      INCLUDE 'emtstor.h'
      INCLUDE 's0.h'
      INCLUDE 's1.h'
      INCLUDE 's2.h'
      INCLUDE 's4.h'
      INCLUDE 'branches.h'
      INCLUDE 'pscadv3.h'
      INCLUDE 'fnames.h'
      INCLUDE 'radiolinks.h'
      INCLUDE 'matlab.h'
      INCLUDE 'rtconfig.h'

!---------------------------------------
! Function/Subroutine Declarations
!---------------------------------------

      REAL    EMTDC_VVDC    ! 
      REAL    VBRANCH       ! 

!---------------------------------------
! Variable Declarations
!---------------------------------------


! Electrical Node Indices
      INTEGER  NT_3(3)

! Control Signals
      REAL     IFlt, Is(3), IBrkA, IBrkB, IBrkC
      REAL     Ia, Ib, Ic, Vload(3)

! Internal Variables
      INTEGER  IVD1_1
      REAL     RVD1_1, RVD1_2

! Indexing variables
      INTEGER ICALL_NO                            ! Module call num
      INTEGER ISTOL, ISTOI, ISTOF, ISTOC, IT_0    ! Storage Indices
      INTEGER IPGB                                ! Control/Monitoring
      INTEGER SS, INODE, IBRCH, IXFMR             ! SS/Node/Branch/Xfmr


!---------------------------------------
! Local Indices
!---------------------------------------

! Dsdyn <-> Dsout transfer index storage

      NTXFR = NTXFR + 1

      ISTOL = TXFR(NTXFR,1)
      ISTOI = TXFR(NTXFR,2)
      ISTOF = TXFR(NTXFR,3)
      ISTOC = TXFR(NTXFR,4)

! Define electric network subsystem number

      SS        = NODE(NNODE+1)

! Increment and assign runtime configuration call indices

      ICALL_NO  = NCALL_NO
      NCALL_NO  = NCALL_NO + 1

! Increment global storage indices

      IPGB      = NPGB
      NPGB      = NPGB + 7
      INODE     = NNODE + 2
      NNODE     = NNODE + 17
      IBRCH     = NBRCH(SS)
      NBRCH(SS) = NBRCH(SS) + 27
      IXFMR     = NXFMR
      NXFMR     = NXFMR + 3
      NCSCS     = NCSCS + 0
      NCSCR     = NCSCR + 0

!---------------------------------------
! Transfers from storage arrays
!---------------------------------------

      IFlt     = STOF(ISTOF + 2)
      IBrkA    = STOF(ISTOF + 9)
      IBrkB    = STOF(ISTOF + 10)
      IBrkC    = STOF(ISTOF + 11)
      Ia       = STOF(ISTOF + 12)
      Ib       = STOF(ISTOF + 13)
      Ic       = STOF(ISTOF + 14)

! Array (1:3) quantities...
      DO IT_0 = 1,3
         Is(IT_0) = STOF(ISTOF + 2 + IT_0)
         Vload(IT_0) = STOF(ISTOF + 14 + IT_0)
      END DO


!---------------------------------------
! Electrical Node Lookup
!---------------------------------------


! Array (1:3) quantities...
      DO IT_0 = 1,3
         NT_3(IT_0) = NODE(INODE + 6 + IT_0)
      END DO

!---------------------------------------
! Configuration of Models
!---------------------------------------

      IF ( TIMEZERO ) THEN
         FILENAME = 'Main.dta'
         CALL EMTDC_OPENFILE
         SECTION = 'DATADSO:'
         CALL EMTDC_GOTOSECTION
      ENDIF
!---------------------------------------
! Generated code from module definition
!---------------------------------------


! 10:[xfmr-3p2w] 3 Phase 2 Winding Transformer 'T1'
      Ia = CDCTR(2,(IXFMR + 1))
      Ib = CDCTR(2,(IXFMR + 2))
      Ic = CDCTR(2,(IXFMR + 3))

! 20:[ammeter] Current Meter 'Is'
      Is(1) = ( CBR((IBRCH+13), SS))
      Is(2) = ( CBR((IBRCH+14), SS))
      Is(3) = ( CBR((IBRCH+15), SS))

! 30:[voltmetergnd] Voltmeter (Line - Ground) 'Vload'
      Vload(1) = EMTDC_VVDC(SS, NT_3(1), 0)
      Vload(2) = EMTDC_VVDC(SS, NT_3(2), 0)
      Vload(3) = EMTDC_VVDC(SS, NT_3(3), 0)

! 60:[pgb] Output Channel '3 Phase Source Current'

      DO IVD1_1 = 1, 3
         PGB(IPGB+1+IVD1_1-1) = Is(IVD1_1)
      ENDDO

! 70:[breaker3] 3 Phase Breaker 'Brk1'
! Three Phase Breaker Currents
      IBrkA = ( CBR((IBRCH+22), SS))
      IBrkB = ( CBR((IBRCH+23), SS))
      IBrkC = ( CBR((IBRCH+24), SS))
      CALL BRK_POWER(SS, (IBRCH+22), (IBRCH+23), (IBRCH+24),0,0,0,IVD1_1&
     &,0.02,RVD1_1,RVD1_2)

! 80:[tpflt] Three Phase Fault 
!
! Multi-phase Fault Currents
!
      IFlt =  ( CBR((IBRCH+21), SS)) + ( CBR((IBRCH+18), SS)) - ( CBR((I&
     &BRCH+17), SS)) 
!

!---------------------------------------
! Feedbacks and transfers to storage
!---------------------------------------

      STOF(ISTOF + 2) = IFlt
      STOF(ISTOF + 9) = IBrkA
      STOF(ISTOF + 10) = IBrkB
      STOF(ISTOF + 11) = IBrkC
      STOF(ISTOF + 12) = Ia
      STOF(ISTOF + 13) = Ib
      STOF(ISTOF + 14) = Ic

! Array (1:3) quantities...
      DO IT_0 = 1,3
         STOF(ISTOF + 2 + IT_0) = Is(IT_0)
         STOF(ISTOF + 14 + IT_0) = Vload(IT_0)
      END DO


!---------------------------------------
! Close Model Data read
!---------------------------------------

      IF ( TIMEZERO ) CALL EMTDC_CLOSEFILE
      RETURN
      END

!=======================================================================

      SUBROUTINE MainDyn_Begin(f)

!---------------------------------------
! Standard includes
!---------------------------------------

      INCLUDE 'nd.h'
      INCLUDE 'emtconst.h'
      INCLUDE 's0.h'
      INCLUDE 's1.h'
      INCLUDE 's4.h'
      INCLUDE 'branches.h'
      INCLUDE 'pscadv3.h'
      INCLUDE 'radiolinks.h'
      INCLUDE 'rtconfig.h'

!---------------------------------------
! Function/Subroutine Declarations
!---------------------------------------


!---------------------------------------
! Variable Declarations
!---------------------------------------


! Subroutine Arguments
      REAL,    INTENT(IN)  :: f

! Electrical Node Indices

! Control Signals

! Internal Variables
      INTEGER  IVD1_1
      REAL     RVD1_1, RVD1_2, RVD1_3, RVD1_4
      REAL     RVD1_5, RVD1_6

! Indexing variables
      INTEGER ICALL_NO                            ! Module call num
      INTEGER IT_0                                ! Storage Indices
      INTEGER SS, INODE, IBRCH, IXFMR             ! SS/Node/Branch/Xfmr


!---------------------------------------
! Local Indices
!---------------------------------------


! Define electric network subsystem number

      SS        = NODE(NNODE+1)

! Increment and assign runtime configuration call indices

      ICALL_NO  = NCALL_NO
      NCALL_NO  = NCALL_NO + 1

! Increment global storage indices

      INODE     = NNODE + 2
      NNODE     = NNODE + 17
      IBRCH     = NBRCH(SS)
      NBRCH(SS) = NBRCH(SS) + 27
      IXFMR     = NXFMR
      NXFMR     = NXFMR + 3
      NCSCS     = NCSCS + 0
      NCSCR     = NCSCR + 0

!---------------------------------------
! Electrical Node Lookup
!---------------------------------------


!---------------------------------------
! Generated code from module definition
!---------------------------------------


! 10:[xfmr-3p2w] 3 Phase 2 Winding Transformer 'T1'
      CALL COMPONENT_ID(ICALL_NO,72575535)
      RVD1_1 = ONE_3RD*100.0
      RVD1_2 = 13.8
      RVD1_3 = 230.0*SQRT_1BY3
      CALL E_TF2W_CFG((IXFMR + 1),0,RVD1_1,60.0,0.1,0.0,RVD1_2,RVD1_3,1.&
     &0)
      CALL E_TF2W_CFG((IXFMR + 2),0,RVD1_1,60.0,0.1,0.0,RVD1_2,RVD1_3,1.&
     &0)
      CALL E_TF2W_CFG((IXFMR + 3),0,RVD1_1,60.0,0.1,0.0,RVD1_2,RVD1_3,1.&
     &0)
      IF (0.0 .LT. 1.0E-6) THEN
        RVD1_5 = 0.0
        RVD1_6 = 0.0
        IVD1_1 = 0
      ELSE
        RVD1_6 = 0.0
        RVD1_4 = 6.0/(100.0*RVD1_6)
        RVD1_5 = RVD1_4*RVD1_2*RVD1_2
        RVD1_6 = RVD1_4*RVD1_3*RVD1_3
        IVD1_1 = 1
      ENDIF
      CALL E_BRANCH_CFG( (IBRCH+4),SS,IVD1_1,0,0,RVD1_5,0.0,0.0)
      CALL E_BRANCH_CFG( (IBRCH+5),SS,IVD1_1,0,0,RVD1_5,0.0,0.0)
      CALL E_BRANCH_CFG( (IBRCH+6),SS,IVD1_1,0,0,RVD1_5,0.0,0.0)
      CALL E_BRANCH_CFG( (IBRCH+7),SS,IVD1_1,0,0,RVD1_6,0.0,0.0)
      CALL E_BRANCH_CFG( (IBRCH+8),SS,IVD1_1,0,0,RVD1_6,0.0,0.0)
      CALL E_BRANCH_CFG( (IBRCH+9),SS,IVD1_1,0,0,RVD1_6,0.0,0.0)
      CALL TSAT1_CFG( (IBRCH+10), (IBRCH+11), (IBRCH+12),SS,RVD1_1,RVD1_&
     &2,0.2,1.25,60.0,1.0,1.0,0.0)

! 70:[breaker3] 3 Phase Breaker 'Brk1'
      CALL COMPONENT_ID(ICALL_NO,761528435)
      RTCF(NRTCF) = ABS(0.0)
      NRTCF = NRTCF + 1

! 80:[tpflt] Three Phase Fault 
      CALL E3PHFLT1_CFG(1000000.0,0.0)

! 90:[datamerge] Merges data signals into an array 

! 100:[pgb] Output Channel 'Two Phase to Ground Fault'

! 110:[pgb] Output Channel '3 Phase Breaker Current'

! 1:[source3] Three Phase Voltage Source Model 1 'Src'
      CALL COMPONENT_ID(ICALL_NO,280552298)
      RVD1_1 = 1.0
      RVD1_2 = 0.1
      CALL ESYS651_CFG(3,4,0,0,0,SS, (IBRCH+1), (IBRCH+2), (IBRCH+3),0,0&
     &,0, 60.0,60.0,0.0,13.8,0.0,0.0,100.0,230.0,230.0, 0.5,80.0,2.0,1.0&
     &,1.0,0.1, 1.0,80.0,RVD1_1,RVD1_2)

      RETURN
      END

!=======================================================================

      SUBROUTINE MainOut_Begin(f)

!---------------------------------------
! Standard includes
!---------------------------------------

      INCLUDE 'nd.h'
      INCLUDE 'emtconst.h'
      INCLUDE 's0.h'
      INCLUDE 's1.h'
      INCLUDE 's4.h'
      INCLUDE 'branches.h'
      INCLUDE 'pscadv3.h'
      INCLUDE 'radiolinks.h'
      INCLUDE 'rtconfig.h'

!---------------------------------------
! Function/Subroutine Declarations
!---------------------------------------


!---------------------------------------
! Variable Declarations
!---------------------------------------


! Subroutine Arguments
      REAL,    INTENT(IN)  :: f

! Electrical Node Indices
      INTEGER  NT_3(3)

! Control Signals

! Internal Variables

! Indexing variables
      INTEGER ICALL_NO                            ! Module call num
      INTEGER IT_0                                ! Storage Indices
      INTEGER SS, INODE, IBRCH, IXFMR             ! SS/Node/Branch/Xfmr


!---------------------------------------
! Local Indices
!---------------------------------------


! Define electric network subsystem number

      SS        = NODE(NNODE+1)

! Increment and assign runtime configuration call indices

      ICALL_NO  = NCALL_NO
      NCALL_NO  = NCALL_NO + 1

! Increment global storage indices

      INODE     = NNODE + 2
      NNODE     = NNODE + 17
      IBRCH     = NBRCH(SS)
      NBRCH(SS) = NBRCH(SS) + 27
      IXFMR     = NXFMR
      NXFMR     = NXFMR + 3
      NCSCS     = NCSCS + 0
      NCSCR     = NCSCR + 0

!---------------------------------------
! Electrical Node Lookup
!---------------------------------------


! Array (1:3) quantities...
      DO IT_0 = 1,3
         NT_3(IT_0) = NODE(INODE + 6 + IT_0)
      END DO

!---------------------------------------
! Generated code from module definition
!---------------------------------------


! 60:[pgb] Output Channel '3 Phase Source Current'

      RETURN
      END
