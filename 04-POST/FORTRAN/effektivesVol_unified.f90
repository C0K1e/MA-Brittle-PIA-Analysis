PROGRAM effektivesVol
!
! Berechnung des effektiven Volumens (Weibull / Weakest-Link) in ANSYS-FEM.
! Aktive Version v5.0 (Stand: 2026-04-29) — vereint Pfade overhead/_GaussNorm/
! _highOrder/_RAWNodes.
!
! Abhaengigkeiten:
!   - overhead_unified.f90  (MODULE overhead_unified)
!   - ANSYS-Makro "x_Effektives_Volumen_NOD.mac" oder "_RAW.mac"
!     erzeugt effVol_Parameter.out / effVol_Elemente.out / effVol_Faces.out
!     / effVol_Nodes.out (NOD) bzw. effVol_NodeCoords.out + effVol_RawStress.out (RAW)
!
! Steuerung via effVol_Parameter.out:
!   stress_mode = 0 (NOD): knotengemittelte Spannungen aus effVol_Nodes.out
!   stress_mode = 1 (RAW): elementlokale Spannungen aus effVol_RawStress.out
!   do_gaussnorm = .TRUE. : domain-spezifische GaussNorm-Korrektur (VOL/SURF/LINE)
!
! Output-Files:
!   {ausgabename}.out             Haupt-Ergebnis (4 Spalten: S1-Veff, PIA-Veff,
!                                                            S1-Seff, PIA-Seff)
!   {ausgabename}_GauszInfo.out   Domain-spezifische Maxima (nur bei do_gaussnorm)
!
! Vollstaendige Versions-Historie + methodische Details:
!   siehe HISTORY.md im selben Verzeichnis (v17.0 ausgelagert)
!
USE overhead_unified

IMPLICIT none

real, dimension(:,:), allocatable :: veff, seff
! 2D Arrays fuer effektives Volumen und Oberflaeche.
logical, dimension(:), allocatable :: elem_nodes, surf_nodes
! Masken Arrays fuer Skalierungsspannung (nur NOD-Modus)
integer :: m, e, i, j, n, symm_fakt, element_type, m_max, error=0
! Zaehler Variablen und Programm Parameter
integer :: nelem, nface, nnode, spalten_elem, spalten_face
integer :: nraw_lines, nraw_expected   ! V9: Fuer RawStress-Validierung
character (len=256) :: ausgabename
character (len=10)  :: error_num1, error_num2
character (len=32)  :: error_loc
character (len=256) :: error_msg
real, parameter :: pi=3.14159265358979
real :: smax_file, smax_elem, smax_surf, breite_2D, rad
! V9: GaussNorm-Variablen
real :: korrektur_faktor
! V9: Hilfsvariablen fuer erweitertes Parameter-Lesen
integer :: io_stat, stress_mode_int, gaussnorm_int
integer :: dummy_int
! v12.2.2: Float-Lesepuffer fuer AVGModeRAW/GaussPNorm (APDL F14.0 schreibt "1." statt "1")
real :: stress_mode_real, gaussnorm_real

!
! #################################################################################################
! Vorbereitung: Einlesen und allokieren von Arrays
! #################################################################################################
!
! Einlesen des Parameter-Files

! v16.1 (B6 + Phase 5 Fatal): IOSTAT-Pattern fuer alle OPEN-Statements
OPEN (101, FILE='effVol_Parameter.out', STATUS='old', IOSTAT=io_stat)
IF (io_stat /= 0) THEN
  WRITE(*,*) ""
  WRITE(*,*) "!!! [Fortran] FATAL: effVol_Parameter.out nicht gefunden oder unlesbar."
  WRITE(*,*) "!!!           Erwartet im aktuellen Arbeitsverzeichnis."
  WRITE(*,*) "!!!           Aktion: APDL-Macro x_Effektives_Volumen_NOD.mac/_RAW.mac"
  WRITE(*,*) "!!!                   muss diese Datei vor dem Fortran-Aufruf erzeugen."
  STOP 1
END IF
READ(101,'(8X,A256)') ausgabename     ! Name des Output-Files
READ(101,'(8X,ES16.8)') smax_file     ! Bekannte Maximalspannung (nodal)
READ(101,'(8X,I16)') gauss_order       ! Grad der Gausz-Legendre Quadratur
READ(101,'(8X,I16)') m_max            ! Weibull Modul bis zu dem berechnet wird
READ(101,'(8X,I16)') symm_fakt        ! Symmetriefaktor fuer gespiegelte Modelle, etc.
READ(101,'(8X,I16)')  element_type    ! Indikator fuer Elementgeometrie
READ(101,'(8X,F16.8)') breite_2D      ! Breite von 2D Elementen, falls zutreffend
! RAWNodes Erweiterung: AVGModeRAW und GaussPNorm (optional, rueckwaertskompatibel)
! v12.2.2 Fix: Lese als FLOAT (F14.0), nicht Integer (I14) —
! APDL *VWRITE schreibt mit (F14.0,' ') Format das Token "            1." (mit Punkt).
! Integer-Parser wuerde am Punkt stolpern und IOSTAT/=0 zurueckgeben.
! Float-Parser akzeptiert beide Formate ("1" und "1.") robust.
READ(101,'(11X,F14.0)', IOSTAT=io_stat) stress_mode_real
IF (io_stat /= 0) stress_mode_real = 0.0  ! Default: NOD
stress_mode_int = NINT(stress_mode_real)
READ(101,'(11X,F14.0)', IOSTAT=io_stat) gaussnorm_real
IF (io_stat /= 0) gaussnorm_real = 0.0    ! Default: kein GaussNorm
gaussnorm_int = NINT(gaussnorm_real)
CLOSE(101)

! RAWNodes: Modul-Variablen setzen
stress_mode = stress_mode_int
do_gaussnorm = (gaussnorm_int == 1)

ausgabename = ADJUSTL(ausgabename)

! DEBUGGING
! v16.1 (A8): Klartext-Konversion fuer stress_mode + do_gaussnorm
WRITE(*,'(A8,ES16.8)') 'SMAX  = ',  smax_file
WRITE(*,'(A8,I16)')    'GAUSS = ', gauss_order
WRITE(*,'(A8,I16)')    'MMAX  = ', m_max
WRITE(*,'(A8,I16)')    'SYMM  = ', symm_fakt
WRITE(*,'(A8,I16)')    'GEOM  = ', element_type
WRITE(*,'(A8,F16.8)')  'BREIT = ', breite_2D
IF (stress_mode == 0) THEN
  WRITE(*,'(A22)') 'AVGMode = NOD (knotengemittelt, ANSYS-Default)'
ELSE
  WRITE(*,'(A22)') 'AVGMode = RAW (elementlokal aus PRESOL)'
END IF
IF (do_gaussnorm) THEN
  WRITE(*,'(A28)') 'GaussPNorm = ON (domain-spezifisch VOL/SURF/LINE)'
ELSE
  WRITE(*,'(A20)') 'GaussPNorm = OFF (keine Korrektur)'
END IF

arraybreiten_zuweisung: SELECT CASE (element_type)
  CASE (23)
    spalten_elem = 4
    spalten_face = 2
  CASE (24)
    spalten_elem = 4
    spalten_face = 2
  CASE (33)
    spalten_elem = 8
    spalten_face = 4
  CASE (34)
    spalten_elem = 8
    spalten_face = 4
  CASE DEFAULT
    error_loc='arraybreiten_zuweisung'
    error_msg='Fehler beim Zuweisen der Arraybreiten. ueberpruefen Sie die Elementtype im Parameter-File.'
    CALL exit_error(error_loc,error_msg)
END SELECT arraybreiten_zuweisung

! DEBUGGING
WRITE(*,'(A15,I1)') 'spalten_elem = ', spalten_elem
WRITE(*,'(A15,I1)') 'spalten_face = ', spalten_face
PRINT *, 'Reading Element-File:'

! File-laenge lesen, Anzahl Reihen bestimmen
OPEN (102, FILE='effVol_Elemente.out', STATUS='old', IOSTAT=io_stat)
IF (io_stat /= 0) THEN
  WRITE(*,*) "!!! [Fortran] FATAL: effVol_Elemente.out nicht gefunden (IOSTAT=", io_stat, ")"
  WRITE(*,*) "!!!           APDL-Macro hat die Element-Liste nicht exportiert."
  STOP 1
END IF
nelem = 0
DO
  READ(102,*, IOSTAT=error)
  IF (error < 0) EXIT
  nelem = nelem + 1
END DO
REWIND(102)

! DEBUGGING
WRITE(*,'(A28,I10)') ',.................nelem = ', nelem
PRINT *, 'Allocating Elements'

ALLOCATE (elements(spalten_elem,nelem))

! DEBUGGING
PRINT *, '........................................ Check'
PRINT *, 'Filling Element-Array'

! V9: Mode-abhaengiges Einlesen
IF (stress_mode == 1) THEN
  ! RAW: 9 Spalten (ElemID + 8 Knoten). Erste Spalte ueberspringen.
  DO i = 1, nelem
    READ(102,*) dummy_int, (elements(j,i), j=1,spalten_elem)
  END DO
ELSE
  ! NOD: 8 Spalten wie bisher (Bulk-Read)
  READ(102,*) elements
END IF
CLOSE(102)

! DEBUGGING
PRINT *, '........................................ Done'
PRINT *, 'Reading Face-File'


OPEN (103, FILE='effVol_Faces.out', STATUS='old', IOSTAT=io_stat)
IF (io_stat /= 0) THEN
  WRITE(*,*) "!!! [Fortran] FATAL: effVol_Faces.out nicht gefunden (IOSTAT=", io_stat, ")"
  WRITE(*,*) "!!!           APDL-Macro hat die Face-Liste nicht exportiert."
  STOP 1
END IF
nface = 0
DO
  READ(103,*, IOSTAT=error)
  IF (error < 0) EXIT
  nface = nface + 1
END DO
REWIND(103)

! DEBUGGING
WRITE(*,'(A28,I10)') ',................nface = ', nface
PRINT *, 'Allocating Faces'

ALLOCATE (faces(spalten_face,nface))

! DEBUGGING
PRINT *, '........................................ Check'
PRINT *, 'Filling Face-Array'

! V9: Mode-abhaengiges Einlesen
IF (stress_mode == 1) THEN
  ! RAW: 5 Spalten (ParentIdx + 4 Knoten). Erste Spalte = Parent-Element Index.
  ALLOCATE(face_parent_elem(nface))
  DO i = 1, nface
    READ(103,*) face_parent_elem(i), (faces(j,i), j=1,spalten_face)
  END DO
  ! Validierung: Parent-Indizes muessen im Bereich [1, nelem] liegen
  DO i = 1, nface
    IF (face_parent_elem(i) < 1 .OR. face_parent_elem(i) > nelem) THEN
      error_loc='face_parent_check'
      WRITE(error_num1,'(I10)') face_parent_elem(i)
      WRITE(error_num2,'(I10)') i
      error_msg='Face '//TRIM(ADJUSTL(error_num2))//' hat unguelitgen Parent-Index '//TRIM(ADJUSTL(error_num1))
      CALL exit_error(error_loc,error_msg)
      STOP
    END IF
  END DO
ELSE
  ! NOD: 4 Spalten wie bisher (Bulk-Read)
  READ(103,*) faces
END IF
CLOSE(103)

! DEBUGGING
PRINT *, '........................................ Done'
PRINT *, 'Checking for Tets, Tris, rearranging'

SELECT CASE (element_type)
  CASE (23)
    elements(4,:) = elements(3,:)
    element_type = 24
  CASE (33)
    elements(8,:) = elements(4,:)
    elements(7,:) = elements(4,:)
    elements(6,:) = elements(4,:)
    elements(5,:) = elements(4,:)
    elements(4,:) = elements(3,:)
    faces(4,:)    = faces(3,:)
    element_type  = 34
END SELECT

! DEBUGGING
PRINT *, '........................................ Check'

! ==================================================================================
! V9: Mode-abhaengiges Node/Stress-File Einlesen
! ==================================================================================

IF (stress_mode == 0) THEN
  ! ----- NOD-Modus: effVol_Nodes.out (6 Spalten: x,y,z,s1,s2,s3) -----
  PRINT *, 'Node-File Zeilen zaehlen (NOD-Modus)'

  OPEN (104, FILE='effVol_Nodes.out', STATUS='old', IOSTAT=io_stat)
  IF (io_stat /= 0) THEN
    WRITE(*,*) "!!! [Fortran] FATAL: effVol_Nodes.out nicht gefunden (NOD-Modus, IOSTAT=", io_stat, ")"
    WRITE(*,*) "!!!           APDL-Macro hat die Knoten-Spannungsliste nicht exportiert."
    STOP 1
  END IF
  nnode = 0
  DO
    READ(104,*, IOSTAT=error)
    IF (error < 0) EXIT
    nnode = nnode + 1
  END DO
  REWIND(104)

  WRITE(*,'(A28,I10)') ',................nnode = ', nnode
  PRINT *, 'Reading Node-File'

  ALLOCATE (nodes(nnode,6))
  DO i=1,nnode
    READ(104,*) (nodes(i,j),j=1,6)
  END DO
  CLOSE(104)

  PRINT *, '........................................ Done'

ELSE
  ! ----- RAW-Modus: Koordinaten + Elementlokale Spannungen getrennt -----
  PRINT *, 'NodeCoords-File Zeilen zaehlen (RAW-Modus)'

  ! 1. Koordinaten-File: effVol_NodeCoords.out (3 Spalten: x,y,z)
  OPEN (104, FILE='effVol_NodeCoords.out', STATUS='old', IOSTAT=io_stat)
  IF (io_stat /= 0) THEN
    WRITE(*,*) "!!! [Fortran] FATAL: effVol_NodeCoords.out nicht gefunden (RAW-Modus, IOSTAT=", io_stat, ")"
    WRITE(*,*) "!!!           Python-Pipeline hat die Node-Koordinaten nicht exportiert."
    STOP 1
  END IF
  nnode = 0
  DO
    READ(104,*, IOSTAT=error)
    IF (error < 0) EXIT
    nnode = nnode + 1
  END DO
  REWIND(104)

  WRITE(*,'(A28,I10)') ',................nnode = ', nnode
  PRINT *, 'Reading NodeCoords-File'

  ALLOCATE (nodes_coords(nnode,3))
  DO i=1,nnode
    READ(104,*) (nodes_coords(i,j),j=1,3)
  END DO
  CLOSE(104)

  PRINT *, '........................................ Done'

  ! 2. RawStress-File: effVol_RawStress.out (3 Spalten: s1,s2,s3)
  !    nelem * spalten_elem Zeilen, gruppiert nach Element
  PRINT *, 'Reading RawStress-File'

  OPEN (105, FILE='effVol_RawStress.out', STATUS='old', IOSTAT=io_stat)
  IF (io_stat /= 0) THEN
    error_loc='raw_stress_open'
    error_msg='effVol_RawStress.out nicht gefunden! Im RAW-Modus (AVGModeRAW=1) ist diese Datei erforderlich.'
    CALL exit_error(error_loc,error_msg)
    STOP
  END IF

  ! Zeilen zaehlen zur Validierung
  nraw_lines = 0
  DO
    READ(105,*, IOSTAT=error)
    IF (error < 0) EXIT
    nraw_lines = nraw_lines + 1
  END DO
  REWIND(105)

  nraw_expected = nelem * spalten_elem
  IF (nraw_lines /= nraw_expected) THEN
    error_loc='raw_stress_lines'
    WRITE(error_num1,'(I10)') nraw_lines
    WRITE(error_num2,'(I10)') nraw_expected
    error_msg='RawStress hat '//TRIM(ADJUSTL(error_num1))//' Zeilen, erwartet '//TRIM(ADJUSTL(error_num2))//' (nelem*spalten_elem).'
    CALL exit_error(error_loc,error_msg)
    STOP
  END IF

  WRITE(*,'(A28,I10)') ',.........nraw_lines = ', nraw_lines

  ALLOCATE (raw_stress(nelem, spalten_elem, 3))
  DO e = 1, nelem
    DO n = 1, spalten_elem
      READ(105,*) raw_stress(e, n, 1), raw_stress(e, n, 2), raw_stress(e, n, 3)
    END DO
  END DO
  CLOSE(105)

  PRINT *, '........................................ Done'

END IF

! ==================================================================================
! Singularity Check + Stress Preparation
! ==================================================================================

PRINT *, 'Priming node mask Arrays / Singularity Check'

IF (stress_mode == 0) THEN
  ! ----- NOD-Modus: Bestehende Logik mit Node-Masken -----
  ALLOCATE (elem_nodes(nnode), surf_nodes(nnode))
  elem_nodes = .FALSE.
  surf_nodes = .FALSE.

  PRINT *, '........................................ Check'
  PRINT *, 'Element SigMax'

  singularity_check_nod: IF (smax_file>0.) THEN

    loop_smax_el1: DO i=1,spalten_elem
      loop_smax_el2: DO j=1,nelem
        IF (elements(i,j)/=0) elem_nodes(elements(i,j))=.TRUE.
      END DO loop_smax_el2
    END DO loop_smax_el1

    smax_elem=MAXVAL(nodes(:,4),elem_nodes)

    IF ((smax_elem<smax_file*0.95) .OR. (smax_elem>smax_file*1.05)) THEN
      error_loc='singularity_check'
      WRITE (error_num1,'(ES10.2)') smax_file
      WRITE (error_num2,'(ES10.2)') smax_elem
      error_msg='Die im File gefundene maximale Element-Spannung ('
      error_msg=TRIM(error_msg)//TRIM(error_num1)
      error_msg=TRIM(error_msg)//') unterscheidet sich um mehr als 5% von der gegebenen Maximalspannung ('
      error_msg=TRIM(error_msg)//TRIM(error_num2)
      error_msg=TRIM(error_msg)//').'
      CALL exit_error(error_loc,error_msg)
    END IF

    PRINT *, '........................................ Check'
    PRINT *, 'Face SigMax'

    loop_smax_sf1: DO i=1,spalten_face
      loop_smax_sf2: DO j=1,nface
        IF (faces(i,j)/=0) surf_nodes(faces(i,j))=.TRUE.
      END DO loop_smax_sf2
    END DO loop_smax_sf1

    smax_surf=MAXVAL(nodes(:,4),surf_nodes)

    IF ((smax_surf<smax_file*0.95) .OR. (smax_surf>smax_file*1.05)) THEN
      error_loc='singularity_check'
      WRITE (error_num1,'(ES10.2)') smax_file
      WRITE (error_num2,'(ES10.2)') smax_surf
      error_msg='Die im File gefundene maximale Oberflaechen-Spannung ('
      error_msg=TRIM(error_msg)//TRIM(error_num1)
      error_msg=TRIM(error_msg)//') unterscheidet sich um mehr als 5% von der gegebenen Maximalspannung ('
      error_msg=TRIM(error_msg)//TRIM(error_num2)//').'
      CALL exit_error(error_loc,error_msg)
    END IF

    PRINT *, '......................................... Check'

  ELSE

    DO i=1,spalten_elem
      DO j=1,nelem
        IF (elements(i,j)/=0) elem_nodes(elements(i,j))=.TRUE.
      END DO
    END DO
    smax_file=MAXVAL(nodes(:,4),elem_nodes)
    PRINT *, '........ Used max stress from node-file'

  END IF singularity_check_nod

  ! Negative Spannungen verwerfen
  PRINT *, 'Discarding negative stresses'
  WHERE (nodes(:,4:6)<0)
        nodes(:,4:6)=0.0
  END WHERE

  ! Normieren
  PRINT *, '........................................ Check'
  PRINT *, 'Scaling stresses'
  nodes(:,4:6)=nodes(:,4:6)/smax_file

ELSE
  ! ----- RAW-Modus: Vereinfachter Singularity Check -----
  PRINT *, 'RAW Singularity Check'

  singularity_check_raw: IF (smax_file > 0.) THEN
    smax_elem = MAXVAL(raw_stress(:,:,1))

    IF ((smax_elem<smax_file*0.95) .OR. (smax_elem>smax_file*1.05)) THEN
      error_loc='singularity_check_raw'
      WRITE (error_num1,'(ES10.2)') smax_file
      WRITE (error_num2,'(ES10.2)') smax_elem
      error_msg='RAW: Max. Elementspannung ('
      error_msg=TRIM(error_msg)//TRIM(error_num2)
      error_msg=TRIM(error_msg)//') weicht um >5% von smax ('
      error_msg=TRIM(error_msg)//TRIM(error_num1)//').'
      CALL exit_error(error_loc,error_msg)
    END IF

    PRINT *, '......................................... Check'

  ELSE
    smax_file = MAXVAL(raw_stress(:,:,1))
    PRINT *, '........ Used max stress from raw-stress-file'

  END IF singularity_check_raw

  ! Negative Spannungen verwerfen
  PRINT *, 'Discarding negative stresses (RAW)'
  WHERE (raw_stress < 0)
        raw_stress = 0.0
  END WHERE

  ! Normieren
  PRINT *, '........................................ Check'
  PRINT *, 'Scaling stresses (RAW)'
  raw_stress = raw_stress / smax_file

END IF

! DEBUGGING
PRINT *, '........................................ Check'
PRINT *, 'Allocating Gausz-Arrays'

! Gausz-Integrationspunkte zuteilen
ALLOCATE (gausz_r(gauss_order),gausz_w(gauss_order))

gausz_zuweisung: SELECT CASE (gauss_order)
  CASE (1)
    gausz_r(1) = 0.0

    gausz_w(1) = 2.0

  CASE (2)
    gausz_r(1) = -0.5773502691896
    gausz_r(2) =  0.5773502691896

    gausz_w(1) =  1.0
    gausz_w(2) =  1.0

  CASE (3)
    gausz_r(1) = -0.7745966692414
    gausz_r(2) =  0.0
    gausz_r(3) =  0.7745966692414

    gausz_w(1) =  0.5555555555555
    gausz_w(2) =  0.8888888888888
    gausz_w(3) =  0.5555555555555

  CASE (4)
    gausz_r(1) = -0.8611363115940
    gausz_r(2) = -0.3399810435848
    gausz_r(3) =  0.3399810435848
    gausz_r(4) =  0.8611363115940

    gausz_w(1) =  0.3478548451374
    gausz_w(2) =  0.6521451548625
    gausz_w(3) =  0.6521451548625
    gausz_w(4) =  0.3478548451374

  CASE (5)
    gausz_r(1) = -0.9061798459386
    gausz_r(2) = -0.5384693101056
    gausz_r(3) =  0.0
    gausz_r(4) =  0.5384693101056
    gausz_r(5) =  0.9061798459386

    gausz_w(1) =  0.2369268850561
    gausz_w(2) =  0.4786286704993
    gausz_w(3) =  0.5688888888888
    gausz_w(4) =  0.4786286704993
    gausz_w(5) =  0.2369268850561

  CASE (6)
    gausz_r(1) = -0.9324695142031
    gausz_r(2) = -0.6612093864662
    gausz_r(3) = -0.2386191860831
    gausz_r(4) =  0.2386191860831
    gausz_r(5) =  0.6612093864662
    gausz_r(6) =  0.9324695142031

    gausz_w(1) =  0.1713244923791
    gausz_w(2) =  0.3607615730481
    gausz_w(3) =  0.4679139345726
    gausz_w(4) =  0.4679139345726
    gausz_w(5) =  0.3607615730481
    gausz_w(6) =  0.1713244923791

  CASE (7)
    gausz_r(1) = -0.9491079123427
    gausz_r(2) = -0.7415311855993
    gausz_r(3) = -0.4058451513773
    gausz_r(4) =  0.0
    gausz_r(5) =  0.4058451513773
    gausz_r(6) =  0.7415311855993
    gausz_r(7) =  0.9491079123427

    gausz_w(1) =  0.1294849661688
    gausz_w(2) =  0.2797053914892
    gausz_w(3) =  0.3818300505051
    gausz_w(4) =  0.4179591836734
    gausz_w(5) =  0.3818300505051
    gausz_w(6) =  0.2797053914892
    gausz_w(7) =  0.1294849661688

  CASE (8)
    gausz_r(1) = -0.9602898564975
    gausz_r(2) = -0.7966664774136
    gausz_r(3) = -0.5255324099163
    gausz_r(4) = -0.1834346424956
    gausz_r(5) =  0.1834346424956
    gausz_r(6) =  0.5255324099163
    gausz_r(7) =  0.7966664774136
    gausz_r(8) =  0.9602898564975

    gausz_w(1) =  0.1012285362903
    gausz_w(2) =  0.2223810344533
    gausz_w(3) =  0.3137066458778
    gausz_w(4) =  0.3626837833783
    gausz_w(5) =  0.3626837833783
    gausz_w(6) =  0.3137066458778
    gausz_w(7) =  0.2223810344533
    gausz_w(8) =  0.1012285362903

  CASE (9)
    gausz_r(1) = -0.9681602395076
    gausz_r(2) = -0.8360311073266
    gausz_r(3) = -0.6133714327005
    gausz_r(4) = -0.3242534234038
    gausz_r(5) =  0.0
    gausz_r(6) =  0.3242534234038
    gausz_r(7) =  0.6133714327005
    gausz_r(8) =  0.8360311073266
    gausz_r(9) =  0.9681602395076

    gausz_w(1) =  0.0812743883615
    gausz_w(2) =  0.1806481606948
    gausz_w(3) =  0.2606106964029
    gausz_w(4) =  0.3123470770400
    gausz_w(5) =  0.3302393550012
    gausz_w(6) =  0.3123470770400
    gausz_w(7) =  0.2606106964029
    gausz_w(8) =  0.1806481606948
    gausz_w(9) =  0.0812743883615

  CASE (26)
    gausz_r(1) = -0.9958857011456169194829613
    gausz_r(2) = -0.9783854459564709227237245
    gausz_r(3) = -0.9471590666617142328931322
    gausz_r(4) = -0.9026378619843070660877515
    gausz_r(5) = -0.8454459427884979394463016
    gausz_r(6) = -0.7763859488206789061237600
    gausz_r(7) = -0.6964272604199572835881327
    gausz_r(8) = -0.6066922930176180672745545
    gausz_r(9) = -0.5084407148245057017632575
    gausz_r(10) = -0.4030517551234863438125444
    gausz_r(11) = -0.2920048394859569018677803
    gausz_r(12) = -0.1768588203568901839890515
    gausz_r(13) = -0.0592300934293132075314503
    gausz_r(14) =  0.0592300934293132075314503
    gausz_r(15) =  0.1768588203568901839890515
    gausz_r(16) =  0.2920048394859569018677803
    gausz_r(17) =  0.4030517551234863438125444
    gausz_r(18) =  0.5084407148245057017632575
    gausz_r(19) =  0.6066922930176180672745545
    gausz_r(20) =  0.6964272604199572835881327
    gausz_r(21) =  0.7763859488206789061237600
    gausz_r(22) =  0.8454459427884979394463016
    gausz_r(23) =  0.9026378619843070660877515
    gausz_r(24) =  0.9471590666617142328931322
    gausz_r(25) =  0.9783854459564709227237245
    gausz_r(26) =  0.9958857011456169194829613

    gausz_w(1) =  0.0105513726173433759064624
    gausz_w(2) =  0.0244178510926318440010796
    gausz_w(3) =  0.0379623832943629946345965
    gausz_w(4) =  0.0509758252971477948678469
    gausz_w(5) =  0.0632740463295750482641822
    gausz_w(6) =  0.0746841497656595687537617
    gausz_w(7) =  0.0850458943134851791390005
    gausz_w(8) =  0.0942138003559141040677005
    gausz_w(9) =  0.1020591610944254074011539
    gausz_w(10) =  0.1084718405285764880607857
    gausz_w(11) =  0.1133618165463196464370910
    gausz_w(12) =  0.1166604434852965138658121
    gausz_w(13) =  0.1183214152792621959298103
    gausz_w(14) =  0.1183214152792621959298103
    gausz_w(15) =  0.1166604434852965138658121
    gausz_w(16) =  0.1133618165463196464370910
    gausz_w(17) =  0.1084718405285764880607857
    gausz_w(18) =  0.1020591610944254074011539
    gausz_w(19) =  0.0942138003559141040677005
    gausz_w(20) =  0.0850458943134851791390005
    gausz_w(21) =  0.0746841497656595687537617
    gausz_w(22) =  0.0632740463295750482641822
    gausz_w(23) =  0.0509758252971477948678469
    gausz_w(24) =  0.0379623832943629946345965
    gausz_w(25) =  0.0244178510926318440010796
    gausz_w(26) =  0.0105513726173433759064624

  CASE DEFAULT
    error_loc='gausz_zuweisung'
    error_msg="Fehler beim Zuweisen der Gausz'schen Integrationspunkte und Gewichtungsfaktoren. ueberpruefen Sie das Parameter-File."
    CALL exit_error(error_loc,error_msg)
END SELECT gausz_zuweisung

! DEBUGGING
PRINT *, '........................................ Check'
PRINT *, ''
PRINT *, '###############################################'
PRINT *, ''
PRINT *, 'MAIN ROUTINE:'
PRINT *, ''

ALLOCATE (veff(2,0:m_max),seff(2,0:m_max))

! V12.x Unified: GaussNorm-Tracking initialisieren (alle Domains)
IF (do_gaussnorm) THEN
  max_sigma_ratio_volume_global  = 0.0
  max_sigma_ratio_surface_global = 0.0
  max_sigma_ratio_line_global    = 0.0
  PRINT *, 'GaussNorm-Tracking aktiviert (domain-spezifisch: VOL/SURF/LINE).'
END IF

! #################################################################################################
! Hauptteil: Effektive Volumen und Flaechen errechnen
! #################################################################################################

SELECT CASE (element_type)
  CASE (24)
    veff=seff_gausz(elements, nelem, m_max, breite_2D)
    seff=leff_gausz(faces,     nface, m_max, breite_2D)
  CASE (34)
    veff=veff_gausz(elements, nelem, m_max)
    seff=seff_gausz(faces,    nface, m_max, 1.)
END SELECT

! #################################################################################################
! V12.x Unified: Domain-spezifische GaussNorm-Korrektur
! #################################################################################################
! Methodik:
!   element_type=24 (2D):  veff-Spalten ← seff_gausz → korrigiere mit SURFACE-Maximum
!                          seff-Spalten ← leff_gausz → korrigiere mit LINE-Maximum
!   element_type=34 (3D):  veff-Spalten ← veff_gausz → korrigiere mit VOLUME-Maximum
!                          seff-Spalten ← seff_gausz → korrigiere mit SURFACE-Maximum
!
! veff und seff werden also unabhaengig korrigiert mit den Maxima ihrer eigenen
! Integrationsdomaene — kein gemeinsames globales Maximum mehr.

IF (do_gaussnorm) THEN
  PRINT *, ''
  PRINT *, '###############################################'
  PRINT *, 'Apply Gauss-Normalization (domain-specific):'
  PRINT *, ''

  WRITE (*,'(A35,ES16.8)') 'Max. Nodal Stress (alt):       ', smax_file
  WRITE (*,'(A35,ES16.8)') 'Ratio Volume:                  ', max_sigma_ratio_volume_global
  WRITE (*,'(A35,ES16.8)') 'Ratio Surface:                 ', max_sigma_ratio_surface_global
  WRITE (*,'(A35,ES16.8)') 'Ratio Line:                    ', max_sigma_ratio_line_global

  DO m=0,m_max
    ! ----- veff-Spalten korrigieren -----
    IF (element_type == 34) THEN
      ! 3D: veff = veff_gausz(elements) → Volume-Maximum
      IF (max_sigma_ratio_volume_global > 0.0) THEN
        korrektur_faktor = (1.0 / max_sigma_ratio_volume_global)**m
      ELSE
        korrektur_faktor = 1.0
      END IF
    ELSE
      ! 2D (element_type=24): veff = seff_gausz(elements) → Surface-Maximum
      IF (max_sigma_ratio_surface_global > 0.0) THEN
        korrektur_faktor = (1.0 / max_sigma_ratio_surface_global)**m
      ELSE
        korrektur_faktor = 1.0
      END IF
    END IF
    veff(:,m) = veff(:,m) * korrektur_faktor

    ! ----- seff-Spalten korrigieren -----
    IF (element_type == 34) THEN
      ! 3D: seff = seff_gausz(faces) → Surface-Maximum
      IF (max_sigma_ratio_surface_global > 0.0) THEN
        korrektur_faktor = (1.0 / max_sigma_ratio_surface_global)**m
      ELSE
        korrektur_faktor = 1.0
      END IF
    ELSE
      ! 2D (element_type=24): seff = leff_gausz(faces) → Line-Maximum
      IF (max_sigma_ratio_line_global > 0.0) THEN
        korrektur_faktor = (1.0 / max_sigma_ratio_line_global)**m
      ELSE
        korrektur_faktor = 1.0
      END IF
    END IF
    seff(:,m) = seff(:,m) * korrektur_faktor
  END DO

  PRINT *, 'Domain-specific Normalization Done.'
END IF

! Console-Output
WRITE (*,'(/,A8,4A16,/)') 'm', 'S1-V_eff', 'PIA-V_eff', 'S1-S_eff', 'PIA-S_eff'
DO m=0,m_max
  WRITE (*,'(I8,4ES16.8)') m, veff(:,m), seff(:,m)
END DO

! #################################################################################################
! Nachbereitung: Ausgabe des Ergebnis-files
! #################################################################################################

ausgabename=ADJUSTL(ausgabename)

! 1. Haupt-Output File (identisches Format wie HighOrder/Standard)
OPEN (201, FILE=TRIM(ausgabename)//'.out', STATUS='replace')
WRITE (201,'(A8,ES16.8)') 'SMAX = ',  smax_file
WRITE (201,'(A8,I16)')    'GAUSS = ', gauss_order
WRITE (201,'(A8,I16)')    'MMAX  = ', m_max
WRITE (201,'(A8,I16)')    'SYMM  = ', symm_fakt
WRITE (201,'(A8,I16)')    'GEOM  = ', element_type
WRITE (201,'(A8,F16.8)')  'BREIT = ', breite_2D

WRITE (201,'(/,A8,4A16,/)') 'm', 'S1-V_eff', 'PIA-V_eff', 'S1-S_eff', 'PIA-S_eff'

DO i=0,m_max
  WRITE (201,'(I8,4ES16.8)') i, veff(:,i), seff(:,i)
END DO
CLOSE(201)

! 2. v13.1 Unified: GaussNorm-Info File — nur domain-spezifische Felder
!    Legacy-Felder SMAX_GAUSS / RATIO_GAUSS (globales Maximum aller drei Domains)
!    wurden in v13.1 entfernt — nicht mehr methodisch sinnvoll, weil die Korrektur
!    seit v13.0 domain-spezifisch erfolgt.
IF (do_gaussnorm) THEN
  OPEN (202, FILE=TRIM(ausgabename)//'_GauszInfo.out', STATUS='replace')
  ! --- Metadaten ---
  WRITE (202,'(A20,ES16.8)') 'SMAX_NODAL      = ',  smax_file
  WRITE (202,'(A20,I16)')    'GAUSS_FAKT      = ',  gauss_order
  WRITE (202,'(A20,I16)')    'MMAX            = ',  m_max
  WRITE (202,'(A20,I16)')    'SYMM            = ',  symm_fakt
  WRITE (202,'(A20,I16)')    'GEOM            = ',  element_type
  WRITE (202,'(A20,F16.8)')  'BREITE          = ',  breite_2D
  ! --- Domain-spezifische Maxima ---
  WRITE (202,'(A20,ES16.8)') 'SMAX_GAUSS_VOL  = ',  max_sigma_ratio_volume_global  * smax_file
  WRITE (202,'(A20,ES16.8)') 'SMAX_GAUSS_SURF = ',  max_sigma_ratio_surface_global * smax_file
  WRITE (202,'(A20,ES16.8)') 'SMAX_GAUSS_LINE = ',  max_sigma_ratio_line_global    * smax_file
  WRITE (202,'(A20,ES16.8)') 'RATIO_GAUSS_VOL = ',  max_sigma_ratio_volume_global
  WRITE (202,'(A20,ES16.8)') 'RATIO_GAUSS_SURF= ',  max_sigma_ratio_surface_global
  WRITE (202,'(A20,ES16.8)') 'RATIO_GAUSS_LINE= ',  max_sigma_ratio_line_global
  CLOSE(202)
END IF

! Aufraeumen
DEALLOCATE (elements, faces, gausz_r, gausz_w)
IF (ALLOCATED(nodes)) DEALLOCATE(nodes)
IF (ALLOCATED(elem_nodes)) DEALLOCATE(elem_nodes)
IF (ALLOCATED(surf_nodes)) DEALLOCATE(surf_nodes)
IF (ALLOCATED(nodes_coords)) DEALLOCATE(nodes_coords)
IF (ALLOCATED(raw_stress)) DEALLOCATE(raw_stress)
IF (ALLOCATED(face_parent_elem)) DEALLOCATE(face_parent_elem)

END PROGRAM effektivesVol
