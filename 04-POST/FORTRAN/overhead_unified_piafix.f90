MODULE overhead_unified_piafix
!
! Module overhead_unified_piafix — v6.0 PIAFix-Variante (Stand: 2026-05-01).
! Identisch zu overhead_unified, AUSSER dass die PIA-Auswertung am Gauss-Punkt
! statt an Knoten passiert (siehe Detail unten).
!
! Features (identisch zu Legacy):
!   - stress_mode: 0=NOD (knotengemittelt), 1=RAW (elementlokal aus PRESOL)
!   - do_gaussnorm: optionales GaussNorm-Tracking
!   - Domain-spezifische GaussNorm-Maxima (VOL/SURF/LINE getrennt)
!   - Gauss-Quadratur-Ordnungen 1-9 und 26 (tabelliert)
!
! v6.0 PIAFix Aenderung gegenueber Legacy overhead_unified v5.0:
!   - Legacy: PIA wird an Knoten aggregiert (sigma_eq_node = (S1^m + S2^m + S3^m)^(1/m)),
!             dann zu Gauss-Punkten interpoliert, dann ^m integriert.
!   - PIAFix: S1, S2, S3 werden separat zu Gauss-Punkten interpoliert,
!             dann lokal am Gauss-Punkt ausgewertet:
!               pia_integrand = MAX(s1,0)^m + MAX(s2,0)^m + MAX(s3,0)^m
!             Tensile cutoff direkt am Gauss-Punkt; keine 1/m-Wurzel-Aggregation mehr.
!   - m=0: beide Pfade integrieren das geometrische Mass (User-Spec).
!   - Limitation: keine volle Tensor-Interpolation (Hauptspannungs-Rotation
!     zwischen Knoten wird nicht beruecksichtigt). v20+ Scope.
!
! Vollstaendige Versions-Historie + Architektur-Details:
!   siehe HISTORY.md im selben Verzeichnis
!
IMPLICIT none
SAVE

! --- Bestehende Arrays (NOD-Modus) ---
integer, dimension(:,:), allocatable :: elements, faces
! 2D Arrays zum Verknuepfen von Knoten zu Elementen bzw. Flaechen
! Laenge des Arrays ist die Anzahl der Nodes pro Element (bzw. Flaeche)
! Tiefe ist die Anzahl der Elemente (bzw. Flaechen), die zu berechnen sind.
! Diese Reihenfolge ermoeglicht leichteres Einlesen und besseren Zugriff,
! da die Array-Elemente mit gleichem 1. Index im Speicher gruppiert sind.
real, dimension(:,:), allocatable :: nodes
! 2D Array zum Halten der Knoteninformationen (NOD-Modus)
! Laenge ist die Anzahl der Knoten
! Tiefe ist 6: 3 Raumkoordinaten + 3 Hauptnormalspannungen
real, dimension(:), allocatable :: gausz_r, gausz_w
! Arrays mit den Abszissen und Gewichten der Gausz-Integration
integer :: gauss_order
! Grad der Gauszintegration

! --- V9 Erweiterungen ---
integer :: stress_mode = 0
! 0=NOD: Spannungen knotengemittelt in nodes(:,4:6) (wie bisher)
! 1=RAW: Spannungen elementlokal in raw_stress(:,:,:)
logical :: do_gaussnorm = .FALSE.
! GaussNorm-Tracking: Speichert max. interpolierte Gausspunkt-Spannung

! --- v22.0 GP_DEBUG Erweiterung ---
logical :: do_gp_debug = .FALSE.
! GP_DEBUG-Modus: Dumpt die Fortran-Shape-Function-interpolierten GP-Werte
! (sigma1/2/3 an jedem Gauss-Punkt jedes Elements) nach 'gp_interp_debug.out'.
! Wird von Python mit ANSYS-GP-Wahrheit (aus -S-GP.out via ERESX,NO) verglichen.
! Aktivierung via Parameter `GPDebug:` in effVol_Parameter.out.
real :: smax_file_module = 1.0
! Normalisierungs-Wert (= smax_file aus dem Driver). Wird vom Driver vor dem
! Aufruf von veff_gausz/seff_gausz gesetzt, damit der GP_DEBUG-Dump die
! unnormalisierten Werte zurueckrechnen kann (sigma_norm * smax_file_module).

! --- V12.x Unified: Domain-spezifische GaussNorm-Maxima ---
real :: max_sigma_ratio_volume_global  = 0.0
! Maximales Verhaeltnis interpolierter S1 / smax_nodal in der Volumen-Integration (veff_gausz)
real :: max_sigma_ratio_surface_global = 0.0
! Maximales Verhaeltnis interpolierter S1 / smax_nodal in der Flaechen-Integration (seff_gausz)
real :: max_sigma_ratio_line_global    = 0.0
! Maximales Verhaeltnis interpolierter S1 / smax_nodal in der Linien-Integration (leff_gausz)

real, dimension(:,:), allocatable :: nodes_coords
! (nnode, 3) — Nur Raumkoordinaten fuer RAW-Modus
real, dimension(:,:,:), allocatable :: raw_stress
! (nelem_vol, eckn_vol, 3) — Elementlokale Hauptspannungen (S1,S2,S3)
! Zeilen korrespondieren 1:1 zu den sequentiellen Elementindizes in elements(:,:)
! Spalten korrespondieren zu den Eckknoten in ANSYS-Konnektivitaetsreihenfolge
integer, dimension(:), allocatable :: face_parent_elem
! (nface) — Sequentieller Index des Parent-Volumenelements fuer jede Face
! Wert i verweist auf Element elements(:,i) und raw_stress(i,:,:)


CONTAINS

FUNCTION leff_gausz(nlist, nelem, m_max, breite)

  integer, INTENT(IN) :: m_max
    ! Maximaler Weibull-Modul bis zu dem hinberechnet werden soll
  integer, parameter :: eckn = 2
    ! Anzahl der Ecknodes fuer ein Linienelement
  integer, INTENT(IN) :: nelem
    ! Tiefe der eingehenden Elementliste nlist
  integer, dimension(eckn,nelem), INTENT(IN) :: nlist
    ! Zu behandelnde Elementliste
  integer :: m,n,e,i
    ! Diverse Zaehler
  integer, dimension(eckn), parameter :: xin = (/-1, 1/)
    ! Transformierte Koordinaten der Eckknoten
  real, dimension(:,:), allocatable :: form_function_h, dh_dxi, dh_deta, dh_dzeta
  real, dimension(:,:,:), allocatable :: element_nodes
  ! v19.0 PIAFix: sigma2, sigma3 NEU; sigmaV (PIA-Working-Buffer) entfaellt
  real, dimension(:,:), allocatable :: gram, sigma1, sigma2, sigma3, intLen1, intLenV, dx_dxi, dy_dxi, dz_dxi
  real, dimension(2,0:m_max) :: leff_gausz
  real, INTENT(IN) :: breite
    ! Breite der Elemente im 2D Fall

  !/DEBUGGING
  PRINT *, '====================================================='
  PRINT *, ''
  PRINT *, 'Linien Gausz-Quadratur (PIAFix-Variante)'
  PRINT *, ''
  PRINT *, 'Array Allocation'

  ALLOCATE( form_function_h(eckn,gauss_order), dh_dxi(eckn,gauss_order), dh_deta(eckn,gauss_order), dh_dzeta(eckn,gauss_order) )
  ALLOCATE( element_nodes(nelem,eckn,7) )
  ! v19.0 PIAFix: sigma2, sigma3 statt sigmaV
  ALLOCATE( gram(nelem,gauss_order),  sigma1(nelem,gauss_order), sigma2(nelem,gauss_order), sigma3(nelem,gauss_order), intLen1(nelem,gauss_order), intLenV(nelem,gauss_order) )
  ALLOCATE( dx_dxi(nelem,gauss_order), dy_dxi(nelem,gauss_order), dz_dxi(nelem,gauss_order) )


  PRINT *, '................................................. Done'
  PRINT *, 'Ansatzfunktion:'

  DO n = 1, eckn
  DO i = 1, gauss_order

        form_function_h(n,i) = (1 + xin(n)*gausz_r(i))/2.
        dh_dxi(n,i)   = (      xin(n)          )/2.

END DO
END DO

!/DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'Element-Node-Array:'

! V9: RAW-Modus fuer Linien nicht unterstuetzt (nicht im Workflow verwendet)
IF (stress_mode == 1) THEN
  PRINT *, '!!! WARNUNG: leff_gausz im RAW-Modus nicht unterstuetzt.'
  PRINT *, '    Verwende Fallback auf NOD-Logik (nodes-Array).'
END IF

DO e = 1, nelem
  IF (nlist(1,e)/=0) THEN
    DO n = 1, eckn
      element_nodes(e,n,:6) = nodes(nlist(n,e),:)
    END DO
  ELSE
    element_nodes(e,:,:) = 0.
  END IF
END DO

!/DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'Ableitungen der Ansatzfunktion:'

dx_dxi = 0
dy_dxi = 0
dz_dxi = 0

DO i = 1, gauss_order
  DO n = 1, eckn

    dx_dxi(:,i) = dx_dxi(:,i) + element_nodes(:,n,1)*dh_dxi(n,i)
    dy_dxi(:,i) = dy_dxi(:,i) + element_nodes(:,n,2)*dh_dxi(n,i)
    dz_dxi(:,i) = dz_dxi(:,i) + element_nodes(:,n,3)*dh_dxi(n,i)

  END DO
END DO

!/DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'Gramsche Determinante:'

gram = dx_dxi**2 + dy_dxi**2 + dz_dxi**2
gram = SQRT(gram)

!/DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'Spannung in Integrationspunkten:'

! v19.0 PIAFix: S1, S2, S3 separat zu Gauss-Punkten interpolieren (vor m-Loop)
sigma1 = 0
sigma2 = 0
sigma3 = 0

DO i = 1, gauss_order
  gram(:,i) = gram(:,i) * gausz_w(i)
  DO n = 1, eckn
    sigma1(:,i) = sigma1(:,i) + element_nodes(:,n,4)*form_function_h(n,i)
    sigma2(:,i) = sigma2(:,i) + element_nodes(:,n,5)*form_function_h(n,i)
    sigma3(:,i) = sigma3(:,i) + element_nodes(:,n,6)*form_function_h(n,i)
  END DO
END DO

! V12.x Unified: GaussNorm-Tracking (domain-spezifisch + global Legacy)
! In leff_gausz (1D Linie) → max_sigma_ratio_line_global
! v19.0 PIAFix: Tracking bleibt auf S1 (User-Spec)
IF (do_gaussnorm) THEN
  max_sigma_ratio_line_global = MAX(max_sigma_ratio_line_global, MAXVAL(sigma1))
END IF

!/DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, '--------------------------------------------------------,'
PRINT *, 'M-Loop (PIAFix: PIA at Gauss-Point):'

DO m = 0, m_max

      WRITE(*,'(A28,I3,A12,I3,A1,A1,$)') 'Linien Gausz-Quadratur , m = ', m,', Progress: ', (m*100)/m_max, '%', achar(13)

      ! v19.0 PIAFix: m=0 Special-Case — geometrisches Mass (User-Spec).
      ! Bei m>0: S1-Pfad unveraendert; PIA-Pfad lokal am Gauss-Punkt mit
      ! tensile cutoff (MAX(...,0)) und Sum-of-Powers (keine 1/m-Wurzel mehr).
      IF (m == 0) THEN
        intLen1 = gram   ! geometrisches Mass
        intLenV = gram   ! geometrisches Mass (User-Spec)
      ELSE
        ! v20.2: S1-Pfad mit tensile cutoff am GP (nach Entfernung Pre-Clipping im Driver)
        intLen1 = MAX(sigma1, 0.0)**m * gram   ! S1-Pfad: tensile cutoff
        intLenV = ( MAX(sigma1, 0.0)**m  &
                  + MAX(sigma2, 0.0)**m  &
                  + MAX(sigma3, 0.0)**m ) * gram   ! PIA-Pfad: PIAFix-Variante
      END IF

      leff_gausz(1,m) = SUM(intLen1)
      leff_gausz(2,m) = SUM(intLenV)

END DO
!/DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'deallocate Arrays:'
DEALLOCATE(form_function_h, dh_dxi, dh_deta, dh_dzeta, element_nodes, gram, sigma1, sigma2, sigma3, intLen1, intLenV, dx_dxi, dy_dxi, dz_dxi)
PRINT *, '...................................................... Done'

END FUNCTION


FUNCTION seff_gausz(nlist, nelem, m_max, breite)

  integer, INTENT(IN) :: m_max
  integer, parameter  :: eckn = 4
  integer, INTENT(IN) :: nelem
  integer, dimension(eckn,nelem), INTENT(IN) :: nlist
  integer :: m,n,e,i,j
  ! V9: Hilfsvariablen fuer Parent-Element-Lookup
  integer :: nn, parent_e
  ! v22.0: GP_DEBUG — Existenz-Check fuer Append-vs-Replace bei gp_interp_debug.out
  logical :: do_gp_debug_dump_exists

  integer, dimension(eckn), parameter :: xin  = (/-1,  1,  1, -1/)
  integer, dimension(eckn), parameter :: etan = (/-1, -1,  1,  1/)

  real, dimension(:,:,:),   allocatable :: form_function_h, dh_dxi, dh_deta, dh_dzeta
  real, dimension(:,:,:),   allocatable :: element_nodes
  ! v19.0 PIAFix: sigma2, sigma3 NEU; sigmaV (PIA-Working-Buffer) entfaellt
  real, dimension(:,:,:),   allocatable :: gram, sigma1, sigma2, sigma3, intSurf1, intSurfV
  real, dimension(:,:,:),   allocatable :: dx_dxi, dx_deta, dy_dxi, dy_deta, dz_dxi, dz_deta

  real, dimension(2,0:m_max) :: seff_gausz
  real, INTENT(IN) :: breite

  !/DEBUGGING
  PRINT *, '====================================================='
  PRINT *, ''
  PRINT *, 'Flaechen Gausz-Quadratur (PIAFix-Variante)'
  PRINT *, ''
  PRINT *, 'Array Allocation'

  ALLOCATE( form_function_h(eckn,gauss_order,gauss_order), dh_dxi (eckn,gauss_order,gauss_order), dh_deta(eckn,gauss_order,gauss_order), dh_dzeta(eckn,gauss_order,gauss_order) )
  ALLOCATE( element_nodes(nelem,eckn,7) )
  ! v19.0 PIAFix: sigma2, sigma3 statt sigmaV
  ALLOCATE( gram    (nelem,gauss_order,gauss_order), sigma1  (nelem,gauss_order,gauss_order), sigma2  (nelem,gauss_order,gauss_order), sigma3  (nelem,gauss_order,gauss_order), intSurf1(nelem,gauss_order,gauss_order), intSurfV(nelem,gauss_order,gauss_order) )
  ALLOCATE( dx_dxi (nelem,gauss_order,gauss_order), dx_deta(nelem,gauss_order,gauss_order), dy_dxi (nelem,gauss_order,gauss_order), dy_deta(nelem,gauss_order,gauss_order) )
  ALLOCATE( dz_dxi (nelem,gauss_order,gauss_order), dz_deta(nelem,gauss_order,gauss_order) )


  PRINT *, '...................................................... Done'
  PRINT *, 'Ansatzfunktion:'

  DO n = 1, eckn
    DO i = 1, gauss_order
      DO j = 1, gauss_order

        form_function_h(n,i,j) = ((1 + xin(n)*gausz_r(i)) * (1 + etan(n)*gausz_r(j)))/4.
        dh_dxi  (n,i,j) = (  xin(n)            * (1 + etan(n)*gausz_r(j))   )/4.
        dh_deta(n,i,j) = ((1 + xin(n)*gausz_r(i)) * etan(n))/4.

END DO
END DO
END DO

!DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'Element-Node-Array:'

! V9: Dual-Mode Element-Node-Array Befuellung
IF (stress_mode == 0) THEN
  ! NOD-Modus: Wie bisher — nodes(:,1:6) enthaelt Coords + Stresses
  DO e = 1, nelem
    IF (nlist(1,e) /= 0) THEN
      DO n = 1, eckn
        element_nodes(e,n,:6) = nodes(nlist(n,e),:)
      END DO
    ELSE
      element_nodes(e,:,:) = 0.
    END IF
  END DO
ELSE
  ! RAW-Modus: Coords aus nodes_coords, Stresses via Parent-Element-Lookup
  DO e = 1, nelem
    IF (nlist(1,e) /= 0) THEN
      parent_e = face_parent_elem(e)
      DO n = 1, eckn
        ! Koordinaten: Global aus nodes_coords
        element_nodes(e,n,1:3) = nodes_coords(nlist(n,e), 1:3)
        ! Spannungen: Finde lokalen Knotenindex im Parent-Volumenelement
        element_nodes(e,n,4:6) = 0.0  ! Fallback falls Knoten nicht gefunden
        DO nn = 1, 8
          IF (elements(nn, parent_e) == nlist(n,e)) THEN
            element_nodes(e,n,4:6) = raw_stress(parent_e, nn, 1:3)
            EXIT
          END IF
        END DO
      END DO
    ELSE
      element_nodes(e,:,:) = 0.
    END IF
  END DO
END IF

!DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'Ableitungen der Ansatzfunktion:'

dx_dxi  = 0.
dx_deta = 0.
dy_dxi  = 0.
dy_deta = 0.
dz_dxi  = 0.
dz_deta = 0.

DO i = 1, gauss_order
  DO j = 1, gauss_order
    DO n = 1, eckn

      dx_dxi(:,i,j)  = dx_dxi(:,i,j)  + element_nodes(:,n,1)*dh_dxi (n,i,j)
      dx_deta(:,i,j) = dx_deta(:,i,j) + element_nodes(:,n,1)*dh_deta(n,i,j)

      dy_dxi(:,i,j)  = dy_dxi(:,i,j)  + element_nodes(:,n,2)*dh_dxi (n,i,j)
      dy_deta(:,i,j) = dy_deta(:,i,j) + element_nodes(:,n,2)*dh_deta(n,i,j)

      dz_dxi(:,i,j)  = dz_dxi(:,i,j)  + element_nodes(:,n,3)*dh_dxi (n,i,j)
      dz_deta(:,i,j) = dz_deta(:,i,j) + element_nodes(:,n,3)*dh_deta(n,i,j)

    END DO
  END DO
END DO

!DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'Gramsche Determinante:'

gram = dx_dxi**2 * (dy_deta**2 + dz_deta**2)
gram = gram + (dy_dxi*dz_deta - dy_deta*dz_dxi)**2
gram = gram - 2*dx_dxi*dx_deta*(dy_dxi*dy_deta + dz_dxi*dz_deta)
gram = gram + dx_deta**2 * (dy_dxi**2 + dz_dxi**2)
gram = SQRT(gram)

!DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'Spannung in Integrationspunkten:'

! v19.0 PIAFix: S1, S2, S3 separat zu Gauss-Punkten interpolieren (vor m-Loop)
sigma1 = 0.
sigma2 = 0.
sigma3 = 0.

DO i = 1, gauss_order
  DO j = 1, gauss_order
    gram(:,i,j) = gram(:,i,j) * gausz_w(i) * gausz_w(j)
    DO n = 1, eckn
      sigma1(:,i,j) = sigma1(:,i,j) + element_nodes(:,n,4)*form_function_h(n,i,j)
      sigma2(:,i,j) = sigma2(:,i,j) + element_nodes(:,n,5)*form_function_h(n,i,j)
      sigma3(:,i,j) = sigma3(:,i,j) + element_nodes(:,n,6)*form_function_h(n,i,j)
    END DO
  END DO
END DO

! V12.x Unified: GaussNorm-Tracking (domain-spezifisch + global Legacy)
! In seff_gausz (2D Flaeche) → max_sigma_ratio_surface_global
! v19.0 PIAFix: Tracking bleibt auf S1 (User-Spec)
IF (do_gaussnorm) THEN
  max_sigma_ratio_surface_global = MAX(max_sigma_ratio_surface_global, MAXVAL(sigma1))
END IF

! ============================================================
! v22.0: GP_DEBUG — Surface-Pfad-Dump (Append-Mode an gp_interp_debug.out)
! ============================================================
! Wird nur ausgeloest wenn veff_gausz vorher NICHT gelaufen ist (sonst doppelt).
! Driver-Konvention: bei element_type=34 ruft veff_gausz auf -> Volume-Dump schreibt;
! bei element_type=24 ruft seff_gausz fuer den Haupt-Veff -> Surface-Dump schreibt.
! Hier: nur dumpen wenn noch keine Datei existiert (gating Driver-seitig sauberer,
! aber dieser Check ist eine Safety-Net).
! ============================================================
IF (do_gp_debug) THEN
  ! Append-only wenn Datei schon vom veff_gausz-Aufruf existiert.
  INQUIRE(FILE='gp_interp_debug.out', EXIST=do_gp_debug_dump_exists)
  IF (do_gp_debug_dump_exists) THEN
    OPEN(302, FILE='gp_interp_debug.out', STATUS='old', POSITION='append')
    WRITE(302,'(A)') '# --- seff_gausz block (Surface 2D Gauss-Points) ---'
  ELSE
    OPEN(302, FILE='gp_interp_debug.out', STATUS='replace')
    WRITE(302,'(A)') '# face gp_i gp_j sigma1 sigma2 sigma3 (Surface 2D, seff_gausz)'
  END IF
  DO e = 1, nelem
    DO j = 1, gauss_order
      DO i = 1, gauss_order
        WRITE(302,'(3I8,3ES20.10)') e, i, j, &
             sigma1(e,i,j) * smax_file_module, &
             sigma2(e,i,j) * smax_file_module, &
             sigma3(e,i,j) * smax_file_module
      END DO
    END DO
  END DO
  CLOSE(302)
  PRINT *, 'GP_DEBUG: Surface-Dump abgeschlossen (', nelem, ' Faces x ', &
           gauss_order**2, ' GPs).'
END IF

!DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, '--------------------------------------------------------'
PRINT *, 'M-Loop (PIAFix: PIA at Gauss-Point):'

DO m = 0, m_max

  WRITE(*,'(A28,I3,A12,I3,A1,A1,$)') 'Flaechen Gausz-Quadratur, m = ', m, ', Progress: ', (m*100)/m_max, '%', achar(13)

  ! v19.0 PIAFix: m=0 -> geometrisches Mass; m>0 -> lokale PIA-Auswertung am GP
  IF (m == 0) THEN
    intSurf1 = gram
    intSurfV = gram
  ELSE
    ! v20.2: S1-Pfad mit tensile cutoff am GP
    intSurf1 = MAX(sigma1, 0.0)**m * gram
    intSurfV = ( MAX(sigma1, 0.0)**m  &
               + MAX(sigma2, 0.0)**m  &
               + MAX(sigma3, 0.0)**m ) * gram
  END IF

  seff_gausz(1,m) = SUM(intSurf1)
  seff_gausz(2,m) = SUM(intSurfV)

END DO
!DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'deallocate Arrays:'
DEALLOCATE(form_function_h, dh_dxi, dh_deta, dh_dzeta, element_nodes, gram, sigma1, sigma2, sigma3, intSurf1, intSurfV, dx_dxi, dx_deta, dy_dxi, dy_deta, dz_dxi, dz_deta)
PRINT *, '...................................................... Done'

END FUNCTION

FUNCTION veff_gausz(nlist, nelem, m_max)

  integer, INTENT(IN)            :: m_max
  integer, parameter             :: eckn = 8
  integer, INTENT(IN)            :: nelem
  integer, dimension(:,:), INTENT(IN) :: nlist
  integer                        :: m, n, e, i, j, k
  integer, dimension(8), parameter :: xin   = (/-1, 1, 1, -1, -1, 1, 1, -1/)
  integer, dimension(8), parameter :: etan  = (/-1,-1, 1,  1, -1,-1, 1,  1/)
  integer, dimension(8), parameter :: zetan = (/-1,-1,-1, -1,  1, 1, 1,  1/)
  real, dimension(:,:,:,:), allocatable :: form_function_h, dh_dxi, dh_deta, dh_dzeta
  real, dimension(:,:,:),   allocatable :: element_nodes
  ! v19.0 PIAFix: sigma2, sigma3 NEU; sigmaV (PIA-Working-Buffer) entfaellt
  real, dimension(:,:,:,:), allocatable :: jacobi, sigma1, sigma2, sigma3, intVol1, intVolV
  real, dimension(:,:,:,:), allocatable :: dx_dxi, dx_deta, dx_dzeta, dy_dxi, dy_deta, dy_dzeta
  real, dimension(:,:,:,:), allocatable :: dz_dxi, dz_deta, dz_dzeta
  real, dimension (2,0:m_max) :: veff_gausz

  !DEBUGGING
  PRINT *, '===================================================='
  PRINT *, ''
  PRINT *, 'Volumen Gausz-Quadratur (PIAFix-Variante)'
  PRINT *, ''
  PRINT *, 'Array Allocation'


  ALLOCATE( form_function_h(eckn,gauss_order,gauss_order,gauss_order), dh_dxi (eckn,gauss_order,gauss_order,gauss_order), dh_deta(eckn,gauss_order,gauss_order,gauss_order), dh_dzeta(eckn,gauss_order,gauss_order,gauss_order) )
  ALLOCATE( element_nodes(nelem,eckn,7) )
  ! v19.0 PIAFix: sigma2, sigma3 statt sigmaV
  ALLOCATE( jacobi (nelem,gauss_order,gauss_order,gauss_order), sigma1 (nelem,gauss_order,gauss_order,gauss_order), sigma2 (nelem,gauss_order,gauss_order,gauss_order), sigma3 (nelem,gauss_order,gauss_order,gauss_order), intVol1(nelem,gauss_order,gauss_order,gauss_order), intVolV(nelem,gauss_order,gauss_order,gauss_order) )
  ALLOCATE( dx_dxi  (nelem,gauss_order,gauss_order,gauss_order), dx_deta (nelem,gauss_order,gauss_order,gauss_order), dx_dzeta(nelem,gauss_order,gauss_order,gauss_order), dy_dxi  (nelem,gauss_order,gauss_order,gauss_order), dy_deta (nelem,gauss_order,gauss_order,gauss_order) )
  ALLOCATE( dy_dzeta(nelem,gauss_order,gauss_order,gauss_order), dz_dxi  (nelem,gauss_order,gauss_order,gauss_order), dz_deta (nelem,gauss_order,gauss_order,gauss_order), dz_dzeta(nelem,gauss_order,gauss_order,gauss_order) )


  PRINT *, '...................................................... Done'
  PRINT *, 'Ansatzfunktion:'

  DO n = 1, eckn
    DO i = 1, gauss_order
      DO j = 1, gauss_order
        DO k = 1, gauss_order

          form_function_h(n,i,j,k) = ( (1 + xin(n)*gausz_r(i)) * (1 + etan(n)*gausz_r(j)) * (1 + zetan(n)*gausz_r(k)) )/8.
          dh_dxi (n,i,j,k)  = (  xin(n)              * (1 + etan(n)*gausz_r(j)) * (1 + zetan(n)*gausz_r(k)) )/8.
          dh_deta(n,i,j,k)  = ( (1 + xin(n)*gausz_r(i)) *  etan(n)              *  (1 + zetan(n)*gausz_r(k)) )/8.
          dh_dzeta(n,i,j,k) = ( (1 + xin(n)*gausz_r(i)) * (1 + etan(n)*gausz_r(j)) * zetan(n) )/8.

        END DO
      END DO
    END DO
  END DO

  !/DEBUGGING
  PRINT *, '...................................................... Done'
  PRINT *, 'Element-Node-Array:'

  ! V9: Dual-Mode Element-Node-Array Befuellung
  IF (stress_mode == 0) THEN
    ! NOD-Modus: Wie bisher — nodes(:,1:6) enthaelt Coords + Stresses
    DO e = 1, nelem
      IF (nlist(1,e)/=0) THEN
        DO n = 1, eckn
          element_nodes(e,n,:6) = nodes(nlist(n,e),:)
        END DO
      ELSE
        element_nodes(e,:,:) = 0.
      END IF
    END DO
  ELSE
    ! RAW-Modus: Coords global, Stresses elementlokal
    DO e = 1, nelem
      IF (nlist(1,e)/=0) THEN
        DO n = 1, eckn
          ! Koordinaten: Global aus nodes_coords (ein Wert pro Knoten)
          element_nodes(e,n,1:3) = nodes_coords(nlist(n,e), 1:3)
          ! Spannungen: Direkt aus raw_stress (elementlokal, gleicher Index)
          element_nodes(e,n,4:6) = raw_stress(e, n, 1:3)
        END DO
      ELSE
        element_nodes(e,:,:) = 0.
      END IF
    END DO
  END IF

  !/DEBUGGING
  PRINT *, '...................................................... Done'
  PRINT *, 'Ableitungen der Ansatzfunktion:'

  dx_dxi   = 0.
  dx_deta  = 0.
  dx_dzeta = 0.
  dy_dxi   = 0.
  dy_deta  = 0.
  dy_dzeta = 0.
  dz_dxi   = 0.
  dz_deta  = 0.
  dz_dzeta = 0.

  DO i = 1, gauss_order
  DO j = 1, gauss_order
    DO k = 1, gauss_order
      DO n = 1, eckn

        dx_dxi  (:,i,j,k) = dx_dxi  (:,i,j,k) + element_nodes(:,n,1)*dh_dxi  (n,i,j,k)
        dx_deta (:,i,j,k) = dx_deta (:,i,j,k) + element_nodes(:,n,1)*dh_deta (n,i,j,k)
        dx_dzeta(:,i,j,k) = dx_dzeta(:,i,j,k) + element_nodes(:,n,1)*dh_dzeta(n,i,j,k)

        dy_dxi  (:,i,j,k) = dy_dxi  (:,i,j,k) + element_nodes(:,n,2)*dh_dxi  (n,i,j,k)
        dy_deta (:,i,j,k) = dy_deta (:,i,j,k) + element_nodes(:,n,2)*dh_deta (n,i,j,k)
        dy_dzeta(:,i,j,k) = dy_dzeta(:,i,j,k) + element_nodes(:,n,2)*dh_dzeta(n,i,j,k)

        dz_dxi  (:,i,j,k) = dz_dxi  (:,i,j,k) + element_nodes(:,n,3)*dh_dxi  (n,i,j,k)
        dz_deta (:,i,j,k) = dz_deta (:,i,j,k) + element_nodes(:,n,3)*dh_deta (n,i,j,k)
        dz_dzeta(:,i,j,k) = dz_dzeta(:,i,j,k) + element_nodes(:,n,3)*dh_dzeta(n,i,j,k)

      END DO
    END DO
  END DO
  END DO

  !/DEBUGGING
  PRINT *, '........................................................ Done'
  PRINT *, ' Jacobi Determinante:'

  jacobi = dx_dxi*dy_deta*dz_dzeta + dx_dzeta*dy_dxi*dz_deta + dx_deta*dy_dzeta*dz_dxi
  jacobi = jacobi - dx_dzeta*dy_deta*dz_dxi - dx_dxi*dy_dzeta*dz_deta - dx_deta*dy_dxi*dz_dzeta

  !/DEBUGGING
  PRINT *, '........................................................ Done'
  PRINT *, ' Spannung in Integrationspunkten:'

  ! v19.0 PIAFix: S1, S2, S3 separat zu Gauss-Punkten interpolieren (vor m-Loop)
  sigma1 = 0
  sigma2 = 0
  sigma3 = 0

  DO i = 1, gauss_order
    DO j = 1, gauss_order
      DO k = 1, gauss_order
        jacobi(:,i,j,k) = jacobi(:,i,j,k) * gausz_w(i) * gausz_w(j) * gausz_w(k)
        DO n = 1, eckn
          sigma1(:,i,j,k) = sigma1(:,i,j,k) + element_nodes(:,n,4)*form_function_h(n,i,j,k)
          sigma2(:,i,j,k) = sigma2(:,i,j,k) + element_nodes(:,n,5)*form_function_h(n,i,j,k)
          sigma3(:,i,j,k) = sigma3(:,i,j,k) + element_nodes(:,n,6)*form_function_h(n,i,j,k)
        END DO
      END DO
    END DO
  END DO

  ! V12.x Unified: GaussNorm-Tracking (domain-spezifisch + global Legacy)
  ! In veff_gausz (3D Volumen) → max_sigma_ratio_volume_global
  ! v19.0 PIAFix: Tracking bleibt auf S1 (User-Spec)
  IF (do_gaussnorm) THEN
    max_sigma_ratio_volume_global = MAX(max_sigma_ratio_volume_global, MAXVAL(sigma1))
  END IF

  ! ============================================================
  ! v22.0: GP_DEBUG — Dump der interpolierten GP-Spannungen
  ! ============================================================
  ! sigma1/2/3 sind hier durch /smax_file normalisiert (siehe Driver Z. 478).
  ! Beim Dump multiplizieren wir zurueck, damit Python direkt mit -S-GP-clean.csv
  ! vergleichen kann (gleiche physikalische Einheit).
  !
  ! Format: 7 Spalten — elem, gp_i, gp_j, gp_k, sigma1, sigma2, sigma3
  ! Mapping (i,j,k) -> ANSYS-GP-Reihenfolge erfolgt Python-seitig in
  ! gp_interp_vs_truth.py via Vorzeichen-Muster gausz_r vs (xin,etan,zetan).
  ! ============================================================
  IF (do_gp_debug) THEN
    PRINT *, ''
    PRINT *, 'GP_DEBUG: Dumping interpolated GP stresses to gp_interp_debug.out'
    OPEN(301, FILE='gp_interp_debug.out', STATUS='replace')
    WRITE(301,'(A)') '# elem gp_i gp_j gp_k sigma1 sigma2 sigma3'
    WRITE(301,'(A)') '# (sigmaN = Fortran-Lagrange-interpolierter GP-Wert, unnormalisiert)'
    DO e = 1, nelem
      DO k = 1, gauss_order
        DO j = 1, gauss_order
          DO i = 1, gauss_order
            WRITE(301,'(4I8,3ES20.10)') e, i, j, k, &
                 sigma1(e,i,j,k) * smax_file_module, &
                 sigma2(e,i,j,k) * smax_file_module, &
                 sigma3(e,i,j,k) * smax_file_module
          END DO
        END DO
      END DO
    END DO
    CLOSE(301)
    PRINT *, 'GP_DEBUG: Dump abgeschlossen (', nelem, ' Elemente x ', &
             gauss_order**3, ' GPs).'
  END IF

  !/DEBUGGING
  PRINT *, '........................................................ Done'
  PRINT *, '------------------------------------------------------------'
  PRINT *, 'M-Loop (PIAFix: PIA at Gauss-Point):'

  DO m = 0, m_max

    WRITE(*,'(A28,I3,A12,I3,A1,A1,$)') 'Volumen Gausz-Quadratur , m = ',m,', Progress: ',(m*100)/m_max,'% ',achar(13)

    ! v19.0 PIAFix: m=0 -> geometrisches Mass; m>0 -> lokale PIA-Auswertung am GP
    IF (m == 0) THEN
      intVol1 = jacobi
      intVolV = jacobi
    ELSE
      ! v20.2: S1-Pfad mit tensile cutoff am GP
      intVol1 = MAX(sigma1, 0.0)**m * jacobi
      intVolV = ( MAX(sigma1, 0.0)**m  &
                + MAX(sigma2, 0.0)**m  &
                + MAX(sigma3, 0.0)**m ) * jacobi
    END IF

    veff_gausz(1,m) = SUM(intVol1)
    veff_gausz(2,m) = SUM(intVolV)

END DO
!/DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, ' deallocate Arrays:'
DEALLOCATE( form_function_h, dh_dxi, dh_deta, dh_dzeta, element_nodes, jacobi, sigma1, sigma2, sigma3, intVol1, intVolV, dx_dxi, dx_deta, dx_dzeta, dy_dxi, dy_deta, dy_dzeta, dz_dxi, dz_deta, dz_dzeta )
PRINT *, '...................................................... Done'

END FUNCTION

SUBROUTINE exit_error(location, message)

  character(len=32),  INTENT(IN) :: location
  character(len=256), INTENT(IN) :: message

  WRITE (*,'(A72)')  'Ein Fehler in folgender Routine hat zum Abbrechen des Programms gefuehrt:'
  WRITE (*,'(X)')
  WRITE (*,'(A32)')  location
  WRITE (*,'(A256)') message
  WRITE (*,'(A25)')  'Bitte bestaetigen mit "y":'
END SUBROUTINE

END MODULE overhead_unified_piafix
