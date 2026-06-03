MODULE overhead_unified
!
! Module overhead_unified — Aktive Version v5.0 (Stand: 2026-04-29).
! Vereint die historischen Pfade overhead / _GaussNorm / _highOrder / _RAWNodes.
!
! Features (Kurzfassung):
!   - stress_mode: 0=NOD (knotengemittelt), 1=RAW (elementlokal aus PRESOL)
!   - do_gaussnorm: optionales GaussNorm-Tracking
!   - Domain-spezifische GaussNorm-Maxima (VOL/SURF/LINE getrennt)
!   - Gauss-Quadratur-Ordnungen 1-9 und 26 (tabelliert)
!
! Vollstaendige Versions-Historie + Architektur-Details:
!   siehe HISTORY.md im selben Verzeichnis (v17.0 ausgelagert)
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
  real, dimension(:,:), allocatable :: gram, sigma1, sigmaV, intLen1, intLenV, dx_dxi, dy_dxi, dz_dxi
  real, dimension(2,0:m_max) :: leff_gausz
  real, INTENT(IN) :: breite
    ! Breite der Elemente im 2D Fall

  !/DEBUGGING
  PRINT *, '====================================================='
  PRINT *, ''
  PRINT *, 'Linien Gausz-Quadratur'
  PRINT *, ''
  PRINT *, 'Array Allocation'

  ALLOCATE( form_function_h(eckn,gauss_order), dh_dxi(eckn,gauss_order), dh_deta(eckn,gauss_order), dh_dzeta(eckn,gauss_order) )
  ALLOCATE( element_nodes(nelem,eckn,7) )
  ALLOCATE( gram(nelem,gauss_order),  sigma1(nelem,gauss_order), sigmaV(nelem,gauss_order), intLen1(nelem,gauss_order), intLenV(nelem,gauss_order) )
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

sigma1 = 0

DO i = 1, gauss_order
  gram(:,i) = gram(:,i) * gausz_w(i)
  DO n = 1, eckn
    sigma1(:,i) = sigma1(:,i) + element_nodes(:,n,4)*form_function_h(n,i)
  END DO
END DO

! V12.x Unified: GaussNorm-Tracking (domain-spezifisch + global Legacy)
! In leff_gausz (1D Linie) → max_sigma_ratio_line_global
IF (do_gaussnorm) THEN
  max_sigma_ratio_line_global = MAX(max_sigma_ratio_line_global, MAXVAL(sigma1))
END IF

!/DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, '--------------------------------------------------------,'
PRINT *, 'M-Loop:'

DO m = 0, m_max

      WRITE(*,'(A28,I3,A12,I3,A1,A1,$)') 'Linien Gausz-Quadratur , m = ', m,', Progress: ', (m*100)/m_max, '%', achar(13)

      ! v16 (Audit-Find #5): m=0 Special-Case — `1./REAL(0)` ist Division durch Null.
      ! Bei m=0 ist die Weibull-Formel degeneriert (Pf = 1-exp(-V), kein Spannungs-Effekt);
      ! der PIA-Mischterm `(σ1^m + σ2^m + σ3^m)^(1/m)` ist nicht definiert.
      ! Setze element_nodes(:,:,7) = element_nodes(:,:,4) (S1) als sicherer Default.
      ! In der Pipeline wird m=0 nur als V_total-Container genutzt (nicht in Pf-Berechnung).
      IF (m == 0) THEN
        element_nodes(:,:,7) = element_nodes(:,:,4)
      ELSE
        element_nodes(:,:,7) = ( element_nodes(:,:,4)**m + element_nodes(:,:,5)**m + element_nodes(:,:,6)**m )**(1./REAL(m))
      END IF
      sigmaV = 0

      DO i = 1, gauss_order
      DO n = 1, eckn
        sigmaV(:,i) = sigmaV(:,i) + element_nodes(:,n,7)*form_function_h(n,i)
      END DO
      END DO

      intLen1 = sigma1**m * gram
      intLenV = sigmaV**m * gram

      leff_gausz(1,m) = SUM(intLen1)
      leff_gausz(2,m) = SUM(intLenV)

END DO
!/DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'deallocate Arrays:'
DEALLOCATE(form_function_h, dh_dxi, dh_deta, dh_dzeta, element_nodes, gram, sigma1, sigmaV, intLen1, intLenV, dx_dxi, dy_dxi, dz_dxi)
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

  integer, dimension(eckn), parameter :: xin  = (/-1,  1,  1, -1/)
  integer, dimension(eckn), parameter :: etan = (/-1, -1,  1,  1/)

  real, dimension(:,:,:),   allocatable :: form_function_h, dh_dxi, dh_deta, dh_dzeta
  real, dimension(:,:,:),   allocatable :: element_nodes
  real, dimension(:,:,:),   allocatable :: gram, sigma1, sigmaV, intSurf1, intSurfV
  real, dimension(:,:,:),   allocatable :: dx_dxi, dx_deta, dy_dxi, dy_deta, dz_dxi, dz_deta

  real, dimension(2,0:m_max) :: seff_gausz
  real, INTENT(IN) :: breite

  !/DEBUGGING
  PRINT *, '====================================================='
  PRINT *, ''
  PRINT *, 'Flaechen Gausz-Quadratur'
  PRINT *, ''
  PRINT *, 'Array Allocation'

  ALLOCATE( form_function_h(eckn,gauss_order,gauss_order), dh_dxi (eckn,gauss_order,gauss_order), dh_deta(eckn,gauss_order,gauss_order), dh_dzeta(eckn,gauss_order,gauss_order) )
  ALLOCATE( element_nodes(nelem,eckn,7) )
  ALLOCATE( gram    (nelem,gauss_order,gauss_order), sigma1  (nelem,gauss_order,gauss_order), sigmaV  (nelem,gauss_order,gauss_order), intSurf1(nelem,gauss_order,gauss_order), intSurfV(nelem,gauss_order,gauss_order) )
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

sigma1 = 0.

DO i = 1, gauss_order
  DO j = 1, gauss_order
    gram(:,i,j) = gram(:,i,j) * gausz_w(i) * gausz_w(j)
    DO n = 1, eckn
      sigma1(:,i,j) = sigma1(:,i,j) + element_nodes(:,n,4)*form_function_h(n,i,j)
    END DO
  END DO
END DO

! V12.x Unified: GaussNorm-Tracking (domain-spezifisch + global Legacy)
! In seff_gausz (2D Flaeche) → max_sigma_ratio_surface_global
IF (do_gaussnorm) THEN
  max_sigma_ratio_surface_global = MAX(max_sigma_ratio_surface_global, MAXVAL(sigma1))
END IF

!DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, '--------------------------------------------------------'
PRINT *, 'M-Loop:'

DO m = 0, m_max

  WRITE(*,'(A28,I3,A12,I3,A1,A1,$)') 'Flaechen Gausz-Quadratur, m = ', m, ', Progress: ', (m*100)/m_max, '%', achar(13)

  ! v16 (Audit-Find #5): m=0 Special-Case — siehe leff_gausz Z. 198 fuer Begruendung
  IF (m == 0) THEN
    element_nodes(:,:,7) = element_nodes(:,:,4)
  ELSE
    element_nodes(:,:,7) = ( element_nodes(:,:,4)**m + element_nodes(:,:,5)**m + element_nodes(:,:,6)**m )**(1./REAL(m))
  END IF

  sigmaV = 0.

  DO i = 1, gauss_order
    DO j = 1, gauss_order
      DO n = 1, eckn
        sigmaV(:,i,j) = sigmaV(:,i,j) + element_nodes(:,n,7)*form_function_h(n,i,j)
      END DO
    END DO
  END DO

  intSurf1 = sigma1**m * gram
  intSurfV = sigmaV**m * gram


  seff_gausz(1,m) = SUM(intSurf1)
  seff_gausz(2,m) = SUM(intSurfV)




END DO
!DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, 'deallocate Arrays:'
DEALLOCATE(form_function_h, dh_dxi, dh_deta, dh_dzeta, element_nodes, gram, sigma1, sigmaV, intSurf1, intSurfV, dx_dxi, dx_deta, dy_dxi, dy_deta, dz_dxi, dz_deta)
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
  real, dimension(:,:,:,:), allocatable :: jacobi, sigma1, sigmaV, intVol1, intVolV
  real, dimension(:,:,:,:), allocatable :: dx_dxi, dx_deta, dx_dzeta, dy_dxi, dy_deta, dy_dzeta
  real, dimension(:,:,:,:), allocatable :: dz_dxi, dz_deta, dz_dzeta
  real, dimension (2,0:m_max) :: veff_gausz

  !DEBUGGING
  PRINT *, '===================================================='
  PRINT *, ''
  PRINT *, 'Volumen Gausz-Quadratur'
  PRINT *, ''
  PRINT *, 'Array Allocation'


  ALLOCATE( form_function_h(eckn,gauss_order,gauss_order,gauss_order), dh_dxi (eckn,gauss_order,gauss_order,gauss_order), dh_deta(eckn,gauss_order,gauss_order,gauss_order), dh_dzeta(eckn,gauss_order,gauss_order,gauss_order) )
  ALLOCATE( element_nodes(nelem,eckn,7) )
  ALLOCATE( jacobi (nelem,gauss_order,gauss_order,gauss_order), sigma1 (nelem,gauss_order,gauss_order,gauss_order), sigmaV (nelem,gauss_order,gauss_order,gauss_order), intVol1(nelem,gauss_order,gauss_order,gauss_order), intVolV(nelem,gauss_order,gauss_order,gauss_order) )
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

  sigma1 = 0

  DO i = 1, gauss_order
    DO j = 1, gauss_order
      DO k = 1, gauss_order
        jacobi(:,i,j,k) = jacobi(:,i,j,k) * gausz_w(i) * gausz_w(j) * gausz_w(k)
        DO n = 1, eckn
          sigma1(:,i,j,k) = sigma1(:,i,j,k) + element_nodes(:,n,4)*form_function_h(n,i,j,k)
        END DO
      END DO
    END DO
  END DO

  ! V12.x Unified: GaussNorm-Tracking (domain-spezifisch + global Legacy)
  ! In veff_gausz (3D Volumen) → max_sigma_ratio_volume_global
  IF (do_gaussnorm) THEN
    max_sigma_ratio_volume_global = MAX(max_sigma_ratio_volume_global, MAXVAL(sigma1))
  END IF

  !/DEBUGGING
  PRINT *, '........................................................ Done'
  PRINT *, '------------------------------------------------------------'
  PRINT *, 'M-Loop:'

  DO m = 0, m_max

    WRITE(*,'(A28,I3,A12,I3,A1,A1,$)') 'Volumen Gausz-Quadratur , m = ',m,', Progress: ',(m*100)/m_max,'% ',achar(13)

    ! v16 (Audit-Find #5): m=0 Special-Case — siehe leff_gausz Z. 198 fuer Begruendung
    IF (m == 0) THEN
      element_nodes(:,:,7) = element_nodes(:,:,4)
    ELSE
      element_nodes(:,:,7) = ( element_nodes(:,:,4)**m + element_nodes(:,:,5)**m + element_nodes(:,:,6)**m )**(1./REAL(m))
    END IF
    sigmaV = 0

    DO i = 1, gauss_order
      DO j = 1, gauss_order
        DO k = 1, gauss_order
          DO n = 1, eckn
            sigmaV(:,i,j,k) = sigmaV(:,i,j,k) + element_nodes(:,n,7)*form_function_h(n,i,j,k)
          END DO
        END DO
      END DO
    END DO

    intVol1 = sigma1**m * jacobi
    intVolV = sigmaV**m * jacobi


    veff_gausz(1,m) = SUM(intVol1)
    veff_gausz(2,m) = SUM(intVolV)


END DO
!/DEBUGGING
PRINT *, '...................................................... Done'
PRINT *, ' deallocate Arrays:'
DEALLOCATE( form_function_h, dh_dxi, dh_deta, dh_dzeta, element_nodes, jacobi, sigma1, sigmaV, intVol1, intVolV, dx_dxi, dx_deta, dx_dzeta, dy_dxi, dy_deta, dy_dzeta, dz_dxi, dz_deta, dz_dzeta )
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

END MODULE overhead_unified
