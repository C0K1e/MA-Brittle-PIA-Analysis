# log_helper.py
"""
v15.3: Zentrale Print-Helfer fuer einheitliche Terminal-Ausgaben.

Pragmatischer Ansatz — KEIN flaechendeckendes Refactor aller bestehenden print()-
Stellen. Stattdessen:
  - Bestehende konsistente Tags ([Analyst], [VTK], [Report]) bleiben unveraendert
  - Diese Helfer werden in (a) neuen Code-Stellen und (b) unbetagged Stellen
    (z.B. 01-run_manager.py Z. 43, 172, 175, ...) eingesetzt
  - Plot-Skripte mit stillem `except: continue` koennen log_warn() nutzen

Verwendung:
    from log_helper import log_info, log_warn, log_error, log_skip, log_ok

    log_info("Manager", "Batch gestartet")
    log_warn("PLOT", f"CSV unlesbar: {path}: {e}")
    log_error("Analyst", "REF=GPX ohne GauszInfo")
    log_skip("VTK", "Gauss-Pfad nicht vorgesehen")
    log_ok("ErrDecomp", "Produktform OK")

Output-Format (3-Zeichen-Indent + Tag):
    [Manager] Batch gestartet
    [PLOT] WARN: CSV unlesbar: ...
    !!! [Analyst] ERROR: REF=GPX ohne GauszInfo
    [VTK] SKIP: Gauss-Pfad nicht vorgesehen
    [ErrDecomp] OK: Produktform OK
"""


def log_info(tag: str, msg: str) -> None:
    """Standard-Info: Status-Meldung im Pipeline-Fluss."""
    print(f"   [{tag}] {msg}")


def log_warn(tag: str, msg: str) -> None:
    """Warnung: Auffaellig, aber nicht fatal."""
    print(f"   [{tag}] WARN: {msg}")


def log_error(tag: str, msg: str) -> None:
    """Fehler: Methodisch oder pipeline-relevant problematisch.

    Wirft NICHT — Aufrufer entscheidet ob raise oder weiter. Fuer hart fail:
    `log_error(...); raise RuntimeError(...)` Pattern.
    """
    print(f"!!! [{tag}] ERROR: {msg}")


def log_skip(tag: str, msg: str) -> None:
    """Skip: Erwartete Skip-Logik (z.B. VTK bei Gauss-Pfad)."""
    print(f"   [{tag}] SKIP: {msg}")


def log_ok(tag: str, msg: str) -> None:
    """Erfolgreich abgeschlossen: positive Status-Meldung."""
    print(f"   [{tag}] OK: {msg}")
