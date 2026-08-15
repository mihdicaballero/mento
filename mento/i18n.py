"""Language of the detailed report output.

``flexure_results_detailed``, ``shear_results_detailed`` and their ``_doc``
counterparts print English by default. Switch the whole package once and every
later report comes out translated::

    import mento

    mento.set_language("es")
    node.shear_results_detailed()        # console tables in Spanish
    node.shear_results_detailed_doc()    # Word document in Spanish

The English text a report builder produces *is* the catalog key, so the design
code modules in ``mento/codes`` stay untouched and monolingual. A catalog only
has to carry the strings that differ; anything missing falls back to English
instead of raising, so a new label always renders, translated or not.

Adding a language is data, not code: write a ``{english: translation}`` mapping
and register it in ``_CATALOGS``.

Scope: report text only. Variable names (``fc``, ``Av``, ``DCR``), units,
numbers, the design code designation and generated file names stay as they are.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Optional, Tuple

if TYPE_CHECKING:
    import pandas as pd

DEFAULT_LANGUAGE = "en"

# --------------------------------------------------------------------------
# Spanish catalog
# --------------------------------------------------------------------------
# Terminology follows CIRSOC 201 usage: "hormigón", "estribo", "cuantía",
# "solicitaciones", "recubrimiento", "altura útil".

ES: Dict[str, str] = {
    # -- console section titles --------------------------------------------
    "===== BEAM FLEXURE DETAILED RESULTS =====": "===== RESULTADOS DETALLADOS DE FLEXIÓN DE VIGA =====",
    "===== BEAM SHEAR DETAILED RESULTS =====": "===== RESULTADOS DETALLADOS DE CORTE DE VIGA =====",
    "===== SLAB FLEXURE DETAILED RESULTS =====": "===== RESULTADOS DETALLADOS DE FLEXIÓN DE LOSA =====",
    "===== SLAB SHEAR DETAILED RESULTS =====": "===== RESULTADOS DETALLADOS DE CORTE DE LOSA =====",
    "===== SHEAR WALL DETAILED RESULTS =====": "===== RESULTADOS DETALLADOS DE TABIQUE =====",
    "MATERIALS": "MATERIALES",
    "GEOMETRY": "GEOMETRÍA",
    "FORCES": "SOLICITACIONES",
    "MAX AND MIN LIMIT CHECKS": "VERIFICACIONES DE LÍMITES MÁXIMOS Y MÍNIMOS",
    "SHEAR STRENGTH": "RESISTENCIA AL CORTE",
    "CONCRETE STRENGTH": "RESISTENCIA DEL HORMIGÓN",
    "FLEXURAL CAPACITY - TOP": "CAPACIDAD A FLEXIÓN - SUPERIOR",
    "FLEXURAL CAPACITY - BOTTOM": "CAPACIDAD A FLEXIÓN - INFERIOR",
    # -- table column headers ----------------------------------------------
    "Materials": "Materiales",
    "Geometry": "Geometría",
    "Design forces": "Solicitaciones",
    "Variable": "Variable",
    "Value": "Valor",
    "Unit": "Unidad",
    "Check": "Verificación",
    "Min.": "Mín.",
    "Max.": "Máx.",
    "Ok?": "¿Ok?",
    "Shear reinforcement strength": "Resistencia de la armadura de corte",
    "Shear strength": "Resistencia al corte",
    "Top reinforcement check": "Verificación de la armadura superior",
    "Bottom reinforcement check": "Verificación de la armadura inferior",
    # -- Word document headings --------------------------------------------
    "Concrete beam flexure check": "Verificación a flexión de viga de hormigón",
    "Concrete beam shear check": "Verificación a corte de viga de hormigón",
    "Concrete slab flexure check": "Verificación a flexión de losa de hormigón",
    "Concrete slab shear check": "Verificación a corte de losa de hormigón",
    "Concrete shear wall check": "Verificación de tabique de hormigón",
    "Beam {label} flexure check": "Verificación a flexión de viga {label}",
    "Beam {label} shear check": "Verificación a corte de viga {label}",
    "Slab {label} flexure check": "Verificación a flexión de losa {label}",
    "Slab {label} shear check": "Verificación a corte de losa {label}",
    "Shear Wall {label} shear check": "Verificación a corte de tabique {label}",
    "Made with mento {version}. Design code: {design_code}": (
        "Generado con mento {version}. Código de diseño: {design_code}"
    ),
    "Limit checks": "Verificaciones de límites",
    "Design checks": "Verificaciones de diseño",
    "Flexural Capacity Top": "Capacidad a flexión superior",
    "Flexural Capacity Bottom": "Capacidad a flexión inferior",
    # -- row labels: materials and geometry --------------------------------
    "Section Label": "Identificación de la sección",
    "Concrete strength": "Resistencia del hormigón",
    "Steel reinforcement yield strength": "Tensión de fluencia del acero",
    "Concrete density": "Densidad del hormigón",
    "Normalweight concrete": "Hormigón de densidad normal",
    "Safety factor for shear": "Coeficiente de seguridad para corte",
    "Safety factor for concrete": "Coeficiente de seguridad del hormigón",
    "Safety factor for steel": "Coeficiente de seguridad del acero",
    "Section height": "Altura de la sección",
    "Section width": "Ancho de la sección",
    "Clear cover": "Recubrimiento libre",
    "Mechanical top cover": "Recubrimiento mecánico superior",
    "Mechanical bottom cover": "Recubrimiento mecánico inferior",
    "Longitudinal tension rebar": "Armadura longitudinal traccionada",
    "Effective height": "Altura útil",
    # -- row labels: forces -------------------------------------------------
    "Axial, positive for compression": "Axial, positivo en compresión",
    "Shear": "Corte",
    "Top max moment": "Momento máximo superior",
    "Bottom max moment": "Momento máximo inferior",
    # -- row labels: limit checks ------------------------------------------
    "Stirrup spacing along length": "Separación de estribos en la dirección longitudinal",
    "Stirrup spacing along width": "Separación de estribos en la dirección transversal",
    "Minimum shear reinforcement": "Armadura mínima de corte",
    "Minimum rebar diameter": "Diámetro mínimo de barra",
    "Min/Max As rebar top": "As mín/máx de la armadura superior",
    "Min/Max As rebar bottom": "As mín/máx de la armadura inferior",
    "Minimum spacing top": "Separación mínima superior",
    "Minimum spacing bottom": "Separación mínima inferior",
    # -- row labels: shear reinforcement -----------------------------------
    "Number of stirrups": "Número de estribos",
    "Stirrup diameter": "Diámetro del estribo",
    "Stirrup spacing": "Separación de estribos",
    "Minimum shear reinforcing": "Armadura mínima de corte",
    "Required shear reinforcing": "Armadura de corte requerida",
    "Defined shear reinforcing": "Armadura de corte adoptada",
    "Shear rebar strength": "Resistencia de la armadura de corte",
    "Steel shear strength": "Resistencia al corte del acero",
    "Concrete shear strength": "Resistencia al corte del hormigón",
    # -- row labels: shear capacity ----------------------------------------
    "Effective shear area": "Área efectiva de corte",
    "Gross shear area": "Área bruta de corte",
    "Longitudinal reinforcement ratio": "Cuantía de armadura longitudinal",
    "Size modification factor": "Factor de modificación por tamaño",
    "Axial stress": "Tensión axial",
    "Concrete effective shear stress": "Tensión de corte efectiva del hormigón",
    "Maximum shear strength": "Resistencia máxima al corte",
    "Maximum shear capacity": "Capacidad máxima de corte",
    "Total shear strength": "Resistencia total al corte",
    "Total shear capacity": "Capacidad total de corte",
    "Max shear check": "Verificación de corte máximo",
    "Demand Capacity Ratio": "Relación demanda-capacidad",
    "Concrete strut angle": "Ángulo de la biela de hormigón",
    "k value": "Valor de k",
    "Coefficient for long term effects and loading effects": (
        "Coeficiente de efectos de larga duración y de aplicación de la carga"
    ),
    # -- row labels: flexural capacity -------------------------------------
    "First layer bars": "Barras de la primera capa",
    "Second layer bars": "Barras de la segunda capa",
    "Minimum rebar reinforcing": "Armadura mínima",
    "Required rebar reinforcing top": "Armadura superior requerida",
    "Required rebar reinforcing bottom": "Armadura inferior requerida",
    "Defined rebar reinforcing top": "Armadura superior adoptada",
    "Defined rebar reinforcing bottom": "Armadura inferior adoptada",
    "Depth of equivalent strength block ratio": "Relación de profundidad del bloque equivalente",
    "Total flexural strength": "Resistencia total a flexión",
    # -- row labels: shear wall --------------------------------------------
    "Wall thickness": "Espesor del tabique",
    "Wall length": "Longitud del tabique",
    "Wall height": "Altura del tabique",
    "Aspect ratio": "Relación de aspecto",
    "Horizontal bar spacing (E.F.)": "Separación de barras horizontales (en cada cara)",
    "Vertical bar spacing (E.F.)": "Separación de barras verticales (en cada cara)",
    "Horizontal reinforcement ratio": "Cuantía de armadura horizontal",
    "Minimum vertical reinf. ratio": "Cuantía vertical mínima",
}

# English is the source language, so its catalog is empty: every lookup falls
# through to the key itself.
_CATALOGS: Dict[str, Dict[str, str]] = {
    "en": {},
    "es": ES,
}

_language: str = DEFAULT_LANGUAGE


def available_languages() -> Tuple[str, ...]:
    """Language codes ``set_language`` accepts."""
    return tuple(sorted(_CATALOGS))


def set_language(language: str) -> None:
    """Set the language of every detailed report produced from now on.

    Parameters
    ----------
    language : str
        ISO 639-1 code, ``"en"`` or ``"es"``.

    Raises
    ------
    ValueError
        If the language has no catalog.
    """
    if language not in _CATALOGS:
        raise ValueError(f"Unknown language {language!r}. Available: {', '.join(available_languages())}.")
    global _language
    _language = language


def get_language() -> str:
    """The language detailed reports are currently rendered in."""
    return _language


def translate(text: str, language: Optional[str] = None, **fields: Any) -> str:
    """Translate one report string, filling ``{placeholders}`` from ``fields``.

    Falls back to ``text`` when the catalog has no entry for it, so an
    untranslated label still renders.
    """
    catalog = _CATALOGS.get(get_language() if language is None else language, {})
    translated = catalog.get(text, text)
    return translated.format(**fields) if fields else translated


def translate_table(data: Mapping[str, List[Any]], language: Optional[str] = None) -> Dict[str, List[Any]]:
    """Translate a report table: its column headers and its label column.

    The first column holds the row labels; the rest are values, units and check
    marks, which stay as they are. Returns a new dict — the caller's table, which
    the section keeps as state, is never modified.
    """
    if not data:
        return dict(data)

    label_column = next(iter(data))
    translated: Dict[str, List[Any]] = {}
    for column, values in data.items():
        if column == label_column:
            values = [translate(v, language) if isinstance(v, str) else v for v in values]
        translated[translate(column, language)] = values
    return translated


def translate_dataframe(df: "pd.DataFrame", language: Optional[str] = None) -> "pd.DataFrame":
    """Same as :func:`translate_table`, for the DataFrames the Word builder takes."""
    if not len(df.columns):
        return df

    out = df.copy()
    label_column = out.columns[0]
    out[label_column] = [translate(v, language) if isinstance(v, str) else v for v in out[label_column]]
    out.columns = [translate(str(c), language) for c in out.columns]
    return out
