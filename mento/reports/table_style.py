"""How a Word table looks -- one style definition, written into the document once.

:mod:`~mento.reports.tables` decides what a table says; this module decides how
it looks. The separation is why this is its own file: a colour choice should
never be a reason to touch the code that builds a results frame.

Word has two ways of making a table look like something, and only one of them
survives being edited. Formatting each cell -- a fill here, a border there --
produces the right picture and nothing else: a reader who adds a row in Word
gets a bare one, because the banding was never a rule, only forty-odd painted
cells. A *style* is the rule. It lives once in the document's ``styles.xml``,
every table points at it by id, and Word applies it to whatever the table
happens to contain, including the row added after we were done.

So mento defines its own style rather than borrowing one of Word's built-ins.
``Light Shading``, which the reports used before, is unreachable in all the
ways that matter: its grey and its rules live in python-docx's template where
no argument can reach them; it names a fill and a *theme* fill at once, and
Word resolves the theme first, so the colour it declares is not the colour on
the page; and it bolds the first column, which is what the reports then had to
undo cell by cell.

The knobs are :class:`TableStyle` and the setter is package level, the same way
the language is set::

    import mento
    from mento.reports.table_style import TableStyle

    mento.set_table_style(TableStyle(band_fill="F4F7FB", header_fill="2C6DC4"))

Two units here are not the ones a user thinks in:

- A **colour** is six hexadecimal digits with no ``#``. Word does not complain
  about ``"#F2F2F2"``; it discards the fill and moves on, so the constructor
  rejects it instead.
- A **thickness** is given in points and written in eighths of a point, which
  is the only thing ``w:sz`` understands: 1.0 pt is ``w:sz="8"``. Word draws
  nothing thinner than 0.25 pt and nothing thicker than 12 pt, so the value is
  clamped to that range. A thickness of exactly ``0`` is not a hairline but the
  absence of a line, written ``w:val="nil"``.

The cell padding is in points too, written in twentieths. It is the field that
decides how tall a long table is, because it is charged once per row: the
detailed annexes are forty rows of it, and 0.4 pt either way is the difference
between one page and two.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional
from xml.sax.saxutils import escape

from docx.oxml.ns import nsdecls

#: A colour in OOXML is six hexadecimal digits, and the ``#`` a stylesheet
#: would carry is not part of it. An empty string is a colour left unset.
_HEX_COLOR = re.compile(r"\A[0-9A-Fa-f]{6}\Z")

#: The thinnest and thickest rule Word will draw, in eighths of a point:
#: 0.25 pt and 12 pt. Asking for less than the first still draws the first, so
#: a thickness meaning "no line" has to be spelled as exactly zero.
MIN_EIGHTHS = 2
MAX_EIGHTHS = 96

#: Where in Word's style gallery the definition sorts. High enough to stay out
#: of the way of the styles a user actually picks from.
_UI_PRIORITY = 59

#: Cell padding left and right, in twentieths of a point: Word's usual indent.
#: The padding above and below is a field, because it is the one dimension that
#: multiplies -- forty rows of it decide whether an annex is one page or two.
_CELL_MARGIN_HORIZONTAL = 108

#: Twentieths of a point in a point, which is what ``w:tblCellMar`` counts in.
_TWIPS_PER_POINT = 20


def points_to_eighths(points: float) -> int:
    """``points`` as the eighths of a point ``w:sz`` is written in.

    Zero passes through as zero, which the border writer reads as "no line".
    Anything else lands in ``[MIN_EIGHTHS, MAX_EIGHTHS]``, because a value
    outside that range is not drawn as asked: Word clamps it silently, and a
    0.05 pt rule requested and a 0.25 pt rule delivered is a difference worth
    resolving here rather than in the page.
    """
    if points < 0:
        raise ValueError(f"A border thickness cannot be negative; got {points}.")
    if points == 0:
        return 0
    return max(MIN_EIGHTHS, min(MAX_EIGHTHS, round(points * 8)))


def _validate_color(value: str, field: str) -> None:
    """Reject a colour Word would discard without saying so."""
    if value and not _HEX_COLOR.match(value):
        raise ValueError(
            f"{field}={value!r} is not a colour. Give six hexadecimal digits and no '#', "
            f"as in 'F2F2F2'; Word drops anything else without reporting it."
        )


def _border(tag: str, points: float, color: str) -> str:
    """One ``w:tblBorders`` or ``w:tcBorders`` child.

    A thickness of zero is written ``nil`` -- not a thin line, but the explicit
    absence of one, which is what stops a border inherited from ``TableNormal``
    from showing through.
    """
    eighths = points_to_eighths(points)
    if eighths == 0:
        return f'<w:{tag} w:val="nil"/>'
    return f'<w:{tag} w:val="single" w:sz="{eighths}" w:space="0" w:color="{color or "auto"}"/>'


def _shading(fill: str) -> str:
    """A cell fill, named and only named.

    Deliberately without ``w:themeFill``: Word resolves a theme fill ahead of
    the literal one beside it, so a definition carrying both shows the theme's
    colour and not the one that was asked for. That is the bug in Word's own
    ``Light Shading``, and repeating it would make every colour in this module
    a suggestion.
    """
    return f'<w:shd w:val="clear" w:color="auto" w:fill="{fill}"/>'


@dataclass(frozen=True)
class TableStyle:
    """The look of every table in a Word report.

    Frozen, because it is written into a document by identity: the same style
    used by forty tables is one definition in ``styles.xml``, and a style that
    could change under the document would make that promise unkeepable.

    Parameters
    ----------
    name : str
        The style's name in Word's gallery. Its ``styleId`` is this name with
        everything that is not a letter or a digit removed, so two styles meant
        to coexist in one document need two names.
    band_fill : str
        Fill of the banded rows, six hex digits. Empty for no banding.
    band_size : int
        How many rows a band spans.
    header_fill : str
        Fill of the header row. Empty for no fill.
    header_color : str
        Text colour of the header row. Empty to use the body's.
    header_bold : bool
        Whether the header row is bold.
    text_color : str
        Text colour of the table. Empty to inherit the document's.
    border_color : str
        Colour of every rule the style draws.
    top_pt : float
        Rule above the header.
    header_rule_pt : float
        Rule between the header and the first data row.
    bottom_pt : float
        Rule closing the table.
    inside_pt : float
        Rule between data rows. Zero by default: the banding already separates
        them, and a line as well makes the table louder than what it says.
    vertical_pt : float
        Rule between columns.
    cell_padding_pt : float
        Air above and below the text in a cell. It is charged once per row, so
        it is the field that decides how tall a long table is. 1.4 pt is
        Word's own comfortable padding; the default is tighter so that a
        detailed annex closes on one page.
    """

    name: str = "Mento Table"
    band_fill: str = "F2F2F2"
    band_size: int = 1
    header_fill: str = ""
    header_color: str = ""
    header_bold: bool = True
    text_color: str = "323232"
    border_color: str = "404040"
    top_pt: float = 1.0
    header_rule_pt: float = 1.0
    bottom_pt: float = 1.0
    inside_pt: float = 0.0
    vertical_pt: float = 0.0
    cell_padding_pt: float = 0.8

    def __post_init__(self) -> None:
        if not self.style_id:
            raise ValueError(f"name={self.name!r} has no letter or digit to make a style id from.")
        for field in ("band_fill", "header_fill", "header_color", "text_color", "border_color"):
            _validate_color(getattr(self, field), field)
        if self.cell_padding_pt < 0:
            raise ValueError(f"cell_padding_pt={self.cell_padding_pt} is not a padding; it cannot be negative.")
        if self.band_size < 1:
            raise ValueError(f"band_size={self.band_size} is not a number of rows; it must be 1 or more.")
        # Validating by conversion: every thickness has to survive the same
        # arithmetic that writes it, so a negative one is refused here rather
        # than reaching the document as a missing border.
        for field in ("top_pt", "header_rule_pt", "bottom_pt", "inside_pt", "vertical_pt"):
            points_to_eighths(getattr(self, field))

    @property
    def style_id(self) -> str:
        """The name reduced to what OOXML accepts as a ``w:styleId``."""
        return "".join(character for character in self.name if character.isalnum())

    def definition_xml(self, style_id: Optional[str] = None) -> str:
        """The ``w:style`` element to append to a document's ``styles.xml``.

        ``style_id`` overrides :attr:`style_id`, which is what lets a document
        already holding a style of that id take this one under another.
        """
        style_id = style_id or self.style_id
        padding = int(round(self.cell_padding_pt * _TWIPS_PER_POINT))
        # A document that already held this id gets the definition under
        # another; the gallery name follows it, so two styles in one document
        # are two entries a user can tell apart rather than two "Mento Table"s.
        name = self.name + style_id[len(self.style_id) :]
        name = escape(name, {'"': "&quot;"})

        header_run = "<w:b/><w:bCs/>" if self.header_bold else ""
        if self.header_color:
            header_run += f'<w:color w:val="{self.header_color}"/>'

        # The header's rules are set on its cells rather than on the table: the
        # line above it is the table's own top border, but the one under it
        # separates two rows, and no table-wide setting names that edge.
        header_cell = (
            "<w:tcBorders>"
            f"{_border('top', self.top_pt, self.border_color)}"
            f"{_border('bottom', self.header_rule_pt, self.border_color)}"
            "</w:tcBorders>"
        ) + (_shading(self.header_fill) if self.header_fill else "")

        conditional = (
            '<w:tblStylePr w:type="firstRow">'
            + (f"<w:rPr>{header_run}</w:rPr>" if header_run else "")
            + f"<w:tcPr>{header_cell}</w:tcPr>"
            + "</w:tblStylePr>"
        )
        if self.band_fill:
            conditional += (
                f'<w:tblStylePr w:type="band1Horz"><w:tcPr>{_shading(self.band_fill)}</w:tcPr></w:tblStylePr>'
            )

        return (
            f'<w:style {nsdecls("w")} w:type="table" w:styleId="{style_id}">'
            f'<w:name w:val="{name}"/>'
            '<w:basedOn w:val="TableNormal"/>'
            f'<w:uiPriority w:val="{_UI_PRIORITY}"/>'
            '<w:pPr><w:spacing w:after="0" w:line="240" w:lineRule="auto"/></w:pPr>'
            + (f'<w:rPr><w:color w:val="{self.text_color}"/></w:rPr>' if self.text_color else "")
            + "<w:tblPr>"
            f'<w:tblStyleRowBandSize w:val="{self.band_size}"/>'
            '<w:tblInd w:w="0" w:type="dxa"/>'
            "<w:tblBorders>"
            f"{_border('top', self.top_pt, self.border_color)}"
            f"{_border('bottom', self.bottom_pt, self.border_color)}"
            f"{_border('insideH', self.inside_pt, self.border_color)}"
            f"{_border('insideV', self.vertical_pt, self.border_color)}"
            "</w:tblBorders>"
            "<w:tblCellMar>"
            f'<w:top w:w="{padding}" w:type="dxa"/>'
            f'<w:left w:w="{_CELL_MARGIN_HORIZONTAL}" w:type="dxa"/>'
            f'<w:bottom w:w="{padding}" w:type="dxa"/>'
            f'<w:right w:w="{_CELL_MARGIN_HORIZONTAL}" w:type="dxa"/>'
            "</w:tblCellMar>"
            "</w:tblPr>" + conditional + "</w:style>"
        )


#: The style every :class:`~mento.results.DocumentBuilder` starts from. Module
#: level and mutable through the setter below, the same shape as the language:
#: a report is configured once, at the top of a script, not per table.
_table_style: TableStyle = TableStyle()


def set_table_style(style: TableStyle) -> None:
    """Set the look of every Word table produced from now on.

    Parameters
    ----------
    style : TableStyle
        The style to write into each new document.

    Examples
    --------
    >>> import mento
    >>> from mento.reports.table_style import TableStyle
    >>> mento.set_table_style(TableStyle(band_fill="F4F7FB", header_fill="2C6DC4"))
    """
    if not isinstance(style, TableStyle):
        raise TypeError(f"set_table_style takes a TableStyle, not {type(style).__name__}.")
    global _table_style
    _table_style = style


def get_table_style() -> TableStyle:
    """The style new Word documents are currently given."""
    return _table_style
