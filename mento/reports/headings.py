"""How a report's headings look: their colour, and their numbering.

Like :mod:`~mento.reports.table_style`, this is appearance rather than content,
and for the same reason it is written into the document's styles rather than
onto the paragraphs: a heading added in Word afterwards is coloured and
numbered because ``Heading 1`` says so, not because the builder painted the
ones it happened to write.

The numbering is the part Word will not do from the object model. A number in
front of a heading is not text: it is a *multilevel list* defined once in
``numbering.xml``, whose levels name the styles they belong to, with each
heading style pointing back at the level it takes its number from. Set up that
way the numbers are Word's own -- they renumber when a section is moved,
deleted or inserted, which is the whole point of not writing "1.1" into the
text.

Word's built-in heading styles are named in terms of the document *theme*
twice over: ``w:themeColor`` beside ``w:color``, and ``w:asciiTheme`` beside
``w:ascii``. Word resolves the theme side first, so a colour or a font named
next to one is a colour or a font that never appears -- the same trap as the
theme fill in :mod:`~mento.reports.table_style`. python-docx's
``font.color.rgb`` replaces the whole ``w:color`` element and so takes the
theme attributes with it; ``font.name`` only adds an attribute to the
``w:rFonts`` already there, so the theme names have to be removed by hand.

That one bites in a place nobody looks: the heading *text* is set run by run
and comes out right, while the paragraph mark keeps the style's theme font --
and a heading's number is drawn in the mark's font. A Lato report with Calibri
numbers in front of its headings is what that looks like.
"""

from __future__ import annotations

from typing import Any

# python-docx is imported inside the functions that use it: this module is on
# the import path of every element, and a calculation never opens a document.

#: The heading levels that are numbered, and the level of the list each takes
#: its number from: "1" for a Heading 1 and "1.1" for a Heading 2. Deeper
#: headings are left unnumbered -- the reports do not go deeper, and a level
#: nobody writes is a definition nobody can check.
NUMBERED_LEVELS = {1: 0, 2: 1}

#: How the number is separated from the heading text. A tab -- Word's default
#: -- jumps to the next tab stop and leaves a gap wide enough to read as an
#: indent; a space puts the number in front of the title, which is what a
#: numbered heading looks like.
_NUMBER_SUFFIX = "space"


#: The heading levels styled here. Word's template defines nine; the reports
#: use two, and the rest are given the same treatment so a document that grows
#: a level does not grow a Calibri heading with it.
_STYLED_LEVELS = range(1, 5)


def _style_id(level: int) -> str:
    return f"Heading{level}"


def _name_font(style: Any, font_name: str) -> None:
    """Name the font of a style, and unname the theme's.

    ``font.name`` writes ``w:ascii`` into the ``w:rFonts`` that is already
    there, next to the ``w:asciiTheme`` the built-in heading styles carry --
    and Word reads the theme one. The attribute has to go, or the name beside
    it is decoration.
    """
    style.font.name = font_name
    from docx.oxml.ns import qn

    fonts = style.element.find(qn("w:rPr")).find(qn("w:rFonts"))
    for attribute in ("asciiTheme", "hAnsiTheme", "eastAsiaTheme", "cstheme"):
        if fonts.get(qn(f"w:{attribute}")) is not None:
            del fonts.attrib[qn(f"w:{attribute}")]
    fonts.set(qn("w:cs"), font_name)


def style_headings(
    document: Any,
    font_name: str,
    heading_color: str,
    subheading_color: str,
    text_color: str,
) -> None:
    """Give the document's text and heading styles their font and colour.

    ``heading_color`` is for ``Heading 1``, ``subheading_color`` for every
    deeper heading, and ``text_color`` for the body -- set on ``Normal``, which
    the rest inherit from.
    """
    normal = document.styles["Normal"]
    from docx.shared import RGBColor

    normal.font.color.rgb = RGBColor.from_string(text_color)
    _name_font(normal, font_name)

    for level in _STYLED_LEVELS:
        style = document.styles[f"Heading {level}"]
        style.font.color.rgb = RGBColor.from_string(heading_color if level == 1 else subheading_color)
        _name_font(style, font_name)


def _abstract_numbering_xml(abstract_id: int, font_name: str) -> str:
    """A two-level list whose levels belong to the heading styles.

    ``w:pStyle`` inside a level is what ties it to the style, and
    ``w:lvlText`` is the number itself: ``%1`` is this level's counter and
    ``%1.%2`` the first level's followed by the second's.

    The level names its font as well. A number is drawn in the paragraph
    mark's font unless the level says otherwise, and saying so here is what
    keeps a heading's number from being the one word of the report in a font
    nobody chose.
    """
    levels = "".join(
        f'<w:lvl w:ilvl="{ilvl}">'
        '<w:start w:val="1"/>'
        '<w:numFmt w:val="decimal"/>'
        f'<w:pStyle w:val="{_style_id(level)}"/>'
        f'<w:suff w:val="{_NUMBER_SUFFIX}"/>'
        f'<w:lvlText w:val="{".".join(f"%{n + 1}" for n in range(ilvl + 1))}"/>'
        '<w:lvlJc w:val="left"/>'
        '<w:pPr><w:ind w:left="0" w:firstLine="0"/></w:pPr>'
        f'<w:rPr><w:rFonts w:ascii="{font_name}" w:hAnsi="{font_name}" w:cs="{font_name}"/></w:rPr>'
        "</w:lvl>"
        for level, ilvl in sorted(NUMBERED_LEVELS.items())
    )
    from docx.oxml.ns import nsdecls

    return (
        f'<w:abstractNum {nsdecls("w")} w:abstractNumId="{abstract_id}">'
        '<w:multiLevelType w:val="multilevel"/>'
        f"{levels}"
        "</w:abstractNum>"
    )


def _next_free_id(numbering: Any, tag: str, attribute: str) -> int:
    """One past the highest id already used, so nothing is redefined.

    The template python-docx builds a document from already carries nine lists
    -- the bulleted and numbered ones of Word's own gallery -- and taking an id
    that is in use would silently renumber one of them.
    """
    from docx.oxml.ns import qn

    used = [int(element.get(qn(attribute))) for element in numbering.findall(qn(tag))]
    return max(used, default=-1) + 1


def number_headings(document: Any, font_name: str) -> int:
    """Number the heading styles ``1`` and ``1.1``, and return the list's id.

    The definition goes into ``numbering.xml`` once; each heading style is then
    pointed at the level it takes its number from. Calling this twice on one
    document would define the list twice, so it is called once, when the
    document is built.
    """
    from docx.oxml import parse_xml
    from docx.oxml.ns import nsdecls

    numbering = document.part.numbering_part.element

    abstract_id = _next_free_id(numbering, "w:abstractNum", "w:abstractNumId")
    numbering.append(parse_xml(_abstract_numbering_xml(abstract_id, font_name)))

    # A style points at a `w:num`, which points at the definition. The
    # indirection is Word's: it is what lets two lists share a shape.
    num_id = max(_next_free_id(numbering, "w:num", "w:numId"), 1)
    numbering.append(
        parse_xml(f'<w:num {nsdecls("w")} w:numId="{num_id}"><w:abstractNumId w:val="{abstract_id}"/></w:num>')
    )

    for level, ilvl in NUMBERED_LEVELS.items():
        properties = document.styles[f"Heading {level}"].element.get_or_add_pPr()
        number_properties = properties.get_or_add_numPr()
        number_properties.get_or_add_ilvl().val = ilvl
        number_properties.get_or_add_numId().val = num_id

    return num_id
