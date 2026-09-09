"""Heading colour and heading numbering, read out of the document's styles.

The numbers are the part worth testing at this level: a number in front of a
heading is not text, it is a multilevel list in ``numbering.xml`` whose levels
name the styles they belong to. Nothing about it is visible in ``paragraph.text``
-- Word renders it -- so the definition and the link from the style to it are
what there is to check.
"""

from typing import Any, List

import pytest
from docx import Document
from docx.oxml.ns import qn

from mento.reports.headings import NUMBERED_LEVELS, color_headings, number_headings
from mento.results import HEADING_COLOR, TEXT_COLOR, DocumentBuilder


@pytest.fixture
def document() -> Any:
    return Document()


@pytest.fixture
def builder() -> DocumentBuilder:
    return DocumentBuilder(title="Test Document")


def _numbering(document: Any) -> Any:
    return document.part.numbering_part.element


def _levels(document: Any, abstract_id: str) -> List[Any]:
    for abstract in _numbering(document).findall(qn("w:abstractNum")):
        if abstract.get(qn("w:abstractNumId")) == abstract_id:
            return abstract.findall(qn("w:lvl"))
    raise AssertionError(f"no abstractNum {abstract_id}")


def _abstract_id_of(document: Any, num_id: str) -> str:
    for num in _numbering(document).findall(qn("w:num")):
        if num.get(qn("w:numId")) == num_id:
            return str(num.find(qn("w:abstractNumId")).get(qn("w:val")))
    raise AssertionError(f"no num {num_id}")


def _style_num_id(document: Any, level: int) -> Any:
    properties = document.styles[f"Heading {level}"].element.find(qn("w:pPr"))
    number = properties.find(qn("w:numPr")) if properties is not None else None
    return None if number is None else number.find(qn("w:numId")).get(qn("w:val"))


# --- Colour ---


def test_the_report_colours_reach_the_styles(builder: DocumentBuilder) -> None:
    """A heading blue, and one grey for everything else."""
    styles = builder.doc.styles

    assert str(styles["Heading 1"].font.color.rgb) == HEADING_COLOR
    assert str(styles["Heading 2"].font.color.rgb) == TEXT_COLOR
    assert str(styles["Normal"].font.color.rgb) == TEXT_COLOR


def test_a_heading_colour_carries_no_theme_attribute(document: Any) -> None:
    """Word's heading styles name a theme colour beside the literal one.

    It resolves the theme first, so a colour set next to one is a colour that
    never appears -- the same trap as the theme fill in a table style.
    """
    color_headings(document, "0A3E81", "323232", "323232")
    run_properties = document.styles["Heading 1"].element.find(qn("w:rPr"))
    color = run_properties.find(qn("w:color"))

    assert color.get(qn("w:val")) == "0A3E81"
    assert color.get(qn("w:themeColor")) is None
    assert color.get(qn("w:themeShade")) is None


def test_deeper_headings_take_the_subheading_colour(document: Any) -> None:
    color_headings(document, "0A3E81", "323232", "323232")

    for level in (2, 3, 4):
        assert str(document.styles[f"Heading {level}"].font.color.rgb) == "323232"


# --- Numbering ---


def test_headings_are_numbered_one_and_one_one(builder: DocumentBuilder) -> None:
    """``%1`` on the title, ``%1.%2`` on the sections under it."""
    num_id = _style_num_id(builder.doc, 1)
    levels = _levels(builder.doc, _abstract_id_of(builder.doc, num_id))

    assert [level.find(qn("w:lvlText")).get(qn("w:val")) for level in levels] == ["%1", "%1.%2"]
    assert [level.find(qn("w:numFmt")).get(qn("w:val")) for level in levels] == ["decimal", "decimal"]


def test_each_level_names_the_style_it_belongs_to(builder: DocumentBuilder) -> None:
    """``w:pStyle`` inside a level is what ties the number to the style."""
    num_id = _style_num_id(builder.doc, 1)
    levels = _levels(builder.doc, _abstract_id_of(builder.doc, num_id))

    assert [level.find(qn("w:pStyle")).get(qn("w:val")) for level in levels] == ["Heading1", "Heading2"]


def test_both_heading_styles_point_at_the_same_list(builder: DocumentBuilder) -> None:
    """One list, two levels -- otherwise 1.1 counts against nothing."""
    assert _style_num_id(builder.doc, 1) == _style_num_id(builder.doc, 2)

    for level, ilvl in NUMBERED_LEVELS.items():
        properties = builder.doc.styles[f"Heading {level}"].element.find(qn("w:pPr"))
        assert properties.find(qn("w:numPr")).find(qn("w:ilvl")).get(qn("w:val")) == str(ilvl)


def test_deeper_headings_are_left_unnumbered(builder: DocumentBuilder) -> None:
    """The reports do not go past level 2, and an unwritten level is one that
    cannot be checked."""
    assert _style_num_id(builder.doc, 3) is None


def test_the_list_does_not_take_an_id_the_template_uses(document: Any) -> None:
    """python-docx's template already carries nine lists.

    Redefining one of their ids would renumber Word's own bulleted and
    numbered styles, and nothing would say so.
    """
    before = {element.get(qn("w:numId")) for element in _numbering(document).findall(qn("w:num"))}
    before_abstract = {
        element.get(qn("w:abstractNumId")) for element in _numbering(document).findall(qn("w:abstractNum"))
    }

    num_id = str(number_headings(document))

    assert num_id not in before
    assert _abstract_id_of(document, num_id) not in before_abstract


def test_the_number_is_separated_by_a_space(builder: DocumentBuilder) -> None:
    """Word's default is a tab, which jumps to a tab stop and reads as an
    indent rather than as a number in front of the title."""
    num_id = _style_num_id(builder.doc, 1)

    for level in _levels(builder.doc, _abstract_id_of(builder.doc, num_id)):
        assert level.find(qn("w:suff")).get(qn("w:val")) == "space"
        assert level.find(qn("w:pPr")).find(qn("w:ind")).get(qn("w:left")) == "0"


def test_a_document_defines_the_list_once(builder: DocumentBuilder) -> None:
    """It is set up when the document is built, not per heading."""
    for section in ("Materials", "Limit checks", "Design checks"):
        builder.add_heading(section, level=2)

    abstract = _numbering(builder.doc).findall(qn("w:abstractNum"))
    mento_lists = [
        element
        for element in abstract
        if any(level.find(qn("w:pStyle")) is not None for level in element.findall(qn("w:lvl")))
        and element.find(qn("w:lvl")).find(qn("w:pStyle")).get(qn("w:val")) == "Heading1"
    ]

    assert len(mento_lists) == 1
