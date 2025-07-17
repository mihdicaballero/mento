import pytest
import pandas as pd
from docx import Document
from typing import List, Any
from docx.shared import Cm
from mento.results import Formatter, TablePrinter, DocumentBuilder


@pytest.fixture
def formatter() -> Formatter:
    return Formatter()


@pytest.fixture
def table_printer() -> TablePrinter:
    return TablePrinter()


@pytest.fixture
def document_builder() -> DocumentBuilder:
    return DocumentBuilder(title="Test Document")


# Tests for Formatter class
def test_formatter_DCR(formatter: Formatter) -> None:
    assert formatter.DCR(0.85) == "$\\color{#439b00}{\\text{DCR}=0.85}$"
    assert formatter.DCR(0.95) == "$\\color{#efc200}{\\text{DCR}=0.95}$"
    assert formatter.DCR(1.2) == "$\\color{#d43e36}{\\text{DCR}=1.2}$"


def test_formatter_DCR_value(formatter: Formatter) -> None:
    assert formatter.DCR_value(0.85) == "$\\color{#439b00}{0.85}$"
    assert formatter.DCR_value(0.95) == "$\\color{#efc200}{0.95}$"
    assert formatter.DCR_value(1.2) == "$\\color{#d43e36}{1.2}$"


def test_formatter_is_lower(formatter: Formatter) -> None:
    assert formatter.is_lower(0.8, 1.0) == "$\\color{#439b00}{\\, \\checkmark}$"
    assert formatter.is_lower(1.2, 1.0) == "$\\color{#d43e36}{\\, \\times}$"


def test_formatter_is_greater(formatter: Formatter) -> None:
    assert formatter.is_greater(1.2, 1.0) == "$\\color{#439b00}{\\, \\checkmark}$"
    assert formatter.is_greater(0.8, 1.0) == "$\\color{#d43e36}{\\, \\times}$"


def test_formatter_color_DCR_df(formatter: Formatter) -> None:
    df = pd.DataFrame({"DCR": [0.85, 0.95, 1.2]})
    styled_df = formatter.color_DCR_df(df, ["DCR"])
    assert (
        styled_df is not None
    )  # Simply check if the styled DataFrame is generated without errors


# Tests for TablePrinter class
def test_table_printer_data_output(table_printer: TablePrinter, capsys: Any) -> None:
    data = {"Column1": ["Value1", "Value2"], "Column2": ["a", "b"], "Column3": [5.5, 6.7], "Column4": ["cm", "mm"]}
    headers = ["Column1", "Column2", "Column3", "Column4"]
    table_printer.print_table_data(data, headers)
    captured = capsys.readouterr()
    assert "Value1" in captured.out
    assert "Value2" in captured.out
    assert "5.5" in captured.out
    assert "6.7" in captured.out


# Tests for DocumentBuilder class
def test_document_builder_initialization(document_builder: DocumentBuilder) -> None:
    assert document_builder.title == "Test Document"
    assert isinstance(document_builder.doc, type(Document()))


def test_document_builder_heading(document_builder: DocumentBuilder) -> None:
    document_builder.add_heading("Heading Level 1", level=1)
    assert len(document_builder.doc.paragraphs) > 0
    assert document_builder.doc.paragraphs[0].text == "Heading Level 1"


def test_document_builder_table(document_builder: DocumentBuilder) -> None:
    df = pd.DataFrame({"Name": ["Alice", "Bob"], "Age": [25, 30]})
    column_widths: List[Cm] = [Cm(5), Cm(2)]
    document_builder.add_table(df, column_widths)
    assert len(document_builder.doc.tables) > 0


def test_document_builder_save(
    document_builder: DocumentBuilder, tmp_path: Any
) -> None:
    filename = tmp_path / "test_document.docx"
    document_builder.save(str(filename))
    assert filename.exists()


if __name__ == "__main__":
    pytest.main()
