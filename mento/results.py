import pandas as pd
from typing import Optional, List, Any
from tabulate import tabulate
from docx import Document
from docx.shared import Pt, Cm, RGBColor
import seaborn as sns
import matplotlib.pyplot as plt
from docx.oxml import parse_xml
from docx.oxml.ns import nsdecls

CUSTOM_COLORS = {
        'blue': '#1f77b4',  # Default Matplotlib blue
        'red': '#d62728',   # Default Matplotlib red
    }

# Default Matplotlib and Seaborn settings
def configure_plot_settings() -> None:
    """
    Configures global settings for Matplotlib and Seaborn plots.
    """
    # Basic style for the plot
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", rc={"lines.linewidth": 1.8})
    sns.set_style(rc={'axes.facecolor': '#F8F8F8'})

    # Matplotlib specific settings
    plt.rcParams.update({
        'text.usetex': True,
        'font.family': 'serif',
        'font.serif': ['Lato'],
        'axes.titlesize': 12,
        'axes.labelsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
    })

class Formatter:
    """
    Class to stlye the results in a Markdown display and a DataFrame output.
    """
    def __init__(self) -> None:
        # Define colors in the constructor
        self.green = '#439b00'
        self.red = '#d43e36'
        self.yellow = '#efc200'
        self.mid_value=0.95
        self.max_value=1
    
    def DCR(self, DCR: float) -> str:
        # Determine color based on DCR value
        if self.mid_value > DCR:
            color = self.green
        elif self.mid_value <= DCR <= self.max_value:
            color = self.yellow
        else:
            color = self.red
        
        return f"$\\color{{{color}}}{{\\text{{DCR}}={round(DCR,2)}}}$"
    
    def DCR_value(self, DCR: float) -> str:
        # Determine color based on DCR value
        if self.mid_value > DCR:
            color = self.green
        elif self.mid_value <= DCR <= 1:
            color = self.yellow
        else:
            color = self.red
        
        return f"$\\color{{{color}}}{{{round(DCR,2)}}}$"
    
    def is_lower(self, value1: float, value2: float) -> str:
        # Compare two values and return the appropriate formatted output
        if value1 < value2:
            return r"$\color{" + self.green + r"}{\, \checkmark}$"  # Green checkmark
        else:
            return r"$\color{" + self.red + r"}{\, \times}$"  # Red cross
    
    def is_greater(self, value1: float, value2: float) -> str:
        # Compare two values and return the appropriate formatted output
        if value1 > value2:
            return r"$\color{" + self.green + r"}{\, \checkmark}$"  # Green checkmark
        else:
            return r"$\color{" + self.red + r"}{\, \times}$"  # Red cross
    # Formatting functions outside the class

    def DCR_value_df(self,DCR: Optional[float]) -> str:
        # If DCR is None or not a number, return an empty string
        if DCR is None or not isinstance(DCR, (int, float)):
            return ""

        # Determine the color based on DCR value
        if self.mid_value > DCR:
            return f"color: {self.green}"
        elif self.mid_value <= DCR <= self.max_value:
            return f"color: {self.yellow}"
        else:
            return f"color: {self.red}"

    def apply_DCR_style(self,value: float) -> str:
        formatted_value = round(value, 2) if isinstance(value, (int, float)) else value
        return self.DCR_value_df(formatted_value)

    def color_DCR_df(self, df: pd.DataFrame, DCR_columns: list) -> pd.io.formats:
        """
        Apply color styling to specified DCR-related columns in the DataFrame.
        
        :param df: DataFrame with DCR values to color.
        :param DCR_columns: List of column names to apply the DCR styling to.
        :return: A styled DataFrame with colored DCR values.
        """
        return df.style.map(self.apply_DCR_style, subset=DCR_columns).format(precision=2)

class TablePrinter:
    """
    A class for printing tables with customizable formatting options.
    
    Attributes:
    -----------
    title : Optional[str]
        Optional title displayed above the printed table.

    Methods:
    --------
    print_table_data(data: List[List[Any]], headers: List[str], tablefmt: str = "fancygrid", numalign: str = "right")
      -> None
        Prints table data with customizable formatting.
        
    print_table_min_max(data: List[List[Any]], headers: List[str], tablefmt: str = "fancygrid", numalign: str = "right")
      -> None
        Prints table data with column alignment for minimum and maximum values.
    """
    
    def __init__(self, title: Optional[str] = None) -> None:
        """
        Initializes the TablePrinter with an optional title.
        
        Parameters:
        -----------
        title : Optional[str], default=None
            Title to be displayed above the table, if provided.
        """
        self.title = title

    def print_table_data(self, data: dict[str, list[Any]], headers: str, tablefmt: str = "fancygrid",
                         numalign: str = "right") -> None:
        """
        Prints table data with customizable formatting options.
        
        Parameters:
        -----------
        data : List[List[Any]]
            The data to be printed in table format, where each inner list represents a row.
        
        headers : List[str]
            Column headers for the table.
        
        tablefmt : str, default="fancygrid"
            Table formatting style supported by tabulate (e.g., "plain", "grid", "fancygrid").
        
        numalign : str, default="right"
            Number alignment in the table. Common options are "right", "center", or "left".
        
        Returns:
        --------
        table:  Return the formatted table string
        """
        # if self.title:
        #     print(f"======= {self.title} ======= \n")
        
        colalign = ("left", "center", "right", "left")
        table = tabulate(
            data,
            headers=headers,
            tablefmt=tablefmt,
            numalign=numalign,
            colalign=colalign
        )
        print(table, '\n')
        return None  
        
    def print_table_min_max(self, data: dict[str, list[Any]], headers: str, tablefmt: str = "fancygrid",
                            numalign: str = "right") -> None:
        """
        Prints table data with column alignment for minimum and maximum values.
        
        Parameters:
        -----------
        data : List[List[Any]]
            The data to be printed in table format, where each inner list represents a row.
        
        headers : List[str]
            Column headers for the table.
        
        tablefmt : str, default="fancygrid"
            Table formatting style supported by tabulate (e.g., "plain", "grid", "fancygrid").
        
        numalign : str, default="right"
            Number alignment in the table. Common options are "right", "center", or "left".
        
        Returns:
        --------
        None
        """
        if self.title:
            print(f"======= {self.title} ======= \n")
        
        colalign = ("left", "center", "center", "center", "center", "center", "left")
        table = tabulate(
            data,
            headers=headers,
            tablefmt=tablefmt,
            numalign=numalign,
            colalign=colalign
        )
        print(table, '\n')

class DocumentBuilder:
    """
    A class to build and style a Word document, including adding headings and tables from data frames.
    
    Attributes:
    -----------
    title : str
        Title of the document, displayed at the beginning.
        
    font_name : str
        Default font name for the document text.
        
    font_size : int
        Default font size for the document text.
        
    doc : Document
        The Document object representing the Word document being built.

    Methods:
    --------
    set_document_style() -> None
        Configures the document's default style with the specified font name and size.
        
    add_heading(text: str, level: int) -> None
        Adds a heading to the document at the specified level.
        
    set_col_widths(table: 'docx.table.Table', column_widths: List[Cm]) -> None
        Sets the width for each column in a given table.
        
    add_table(df: pd.DataFrame, column_widths: List[Cm]) -> None
        Adds a table to the document, with data from a DataFrame and custom column widths.
        
    save(filename: str) -> None
        Saves the document to a specified filename.
    """
    
    def __init__(self, title: str, font_name: str = 'Lato', font_size: int = 9) -> None:
        """
        Initializes the DocumentBuilder with a title, font name, and font size.
        
        Parameters:
        -----------
        title : str
            Title to be displayed in the document.
            
        font_name : str, default='Lato'
            Font name to be used for the document text.
            
        font_size : int, default=9
            Font size for the document text.
        """
        self.doc = Document()
        self.title = title
        self.font_name = font_name
        self.font_size = font_size
        self.set_document_style()

    def set_document_style(self) -> None:
        """
        Sets the default style of the document, applying the font name and size.
        
        Returns:
        --------
        None
        """
        style = self.doc.styles['Normal']
        font = style.font
        font.name = self.font_name
        font.size = Pt(self.font_size)

    def add_heading(self, text: str, level: int) -> None:
        """
        Adds a heading to the document at the specified level.
        
        Parameters:
        -----------
        text : str
            The text for the heading.
            
        level : int
            The heading level (e.g., 0 for title, 1 for main headings, etc.).
        
        Returns:
        --------
        None
        """
        heading = self.doc.add_heading(text, level=level)
        heading.paragraph_format.space_before = Pt(0)  # Remove space after the paragraph

    def add_text(self, text: str) -> None:
        """ Adds a paragraph to the document """
        self.doc.add_paragraph(text)
        
    def set_col_widths(self, table: Any, column_widths: List[Cm]) -> None:
        """
        Sets the width of each column in the table.
        
        Parameters:
        -----------
        table : docx.table.Table
            The table object whose column widths need to be set.
            
        column_widths : List[Cm]
            List of widths for each column, in centimeters.
        
        Returns:
        --------
        None
        """
        for row in table.rows:
            for idx, width in enumerate(column_widths):
                row.cells[idx].width = width

    def add_table(self, df: pd.DataFrame, column_widths: List[Cm]) -> None:
        """
        Adds a table to the document, populated with data from a DataFrame.
        
        Parameters:
        -----------
        df : pd.DataFrame
            The DataFrame containing the data for the table.
            
        column_widths : List[Cm]
            List of column widths for the table.
        
        Returns:
        --------
        Table
        """
        # Initialize the table with the DataFrame dimensions
        table = self.doc.add_table(rows=df.shape[0] + 1, cols=df.shape[1])
        # Apply table style
        table.style = "Light Shading"
        # Set column widths
        self.set_col_widths(table, column_widths)
        # Populate column headers
        for j in range(df.shape[1]):
            table.cell(0, j).text = str(df.columns[j])
        # Populate table data
        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                cell = df.iat[i, j]
                table.cell(i + 1, j).text = str(cell)
        # Customize style for first column
        for row in table.rows[1:]:
            cell = row.cells[0]
            run = cell.paragraphs[0].runs[0]
            run.font.bold = False
        # Add a spacer paragraph after the table
        spacer = self.doc.add_paragraph()
        spacer.paragraph_format.space_after = Pt(0)  # Remove space after the paragraph

        return None

    def add_table_data(self, df: pd.DataFrame) -> None:
        column_widths = [Cm(12), Cm(2), Cm(2), Cm(2)]
        self.add_table(df, column_widths)

    def add_table_dcr(self, df: pd.DataFrame) -> None:
        """
        Adds a table using `add_table` and modifies the last row's DCR cell color.

        - Green if DCR < 1
        - Red if DCR >= 1
        """
        column_widths = [Cm(12), Cm(2), Cm(2), Cm(2)]
        self.add_table(df, column_widths)
        # Get the last table added to the document
        table = self.doc.tables[-1]  # Retrieve the most recent table

        # Apply color to DCR cell (third column of the last row)
        last_row_idx = df.shape[0]  # Last row index
        dcr_column_idx = 2  # Third column index (zero-based)

        dcr_value = float(df.iat[-1, dcr_column_idx])  # Extract DCR value
        if dcr_value < 1:
            shading_color = "C6EFCE"  # Green
            font_color = "006100"  # Green font
        else:
            shading_color = "FFC7CE"  # Red
            font_color = "9C0006"  # Red font

        # Apply shading and font color to the entire last row
        for cell in table.rows[last_row_idx].cells:
            # Apply background color
            cell._element.get_or_add_tcPr().append(
                parse_xml(f'<w:shd {nsdecls("w")} w:fill="{shading_color}"/>')
            )
            # Apply font color
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.color.rgb = RGBColor.from_string(font_color)
    
    def add_table_min_max(self, df: pd.DataFrame) -> None:
        column_widths = [Cm(10), Cm(2), Cm(2), Cm(2), Cm(2), Cm(1)]
        self.add_table(df, column_widths)


    def save(self, filename: str) -> None:
        """
        Saves the document to a specified file.
        
        Parameters:
        -----------
        filename : str
            The filename (including path) where the document will be saved.
        
        Returns:
        --------
        None
        """
        self.doc.save(filename)



    #  Examples to run in a Jupyter Notebook
    # formatter = Formatter()
    # display(Markdown(formatter.DCR(0.85)))
    # display(Markdown(formatter.DCR_value(0.85)))
    # display(Markdown(formatter.is_lower(0.85,1)))
    # display(Markdown(formatter.is_greater(0.85,1)))
