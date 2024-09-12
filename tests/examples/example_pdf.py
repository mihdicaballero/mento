import sys
import os

# Add the project root (where structurelab is located) to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from structurelab.PDF.pdf_tools import *

def main():
    ## EXAMPLE MERGE FOLDER
    # Input folder containing PDF files
    folder = r"D:\Desktop\Test"
    
    # Merging PDFs with bookmarks enabled
    merge_pdfs_by_folder(folder, bookmark=True)

    ## EXAMPLE MARGE BY NAME
    # Input folder containing PDF files
    folder = r"D:\Desktop\Test"
    # List of input PDFs
    pdfs = ['pdf1', 'pdf2', 'pdf3']

    # Output PDF file
    output = 'merged_pdfs.pdf'

    # Merging PDFs
    merge_pdfs_by_name(folder, pdfs, output, bookmark=True) 

    ## EXAMPLE EXTRACT PAGES
    directory = r"D:\Desktop\Test"
    input_pdf = "Pages from XXXXXXX"
    output_pdf = input_pdf + "_extract"

    # Define the range of pages to extract
    page_start = 1
    page_end = 3
    # Also define specific page numbers to extract
    pages_numbers = [5, 6]

    # Create an array of pages to extract
    page_range = np.arange(page_start - 1, page_end)
    pages_to_extract = np.unique(np.append(page_range, np.add(pages_numbers, -1)))

    # Add this ranges inside the method

    # Extract the specified pages from the PDF
    extract_pdf_pages(directory, input_pdf + ".pdf", output_pdf + ".pdf", pages_to_extract)

    ## EXAMPLE MERGE BACKGROUND TEMPLATE
        # Template PDF file with background
    template_pdf = Path(
        r"D:\Desktop\Engineering\background.pdf")   
    # Working folder containing PDF file from Calcpad
    calcpad_folder = r"D:\Desktop\Engineering"
    # Folder where the annexes are stored
    annex_folder = r"D:\Desktop\Engineering" 
    # Input PDF name. "Example.pdf"
    input_name = 'Input.pdf'
    # Output PDF name. "Example.pdf"
    output_name = "A-Nombre anexo.pdf"
    # Document information
    date = "16 de julio de 2023"
    project_name = "Nombre de proyecto"
    annex_name = "Nombre de anexo"

    # Add the following inside the method
    # Text definition, don't change.
    font_color = (50, 50, 50)  # Typical gray text color. 
    number_coords = [182, 10] # Text location for page numbers. Typical, D
    font_size = 12 # Text size, typical.
    coordinates_list_top = [(12, 285), (12, 279)] # First for project and second for annex name.
    coordinates_list_bot = [(12, 16), (12, 10)] # First for project and second for annex name.

    # Create the PDF with background, number pages and text.
    output_path = merge_template_to_pdf(template_pdf, calcpad_folder, annex_folder, input_name, output_name)
    output_path2 = add_page_numbers(output_path, number_coords[0], number_coords[1],"Helvetica", font_size, font_color)
    output_path3 = insert_text_into_pdf(output_path2, [project_name,annex_name], coordinates_list_top, "Helvetica", font_size, "white")
    output_path4 = insert_text_into_pdf(output_path3, ["Doc.: "+input_name,date], coordinates_list_bot, "Helvetica", font_size, font_color)    
    # remove temporary PDF files
    remove_pdf_file(output_path)
    remove_pdf_file(output_path2)
    remove_pdf_file(output_path3)
    rename_pdf_file(output_path4,output_name)

def calcpad():
        # Template PDF file with background
    template_pdf = Path(
        r"O:\Ingenieria\Memorias y Anexos\Anexos\background\background.pdf")   
    # Working folder containing PDF file from Calcpad
    calcpad_folder = r"P:\24016 SACEEM_Grua Bellevue\2 Calculos\01 Base grúa\Planillas"
    # Folder where the annexes are stored
    annex_folder = r"P:\24016 SACEEM_Grua Bellevue\4 Memoria y anexos\02 Memoria de cálculo\anexos" 
    # Input PDF name. "Example.pdf"
    input_name = '02-A-Diseño cabezal grua.pdf'
    # Output PDF name. "Example.pdf"
    output_name = "A-Cabezal_diseño.pdf"
    # Document information
    date = "05 de setiembre de 2024"
    project_name = "Edificio Bellevue"
    annex_name = "Cabezal de grúa"

    # Add the following inside the method
    # Text definition, don't change.
    font_color = (50, 50, 50)  # Typical gray text color. 
    number_coords = [182, 10] # Text location for page numbers. Typical, D
    font_size = 12 # Text size, typical.
    coordinates_list_top = [(12, 285), (12, 279)] # First for project and second for annex name.
    coordinates_list_bot = [(12, 16), (12, 10)] # First for project and second for annex name.

    # Create the PDF with background, number pages and text.
    output_path = merge_template_to_pdf(template_pdf, calcpad_folder, annex_folder, input_name, output_name)
    output_path2 = add_page_numbers(output_path, number_coords[0], number_coords[1],"Helvetica", font_size, font_color)
    output_path3 = insert_text_into_pdf(output_path2, [project_name,annex_name], coordinates_list_top, "Helvetica", font_size, "white")
    output_path4 = insert_text_into_pdf(output_path3, ["Doc.: "+input_name,date], coordinates_list_bot, "Helvetica", font_size, font_color)    
    # remove temporary PDF files
    remove_pdf_file(output_path)
    remove_pdf_file(output_path2)
    remove_pdf_file(output_path3)
    rename_pdf_file(output_path4,output_name)

if __name__ == "__main__":
    # main()
    calcpad()