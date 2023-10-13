import os
import PyPDF3
from pathlib import Path
from PDF_page_numbers import add_page_numbers, insert_text_into_pdf, remove_pdf_file, rename_pdf_file

'''
A MEJORAR:
- Rellenar info de header y footer desde el script.
'''


def merge_template_to_pdf(template_pdf, calcpad_folder, annex_folder, input_name, output_name):
    existing_pdf = os.path.join(calcpad_folder, input_name)
    output_pdf = os.path.join(annex_folder, output_name)

    with open(existing_pdf, "rb") as existing_file, open(template_pdf, "rb") as template_file:
        existing_reader = PyPDF3.PdfFileReader(existing_file)
        template_reader = PyPDF3.PdfFileReader(template_file)
        template_page = template_reader.getPage(0)

        output_writer = PyPDF3.PdfFileWriter()

        for page_num in range(existing_reader.getNumPages()):
            page = existing_reader.getPage(page_num)
            page.mergePage(template_page)
            output_writer.addPage(page)

        with open(output_pdf, "wb") as output_file:
            output_writer.write(output_file)
    return os.path.join(annex_folder, output_name)

def main():
    # Template PDF file with background
    template_pdf = Path(
        r"O:\Ingeniería\Memorias y Anexos\Template\background\background.pdf")   
    # Working folder containing PDF file from Calcpad
    calcpad_folder = r"P:\23029 Agustin Murgia_Verificación andamios\2 Cálculos\01 Andamio\Calcpad"
    # Folder where the annexes are stored
    annex_folder = r"P:\23029 Agustin Murgia_Verificación andamios\4 Memoria y anexos\02 Memoria de cálculo\anexos" 
    # Input PDF name. "Example.pdf"
    input_name = 'ACC_Viento_UNIT 50-84 General.pdf'
    # Output PDF name. "Example.pdf"
    output_name = "A-Nombre anexo.pdf"
    # Document information
    date = "16 de julio de 2023"
    project_name = "Nombre de proyecto"
    annex_name = "Nombre de anexo"

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
    main()
