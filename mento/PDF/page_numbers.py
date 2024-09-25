import os
from PyPDF3 import PdfFileReader, PdfFileWriter
from pathlib import Path
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from tqdm import tqdm

def remove_pdf_file(file_path):
    """
    Remove a specific PDF file.
    """
    if os.path.exists(file_path):
        os.remove(file_path)
        print(f"File '{file_path}' has been successfully removed.")
    else:
        print(f"File '{file_path}' does not exist.")

def rename_pdf_file(file_path, new_name):
    """
    Rename an existing PDF file.
    """
    if os.path.exists(file_path):
        new_file_path = os.path.join(os.path.dirname(file_path), new_name)
        os.rename(file_path, new_file_path)
        print(f"File '{file_path}' has been successfully renamed to '{new_name}'.")
    else:
        print(f"File '{file_path}' does not exist.")

def create_page_pdf(num, tmp, x=180, y=8, font_name="Helvetica", font_size=11, font_color="black"):
    c = canvas.Canvas(tmp, pagesize=A4)
    
    # Create a custom color using RGB values
    if isinstance(font_color, tuple) and len(font_color) == 3:
        r, g, b = font_color
        font_color = colors.Color(r/255, g/255, b/255)
    
    for i in range(1, num + 1):
        c.setFont(font_name, font_size)
        c.setFillColor(font_color)
        c.drawString(x * mm, y * mm, f"{i}-{num}")
        c.showPage()
    
    c.save()

def create_text_pdf(num, tmp, text_list, coordinates_list, font_name="Helvetica", font_size=11, font_color="black"):
    c = canvas.Canvas(tmp, pagesize=A4)
    
    # Create a custom color using RGB values
    if isinstance(font_color, tuple) and len(font_color) == 3:
        r, g, b = font_color
        font_color = colors.Color(r/255, g/255, b/255)
    
    for i in range(1, num + 1):
        c.setFont(font_name, font_size)
        
        for text, (x, y) in zip(text_list, coordinates_list):
            c.setFillColor(font_color)
            c.drawString(x * mm, y * mm, text)
        
        c.showPage()
    
    c.save()

def createPagePdfWithTotalFirst(num, tmp):
    c = canvas.Canvas(tmp)
    for i in range(1, num + 1):
        if i == 1:
            c.setFont('Helvetica',12)
            c.drawString(180 * mm, 33 * mm, str(num))
        else:
            c.setFont('Helvetica', 9)
            c.drawString(181 * mm, 13 * mm, 'Pág. ' + str(i) + '/' + str(num))
        c.showPage()
    c.save()

def insert_text_into_pdf(pdf_path, text_list, coordinates_list, font_name="Helvetica", font_size=11, font_color="black"):
    """
    Insert a list of texts at specified coordinates into an input PDF.
    """
    tmp_pdf = "__tmp.pdf"
    output_pdf = os.path.splitext(pdf_path)[0] + "_with_text.pdf"

    output = PdfFileWriter()
    with open(pdf_path, 'rb') as file:
        pdf = PdfFileReader(file, strict=False)
        num_pages = pdf.getNumPages()

        # Create new PDF with the inserted text
        create_text_pdf(num_pages, tmp_pdf, text_list, coordinates_list, font_name, font_size, font_color)

        with open(tmp_pdf, 'rb') as tmp_file:
            text_pdf = PdfFileReader(tmp_file)

            for page_num in range(num_pages):
                page = pdf.getPage(page_num)
                text_layer = text_pdf.getPage(page_num)

                # Merge the text layer with the actual page
                page.mergePage(text_layer)
                output.addPage(page)

        # Write the result
            if output.getNumPages() > 0:
                with open(output_pdf, 'wb') as new_file:
                    output.write(new_file)

        os.remove(tmp_pdf)
    return output_pdf

def add_page_numbers(pdf_path, x, y, font_name, font_size, font_color):
    """
    Add page numbers to a PDF and save the result as a new PDF.
    """
    tmp_pdf = "__tmp.pdf"
    output_pdf = os.path.splitext(pdf_path)[0] + "_with_page_numbers.pdf"

    output = PdfFileWriter()
    with open(pdf_path, 'rb') as file:
        pdf = PdfFileReader(file, strict=False)
        num_pages = pdf.getNumPages()

        # Create new PDF with page numbers
        create_page_pdf(num_pages, tmp_pdf, x, y, font_name, font_size, font_color)

        with open(tmp_pdf, 'rb') as tmp_file:
            number_pdf = PdfFileReader(tmp_file)

            for page_num in tqdm(range(num_pages), desc="Adding page numbers"):
                page = pdf.getPage(page_num)
                number_layer = number_pdf.getPage(page_num)

                # Merge number page with the actual page
                page.mergePage(number_layer)
                output.addPage(page)

        # Write the result
            if output.getNumPages() > 0:
                with open(output_pdf, 'wb') as new_file:
                    output.write(new_file)

        os.remove(tmp_pdf)
    return output_pdf

def add_page_numbers_but_first(pdf_path):
    """
    Add page numbers to a pdf, save the result as a new pdf
    @param pdf_path: path to pdf
    """
    tmp = "__tmp.pdf"

    output = PdfFileWriter()
    with open(pdf_path, 'rb') as f:
        pdf = PdfFileReader(f, strict=False)
        n = pdf.getNumPages()

        # create new PDF with page numbers
        create_page_pdf(n, tmp)

        with open(tmp, 'rb') as ftmp:
            numberPdf = PdfFileReader(ftmp)
            # iterarte pages
            for p in range(n):
                page = pdf.getPage(p)
                numberLayer = numberPdf.getPage(p)
                if p == 0:
                    output.addPage(page)
                    continue
                # merge number page with actual page
                page.mergePage(numberLayer)
                output.addPage(page)
            # write result
            if output.getNumPages():
                newpath = pdf_path[:-4] + "_numbered.pdf"
                with open(newpath, 'wb') as f:
                    output.write(f)
        os.remove(tmp)


def add_page_numbers_with_total_first(pdf_path):
    """
    Add page numbers to a pdf from second ot last page, with total in first page, save the result as a new pdf
    @param pdf_path: path to pdf
    """
    tmp = "__tmp.pdf"

    output = PdfFileWriter()
    with open(pdf_path, 'rb') as f:
        pdf = PdfFileReader(f, strict=False)
        n = pdf.getNumPages()
        # Get bookmarks from original PDF
        bookmarks = pdf.getOutlines()
        # create new PDF with page numbers
        createPagePdfWithTotalFirst(n, tmp)
        with open(tmp, 'rb') as ftmp:
            numberPdf = PdfFileReader(ftmp)
            # iterarte pages
            for p in range(n):
                page = pdf.getPage(p)
                numberLayer = numberPdf.getPage(p)
                # merge number page with actual page
                page.mergePage(numberLayer)
                output.addPage(page)
            # Add bookmarks
            for b in bookmarks:
                title = b['/Title']
                page = b['/Page'] + 1
                output.addBookmark(title, page, parent=None)  # add  bookmark
            # write result
            if output.getNumPages():
                newpath = pdf_path[:-4] + "_numbered.pdf"
                with open(newpath, 'wb') as f:
                    output.write(f)
        os.remove(tmp)


if __name__ == "__main__":
    # # Indicar el directorio de la carpeta con el PDF (va una r al principio para que funcione)
    source_dir = Path(
        r'C:\Users\Equipo 32\Desktop\Engineering')
    # Lista con todos los PDFs del directorio
    input_file = (source_dir / 'Test.pdf')
    # Leo el archivo PDF
    pdf = PdfFileReader(str(input_file), 'rb')
    # Cantidad total de páginas del PDF
    num_pages = pdf.getNumPages()
    # Example usage
    output_pdf = "numbered_pages.pdf"
    rgb_color = (50, 50, 50)  # Custom color in RGB values
    add_page_numbers(str(input_file),x=180, y=8, font_name="Helvetica", font_size=12, font_color=rgb_color)
