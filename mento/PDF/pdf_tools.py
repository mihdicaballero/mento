import os
import PyPDF3
import numpy as np
from PyPDF3 import PdfFileMerger
from pathlib import Path
from structurelab.PDF.page_numbers import add_page_numbers, insert_text_into_pdf, remove_pdf_file, rename_pdf_file


def extract_pdf_pages(directory, input_file, output_file, pages):
    """
    Extract specific pages from an input PDF and save them as a new PDF.

    Args:
        directory (str): The path to the directory containing the input and output PDF files.
        input_file (str): The name of the input PDF file to extract pages from.
        output_file (str): The desired name of the output PDF file.
        pages (list): A list of page numbers (0-based index) to be extracted from the input PDF.

    Returns:
        None

    Raises:
        FileNotFoundError: If the input file is not found.
    
    Usage:
        extract_pdf_pages('path/to/directory', 'input.pdf', 'output.pdf', [0, 2, 4])

    This function takes an input PDF, specified by its file name, and extracts the pages listed in the 'pages' parameter. It then saves the extracted pages as a new PDF with the specified 'output_file' in the same directory. The function is a practical tool for customizing PDFs by selecting specific pages to create tailored documents.
    """
    # Join the directory and file paths to create the full input file path
    input_path = os.path.join(directory, input_file)

    # Open the PDF file for reading
    with open(input_path, "rb") as file:
        # Create a PDF object
        pdf = PyPDF3.PdfFileReader(file)

        # Create a new PDF file for writing
        output = PyPDF3.PdfFileWriter()

        # Loop through the pages and add them to the output file
        for page_num in pages:
            page = pdf.getPage(page_num)
            output.addPage(page)

        # Join the directory and file paths to create the full output file path
        output_path = os.path.join(directory, output_file)

        # Save the output file
        with open(output_path, "wb") as outfile:
            output.write(outfile)

def merge_pdfs_by_folder(folder, bookmark=True):
    """
    Merge multiple PDFs from a folder into a single PDF file.

    Args:
        folder (str): The path to the folder containing the input PDF files.
        bookmark (bool): Optional. If True, bookmarks are added for each merged PDF using the PDF filename (default is True).

    Returns:
        None

    Usage:
        merge_pdfs_by_folder('input_folder', bookmark=True)

    This function combines multiple PDF files located in the specified 'folder' into a single output PDF. If 'bookmark' is set to True, bookmarks are generated for each merged PDF using the PDF filename. It's a practical utility for consolidating and indexing multiple PDF documents in a folder.
    """
    # Getting a list of all PDF files in the folder
    pdfs = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith('.pdf')]

    # Creating a pdf writer object
    pdf_writer = PyPDF3.PdfFileWriter()

    # Adding pages to the writer object
    for pdf in pdfs:
        pdf_reader = PyPDF3.PdfFileReader(pdf)
        for page in range(pdf_reader.numPages):
            pdf_writer.addPage(pdf_reader.getPage(page))

    # Setting bookmark names if enabled
    if bookmark:
        for i, pdf in enumerate(pdfs):
            bookmark_name = os.path.splitext(os.path.basename(pdf))[0]
            pdf_writer.addBookmark(title=bookmark_name, pagenum=i, parent=None)

    # Writing combined PDF to output file
    output = os.path.join(folder, 'merged_pdfs.pdf')
    with open(output, 'wb') as f:
        pdf_writer.write(f)

def merge_pdfs_by_name(folder, pdfs, output, bookmark=True):
    """
    Merge multiple PDFs by filename in a specified folder into a single PDF file.

    Args:
        folder (str): The path to the folder containing the input PDF files.
        pdfs (list): A list of PDF filenames (without the ".pdf" extension) to merge.
        output (str): The name of the output PDF file to be created.
        bookmark (bool): Optional. If True, bookmarks are added for each merged PDF using the PDF filename (default is True).

    Returns:
        None

    Usage:
        merge_pdfs_by_name('input_folder', ['document1', 'document2', 'document3'], 'merged.pdf', bookmark=True)

    This function combines multiple PDF files located in the specified 'folder' based on their filenames and creates a single output PDF in the same folder. If 'bookmark' is set to True, bookmarks will be generated for each merged PDF with the PDF filename. It's a handy tool for merging and organizing multiple PDF documents into a single, indexed file.
    """
    merger = PdfFileMerger(strict=False)

    for pdf in pdfs:
        input_file = os.path.join(folder, pdf + ".pdf")
        bookmark_name = os.path.splitext(os.path.basename(input_file))[0] if bookmark else None
        merger.append(fileobj=open(input_file, 'rb'), import_bookmarks=False, bookmark=bookmark_name)

    output_path = os.path.join(folder, output)
    merger.write(fileobj=open(output_path, 'wb'))
    merger.close()

def merge_template_to_pdf(template_pdf, input_folder, annex_folder, input_name, output_name):
    """
    Merge a template PDF with an existing PDF to create a customized PDF document.

    Args:
        template_pdf (str): The path to the template PDF file that contains content to be merged.
        input_folder (str): The folder containing the existing PDF file.
        annex_folder (str): The folder where the output PDF will be saved.
        input_name (str): The name of the existing PDF file within the 'input_folder'.
        output_name (str): The name of the output PDF file to be saved within the 'annex_folder'.

    Returns:
        str: The path to the newly created PDF with merged content.

    Usage:
        merge_template_to_pdf('template.pdf', 'input_folder', 'output_folder', 'input.pdf', 'output.pdf')

    This function combines a template PDF with an existing PDF, creating a new PDF document with customized content. It opens and reads the existing PDF and template PDF, iterates through the pages of the existing PDF, and merges each page with the template page. The resulting PDF is saved in the 'annex_folder' with the specified 'output_name'. The function is a useful tool for generating personalized PDFs by integrating template content into each page of the original document.
    """
    existing_pdf = os.path.join(input_folder, input_name)
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