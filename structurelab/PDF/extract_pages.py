import os
import PyPDF3
import numpy as np


def extract_pdf_pages(directory, input_file, output_file, pages):
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


def main():
    # Example usage
    directory = r"C:\Users\Equipo 32\OneDrive - CHM Structural Engineers\Peachtree - Society Atlanta\RFI's\RFI 118 Bathroom and Exterior Signange"
    input_pdf = "Pages from Society Atlanta - CB Specialties Submittals 11212022"
    output_pdf = input_pdf + "_extract"

    # Define the range of pages to extract
    page_start = 1
    page_end = 3
    # Also define specific page numbers to extract
    pages_numbers = [5, 6]

    # Create an array of pages to extract
    page_range = np.arange(page_start - 1, page_end)
    pages_to_extract = np.unique(np.append(page_range, np.add(pages_numbers, -1)))

    # Extract the specified pages from the PDF
    extract_pdf_pages(directory, input_pdf + ".pdf", output_pdf + ".pdf", pages_to_extract)


if __name__ == "__main__":
    main()
