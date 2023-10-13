import os
from PyPDF3 import PdfFileMerger


def merge_pdfs_by_name(folder, pdfs, output, bookmark=True):
    """
    Merge a list of PDF files specified by their names within the `folder` and save the combined result into the `output` file.
    If `bookmark` is True, bookmarks will be added to the output file to navigate directly to the input file sections.
    """
    merger = PdfFileMerger(strict=False)

    for pdf in pdfs:
        input_file = os.path.join(folder, pdf + ".pdf")
        bookmark_name = os.path.splitext(os.path.basename(input_file))[0] if bookmark else None
        merger.append(fileobj=open(input_file, 'rb'), import_bookmarks=False, bookmark=bookmark_name)

    output_path = os.path.join(folder, output)
    merger.write(fileobj=open(output_path, 'wb'))
    merger.close()

def main():
    # Input folder containing PDF files
    folder = r"D:\Desktop\PDF"
    # List of input PDFs
    pdfs = ['pdf1', 'pdf2', 'pdf3']

    # Output PDF file
    output = 'merged_pdfs.pdf'

    # Merging PDFs
    merge_pdfs_by_name(folder, pdfs, output, bookmark=True)

if __name__ == "__main__":
    main()
