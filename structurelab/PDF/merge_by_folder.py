import os
import PyPDF3


def merge_pdfs_by_folder(folder, bookmark=True):
    """
    Merge all PDF files within the `folder` and save the combined result into the 'merged_pdfs.pdf' file.
    If `bookmark` is True, bookmarks will be added to the output file to navigate directly to the input file sections.
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


def main():
    # Input folder containing PDF files
    folder = r"D:\Desktop\Engineering"
    
    # Merging PDFs with bookmarks enabled
    merge_pdfs_by_folder(folder, bookmark=True)


if __name__ == "__main__":
    main()
