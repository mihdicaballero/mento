from structurelab.PDF.merge_by_folder import merge_pdfs_by_folder

def main():
    # Input folder containing PDF files
    folder = r"D:\Desktop\Engineering"
    
    # Merging PDFs with bookmarks enabled
    merge_pdfs_by_folder(folder, bookmark=True)