from structurelab.PDF.merge_by_folder import merge_pdfs_by_folder

def main():
    # Input folder containing PDF files
    folder = r"D:\Box\Proyectos\23031 Ingener_CMQX\6 Correspondencia\Salida\231012 APC MVQ Rev. 3 adelanto\PDF"
    
    # Merging PDFs with bookmarks enabled
    merge_pdfs_by_folder(folder, bookmark=True)