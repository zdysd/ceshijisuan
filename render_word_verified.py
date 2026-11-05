"""Use the skill renderer with native Word PDF export on this Windows host.

The packaged LibreOffice path was attempted and diagnosed as unavailable.
All OOXML preparation, PDF rasterization and page naming remain in render_docx.
"""
import importlib.util,sys
from pathlib import Path
import win32com.client

SKILL=Path('E:/codex/Codex_Data/plugins/cache/openai-primary-runtime/documents/26.1011.11113/skills/documents/render_docx.py')
spec=importlib.util.spec_from_file_location('packaged_docx_renderer',SKILL)
renderer=importlib.util.module_from_spec(spec);spec.loader.exec_module(renderer)


def word_pdf(doc_path,user_profile,convert_tmp_dir,stem,verbose):
    target=Path(convert_tmp_dir)/(stem+'.pdf')
    app=win32com.client.DispatchEx('Word.Application');document=None
    try:
        app.Visible=False;app.DisplayAlerts=0
        document=app.Documents.Open(str(Path(doc_path).resolve()),ReadOnly=True,AddToRecentFiles=False,Visible=False)
        document.Repaginate();document.Fields.Update()
        document.ExportAsFixedFormat(str(target),17)
    finally:
        if document is not None:document.Close(False)
        app.Quit()
    return str(target),'Native Microsoft Word PDF export; LibreOffice unavailable on this host.'


renderer.convert_to_pdf=word_pdf
if __name__=='__main__':renderer.main()
