from pathlib import Path
import sys,json,subprocess,csv,hashlib,gc
from concurrent.futures import ThreadPoolExecutor
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(__file__).resolve().parent/'_result'
skill=Path(sys.argv[1]);qa=out/'07_qa';qa.mkdir(exist_ok=True)
def run_audit(pdf):
    target=qa/'collision'/pdf.parent.name;target.mkdir(parents=True,exist_ok=True)
    p=subprocess.run([sys.executable,str(skill/'scripts/audit_figure_collisions.py'),str(pdf),'--json-out',str(target/(pdf.stem+'.json'))],capture_output=True,text=True,encoding='utf-8',errors='replace')
    t=subprocess.run([sys.executable,str(skill/'scripts/audit_pdf_text.py'),str(pdf),'--min-pt','5','--json'],capture_output=True,text=True,encoding='utf-8',errors='replace')
    (target/(pdf.stem+'.text.json')).write_text(t.stdout,encoding='utf-8')
    return {'figure':str(pdf.relative_to(out)),'collision_exit':p.returncode,'font_exit':t.returncode,'collision_summary':p.stdout[-1600:]}
pdfs=sorted((out/'04_figures').rglob('*.pdf'))
with ThreadPoolExecutor(max_workers=4) as pool: reports=list(pool.map(run_audit,pdfs))
(qa/'figure_audit.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2),encoding='utf-8')
s=subprocess.run([sys.executable,str(skill/'scripts/validate_figure.py'),str(out/'01_model/export_result.py'),'--json'],capture_output=True,text=True,encoding='utf-8',errors='replace')
(qa/'source_audit.json').write_text(s.stdout,encoding='utf-8')
images=[p.with_suffix('.png') for p in pdfs]
for start in range(0,len(images),6):
    fig,axes=plt.subplots(3,2,figsize=(15,11));fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.98,hspace=.06,wspace=.03)
    for ax in axes.flat:ax.axis('off')
    for ax,p in zip(axes.flat,images[start:start+6]):
        ax.imshow(plt.imread(p));ax.set_title(p.stem,fontsize=8,pad=1)
    fig.savefig(qa/f'contact_{start//6+1:02d}.png',dpi=130);plt.close(fig);gc.collect()
print(json.dumps({'figures':len(reports),'collision_failed':[r['figure'] for r in reports if r['collision_exit']==1],
 'audit_errors':[r['figure'] for r in reports if r['collision_exit'] not in [0,1]],
 'font_failed':[r['figure'] for r in reports if r['font_exit']!=0]},ensure_ascii=False))
