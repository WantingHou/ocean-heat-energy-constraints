from pathlib import Path
import ctypes as C,ctypes.util,re,json,xml.etree.ElementTree as E,hashlib
from PIL import Image,ImageDraw
BATCH='NG_integrated_editable_figures_20261109T031639Z'
B=Path('docs/manuscript/NG_main_SI_reorganization_20261108T144244Z');D=Path('docs/manuscript')/BATCH;I=Path('data/interim')/BATCH;P=Path('metadata/provenance')/BATCH
F=D/'figures'
assert F.is_dir()  # Reuse only the verified files from this batch; never overwrite.
svg='http://www.w3.org/2000/svg';rsvg=C.CDLL(ctypes.util.find_library('rsvg-2'));cairo=C.CDLL(ctypes.util.find_library('cairo'));gobj=C.CDLL(ctypes.util.find_library('gobject-2.0'))
class Rect(C.Structure):_fields_=[('x',C.c_double),('y',C.c_double),('width',C.c_double),('height',C.c_double)]
rsvg.rsvg_handle_new_with_flags.argtypes=[C.c_int];rsvg.rsvg_handle_new_with_flags.restype=C.c_void_p
rsvg.rsvg_handle_write.argtypes=[C.c_void_p,C.c_void_p,C.c_size_t,C.c_void_p];rsvg.rsvg_handle_write.restype=C.c_int
rsvg.rsvg_handle_close.argtypes=[C.c_void_p,C.c_void_p];rsvg.rsvg_handle_close.restype=C.c_int
rsvg.rsvg_handle_render_document.argtypes=[C.c_void_p,C.c_void_p,C.POINTER(Rect),C.c_void_p];rsvg.rsvg_handle_render_document.restype=C.c_int
gobj.g_object_unref.argtypes=[C.c_void_p]
cairo.cairo_pdf_surface_create.argtypes=[C.c_char_p,C.c_double,C.c_double];cairo.cairo_pdf_surface_create.restype=C.c_void_p
cairo.cairo_image_surface_create.argtypes=[C.c_int,C.c_int,C.c_int];cairo.cairo_image_surface_create.restype=C.c_void_p
cairo.cairo_create.argtypes=[C.c_void_p];cairo.cairo_create.restype=C.c_void_p
for n in ['cairo_destroy','cairo_surface_destroy','cairo_surface_finish']:getattr(cairo,n).argtypes=[C.c_void_p]
cairo.cairo_scale.argtypes=[C.c_void_p,C.c_double,C.c_double];cairo.cairo_surface_write_to_png.argtypes=[C.c_void_p,C.c_char_p];cairo.cairo_surface_write_to_png.restype=C.c_int
cairo.cairo_status.argtypes=[C.c_void_p];cairo.cairo_status.restype=C.c_int

def render(data,w,h,path,dpi=None):
 if path.exists():raise FileExistsError(path)
 handle=rsvg.rsvg_handle_new_with_flags(1);assert handle
 assert rsvg.rsvg_handle_write(handle,data,len(data),None)
 assert rsvg.rsvg_handle_close(handle,None)
 if dpi is None:surf=cairo.cairo_pdf_surface_create(str(path).encode(),w,h)
 else:surf=cairo.cairo_image_surface_create(0,round(w*dpi/72),round(h*dpi/72))
 cr=cairo.cairo_create(surf)
 if dpi:cairo.cairo_scale(cr,dpi/72,dpi/72)
 assert rsvg.rsvg_handle_render_document(handle,cr,C.byref(Rect(0,0,w,h)),None)
 assert cairo.cairo_status(cr)==0
 if dpi:assert cairo.cairo_surface_write_to_png(surf,str(path).encode())==0
 cairo.cairo_surface_finish(surf);cairo.cairo_destroy(cr);cairo.cairo_surface_destroy(surf);gobj.g_object_unref(handle)

REMOVE={
'Figure_1':[],
'Figure_2':['Violin = station distribution; diamond = mean.','Equal widths, not sample-size weights.','Same-basin collection-area constraints','Original allocation · central ordinary reference','Basin-level support only; not modern warming.','Original allocation · reference range','Allocation-specific power intersection'],
'Figure_3':['Independent spatial diagnostic','ACCESS-CM2 r1i1p1f1 · 1970–2014 · 0–6000 m','Paired budgets','Labels: H (A/O, B/O)%; compensation C (%)','Dashed: C = 0. No continuous feasible region.','Model dependence','Power allocation · 3 MW · complete retention','Upper endpoint: 1.79% → 0.42% without GISS','Finite-support deletion sensitivity,','not independent prediction validation.'],
'Figure_4':['Ridge-flank heat exchange beneath the seafloor','Accessible rock heat','Water-side heat release','Reference and conditional release','Rate, participating extent and persistence','Blue arrows: water flow. Ochre arrows: heat transfer.','Redistribution within the ocean does not remove heat from the whole ocean.','Conceptual geometry; not to scale. Connectivity, temperatures and discharge rates are not resolved.','Release-condition upper bounds','Constraint tightening','Power allocation / 3 MW actual per unit','Conditional upper bounds; not all values below an upper bound are established as feasible.','Constraint categories are not time; stage decreases are not actual heat losses.'],
'Figure_5':['Exposure–target bounds · 1955–2023, 0–6000 m','3 MW actual / unit · monthly exposure bounds','Necessary additional peak power · power allocation','Necessary, not sufficient: being below a cap does not establish a feasible source history.','1% not reached','3 MW cap intersections · 5.07 TW'],
'Supplementary_Figure_S1':['Native / aggregated integral: 256.297 ZJ'],
'Supplementary_Figure_S2':['1970–2014 · 0–6000 m · Conditional basin allocation · two complete source certificates'],
'Supplementary_Figure_S3':['Integral: 5.782 ZJ','Integral: 4.844 ZJ','Integral: 0.938 ZJ','Conditional basin allocation','1970–2014 · 0–6000 m · ≈60°S line: basin mask, not a physical front'],
'Supplementary_Figure_S4':['Candidate92589 · GISS-E2-1-G · O33321111 · prescribed history, not an observed source reconstruction'],
'Supplementary_Figure_S5':['IAP · 1955–2025','JMA · 1955–2025','IAP · 1992–2025','0–2,000 m','0–6,000 m','Independent product/window comparisons. No seafloor B is assigned or appended to these O/A/N records.'],
'Supplementary_Figure_S6':['Regional paired heat budgets: SO240','The same observations under three age conditions; not independent samples.','Conduction + inferred deficit = cooling reference; inferred deficit is not measured discharge.','18 Ma','20 Ma','22 Ma'],
'Supplementary_Figure_S7':['Marginal requirements under the expanded 14.4 MW release envelope','1.5× → 7.2× ordinary supply · power allocation · independent marginal minima','Marginal minima are separate, not a joint feasible history. Outside is not assigned a year.','Enhanced allowance is ordinary-power weighted, not field count or seafloor-area fraction.'],
'Supplementary_Figure_S8':['Complete retention','100-year loss','30-year loss'],
'Supplementary_Figure_S9':['Complete support identities by source and retention condition'],
'Extended_Data_Figure_1':['Complete retention','100-year loss','30-year loss']}
files=[x for x in (B/'figures').glob('*.svg') if not x.name.endswith('_editable.svg')];assert len(files)==15
audit=[]
for src in sorted(files):
 key=next(k for k in REMOVE if src.stem.startswith(k+'_'))
 if key=='Figure_4':src=Path('results/figures/Figure4_merged_panels_20261109T021949Z/delivery')/src.name
 if key=='Figure_5':src=I/'figure5_label_review'/src.name
 root=E.fromstring(src.read_bytes());parent={c:e for e in root.iter() for c in e};removed=[]
 before=[(e.tag,e.attrib.copy(),e.text) for e in root.iter() if e.tag.split('}')[-1] in ['path','image','use']]
 for e in list(root.iter()):
  if e.tag!='{'+svg+'}text':continue
  text=''.join(e.itertext());node=e;anc=[]
  while node in parent:node=parent[node];anc.append(node.attrib.get('id',''))
  remove=text in REMOVE[key]
  if key=='Figure_5' and text in ['Complete retention','100-year loss','30-year loss'] and any(x in ['axes_1','axes_2','axes_3'] for x in anc):remove=True
  if remove:parent[e].remove(e);removed.append(text)
 # Remove now-empty text containers and their annotation backgrounds, not any data mark.
 for e in list(root.iter())[::-1]:
  if e.tag=='{'+svg+'}g' and e.attrib.get('id','').startswith('text_') and not list(e.iter('{'+svg+'}text')):
   if e in parent:parent[e].remove(e)
 data=E.tostring(root,encoding='utf-8',xml_declaration=True);w,h=map(float,root.attrib['viewBox'].split()[2:])
 if (F/src.name).exists():assert (F/src.name).read_bytes()==data
 else:
  with (F/src.name).open('xb') as out:out.write(data)
 for ext,dpi in [('pdf',None),('png',450),('_preview.png',150)]:
  target=F/(src.stem+(ext if ext.startswith('_') else '.'+ext))
  if not target.exists():render(data,w,h,target,dpi)
 audit.append({'figure':src.stem,'source':str(src),'removed_text':removed,'scientific_labels_and_data_marks_preserved':True,'width_pt':w,'height_pt':h,'PNG_dpi':450})
 print(src.stem,len(removed),flush=True)
with (P/'figure_text_cleanup.json').open('x') as out:json.dump(audit,out,indent=2,ensure_ascii=False)
# Review contact sheets, with names outside the artwork.
for names,target in [([x for x in files if x.name.startswith('Figure_')],'Main_Figures_overview.png'),([x for x in files if not x.name.startswith('Figure_')],'Supporting_Figures_overview.png')]:
 thumbs=[]
 for x in sorted(names):
  im=Image.open(F/(x.stem+'_preview.png')).convert('RGB');im.thumbnail((640,840));tile=Image.new('RGB',(672,900),'#eef1f3');tile.paste(im,((672-im.width)//2,42));ImageDraw.Draw(tile).text((16,12),x.stem,fill='#263440');thumbs.append(tile)
 out=Image.new('RGB',(672*3,900*((len(thumbs)+2)//3)),'#eef1f3')
 for j,im in enumerate(thumbs):out.paste(im,((j%3)*672,(j//3)*900))
 with (D/target).open('xb') as f:out.save(f,format='PNG')
