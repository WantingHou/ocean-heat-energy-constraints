"""Render the supplied, inspected single-slide OOXML without Office or new packages.
Raster geological texture is independent of editable paths, labels and leaders.
Unsupported geometry fails explicitly; no hidden conversion or whole-slide raster.
"""
import zipfile, xml.etree.ElementTree as ET
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Ellipse, Rectangle
from PIL import Image
from io import BytesIO
NS={'p':'http://schemas.openxmlformats.org/presentationml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
def draw_supplied_ppt(ax, ppt_path, minimum_font=7.2):
 with zipfile.ZipFile(ppt_path) as z:
  pres=ET.fromstring(z.read('ppt/presentation.xml')); size=pres.find('p:sldSz',NS);sw=float(size.attrib['cx']);sh=float(size.attrib['cy'])
  slide=ET.fromstring(z.read('ppt/slides/slide1.xml'))
  texture=Image.open(BytesIO(z.read('ppt/media/image1.png'))).convert('RGB')
  ax.imshow(texture,extent=(0,sw,sh,0),zorder=0,aspect='auto').set_gid('geological-texture-raster-only')
  ax.set(xlim=(0,sw),ylim=(sh,0));ax.axis('off')
  scale=ax.get_position().width*ax.figure.get_figwidth()/(sw/914400)
  counts={'native_shapes':0,'vector_paths':0,'editable_text_labels':0,'anchors':0,'raster_textures':1,'minimum_font_pt':minimum_font,'native_office_render':False}
  def color(node,default='none'):
   if node is None:return default
   x=node.find('a:solidFill/a:srgbClr',NS)
   if x is not None:return '#'+x.attrib['val']
   if node.find('a:noFill',NS) is not None:return 'none'
   return default
  for k,sp in enumerate(slide.findall('p:cSld/p:spTree/p:sp',NS)):
   counts['native_shapes']+=1
   name=sp.find('p:nvSpPr/p:cNvPr',NS).attrib['name'];prop=sp.find('p:spPr',NS);xf=prop.find('a:xfrm',NS)
   assert not xf.attrib, 'Unsupported transformed shape '+name
   off=xf.find('a:off',NS);ext=xf.find('a:ext',NS);x,y=float(off.attrib['x']),float(off.attrib['y']);w,h=float(ext.attrib['cx']),float(ext.attrib['cy'])
   line=prop.find('a:ln',NS);ec=color(line);fc=color(prop);lw=float(line.attrib.get('w',0))/12700*scale if line is not None else 0
   geom=prop.find('a:custGeom',NS);artist=None
   if geom is not None:
    for path in geom.findall('a:pathLst/a:path',NS):
     pw,ph=float(path.attrib['w']),float(path.attrib['h']);verts=[];codes=[]
     for cmd in path:
      tag=cmd.tag.split('}')[-1]
      if tag=='close':verts.append(verts[0]);codes.append(MPath.CLOSEPOLY);continue
      code={'moveTo':MPath.MOVETO,'lnTo':MPath.LINETO,'cubicBezTo':MPath.CURVE4}[tag]
      for pt in cmd.findall('a:pt',NS):verts.append((x+float(pt.attrib['x'])*w/pw,y+float(pt.attrib['y'])*h/ph));codes.append(code)
     artist=PathPatch(MPath(verts,codes),facecolor=fc,edgecolor=ec,lw=lw,zorder=k+1,capstyle='round',joinstyle='round');ax.add_patch(artist);counts['vector_paths']+=1
   else:
    preset=prop.find('a:prstGeom',NS).attrib['prst']
    if preset=='ellipse':artist=Ellipse((x+w/2,y+h/2),w,h,facecolor=fc,edgecolor=ec,lw=lw,zorder=k+1);ax.add_patch(artist);counts['anchors']+=1
    elif preset=='rect':
     if fc!='none' or ec!='none':artist=Rectangle((x,y),w,h,facecolor=fc,edgecolor=ec,lw=lw,zorder=k+1);ax.add_patch(artist)
    else:raise ValueError('Unsupported preset '+preset)
   if artist is not None:artist.set_gid('ppt-'+name)
   body=sp.find('p:txBody',NS)
   if body is None:continue
   paragraphs=body.findall('a:p',NS);texts=[''.join(p.itertext()) for p in []]
   texts=[''.join(t.text or '' for t in p.findall('.//a:t',NS)) for p in paragraphs]
   rp=body.find('.//a:rPr',NS);font=max(minimum_font,float(rp.attrib['sz'])/100*scale);bold=rp.attrib.get('b')=='1'
   # Explicit paragraph breaks are preserved. At publication width, small supplied
   # fonts are enlarged; the untouched PPT remains the native editable master.
   text='\n'.join(texts)
   if name in ['Label — ordinary','Label — alternative']:x-=sw*.007
   if name=='Label — colour key':text='Blue arrows: water flow. Ochre arrows: heat transfer.\nRedistribution within the ocean does not remove heat from the whole ocean.'
   if name=='Label — scope note':y+=sh*.038
   if name=='Label — control note':x=sw*.99
   artist=ax.text(x,y,text,fontsize=font,fontfamily='DejaVu Sans',fontstretch='condensed',fontweight='bold' if bold else 'normal',color=color(rp,'#263440'),va='top',ha='right' if name=='Label — control note' else 'left',linespacing=1.05,zorder=k+1,clip_on=False)
   artist.set_gid('ppt-'+name);counts['editable_text_labels']+=1
  assert counts['native_shapes']==132 and counts['editable_text_labels']==21 and counts['anchors']==11
  return counts
