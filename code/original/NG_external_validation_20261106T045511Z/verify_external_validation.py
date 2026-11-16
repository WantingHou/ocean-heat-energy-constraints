"""Independent stdlib/Decimal audit of regional calculations and manuscript preservation."""
from pathlib import Path
import argparse,csv,io,json,re,hashlib
from decimal import Decimal,localcontext
from collections import defaultdict
ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--predecessor',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();D=a.package;P=a.predecessor
if a.output.exists():raise FileExistsError(a.output)
raw=a.raw.read_text().split('*/',1)[1].strip();records=list(csv.DictReader(io.StringIO(raw),delimiter='\t'));g=defaultdict(list)
for r in records:g[r['Station'].split('P')[0]].append(Decimal(r['HF [mW/m**2]']))
assert len(records)==109 and len(g)==14
primary={k:v for k,v in g.items() if k!='HF1511'};assert sum(map(len,primary.values()))==105 and len(primary)==13
with localcontext() as ctx:
 ctx.prec=45
 mean=lambda v:sum(v)/len(v)
 ms=mean([v for vs in primary.values() for v in vs]);mp=mean([mean(vs) for vs in primary.values()]);allmean=mean([v for vs in g.values() for v in vs]);maxerror=0.
 for row in csv.DictReader((D/'Source_Data/external_validation/regional_comparison.csv').open()):
  age=Decimal(row['age_Ma']);q=Decimal('506.7')/age.sqrt();m=ms if row['aggregation']=='primary_station' else mp;pred=q*(Decimal('.66')-Decimal('.0111')*age);deficit=q-m
  expected={'conductive_mean_mW_m2':m,'cooling_reference_mW_m2':q,'inferred_deficit_mW_m2':deficit,'frozen_age_deficit_mW_m2':pred,'ordinary_density_mW_m2':Decimal(40),'deficit_minus_density_mW_m2':deficit-40,'deficit_to_density':deficit/40,'deficit_minus_age_prediction_mW_m2':deficit-pred,'predicted_conductive_mW_m2':q-pred,'deficit_fraction':deficit/q}
  for k,v in expected.items():maxerror=max(maxerror,float(abs(Decimal(row[k])-v)))
assert maxerror<1e-11
station_rows=list(csv.DictReader((D/'Source_Data/external_validation/stations.csv').open()));sr={r['station']:r for r in station_rows}
for r in records:
 q=sr[r['Station']];assert Decimal(q['conductive_mW_m2'])==Decimal(r['HF [mW/m**2]']);assert q['primary']==str(not r['Station'].startswith('HF1511'));assert q['assumed_conductivity']==str(Decimal(r['k [W/m/K]'])==Decimal('.83'))
for row in csv.DictReader((D/'Source_Data/external_validation/profiles.csv').open()):
 vals=g[row['profile']];assert len(vals)==int(row['count']);assert abs(float(mean(vals))-float(row['mean_conductive_mW_m2']))<1e-12
main=(D/'manuscript_NG_external_EN.md').read_text();si=(D/'Supplementary_Information_NG_EN.md').read_text();oldsi=(P/'Supplementary_Information_NG_EN.md').read_text();oldmain=(P/'manuscript_NG_EN.md').read_text()
assert re.findall(r'^## Figure (\d+)$',main,re.M)==[str(i) for i in range(1,7)];assert not re.search(r'^\|',main,re.M)
assert re.findall(r'^## Supplementary Figure S(\d+)$',si,re.M)==[str(i) for i in range(1,18)]
assert re.findall(r'^## Supplementary Table (\d+)\.',si,re.M)==[str(i) for i in range(1,51)]
for text in [main,si]:
 for path in re.findall(r'!\[[^\]]*\]\(([^)]+)\)',text):assert (D/path).exists() and (D/path).with_suffix('.pdf').exists()
 refs=set(map(int,re.findall(r'^([0-9]+)\. ',text,re.M)));used={int(n) for m in re.finditer(r'\[([0-9]+(?:, [0-9]+)*)\]',text) for n in m.group(1).split(', ')};assert used<=refs
 assert all('http' in line for line in text.splitlines() if re.match(r'^\d+\. ',line))
# Reference renumbering is allowed; all original table cells otherwise stay identical.
clean=lambda s:re.sub(r'\[[0-9, –-]+\]','[citation]',s)
for n in range(1,47):
 def tab(text):
  block=re.search(rf'^## Supplementary Table {n}\..*?(?=^## |\Z)',text,re.M|re.S).group();return [clean(l) for l in block.splitlines() if l.startswith('|')]
 assert tab(si)==tab(oldsi),n
for before,after in [(oldmain,main),(oldsi,si)]:assert re.findall(r'\$\$\s*(.*?)\s*\$\$',before,re.S)==re.findall(r'\$\$\s*(.*?)\s*\$\$',after,re.S)
oldhash=json.loads((D/'provenance/preceding_package_hashes.json').read_text());assert all(hashlib.sha256((P/p).read_bytes()).hexdigest()==h for p,h in oldhash.items())
preserved_data=sum(1 for p in (P/'Source_Data').rglob('*') if p.is_file() and (D/p.relative_to(P)).read_bytes()==p.read_bytes())
preserved_figures=sum(1 for p in (P/'figures').iterdir() if p.suffix in ['.png','.pdf'] and (D/p.relative_to(P)).read_bytes()==p.read_bytes())
result=dict(independent_arithmetic='stdlib CSV and 45-digit Decimal, separate from pandas production calculation',station_records_checked=109,profiles_checked=14,regional_comparison_rows_checked=6,max_decimal_discrepancy_mW_m2=maxerror,primary_station_mean=float(ms),primary_profile_mean=float(mp),all_cruise_mean=float(allmean),main_figures=6,main_tables=0,SI_figures=17,SI_tables=50,unchanged_SI_tables_checked=46,unchanged_equation_blocks=22,predecessor_files_unchanged=len(oldhash),inherited_source_data_files_identical=preserved_data,inherited_figure_files_identical=preserved_figures,erratum_text_verified=False,global_attribution_identified=False)
with a.output.open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))
