"""Prospective zero-exceedance design, not an observed confidence result."""
from pathlib import Path
import argparse,csv,math
ap=argparse.ArgumentParser();ap.add_argument('--thresholds',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();rows=[]
with a.thresholds.open() as f:
 for x in csv.DictReader(f):
  if float(x['target_percent'])!=5:continue
  threshold=float(x['critical_kernel_weighted_excess_fraction']);n=math.ceil(math.log(.05)/math.log1p(-threshold));assert (1-threshold)**n<=.05 and (1-threshold)**(n-1)>.05
  rows.append({'allocation':x['allocation'],'loss':x['loss'],'hypothetical_target_fraction':threshold,'independent_zero_exceedance_elements_needed':n,'one_sided_design_alpha':.05,'observed_admissible_sample_count':None,'scope':'prospective only: independent ordinary-power × positive-response-kernel weighted historical elements, known baseline and correct classification; NOT current vent observations or an achieved confidence bound'})
with a.output.open('x',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
