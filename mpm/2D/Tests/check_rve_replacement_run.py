"""Short end-to-end replacement, threaded parity, conf round-trip and invalid-model check."""
from pathlib import Path
import argparse,ast,json,subprocess,shutil
import numpy as np
root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--simulation',type=Path,required=True,help='75-MP case directory containing input_RVE_replacement.txt and VER1000.txt (external simulation data)')
parser.add_argument('--executable',type=Path,default=root/'BUILD_RVE_REPLACEMENT/mpmbox')
args=parser.parse_args()
exe=args.executable.resolve();simulation=args.simulation.resolve()
work=root/'Tests/_work/rve_replacement';work.mkdir(parents=True,exist_ok=True)
base=(simulation/'input_RVE_replacement.txt').read_text()
def snapshot(f):
 lines=f.read_text().splitlines()
 t=float(next(l.split()[1] for l in lines if l.startswith('t ')))
 k=next(i for i,l in enumerate(lines) if l.startswith('MPs '));n=int(lines[k].split()[1])
 rows=np.array([[float(x) for x in line.split()[1:]] for line in lines[k+1:k+1+n]])
 return t,rows[np.argsort(rows[:,0])]

outputs=[]
for threads in [1,4]:
 folder=work/f'threads{threads}';folder.mkdir(exist_ok=True);shutil.copy2(simulation/'VER1000.txt',folder/'VER1000.txt')
 lines=[]
 for line in base.splitlines():
  if line.startswith('finalTime '):line='finalTime 0.000025'
  if line.startswith('confPeriod '):line='confPeriod 1'
  if line.startswith('dt '):line='dt 1e-5'
  if line.startswith('CHCL.replaceRVE '):line='CHCL.replaceRVE 10 MC_RVE' # deliberately force replacement for a short test
  lines.append(line)
 (folder/'input.txt').write_text('\n'.join(lines)+'\n')
 with (folder/'run.log').open('w') as log:subprocess.run([str(exe),'input.txt','-j',str(threads)],cwd=folder,stdout=log,stderr=log,check=True,timeout=90)
 messages=[l for l in (folder/'run.log').read_text().splitlines() if '[RVE replacement] DEM -> Mohr-Coulomb:' in l]
 assert len(messages)==75, ('missing/duplicate transition notifications',len(messages))
 assert all('MP=' in l and 't=' in l and 'threshold=10' in l for l in messages)
 last=max(folder.glob('conf[0-9]*.txt'),key=lambda p:int(p.stem[4:]));text=last.read_text();k=next(i for i,l in enumerate(text.splitlines()) if l.startswith('MPs '));rows=text.splitlines()[k+1:k+76];assert len(rows)==75 and all(r.split()[0]=='MC_RVE' for r in rows)
 assert 'CHCL.replaceRVE 10 MC_RVE' in text
 t,data=snapshot(last);outputs.append(data)
 # Continue after the converted model has been serialized and read back.
 restart=text.replace('finalTime 2.5e-05','finalTime 0.000045')
 (folder/'restart.txt').write_text(restart)
 with (folder/'restart.log').open('w') as log:subprocess.run([str(exe),'restart.txt','-j',str(threads)],cwd=folder,stdout=log,stderr=log,check=True,timeout=90)
np.testing.assert_allclose(outputs[0],outputs[1],rtol=1e-10,atol=1e-10)
folder=work/'invalid';folder.mkdir(exist_ok=True);shutil.copy2(simulation/'VER1000.txt',folder/'VER1000.txt')
text=base.replace('CHCL.replaceRVE 1.5 MC_RVE','CHCL.replaceRVE 1.5 UNKNOWN');(folder/'input.txt').write_text(text)
with (folder/'run.log').open('w') as log:r=subprocess.run([str(exe),'input.txt'],cwd=folder,stdout=log,stderr=log,timeout=90)
assert r.returncode>0 and 'declared valid MohrCoulomb' in (folder/'run.log').read_text()
print('PASS: 75 DEM->MC transitions; 1/4-thread parity; conf serialization and restart; invalid target rejected')
