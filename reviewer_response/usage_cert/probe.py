# timing / hit-rate probe for choosing n_iter (not part of the reported results)
import sys, time, numpy as np
sys.argv=['x']
import usage_cert as U
from multiprocessing import get_context
vecs=U.loadClassVectors(); chem=U.distMatrix(U.CLASS_ORDER,vecs)
c64,nbr64=U.buildSenseGraph(3); fU=U.uniformUsage(10,list(range(10)))
comp,_=U.loadProteome(); fP=U.classUsage(comp,set(U.CLASS_ORDER))
lab=U.randomStart(64,10,np.random.default_rng(0))
t=time.time(); U.anneal(lab.copy(),nbr64,chem,fU,10,1000,1.0,0); print('jit',time.time()-t)
def one(args):
    mult,seed,f=args
    lab=U.randomStart(64,10,np.random.default_rng(seed)); t=time.time()
    U.anneal(lab,nbr64,chem,f,10,mult*64*64,1.0,seed)
    return U.fullDw(lab,nbr64,chem,f,10), time.time()-t
for name,f in (('U',fU),('P',fP)):
  for mult in (1000,4000):
    with get_context('fork').Pool(60) as p: out=p.map(one,[(mult,s,f) for s in range(120)],chunksize=1)
    v=np.array([o[0] for o in out]); tt=[o[1] for o in out]
    print(name,mult,'best',v.min(),'hits',np.sum(v-v.min()<=1e-9),'/120 p50',np.median(v),'sec/restart',np.mean(tt),flush=True)
