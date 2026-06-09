import numpy as np, csv, os, time, math
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
from math import isqrt
t0=time.time()
N=1_000_000_000
s=np.ones(N+1,dtype=bool); s[:2]=False
for p in range(2,isqrt(N)+1):
    if s[p]: s[p*p::p]=False
pr=np.nonzero(s)[0]; del s
pr=pr[pr>3]
w=(pr%6==1).astype(np.int8)
sig=(w[1:]^w[:-1]).astype(np.int8)               # 1=switch
left=pr[:-1]                                      # int64 view, prime at left of each step
m=len(sig)
print("primes>3=%d steps=%d (%.1fs)"%(len(pr),m,time.time()-t0))

# ---------- (1) step autocorrelation + AR(1) / PACF(2) ----------
ps=float(sig.mean()); sf=sig.astype(np.float64)-ps; den=float(np.sum(sf*sf))
rho=[float(np.sum(sf[:-k]*sf[k:])/den) for k in range(1,9)]
r1=rho[0]; r2=rho[1]
ar1_r2=r1*r1; pacf2=(r2-r1*r1)/(1-r1*r1)
with open(os.path.join(_DATA,"fine_autocorr.csv"),"w",newline="") as fh:
    wc=csv.writer(fh); wc.writerow(["lag","rho","sigma","AR1_prediction_rho1_pow_lag"])
    for k in range(1,9):
        wc.writerow([k,"%.6f"%rho[k-1],"%.1f"%(abs(rho[k-1])*(m-k)**0.5),"%.6f"%(r1**k)])
with open(os.path.join(_DATA,"fine_lag2_test.csv"),"w",newline="") as fh:
    wc=csv.writer(fh); wc.writerow(["quantity","value"])
    for k,v in [("P_switch","%.5f"%ps),("rho1","%.6f"%r1),("rho2_measured","%.6f"%r2),
        ("rho2_AR1_pred_r1sq","%.6f"%ar1_r2),("rho2_excess","%.6f"%(r2-ar1_r2)),
        ("rho2_ratio_to_AR1","%.2f"%(r2/ar1_r2)),("PACF2","%.6f"%pacf2),
        ("PACF2_sigma","%.0f"%(abs(pacf2)*m**0.5)),("rho3","%.6f"%rho[2])]:
        wc.writerow([k,v])
print("rho1=%.5f rho2=%.5f (AR1 %.5f, x%.2f) PACF2=%.5f (%.0f sig)"%(r1,r2,ar1_r2,r2/ar1_r2,pacf2,abs(pacf2)*m**0.5))

# ---------- (2) decay law, fine log windows ----------
rows=[]; lo=4.0
while lo<9.01:
    hi=lo+0.25; L=int(10**lo); H=int(10**hi)
    i0=int(np.searchsorted(left,L)); i1=int(np.searchsorted(left,H)); n=i1-i0
    if n>20000:
        sw=sig[i0:i1]; psame=1.0-float(sw.mean()); d1=0.5-psame
        x=sw[:-1].astype(np.float64); y=sw[1:].astype(np.float64); xm=x-x.mean()
        r1w=float(np.sum(xm*(y-y.mean()))/np.sum(xm*xm))
        lnN=math.log(math.sqrt(L*H))
        rows.append((lnN,psame,d1,abs(r1w),n))
    lo+=0.25
LN=np.array([r[0] for r in rows]); D1=np.array([r[2] for r in rows]); R1=np.array([r[3] for r in rows])
def fitpow(y):
    Y=np.log(y); X=np.log(LN); A=np.vstack([X,np.ones_like(X)]).T
    c,_,_,_=np.linalg.lstsq(A,Y,rcond=None); return -c[0],math.exp(c[1])
a1,c1=fitpow(D1); a2,c2=fitpow(R1)
with open(os.path.join(_DATA,"fine_decay.csv"),"w",newline="") as fh:
    wc=csv.writer(fh); wc.writerow(["lnN","P_same","delta1","abs_rho1","n","delta1_times_lnN","ratio_d1_over_rho1"])
    for (lnN,psame,d1,ar,n) in rows:
        wc.writerow(["%.4f"%lnN,"%.5f"%psame,"%.6f"%d1,"%.6f"%ar,n,"%.4f"%(d1*lnN),"%.3f"%(d1/ar)])
with open(os.path.join(_DATA,"fine_decay_fit.csv"),"w",newline="") as fh:
    wc=csv.writer(fh); wc.writerow(["signal","power_alpha","C","note"])
    wc.writerow(["delta1_first_order","%.3f"%a1,"%.3f"%c1,"delta1 ~ C/(logN)^alpha"])
    wc.writerow(["absrho1_second_order","%.3f"%a2,"%.3f"%c2,"|rho1| ~ C/(logN)^alpha"])
    wc.writerow(["exponent_ratio_2nd_over_1st","%.3f"%(a2/a1),"-","1.0 => same law"])
    wc.writerow(["delta1_times_logN_mean","%.3f"%float((D1*LN).mean()),"-","~const => leading 1/logN"])
    wc.writerow(["ratio_delta1_over_rho1_mean","%.3f"%float((D1/R1).mean()),"-","~const => locked"])
print("decay: alpha(d1)=%.3f alpha(rho1)=%.3f  d1*logN_mean=%.3f  d1/rho1_mean=%.3f"%(
    a1,a2,float((D1*LN).mean()),float((D1/R1).mean())))
print("TOTAL %.1fs"%(time.time()-t0))
