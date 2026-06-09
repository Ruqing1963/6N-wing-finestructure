import numpy as np, time, math
from math import isqrt
t0=time.time()
N=1_000_000_000
s=np.ones(N+1,dtype=bool); s[:2]=False
for p in range(2,isqrt(N)+1):
    if s[p]: s[p*p::p]=False
pr=np.nonzero(s)[0]
del s                                        # free 1GB immediately
pr=pr[pr>3]
w=(pr%6==1).astype(np.int8)
print("primes>3=%d (%.1fs)"%(len(pr),time.time()-t0))
sig=(w[1:]^w[:-1]).astype(np.int8)           # 1=switch  (xor of adjacent wings)
left=pr[:-1]                                  # prime at left of each step (int64 view)
same=(1-sig).astype(np.int8)                  # 1=same wing
m=len(sig)

centers=[]; D1=[]; R1=[]; LN=[]
lo=4.0
print("\n logN_c  P(same)  delta1   rho1(win)     n")
while lo<9.01:
    hi=lo+0.25; L=int(10**lo); H=int(10**hi)
    i0=np.searchsorted(left,L); i1=np.searchsorted(left,H)
    n=i1-i0
    if n>20000:
        sw=sig[i0:i1]
        psame=1.0-sw.mean(); d1=0.5-psame
        # lag-1 autocorr of sw on this window
        x=sw[:-1].astype(np.float64); y=sw[1:].astype(np.float64)
        xm=x-x.mean()
        r1=float(np.sum(xm*(y-y.mean()))/np.sum(xm*xm))
        lnN=math.log(math.sqrt(L*H))
        centers.append((lo+hi)/2); D1.append(d1); R1.append(abs(r1)); LN.append(lnN)
        print("  %5.2f  %.5f  %+.5f  %+.5f   %d"%((lo+hi)/2,psame,d1,r1,n))
    lo+=0.25
D1=np.array(D1); R1=np.array(R1); LN=np.array(LN)

def fitpow(y):
    Y=np.log(y); X=np.log(LN); A=np.vstack([X,np.ones_like(X)]).T
    c,_,_,_=np.linalg.lstsq(A,Y,rcond=None); return -c[0],math.exp(c[1])
a1,c1=fitpow(D1); a2,c2=fitpow(R1)
print("\nfit delta1 ~ C/(logN)^alpha :  alpha=%.3f  C=%.3f"%(a1,c1))
print("fit |rho1| ~ C/(logN)^alpha :  alpha=%.3f  C=%.3f"%(a2,c2))
print("ratio of exponents (2nd/1st) = %.3f   [1.0 => same decay law]"%(a2/a1))
# O-S leading-form check: delta1 ~ k*loglogN/logN ?
print("\nO-S form  delta1*logN  and  delta1*logN/loglogN :")
for i in range(0,len(LN),2):
    print("  logN=%.2f  d1*logN=%.4f  d1*logN/loglogN=%.4f  rho1*logN=%.4f"%(
        LN[i],D1[i]*LN[i],D1[i]*LN[i]/math.log(LN[i]),R1[i]*LN[i]))
# also the ratio delta1/|rho1| across scale (should be ~const if same law)
print("\nratio delta1/|rho1| across scale (const => locked together):")
print("  ",np.round(D1/R1,3))
print("\nsave for figure")
np.savetxt("/tmp/decay.csv",np.column_stack([LN,D1,R1]),delimiter=",",header="lnN,delta1,absrho1",comments="")
print("TOTAL %.1fs"%(time.time()-t0))
