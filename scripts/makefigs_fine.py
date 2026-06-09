import csv, numpy as np, matplotlib, os, math
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":0.3,"figure.dpi":150,"savefig.bbox":"tight"})

def rd(n): return list(csv.DictReader(open(os.path.join(_DATA,n))))

# ============ FIGURE 1 : lag-2 echo vs AR(1) ============
fig,ax=plt.subplots(1,2,figsize=(9.4,4.0))
ac=rd("fine_autocorr.csv")
lag=np.array([int(r["lag"]) for r in ac]); rho=np.array([float(r["rho"]) for r in ac])
ar1=np.array([float(r["AR1_prediction_rho1_pow_lag"]) for r in ac])
W=0.38
ax[0].axhline(0,color="0.5",lw=0.8)
ax[0].bar(lag-W/2,rho,W,color="#1f77b4",label="measured $\\rho_k$")
ax[0].bar(lag+W/2,ar1,W,color="#d62728",alpha=0.8,label="AR(1) prediction $\\rho_1^{\\,k}$")
ax[0].set_xlabel("lag $k$ (steps)"); ax[0].set_ylabel("switch-sequence autocorrelation")
ax[0].set_title("(A) lag-1 dominates; lag-2 exceeds AR(1)",fontsize=9)
ax[0].set_xticks(lag); ax[0].legend(fontsize=7.5,loc="lower right")
# inset-like zoom on lag2-4
axin=ax[0].inset_axes([0.45,0.12,0.5,0.42])
sel=(lag>=2)&(lag<=4)
axin.axhline(0,color="0.5",lw=0.6)
axin.bar(lag[sel]-W/2,rho[sel],W,color="#1f77b4")
axin.bar(lag[sel]+W/2,ar1[sel],W,color="#d62728",alpha=0.8)
axin.set_xticks([2,3,4]); axin.set_title("zoom lag 2-4",fontsize=6.5); axin.tick_params(labelsize=6)
# (B) PACF2 bar summary
t2={r["quantity"]:r["value"] for r in rd("fine_lag2_test.csv")}
cats=["$\\rho_1^2$\n(AR1 derived)","measured\n$\\rho_2$","excess\n=PACF(2)"]
vals=[float(t2["rho2_AR1_pred_r1sq"]),float(t2["rho2_measured"]),float(t2["PACF2"])]
cols=["#d62728","#1f77b4","#2ca02c"]
ax[1].bar(range(3),vals,color=cols,width=0.6)
for i,v in enumerate(vals): ax[1].annotate("%.5f"%v,(i,v),fontsize=8,ha="center",va="bottom")
ax[1].set_xticks(range(3)); ax[1].set_xticklabels(cats,fontsize=8)
ax[1].set_ylabel("lag-2 autocorrelation")
ax[1].set_title("(B) lag-2 is $4.3\\times$ the AR(1) value (23$\\sigma$)",fontsize=9)
ax[1].set_ylim(0,0.0052)
fig.suptitle(r"Fine structure of the wing-transition memory on the $6N$ skeleton",fontsize=10)
fig.savefig(os.path.join(_FIG,"p35_fig1.pdf")); print("fig1 done")

# ============ FIGURE 2 : decay law ============
fig2,ax2=plt.subplots(1,2,figsize=(9.4,4.0))
dec=rd("fine_decay.csv")
lnN=np.array([float(r["lnN"]) for r in dec]); d1=np.array([float(r["delta1"]) for r in dec])
ar=np.array([float(r["abs_rho1"]) for r in dec]); d1ln=np.array([float(r["delta1_times_lnN"]) for r in dec])
ratio=np.array([float(r["ratio_d1_over_rho1"]) for r in dec])
# (A) log-log decay with fit
fit=rd("fine_decay_fit.csv"); a1=float(fit[0]["power_alpha"]); c1=float(fit[0]["C"])
ax2[0].loglog(lnN,d1,"o",ms=5,color="#1f77b4",label=r"$\delta_1=\frac{1}{2}-P(\mathrm{same})$")
ax2[0].loglog(lnN,ar,"s",ms=5,color="#2ca02c",label=r"$|\rho_1|$ (2nd order)")
xx=np.linspace(lnN.min(),lnN.max(),50)
ax2[0].loglog(xx,c1/xx**a1,"--",color="#1f77b4",lw=1,label=r"$\sim(\log N)^{-%.2f}$"%a1)
ax2[0].set_xlabel(r"$\log N$"); ax2[0].set_ylabel("signal strength")
ax2[0].set_title("(A) the bias fades as a power of $\\log N$",fontsize=9)
from matplotlib.ticker import ScalarFormatter, NullFormatter
ax2[0].set_xticks([14,16,18,20]); ax2[0].xaxis.set_major_formatter(ScalarFormatter())
ax2[0].xaxis.set_minor_formatter(NullFormatter())
ax2[0].legend(fontsize=7.5,loc="lower left")
# (B) delta1*logN ~ const  and ratio locked
ax2b=ax2[1]
ax2b.plot(lnN,d1ln,"o-",ms=5,color="#1f77b4",label=r"$\delta_1\cdot\log N$ (leading $1/\log N$)")
ax2b.axhline(d1ln.mean(),color="#1f77b4",lw=0.7,ls=":")
ax2b.plot(lnN,ratio,"D-",ms=4,color="#d62728",label=r"$\delta_1/|\rho_1|$ (1st/2nd locked)")
ax2b.axhline(ratio.mean(),color="#d62728",lw=0.7,ls=":")
ax2b.set_xlabel(r"$\log N$"); ax2b.set_ylabel("near-constant ratios")
ax2b.set_title("(B) $\\delta_1\\log N\\approx1$; 1st/2nd locked $\\approx1.7$",fontsize=9)
ax2b.legend(fontsize=7.5,loc="center right"); ax2b.set_ylim(0,2.0)
fig2.savefig(os.path.join(_FIG,"p35_fig2.pdf")); print("fig2 done")
