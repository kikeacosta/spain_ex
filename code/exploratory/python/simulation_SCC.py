"""Simulation of the cohort-indexed decomposition (Step 3 of the manuscript), as described in the text.
Gap defined as e(REF) - e(ESP): negative change = narrowing (Fig. 5 convention).
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# Times New Roman (falls back to a metric-compatible serif if not installed)
plt.rcParams.update({"font.family": "serif",
                     "font.serif": ["Times New Roman", "Liberation Serif", "Times", "DejaVu Serif"],
                     "mathtext.fontset": "stix",
                     "pdf.fonttype": 42, "ps.fonttype": 42})  # embed TrueType fonts in PDF

ages = np.arange(0, 111); years = np.arange(1900, 2020)
A, Y = np.meshgrid(ages, years, indexing="ij")
C = Y - A

def base(x):  # Siler-type baseline hazard (REF-like, c.1900)
    return 0.12*np.exp(-1.2*x) + 0.004 + 0.00008*np.exp(0.09*x)

m_ref = base(A)*np.exp(-0.02*(Y-1900))
# Spain: excess log-mortality, flat until 1940, rapid (nonlinear) catch-up 1940-60, gone by ~1980
G = 0.9 - 0.85/(1+np.exp(-(Y-1950)/4.0)) ; G = np.where(Y > 1980, G*np.exp(-(Y-1980)/5), G)
lin = np.clip(0.9*(1-(Y-1940)/60), 0, 0.9)  # homogeneous (linear) convergence alternative

scen = {
 "A. Period only (nonlinear catch-up 1940-60)": m_ref*np.exp(G),
 "B. A + cohorts 1945-54 with -15% mortality (Spain only)": m_ref*np.exp(G)*np.where((C>=1945)&(C<=1954),0.85,1.0),
 "C. Period only, linear convergence": m_ref*np.exp(lin),
}

def temp_lt(m, a, n):
    q = 1-np.exp(-m[a:a+n]); l = np.concatenate([[1.0], np.cumprod(1-q)])
    L = (l[:-1]+l[1:])/2; T = np.cumsum(L[::-1])[::-1]; T = np.append(T, 0.0)
    return l, L, T

def arriaga(m1, m2, a, n):
    """Contribution of each age in [a,a+n) to e2 - e1 (temporary LE, conditional on survival to a)."""
    l1, L1, T1 = temp_lt(m1, a, n); l2, L2, T2 = temp_lt(m2, a, n)
    d = l1[:-1]*(L2/l2[:-1]-L1/l1[:-1]) + T2[1:]*(l1[:-1]/l2[:-1] - l1[1:]/l2[1:])
    return d

def cohort_contrib(mE, mR, a, n=10, y0=1940, y1=2019):
    yi = {y:i for i,y in enumerate(years)}
    D = {y: arriaga(mE[:, yi[y]], mR[:, yi[y]], a, n) for y in range(y0, y1+1)}  # e_REF - e_ESP
    out = {}
    for c in range(y0-a-n+1, y1-a):
        s, ok = 0.0, True
        for t in range(n):
            y = c+a+t
            if y < y0 or y+1 > y1: ok=False; break
            s += D[y+1][t]-D[y][t]
        if ok: out[c]=s
    cs = np.array(sorted(out)); v = np.array([out[c] for c in cs])
    return cs, v

# sanity check: Arriaga sums to difference
l1,L1,T1 = temp_lt(scen["A. Period only (nonlinear catch-up 1940-60)"][:,45],0,10); l2,L2,T2=temp_lt(m_ref[:,45],0,10)
print("check", arriaga(scen["A. Period only (nonlinear catch-up 1940-60)"][:,45], m_ref[:,45],0,10).sum(), T2[0]-T1[0])

# ---- Figure: all ten 10-year age intervals, laid out as Figure 5 -------------
# rows = age intervals (90-99 at the top, 0-9 at the bottom, as in Figure 5)
# columns = scenarios A, B, C
intervals = list(range(90, -1, -10))
col_titles = ["A. Period only\n(nonlinear catch-up 1940-60)",
              "B. A + cohort advantage\n(cohorts 1945-54, -15%)",
              "C. Period only\n(linear convergence 1940-2000)"]
fig, axes = plt.subplots(len(intervals), len(scen), figsize=(7.5, 10.5),
                         sharex=True, sharey=True)
for j, (name, mE) in enumerate(scen.items()):
    for i, a in enumerate(intervals):
        ax = axes[i, j]
        # cohorts that crossed this age interval during 1940-1960 (dark bars in Fig. 5)
        ax.axvspan(1940 - a - 9, 1960 - a, color="grey", alpha=.22, lw=0)
        if name.startswith("B."):
            ax.axvspan(1945, 1955, color="gold", alpha=.40, lw=0)
        cs, v = cohort_contrib(mE, m_ref, a)
        if len(cs):
            scc = v / np.abs(v).max()
            ax.bar(cs, scc, width=1.0, color=np.where(scc < 0, "#2a6f97", "#c1121f"))
        ax.axhline(0, color="k", lw=.4)
        ax.set_ylim(-1.1, 1.1); ax.set_xlim(1850, 2010)
        ax.set_yticks([-1, 0, 1]); ax.tick_params(labelsize=7, length=2)
        ax.set_xticks(range(1860, 2011, 30))
        if j == 0:
            ax.set_ylabel(f"{a}-{a+9}", fontsize=8, rotation=0, ha="right", va="center")
        if i == 0:
            ax.set_title(col_titles[j], fontsize=8.5)
fig.supxlabel("Birth cohort", fontsize=9, y=0.056)
fig.supylabel("Scaled cohort contribution (SCC) by age interval", fontsize=9, x=0.01)
fig.text(0.5, 0.024,
         "Negative values (blue) = narrowing of the REF-ESP gap; positive (red) = widening. "
         "Grey: cohorts that crossed each age interval during 1940-60.\n"
         "Yellow: cohorts 1945-54 with 15% lower mortality at all ages (scenario B); with data ending in 2019, "
         "these cohorts are not observed over a complete interval beyond ages 60-69.",
         ha="center", fontsize=6.8)
plt.tight_layout(rect=[0.02, 0.035, 1, 1], h_pad=0.3, w_pad=0.6)
out = "simulation_SCC_EN.png"
plt.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.08)
plt.savefig(out.replace(".png", ".pdf"), bbox_inches="tight", pad_inches=0.08)

# numeric summary for scenario B, ages 30-39
cs, vA = cohort_contrib(scen[list(scen)[0]], m_ref, 30); _, vB = cohort_contrib(scen[list(scen)[1]], m_ref, 30)
for c,a_,b_ in zip(cs, vA, vB):
    if 1940<=c<=1958: print(c, round(a_,5), round(b_,5), round(b_-a_,5))
