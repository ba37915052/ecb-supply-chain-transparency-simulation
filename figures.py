import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simulation import build_graph, run_counterfeit, run_cascade

# regenerate per-seed arrays quickly (reuse seeds from main run for CPR only, 30 seeds)
cpr0, cpr1 = [], []
for s in range(30):
    G, nodes, oem = build_graph(1000+s)
    c0,_,_ = run_counterfeit(G, nodes, oem, 20000, 2000+s, False)
    c1,_,_ = run_counterfeit(G, nodes, oem, 20000, 2000+s, True)
    cpr0.append(c0); cpr1.append(c1)
rel = [(a-b)/a*100 for a,b in zip(cpr0,cpr1)]

fig, ax = plt.subplots(1,2, figsize=(10,4))
bp = ax[0].boxplot([np.array(cpr0)*1000, np.array(cpr1)*1000], tick_labels=["Baseline\n(As-Is)","Platform\n(To-Be)"], patch_artist=True)
for p,c in zip(bp['boxes'], ['#c44','#284']): p.set_facecolor(c); p.set_alpha(.6)
ax[0].set_ylabel("Counterfeit penetration rate, ×10$^{-3}$")
ax[0].set_title("(a) CPR across 30 independent runs")
ax[1].hist(rel, bins=10, color='#468', alpha=.8, edgecolor='k')
ax[1].axvline(np.mean(rel), color='k', ls='--', label=f"mean {np.mean(rel):.1f}%")
ax[1].set_xlabel("Relative CPR reduction, %"); ax[1].set_ylabel("Frequency")
ax[1].set_title("(b) Distribution of relative reduction"); ax[1].legend()
plt.tight_layout(); plt.savefig("fig_cpr.png", dpi=200); plt.close()

# functionality curves (mean over 10 seeds) - patch run_cascade to return curve
import networkx as nx
from simulation import TIER_SIZES
def cascade_curve(G, nodes, oem, seed, platform, theta=.30, shock=.10, rho=.10, t_max=60):
    rng = np.random.default_rng(seed)
    U = G.to_undirected(); nbrs = {v:list(U.neighbors(v)) for v in U.nodes}
    failed = {v:False for v in U.nodes}
    pool = nodes[2]+nodes[3]; k=int(round(shock*len(pool)))
    for v in rng.choice(pool, size=k, replace=False): failed[int(v)]=True
    p_rec = .20 if platform else .05; er = rho if platform else 0.0
    def conn():
        H=G.subgraph([v for v in G.nodes if not failed[v]])
        if oem not in H: return 0.0
        return sum(1 for v in nodes[1] if v in H and nx.has_path(H,v,oem))/len(nodes[1])
    curve=[]
    for t in range(t_max):
        newly=[v for v in U.nodes if not failed[v] and nbrs[v] and
               max(0.0, sum(failed[u] for u in nbrs[v])/len(nbrs[v])-er)>theta]
        for v in newly: failed[v]=True
        for v in U.nodes:
            if failed[v] and rng.random()<p_rec: failed[v]=False
        curve.append(conn())
    return curve

c0s, c1s = [], []
for s in range(10):
    G, nodes, oem = build_graph(1000+s)
    c0s.append(cascade_curve(G,nodes,oem,3000+s,False))
    c1s.append(cascade_curve(G,nodes,oem,3000+s,True))
c0m, c1m = np.mean(c0s,0), np.mean(c1s,0)
c0sd, c1sd = np.std(c0s,0), np.std(c1s,0)
t = np.arange(60)
plt.figure(figsize=(7,4))
plt.plot(t,c0m,color='#c44',label='Baseline (As-Is)')
plt.fill_between(t,c0m-c0sd,np.minimum(1,c0m+c0sd),color='#c44',alpha=.2)
plt.plot(t,c1m,color='#284',label='Platform (To-Be)')
plt.fill_between(t,np.maximum(0,c1m-c1sd),np.minimum(1,c1m+c1sd),color='#284',alpha=.2)
plt.xlabel("Simulation step"); plt.ylabel("Tier-1 → OEM connectivity")
plt.title("Functionality curves after Tier-2/3 shock (mean ± SD, 10 runs)")
plt.legend(); plt.tight_layout(); plt.savefig("fig_resilience.png", dpi=200); plt.close()

# sensitivity plot
sens = json.load(open("results_sensitivity.json"))
fig, ax = plt.subplots(1,2, figsize=(10,4))
ds=[0.5,1.0,1.5]; y=[sens[f"delta_scale_{d}"]["mean"]*100 for d in ds]
lo=[sens[f"delta_scale_{d}"]["ci95"][0]*100 for d in ds]; hi=[sens[f"delta_scale_{d}"]["ci95"][1]*100 for d in ds]
ax[0].errorbar(ds,y,yerr=[np.array(y)-lo,np.array(hi)-y],marker='o',color='#284',capsize=4,label='Δd scale')
ps=[0.5,1.0,2.0]; y2=[sens[f"p_scale_{p}"]["mean"]*100 for p in ps]
lo2=[sens[f"p_scale_{p}"]["ci95"][0]*100 for p in ps]; hi2=[sens[f"p_scale_{p}"]["ci95"][1]*100 for p in ps]
ax[0].errorbar(ps,y2,yerr=[np.array(y2)-lo2,np.array(hi2)-y2],marker='s',color='#468',capsize=4,label='p scale')
ax[0].axhspan(40,60,color='gray',alpha=.15,label='hypothesized 40–60% band')
ax[0].set_xlabel("Parameter multiplier"); ax[0].set_ylabel("Relative CPR reduction, %")
ax[0].set_title("(a) Sensitivity of counterfeit-risk reduction"); ax[0].legend(fontsize=8)
th=[0.2,0.3,0.4]; r0=[sens[f"theta_{x}"]["ri0"] for x in th]; r1=[sens[f"theta_{x}"]["ri1"] for x in th]
sf=[0.05,0.1,0.2]; s0=[sens[f"shock_{x}"]["ri0"] for x in sf]; s1=[sens[f"shock_{x}"]["ri1"] for x in sf]
ax[1].plot(th,r0,'o--',color='#c44',label='RI baseline vs θ'); ax[1].plot(th,r1,'o-',color='#284',label='RI platform vs θ')
ax[1].plot(sf,s0,'s--',color='#a60',label='RI baseline vs shock'); ax[1].plot(sf,s1,'s-',color='#068',label='RI platform vs shock')
ax[1].set_xlabel("θ / shock fraction"); ax[1].set_ylabel("Resilience index RI")
ax[1].set_title("(b) Sensitivity of resilience index"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig("fig_sensitivity.png", dpi=200); plt.close()
print("figures done")
