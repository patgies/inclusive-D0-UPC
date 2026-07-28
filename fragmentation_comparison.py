import os

import numpy as np
import matplotlib.pyplot as plt
import lhapdf

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 8,
    "ytick.major.size": 8,
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
})

# BCFY fragmentation 

def BCFY_DP(z, r):
    if z <= 0.0 or z >= 1.0:
        return 0.0
    denom = 1.0 - (1.0 - r) * z
    factor = r * z * (1.0 - z)**2 / denom**6
    poly = (6.0
            - 18.0 * (1.0 - 2.0*r) * z
            + (21.0 - 74.0*r + 68.0*r**2) * z**2
            - 2.0*(1.0 - r)*(6.0 - 19.0*r + 18.0*r**2) * z**3
            + 3.0*(1.0 - r)**2 * (1.0 - 2.0*r + 2.0*r**2) * z**4)
    return factor * poly

def BCFY_DV(z, r):
    if z <= 0.0 or z >= 1.0:
        return 0.0
    denom = 1.0 - (1.0 - r) * z
    factor = 3.0 * r * z * (1.0 - z)**2 / denom**6
    poly = (2.0
            - 2.0*(3.0 - 2.0*r) * z
            + 3.0*(3.0 - 2.0*r + 4.0*r**2) * z**2
            - 2.0*(1.0 - r)*(4.0 - r + 2.0*r**2) * z**3
            + (1.0 - r)**2 * (3.0 - 2.0*r + 2.0*r**2) * z**4)
    return factor * poly

def Dc_to_D0(z, r):
    mD     = 1.8648
    mDstar = 2.0067
    ratio  = mDstar / mD

    DP_part = 0.168 * BCFY_DP(z, r)

    z_scaled = ratio * z
    step     = 1.0 if z < (mD / mDstar) else 0.0
    DV_part  = 0.39 * step * BCFY_DV(z_scaled, r) * ratio

    return DP_part + DV_part

# Kniehl & Kramer fragmentation 

def D_kniehl_kramer(z, N=0.694, eps=0.101):
    if z <= 0.0 or z >= 1.0:
        return 0.0
    return N * z * (1.0 - z)**2 / ((1.0 - z)**2 + eps * z)**2

# LHAPDF: prompt D0 fragmentation (central member) 

def D_lhapdf(z, Q=1.5, setname="prompt-D0-1-109", pid=4):
    pdf = lhapdf.mkPDF(setname, 0)
    if np.isscalar(z):
        if z <= 0.0 or z >= 1.0:
            return 0.0
        return pdf.xfxQ(pid, z, Q) / z
    return np.array([pdf.xfxQ(pid, zi, Q) / zi if 0.0 < zi < 1.0 else 0.0 for zi in z])

# Create the LHAPDF object once and reuse it for the grid evaluation.
central_pdf = lhapdf.mkPDF("prompt-D0-1-109", 0)
all_members = lhapdf.mkPDFs("prompt-D0-1-109")
replicas = all_members[1:]  # members 1-100 are the replica sets


def D_lhapdf_cached(z, Q=1.5, pid=4, pdf_obj=None):
    if pdf_obj is None:
        pdf_obj = central_pdf
    if np.isscalar(z):
        if z <= 0.0 or z >= 1.0:
            return 0.0
        return pdf_obj.xfxQ(pid, z, Q) / z
    return np.array([pdf_obj.xfxQ(pid, zi, Q) / zi if 0.0 < zi < 1.0 else 0.0 for zi in z])

#  Evaluate on a grid

z = np.linspace(1e-3, 0.9999, 2000)

r = 0.1
bcfy         = np.array([Dc_to_D0(zi, r) for zi in z])
kniehl_kramer = np.array([D_kniehl_kramer(zi) for zi in z])
central_member_ff = np.array([D_lhapdf_cached(zi, Q=1.5, pdf_obj=central_pdf) for zi in z])
lhapdf_ff = np.array([D_lhapdf_cached(zi, Q=1.5, pdf_obj=central_pdf) for zi in z])

# Replica band for the LHAPDF FF
vals = np.array([[D_lhapdf_cached(zi, Q=1.5, pdf_obj=rep) for zi in z] for rep in replicas])
median = np.median(vals, axis=0)
lo, hi = np.percentile(vals, [16, 84], axis=0)

#Plot 

fig, ax = plt.subplots(figsize=(6.5, 5.5))

ax.plot(z, bcfy,    color='royalblue', lw=2,
        label=r'BCFY ($r=0.1$) [NLO, Braaten et al., $Q^2=1.5^2$ GeV$^2$]')
ax.plot(z, kniehl_kramer, color='crimson',  lw=2, linestyle='--',
        label=r'Kniehl \& Kramer ($N=0.694,\,\varepsilon=0.101$)  [LO, $Q^2=1.5^2$ GeV$^2$]')
ax.plot(z, central_member_ff, color='darkgreen', lw=2, linestyle='-.',
        label=r'LHAPDF central member  [NLO, $Q^2=1.5^2$ GeV$^2$]')
#ax.plot(z, lhapdf_ff, color='forestgreen', lw=2, linestyle=':',
#        label=r'LHAPDF median over replicas')
ax.fill_between(z, lo, hi, color='forestgreen', alpha=0.2, label=r'LHAPDF 68\% band')
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$D_{c \to D^0}(z)$', labelpad=15)
ax.set_title(r'Fragmentation functions')
ax.legend(fontsize=10)
ax.set_xlim(0, 1)
ax.set_ylim(bottom=0)

# The x-axis "0.0" and y-axis "0" labels both sit right at the
# bottom-left corner and overlap there -- hide the y-axis zero label
# instead (the x-axis one stays, so the origin is still labeled once).
y_ticks = ax.get_yticks()
y_labels = []
for tick in y_ticks:
    if tick == 0:
        y_labels.append('')
    else:
        y_labels.append(f'{tick:g}')
ax.set_yticks(y_ticks)
ax.set_yticklabels(y_labels)

plt.tight_layout()

script_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(script_dir, 'plots', 'fragmentation_comparison.pdf')
plt.savefig(output_path, dpi=150, bbox_inches='tight')
print('Saved ' + output_path)
