import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
'''
Makes a plot like Fig 4 in Verbunt & Phinney (1995)
'''
infile = '../../RGEB_tideinfo.txt'
colors = ['#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026']

systems, eccs, tides = np.loadtxt(infile, usecols=(0,1,2), comments='#', unpack=True)

fig = plt.figure()
cm = plt.cm.get_cmap('YlOrRd')

ax = fig.add_subplot(1,1,1)
plt.axis([5, -8, 0, 0.49])
plt.xlabel(r'$\log [-(\Delta \ln e)]$', size=30)
plt.ylabel(r'$e$', size=30)
plt.axvline(x=0, ls=':', color='k')
for system, ecc, tide, color in zip(systems, eccs, tides, colors):
    points = ax.scatter(tide, ecc, s=200, c=color)
    if system == 3955867:
        ax.annotate(str(int(system)), xy=(tide, ecc), xytext=(tide+0.3, ecc+0.01), size=20)
    else:
        ax.annotate(str(int(system)), xy=(tide, ecc), xytext=(tide+0.8, ecc+0.01), size=20)

axnew = fig.add_subplot(15,1,1)
cmap = mpl.colors.ListedColormap(colors)
bounds = [19.38446, 20.6864, 33.65685, 120.3903, 197.9182, 207.1082, 235.29852, 250]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb = mpl.colorbar.ColorbarBase(axnew, cmap=cmap, norm=norm, orientation='horizontal')
axnew.get_yaxis().set_ticks([])
axnew.get_xaxis().set_ticks([])
#cb.set_label(r'$P_{orb}$', size=30)

ax.text(4.7, 0.44, 'Short $P_{orb}$')
ax.text(-6.3, 0.44, 'Long $P_{orb}$')

plt.show()