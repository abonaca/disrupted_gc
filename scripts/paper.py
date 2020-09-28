from orbit_fits import *

from matplotlib.legend_handler import HandlerBase
from itertools import cycle

class HandlerTupleVert(HandlerBase):
    """
    Handler for Tuple.

    Additional kwargs are passed through to `HandlerBase`.

    Parameters
    ----------
    ndivide : int, optional
        The number of sections to divide the vertical legend area into. If None,
        use the length of the input tuple. Default is 1.


    pad : float, optional
        If None, fall back to ``legend.borderpad`` as the default.
        In units of fraction of font size. Default is None.
    """
    def __init__(self, ndivide=1, pad=None, **kwargs):

        self._ndivide = ndivide
        self._pad = pad
        HandlerBase.__init__(self, **kwargs)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):

        handler_map = legend.get_legend_handler_map()

        if self._ndivide is None:
            ndivide = len(orig_handle)
        else:
            ndivide = self._ndivide

        if self._pad is None:
            pad = legend.borderpad * fontsize
        else:
            pad = self._pad * fontsize

        if ndivide > 1:
            height = (height - pad*(ndivide - 1)) / ndivide
#             width = (width - pad * (ndivide - 1)) / ndivide

#         xds_cycle = cycle(xdescent - (width + pad) * np.arange(ndivide))
        yds_cycle = cycle(ydescent - (height + pad) * (np.arange(ndivide)+0.5*(ndivide-1)))

        a_list = []
        for handle1 in orig_handle:
            handler = legend.get_legend_handler(handler_map, handle1)
            _a_list = handler.create_artists(
                legend, handle1,
                xdescent, next(yds_cycle), width, height, fontsize, trans)
#                 next(xds_cycle), ydescent, width, height, fontsize, trans)
            a_list.extend(_a_list)

        return a_list

def summarize_orbits():
    """"""
    
    names = get_names()
    N = len(names)
    
    ecc = np.empty((N,5))
    rperi = np.empty((N,5))
    rapo = np.empty((N,5))
    
    tout = Table([names, ecc, rperi, rapo], names=('name', 'ecc', 'rperi', 'rapo'))
    
    for i in range(N):
        t = Table.read('../data/orbit_props_{:s}_combined.fits'.format(names[i]))
        ind_fiducial = t['potential']=='fiducial'
        for k in ['ecc', 'rperi', 'rapo']:
            tout[i][k][0] = np.nanmedian(t[k][ind_fiducial])
            tout[i][k][1] = tout[i][k][0] - np.nanpercentile(t[k][ind_fiducial], 16)
            tout[i][k][2] = np.nanpercentile(t[k][ind_fiducial], 84) - tout[i][k][0]
            tout[i][k][3] = np.nanstd(t[k][ind_fiducial])
            tout[i][k][4] = np.nanstd(t[k])
    
    tout.pprint()
    tout.write('../data/orbital_summary.fits', overwrite=True)

def table_orbits():
    """Create a latex table with streams' orbital properties"""
    
    tin = Table.read('../data/orbital_summary.fits')
    N = len(tin)
    
    f = open('../paper/table_properties.tex', 'w')
    for i in range(N):
        label = get_properties(tin['name'][i])['label']
        if label=='Ophiuchus':
            f.write('{:s} & {:4.2f}$^{{+{:.3f}}}_{{-{:.3f}}}$ & {:4.1f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:3.2f}\pm{:3.2f}$ \\\\ \n'.format(label, tin['rperi'][i][0], tin['rperi'][i][2], tin['rperi'][i][1], tin['rapo'][i][0], tin['rapo'][i][2], tin['rapo'][i][1], tin['rapo'][i][0], tin['rapo'][i][2]))
        elif label=='Fimbulthul':
            f.write('{:s} & {:4.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:4.1f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:3.2f}\pm{:3.2f}$ \\\\ \n'.format(label, tin['rperi'][i][0], tin['rperi'][i][2], tin['rperi'][i][1], tin['rapo'][i][0], tin['rapo'][i][2], tin['rapo'][i][1], tin['rapo'][i][0], tin['rapo'][i][2]))
        else:
            f.write('{:s} & {:4.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:4.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & ${:3.2f}\pm{:3.2f}$ \\\\ \n'.format(label, tin['rperi'][i][0], tin['rperi'][i][2], tin['rperi'][i][1], tin['rapo'][i][0], tin['rapo'][i][2], tin['rapo'][i][1], tin['rapo'][i][0], tin['rapo'][i][2]))
        
        
        #f.write('{:s} & ${:.1f}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & '.format(tall['name'][i], tall['ra0'][i], r_fit[i,0,0], r_fit[i,1,0], r_fit[i,2,0], r_fit[i,0,1], r_fit[i,1,1], r_fit[i,2,1], r_fit[i,0,2], r_fit[i,1,2], r_fit[i,2,2], r_fit[i,0,3], r_fit[i,1,3], r_fit[i,2,3], r_fit[i,0,4], r_fit[i,1,4], r_fit[i,2,4]))
        
        #f.write('{:s} & {:4.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:4.2f} & {:4.2f} & {:4.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:4.2f} & {:4.2f} & ${:3.2f}\pm{:3.2f}$ \\\\ \n'.format(label, tin['rperi'][i][0], tin['rperi'][i][2], tin['rperi'][i][1], tin['rperi'][i][3], tin['rperi'][i][4], tin['rapo'][i][0], tin['rapo'][i][2], tin['rapo'][i][1], tin['rapo'][i][3], tin['rapo'][i][4], tin['rapo'][i][0], tin['rapo'][i][2]))

        
        #f.write('{:s} & {:4.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:4.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & ${:3.2f}\pm{:3.2f}$ \\\\ \n'.format(r_rp[i][0], r_rp[i][1], r_rp[i][2], r_ra[i][0], r_ra[i][1], r_ra[i][2], tm['logm'][i], tm['logm_err'][i]))

    f.close()

def plot_stream_fit_orbit(name):
    """"""
    stream = Stream(name)
    sampler = pickle.load(open('../data/fits/mcmc_{:s}.pkl'.format(stream.savename), 'rb'))
    chain = sampler.chain[:,256:,:]
    flatchain = np.reshape(chain,(-1,5))
    np.random.seed(391)
    flatchain_short = np.random.permutation(flatchain)[:1000,:]
    
    # determine orientation
    if stream.fra:
        ix = 0
        iy = 1
        xlabel = 'R.A. [deg]'
        ylabel = 'Dec [deg]'
    else:
        ix = -1
        #iy = 0
        tmp = stream.data['dec'][1]
        stream.data['dec'][1] = stream.data['dec'][0]
        stream.data['dec'][0] = tmp
        xlabel = 'Dec [deg]'
        ylabel = 'R.A. [deg]'
    
    fields = ['dec', 'dist', 'pmra', 'pmdec', 'vr']
    labels = [ylabel, 'Distance [kpc]', '$\mu_\\alpha$ [mas yr$^{-1}$]', '$\mu_\delta$ [mas yr$^{-1}$]', '$V_r$ [km s$^{-1}$]']
    istart, iend = 0, -1
    
    nplot = 20
    dra = 2
    f = 1.1
    
    plt.close()
    fig = plt.figure(figsize=(11*f,8*f))
    gs1 = mpl.gridspec.GridSpec(5,1)
    gs1.update(left=0.08, right=0.6, top=0.95, bottom=0.08, hspace=0.05)

    gs2 = mpl.gridspec.GridSpec(2,1)
    gs2.update(left=0.7, right=0.975, top=0.95, bottom=0.08, hspace=0.25)

    ax0 = fig.add_subplot(gs1[0])
    ax1 = fig.add_subplot(gs1[1], sharex=ax0)
    ax2 = fig.add_subplot(gs1[2], sharex=ax0)
    ax3 = fig.add_subplot(gs1[3], sharex=ax0)
    ax4 = fig.add_subplot(gs1[4], sharex=ax0)
    ax5 = fig.add_subplot(gs2[0])
    ax6 = fig.add_subplot(gs2[1])
    ax = [ax0, ax1, ax2, ax3, ax4, ax5, ax6]

    # Observed coordinates -- short orbits
    # left side

    np.random.seed(748)
    colors = mpl.cm.Blues(np.random.rand(nplot))
    lw = 0.5
    
    # plot models
    for j in range(nplot):
        p0 = [x*y.unit for x, y in zip(flatchain_short[j], stream.p0)]
        dec, dist, pmra, pmdec, vr = p0
        
        c = coord.SkyCoord(ra=stream.ra0*u.deg, dec=dec, distance=dist, pm_ra_cosdec=pmra, pm_dec=pmdec, radial_velocity=vr, frame='icrs')
        w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)

        orbit = stream.ham.integrate_orbit(w0, dt=stream.dt, n_steps=stream.nstep)
        model = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=stream.gc_frame)
        
        if stream.fra:
            model_x = model.ra.wrap_at(stream.wangle)
            model_y = model.dec
        else:
            model_x = model.dec
            model_y = model.ra.wrap_at(stream.wangle)

        model_fields = [model_y, model.distance, model.pm_ra_cosdec, model.pm_dec, model.radial_velocity]
        
        for i in range(5):
            plt.sca(ax[i])

            plt.plot(model_x[istart:iend], model_fields[i][istart:iend], '-', color=colors[j], label='Orbit samples', zorder=0, lw=lw)

            #if (i==0) & (j==0):
                #plt.legend(fontsize='small', handlelength=1)

    color_data = 'orangered'
    data_cols = stream.data.keys()

    for i in range(5):
        plt.sca(ax[i])
        
        if fields[i] in stream.data.keys():
            plt.plot(stream.data[fields[i]][ix], stream.data[fields[i]][1], '.', color=color_data, label='Observed stream')
            plt.errorbar(stream.data[fields[i]][ix].value, stream.data[fields[i]][1].value, yerr=stream.data[fields[i]][2].value, fmt='none', color=color_data, alpha=0.7, label='')
        
        plt.ylabel(labels[i])
        if i<4:
            plt.gca().tick_params(labelbottom=False)
        #plt.minorticks_on()
        
    plt.xlim(np.min(stream.data['dec'][0].to(u.deg).value)-dra, np.max(stream.data['dec'][0].to(u.deg).value)+dra)
    plt.xlabel(xlabel)

    plt.sca(ax[0])
    plt.ylabel('Dec [deg]')
    #plt.xlim(185, 255)

    handles, labels = plt.gca().get_legend_handles_labels()
    h_ = [handles[-1], tuple([handles[1], handles[0], handles[2]])]
    l_ = [labels[-1], labels[0]]
    plt.legend(h_, l_, fontsize='small', frameon=False, loc=0, handler_map={tuple: HandlerTupleVert(ndivide=None, pad=0.5)})

    plt.sca(ax[1])
    plt.ylabel('Distance [kpc]')
    #plt.ylim(3.5, 6.5)

    plt.sca(ax[2])
    plt.ylabel('$\mu_\\alpha$ [mas yr$^{-1}$]')
    #plt.ylim(-4.5, 6)

    plt.sca(ax[3])
    plt.xlabel('R.A. [deg]')
    plt.ylabel('$\mu_\delta$ [mas yr$^{-1}$]')
    #plt.ylim(2.1, 8.5)

    # Galactocentric coordinates -- long orbits
    # right side
    
    # plot models
    n_steps_long = np.abs(int((5*u.Gyr/stream.dt).decompose()))
    
    for j in range(nplot):
        p0 = [x*y.unit for x, y in zip(flatchain_short[j], stream.p0)]
        dec, dist, pmra, pmdec, vr = p0
        
        c = coord.SkyCoord(ra=stream.ra0*u.deg, dec=dec, distance=dist, pm_ra_cosdec=pmra, pm_dec=pmdec, radial_velocity=vr, frame='icrs')
        w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)

        orbit = stream.ham.integrate_orbit(w0, dt=stream.dt, n_steps=n_steps_long)

        plt.sca(ax[5])
        plt.plot(orbit.cartesian.x, orbit.cartesian.y, '-', color=colors[j], lw=lw)
        
        plt.sca(ax[6])
        plt.plot(orbit.cylindrical.rho, orbit.cartesian.z, '-', color=colors[j], lw=lw)
        
    plt.sca(ax[5])
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')

    plt.sca(ax[6])
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.xlabel('$\\rho$ [kpc]')
    plt.ylabel('z [kpc]')

    plt.savefig('../paper/figures/fig1_all/{:s}.png'.format(name))
    if name=='fjorm':
        plt.savefig('../paper/figures/stream_fitting.pdf')


