from orbit_fits import *

from matplotlib.legend_handler import HandlerBase
from itertools import cycle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#from matplotlib import rc, rcParams
#rc('font',**{'family':'serif','serif':['Times New Roman'],'size':14})
#rc('text',usetex=True)
#rc('patch',antialiased=False)
#rcParams['font.size'] = 15.5
#rcParams['xtick.direction'] = 'in'
#rcParams['ytick.direction'] = 'in'
#rcParams['xtick.major.size'] = 6
#rcParams['ytick.major.size'] = 6
#rcParams['xtick.minor.size'] = 3
#rcParams['ytick.minor.size'] = 3

plt.style.use('tex')

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


def summarize_orbits(pot='fiducial'):
    """"""
    
    names = get_names()
    N = len(names)
    
    ecc = np.empty((N,5))
    rperi = np.empty((N,5))
    rapo = np.empty((N,5))
    
    tout = Table([names, ecc, rperi, rapo], names=('name', 'ecc', 'rperi', 'rapo'))
    
    for i in range(N):
        t = Table.read('../data/output/orbit_props_{:s}_combined.fits'.format(names[i]))
        ind_fiducial = t['potential']==pot
        for k in ['ecc', 'rperi', 'rapo']:
            tout[i][k][0] = np.nanmedian(t[k][ind_fiducial])
            tout[i][k][1] = tout[i][k][0] - np.nanpercentile(t[k][ind_fiducial], 16)
            tout[i][k][2] = np.nanpercentile(t[k][ind_fiducial], 84) - tout[i][k][0]
            tout[i][k][3] = np.nanstd(t[k][ind_fiducial])
            tout[i][k][4] = np.nanstd(t[k])
    
    tout.pprint()
    tout.write('../data/orbital_summary_{:s}.fits'.format(pot), overwrite=True)

def table_orbits():
    """Create a latex table with streams' orbital properties"""
    
    tin = Table.read('../data/orbital_summary.fits')
    N = len(tin)
    
    tm = Table.read('../data/gc_masses.txt', format='ascii.commented_header')
    
    f = open('../paper/table_properties.tex', 'w')
    for i in range(N):
        label = get_properties(tin['name'][i])['label']
        if label=='Ophiuchus':
            f.write('{:s} & {:4.2f}$^{{+{:.3f}}}_{{-{:.3f}}}$ & {:4.1f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:3.2f}\pm{:3.2f}\pm{:3.2f}$ \\\\ \n'.format(label, tin['rperi'][i][0], tin['rperi'][i][2], tin['rperi'][i][1], tin['rapo'][i][0], tin['rapo'][i][2], tin['rapo'][i][1], tm['logM0'][i], tm['error_stat'][i], tm['error_sys'][i]))
        elif label=='Fimbulthul':
            f.write('{:s} & {:4.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:4.1f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:3.2f}\pm{:3.2f}\pm{:3.2f}$ \\\\ \n'.format(label, tin['rperi'][i][0], tin['rperi'][i][2], tin['rperi'][i][1], tin['rapo'][i][0], tin['rapo'][i][2], tin['rapo'][i][1], tm['logM0'][i], tm['error_stat'][i], tm['error_sys'][i]))
        else:
            f.write('{:s} & {:4.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:4.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & ${:3.2f}\pm{:3.2f}\pm{:3.2f}$ \\\\ \n'.format(label, tin['rperi'][i][0], tin['rperi'][i][2], tin['rperi'][i][1], tin['rapo'][i][0], tin['rapo'][i][2], tin['rapo'][i][1], tm['logM0'][i], tm['error_stat'][i], tm['error_sys'][i]))
        
        
        #f.write('{:s} & ${:.1f}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & '.format(tall['name'][i], tall['ra0'][i], r_fit[i,0,0], r_fit[i,1,0], r_fit[i,2,0], r_fit[i,0,1], r_fit[i,1,1], r_fit[i,2,1], r_fit[i,0,2], r_fit[i,1,2], r_fit[i,2,2], r_fit[i,0,3], r_fit[i,1,3], r_fit[i,2,3], r_fit[i,0,4], r_fit[i,1,4], r_fit[i,2,4]))
        
        #f.write('{:s} & {:4.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:4.2f} & {:4.2f} & {:4.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:4.2f} & {:4.2f} & ${:3.2f}\pm{:3.2f}$ \\\\ \n'.format(label, tin['rperi'][i][0], tin['rperi'][i][2], tin['rperi'][i][1], tin['rperi'][i][3], tin['rperi'][i][4], tin['rapo'][i][0], tin['rapo'][i][2], tin['rapo'][i][1], tin['rapo'][i][3], tin['rapo'][i][4], tin['rapo'][i][0], tin['rapo'][i][2]))

        
        #f.write('{:s} & {:4.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:4.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & ${:3.2f}\pm{:3.2f}$ \\\\ \n'.format(r_rp[i][0], r_rp[i][1], r_rp[i][2], r_ra[i][0], r_ra[i][1], r_ra[i][2], tm['logm'][i], tm['logm_err'][i]))

    f.close()


def summarize_fits(pot='fiducial'):
    """"""
    names = get_names()
    N = len(names)
    
    ra = np.empty(N)
    dec = np.empty((N,4))
    dist = np.empty((N,4))
    pmra = np.empty((N,4))
    pmdec = np.empty((N,4))
    vr = np.empty((N,4))
    
    tout = Table([names, ra, dec, dist, pmra, pmdec, vr], names=('name', 'ra', 'dec', 'dist', 'pmra', 'pmdec', 'vr'))
    
    for i in range(N):
        if pot=='fiducial':
            t = Table.read('../data/output/samples_{:s}.fits'.format(names[i]))
        else:
            t = Table.read('../data/output/samples_{:s}_{:s}.fits'.format(names[i], pot))
        
        tout[i]['ra'] = np.median(t['ra'])
        
        for k in ['dec', 'dist', 'pmra', 'pmdec', 'vr']:
            tout[i][k][0] = np.nanmedian(t[k])
            tout[i][k][1] = tout[i][k][0] - np.nanpercentile(t[k], 16)
            tout[i][k][2] = np.nanpercentile(t[k], 84) - tout[i][k][0]
            tout[i][k][3] = np.nanstd(t[k])
    
    tout.pprint()
    tout.write('../data/fit_summary_{:s}.fits'.format(pot), overwrite=True)

def table_fits():
    """Create a latex table with streams' orbital properties"""
    
    tin = Table.read('../data/fit_summary_fiducial.fits')
    N = len(tin)
    
    f = open('../paper/table_fits.tex', 'w')
    for i in range(N):
        label = get_properties(tin['name'][i])['label']
        #if label=='Ophiuchus':
            #f.write('{:s} & {:4.2f}$^{{+{:.3f}}}_{{-{:.3f}}}$ & {:4.1f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:3.2f}\pm{:3.2f}\pm{:3.2f}$ \\\\ \n'.format(label, tin['rperi'][i][0], tin['rperi'][i][2], tin['rperi'][i][1], tin['rapo'][i][0], tin['rapo'][i][2], tin['rapo'][i][1], tm['logM0'][i], tm['error_stat'][i], tm['error_sys'][i]))
        #elif label=='Fimbulthul':
            #f.write('{:s} & {:4.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:4.1f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & ${:3.2f}\pm{:3.2f}\pm{:3.2f}$ \\\\ \n'.format(label, tin['rperi'][i][0], tin['rperi'][i][2], tin['rperi'][i][1], tin['rapo'][i][0], tin['rapo'][i][2], tin['rapo'][i][1], tm['logM0'][i], tm['error_stat'][i], tm['error_sys'][i]))
        #else:
            #f.write('{:s} & {:4.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:4.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & ${:3.2f}\pm{:3.2f}\pm{:3.2f}$ \\\\ \n'.format(label, tin['rperi'][i][0], tin['rperi'][i][2], tin['rperi'][i][1], tin['rapo'][i][0], tin['rapo'][i][2], tin['rapo'][i][1], tm['logM0'][i], tm['error_stat'][i], tm['error_sys'][i]))
        
        
        #f.write('{:s} & ${:.1f}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ & {:.1f}$^{{+{:.1f}}}_{{-{:.1f}}}$ \\\\ \n'.format(tin['name'][i], tin['ra'][i], tin['dec'][i][0], tin['dec'][i][2], tin['dec'][i][1], tin['dist'][i][0], tin['dist'][i][2], tin['dist'][i][1], tin['pmra'][i][0], tin['pmra'][i][2], tin['pmra'][i][1], tin['pmdec'][i][0], tin['pmdec'][i][2], tin['pmdec'][i][1], tin['vr'][i][0], tin['vr'][i][2], tin['vr'][i][1]))

        f.write('{:s} & ${:.1f}$ & {:.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & {:.2f}$^{{+{:.2f}}}_{{-{:.2f}}}$ & '.format(label, tin['ra'][i], tin['dec'][i][0], tin['dec'][i][2], tin['dec'][i][1], tin['dist'][i][0], tin['dist'][i][2], tin['dist'][i][1], tin['pmra'][i][0], tin['pmra'][i][2], tin['pmra'][i][1], tin['pmdec'][i][0], tin['pmdec'][i][2], tin['pmdec'][i][1], tin['vr'][i][0], tin['vr'][i][2], tin['vr'][i][1]))
        
        provenance = get_properties(tin['name'][i])['provenance']
        if np.isfinite(provenance[-1]):
            f.write('[{:d}, {:d}, {:d}]\\\\ \n'.format(*provenance))
        else:
            f.write('[{:d}, {:d}, N/A]\\\\ \n'.format(*provenance[:2]))
        
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
        if stream.fra:
            ra = stream.x0*u.deg
            dec, dist, pmra, pmdec, vr = p0
        else:
            dec = stream.x0*u.deg
            ra, dist, pmra, pmdec, vr = p0
        
        c = coord.SkyCoord(ra=ra, dec=dec, distance=dist, pm_ra_cosdec=pmra, pm_dec=pmdec, radial_velocity=vr, frame='icrs')
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
    #plt.ylabel('Dec [deg]')
    #plt.xlim(185, 255)

    handles, labels = plt.gca().get_legend_handles_labels()
    h_ = [handles[-1], tuple([handles[1], handles[0], handles[2]])]
    l_ = [labels[-1], labels[0]]
    plt.legend(h_, l_, fontsize='small', frameon=False, loc=0, handler_map={tuple: HandlerTupleVert(ndivide=None, pad=0.5)})

    plt.sca(ax[1])
    #plt.ylabel('Distance [kpc]')
    #plt.ylim(3.5, 6.5)

    plt.sca(ax[2])
    #plt.ylabel('$\mu_\\alpha$ [mas yr$^{-1}$]')
    #plt.ylim(-4.5, 6)

    plt.sca(ax[3])
    #plt.xlabel(xlabel)
    #plt.ylabel('$\mu_\delta$ [mas yr$^{-1}$]')
    #plt.ylim(2.1, 8.5)

    # Galactocentric coordinates -- long orbits
    # right side
    
    # plot models
    n_steps_long = np.abs(int((5*u.Gyr/stream.dt).decompose()))
    
    for j in range(nplot):
        p0 = [x*y.unit for x, y in zip(flatchain_short[j], stream.p0)]
        #dec, dist, pmra, pmdec, vr = p0
        
        #c = coord.SkyCoord(ra=stream.ra0*u.deg, dec=dec, distance=dist, pm_ra_cosdec=pmra, pm_dec=pmdec, radial_velocity=vr, frame='icrs')
        #w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)

        if stream.fra:
            ra = stream.x0*u.deg
            dec, dist, pmra, pmdec, vr = p0
        else:
            dec = stream.x0*u.deg
            ra, dist, pmra, pmdec, vr = p0
        
        c = coord.SkyCoord(ra=ra, dec=dec, distance=dist, pm_ra_cosdec=pmra, pm_dec=pmdec, radial_velocity=vr, frame='icrs')
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

    rmax = 90
    xmin, xmax = plt.gca().get_xlim()
    if xmin<-rmax: xmin = -rmax
    if xmax>rmax: xmax = rmax
    plt.gca().set_xlim(xmin, xmax)
    
    ymin, ymax = plt.gca().get_ylim()
    if ymin<-rmax: ymin = -rmax
    if ymax>rmax: ymax = rmax
    plt.gca().set_ylim(ymin, ymax)

    plt.sca(ax[6])
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.xlabel('$\\rho$ [kpc]')
    plt.ylabel('z [kpc]')

    rmax = 100
    ymin, ymax = plt.gca().get_ylim()
    if (ymin<-rmax) | (ymax>rmax):
        ymin = -rmax
        ymax = rmax
    plt.gca().set_ylim(ymin, ymax)
    
    xmin, xmax = plt.gca().get_xlim()
    if xmin<0: xmin = 0
    if xmax>rmax: xmax = rmax
    #print(ymin, ymax)
    plt.gca().set_xlim(xmin, xmax)

    plt.savefig('../paper/figures/fig1_all/{:s}.png'.format(name), bbox_inches='tight', pad_inches=0.7)
    if name=='fjorm':
        plt.savefig('../paper/figures/stream_fitting.pdf')

def save_xml():
    """"""
    f = open('../paper/streams.xml', 'w')
    
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n<catalog>\n')

    for name in get_names():
        f.write('<stream>\n')
        f.write('<name>{:s}</name>\n'.format(name))
        label = get_properties(name)['label']
        label = label.replace('\\"{o}', 'o')
        #label = label.replace('\\"{o}', '&ouml;')
        print(label)
        f.write('<label>{:s}</label>\n'.format(label))
        f.write('</stream>\n')

    f.write('</catalog>\n')


def orbital_precision():
    """"""
    tin = Table.read('../data/orbital_summary_fiducial.fits')
    tin.pprint()
    
    print(np.median(tin['rperi'][:,3]), np.median(tin['rperi'][:,3]/tin['rperi'][:,0]))
    print(np.median(tin['rapo'][:,3]), np.median(tin['rapo'][:,3]/tin['rapo'][:,0]))

def get_elz(names):
    """"""
    
    lx = np.zeros(len(names)) * u.kpc**2 / u.Myr
    ly = np.zeros(len(names)) * u.kpc**2 / u.Myr
    lz = np.zeros(len(names)) * u.kpc**2 / u.Myr
    etot = np.zeros(len(names)) * u.kpc**2 * u.Myr**-2
    
    for i, name in enumerate(names[:]):
        stream = Stream(name, ham=ham, save_ext='')

        # best-fit orbit
        res = pickle.load(open('../data/fits/minimization_{:s}.pkl'.format(stream.savename), 'rb'))
        orbit = stream.orbital_properties(pbest=res.x, t=stream.tstream)

        # long orbit
        orbit = stream.orbital_properties(pbest=res.x, t=5*u.Gyr)
        l = orbit.angular_momentum()
        lx[i] = np.nanmedian(l[0])
        ly[i] = np.nanmedian(l[1])
        lz[i] = np.nanmedian(l[2])
        etot[i] = np.nanmedian(orbit.energy())
    
    return (lx, ly, lz, etot)

def overall_summary():
    """"""
    to = Table.read('../data/orbital_summary_fiducial.fits')
    tm = Table.read('../data/gc_masses.txt', format='ascii.commented_header')
    
    if len(to)>len(tm):
        tm = Table.read('../data/gc_masses_tri.txt', format='ascii.commented_header')
    
    lx, ly, lz, etot = get_elz(to['name'])
    label = [get_properties(name)['label'] for name in to['name']]
    
    # origin
    prog_dwarf = ['elqui', 'indus', 'jhelum']
    prog_gc = ['triangulum', 'gd1', 'gjoll', 'atlas', 'aliqa_uma', 'phoenix']
    progenitor = ['n/a'] * len(to)
    
    for e, name in enumerate(to['name']):
        if name in prog_dwarf: progenitor[e] = 'dwarf'
        if name in prog_gc: progenitor[e] = 'gc'
    
    tout = hstack([to, tm])
    tout['lx'] = lx
    tout['ly'] = ly
    tout['lz'] = lz
    tout['etot'] = etot
    tout['label'] = label
    tout['progenitor'] = progenitor
    tout.pprint()
    
    tout.write('../data/overall_summary.fits', overwrite=True)


def print_masses():
    """"""
    t = Table.read('../data/overall_summary.fits')
    mass = np.array(10**t['logM0'])
    
    for i in range(len(t)):
        print('{:s} {:.2e}'.format(t['name'][i], mass[i]))
    
    print('median {:.2e}'.format(np.median(mass)))


def gc_ascii():
    """"""
    t = Table.read('../data/baumgardt_positions.fits')
    t.pprint()
    print(t.colnames)
    t.write('../data/baumgardt_positions.txt', format='ascii.commented_header')

def gc_orbits():
    """"""
    
    t = Table.read('../data/baumgardt_positions.fits')
    t.pprint()
    
    N = len(t)
    
    c = coord.SkyCoord(ra=t['RAJ2000'], dec=t['DEJ2000'], distance=t['Rsun'], pm_ra_cosdec=t['pmRA_'], pm_dec=t['pmDE'], radial_velocity=t['RV'], frame='icrs')
    w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)
    dt = -0.5*u.Myr
    T = 5*u.Gyr
    nstep = abs(int((T/dt).decompose()))

    rperi = np.empty(N)*u.kpc
    rapo = np.empty(N)*u.kpc
    lx = np.empty(N)*u.kpc**2*u.Myr**-1
    ly = np.empty(N)*u.kpc**2*u.Myr**-1
    lz = np.empty(N)*u.kpc**2*u.Myr**-1
    etot = np.empty(N)*u.kpc**2*u.Myr**-2
    
    for i in range(N):
        orbit = ham.integrate_orbit(w0[i], dt=dt, n_steps=nstep)
        l = orbit.angular_momentum()
        lx[i] = np.nanmedian(l[0])
        ly[i] = np.nanmedian(l[1])
        lz[i] = np.nanmedian(l[2])
        
        rperi[i] = orbit.pericenter()
        rapo[i] = orbit.apocenter()
        etot[i] = np.nanmedian(orbit.energy())
    
    tout = Table([t['Name'], rperi, rapo, lx, ly, lz, etot], names=('name', 'rperi', 'rapo', 'lx', 'ly', 'lz', 'etot'))
    tout.pprint()
    tout.write('../data/gc_orbits.fits', overwrite=True)

def conste_rapo(rperi, e):
    """"""
    return rperi * (1+e)/(1-e)

def associations():
    """"""
    t = Table.read('../data/overall_summary.fits')
    tgc = Table.read('../data/gc_orbits.fits')
    
    plt.close()
    fig, ax = plt.subplots(1,1,figsize=(5.5,6.5))
    
    s = 0.15*t['logM0']**5
    im = plt.scatter(t['rperi'][:,0], t['rapo'][:,0], c=t['Lz'], vmin=-2, vmax=2, zorder=1, cmap='PuOr', s=s, edgecolors='k', linewidths=0.5, label='')
    plt.errorbar(t['rperi'][:,0], t['rapo'][:,0], xerr=(t['rperi'][:,1], t['rperi'][:,2]), yerr=(t['rapo'][:,1], t['rapo'][:,2]), fmt='none', zorder=0, color='k', lw=0.5, label='')
    
    
    # constant eccentricity
    x_ = np.linspace(2,100,10)
    xlabel = np.array([7.6, 3.15, 2.5])
    
    for i, e in enumerate([0, 0.5, 0.75]):
        ra_ = conste_rapo(x_, e)
        plt.plot(x_, ra_, 'k--', lw=0.5, alpha=0.3, zorder=0, label='')
        plt.text(xlabel[i], conste_rapo(xlabel[i], e), 'e = {:g}'.format(e), fontsize='xx-small', alpha=0.3, bbox=dict(fc='w', ec='none'), rotation=50, va='center', ha='center')
    
    f = 0.05
    N = len(t)
    #np.random.seed(12846)
    #phase = (np.random.rand(N) + -1) * 0.5*np.pi
    phase = np.array([-0.4, -0.1, -0.3, 0, -0.1, 0.2, -0.2, -0.1, -0.1, 0.2, 0.4, -0.1, 0.1, 0.2, 0.45, -0.2, -0.1, 0.1, 0, -0.3, -0.25, 0.2]) * np.pi
    for i in range(N):
        if t['name'][i]=='willka_yaku':
            label = 'Willka\nYaku'
            va = 'top'
        else:
            label = t['label'][i]
            va = 'center'
        plt.text(t['rperi'][i,0] * (1.03 + f*np.cos(phase[i])), t['rapo'][i,0] * (1 + f*np.sin(phase[i])), '{:s}'.format(label), fontsize='xx-small', va=va) #, bbox=dict(fc='w', ec='none', alpha=0.3))
    
    # legend
    logm = np.array([3.5, 4, 4.5])
    Nm = len(logm)
    for i in range(Nm):
        sz = 0.15 * logm[i]**5
        plt.scatter(16, 200, c='w', ec='k', linewidths=0.5, s=sz, label='log M$_0$/M$_\odot$ = {:.1f}'.format(logm[i]))
    
    plt.legend(loc=4, frameon=False, fontsize='x-small', handlelength=0.6, scatteryoffsets=[0.75])
    
    plt.xlim(2,33)
    plt.ylim(6, 100)
    #plt.gca().set_aspect('equal')
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    
    plt.xticks([5,10,20])
    plt.yticks([10,50,100])
    plt.gca().xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
    plt.gca().yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))

    plt.xlabel('Pericenter [kpc]')
    plt.ylabel('Apocenter [kpc]', labelpad=-15)
    plt.tight_layout()
    
    pos = ax.get_position()
    #cax = plt.axes([pos.x1+0.02, pos.y0, 0.03, pos.y1-pos.y0])
    #plt.colorbar(im, cax=cax, label='z angular momentum [kpc$^2$ Myr$^{-1}$]')
    
    cax = plt.axes([pos.x0, pos.y1+0.015, pos.x1-pos.x0, 0.025])
    plt.colorbar(im, cax=cax, label='$L_z$ [kpc$^2$ Myr$^{-1}$]', orientation='horizontal', ticklocation='top')
    #cax.xaxis.set_label_position('top')
    #cax.xaxis.set_tick_params(labeltop='on')
    
    
    plt.savefig('../plots/streams_peri_apo.png')
    plt.savefig('../paper/figures/streams_peri_apo.pdf', bbox_incehs='tight')

def associations_2panel(label_all=False):
    """"""
    t = Table.read('../data/overall_summary.fits')
    tgc = Table.read('../data/gc_orbits.fits')
    tb = Table.read('../data/baumgardt_positions.fits')

    
    plt.close()
    fig, ax = plt.subplots(1,2,figsize=(10,6), sharex=True, sharey=True)
    
    for j in range(2):
        plt.sca(ax[j])
        
        # plot streams
        s = 0.15*t['logM0']**5
        im = plt.scatter(t['rperi'][:,0], t['rapo'][:,0], c=t['Lz'], vmin=-2, vmax=2, zorder=2, cmap='PuOr', s=s, edgecolors='k', linewidths=0.5, label='')
        plt.errorbar(t['rperi'][:,0], t['rapo'][:,0], xerr=(t['rperi'][:,1], t['rperi'][:,2]), yerr=(t['rapo'][:,1], t['rapo'][:,2]), fmt='none', zorder=0, color='k', lw=0.5, label='')
        
        # constant eccentricity
        x_ = np.linspace(0.1,100,10)
        xlabel = np.array([6, 3.15, 2.05])
        
        for i, e in enumerate([0, 0.5, 0.75]):
            ra_ = conste_rapo(x_, e)
            plt.plot(x_, ra_, 'k--', lw=0.5, alpha=0.3, zorder=0, label='')
            plt.text(xlabel[i], conste_rapo(xlabel[i], e), 'e = {:g}'.format(e), fontsize='xx-small', alpha=0.3, bbox=dict(fc='w', ec='none'), rotation=44, va='center', ha='center')
    
        plt.xlim(1.5, 33)
        plt.ylim(4.5, 100)
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        
        plt.xticks([5,10,20])
        plt.yticks([10,50,100])
        plt.gca().xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        plt.gca().yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))

        plt.xlabel('Pericenter [kpc]')
        if j==0:
            plt.ylabel('Apocenter [kpc]', labelpad=-15)
        
    plt.tight_layout(w_pad=0.5)
    
    # add colorbar
    pos = ax[1].get_position()
    cax = plt.axes([pos.x1+0.01, pos.y0, 0.02, pos.y1-pos.y0])
    plt.colorbar(im, cax=cax, label='z angular momentum [kpc$^2$ Myr$^{-1}$]')
    
    #cax = plt.axes([pos.x0, pos.y1+0.015, pos.x1-pos.x0, 0.025])
    #plt.colorbar(im, cax=cax, label='$L_z$ [kpc$^2$ Myr$^{-1}$]', orientation='horizontal', ticklocation='top')


    plt.sca(ax[0])
    # add stream labels
    f = 0.05
    N = len(t)
    #np.random.seed(12846)
    #phase = (np.random.rand(N) + -1) * 0.5*np.pi
    #print(np.array(t['label']))
    phase = np.array([-0.4, -0.1, -0.3, 0, -0.1, 0.2, -0.2, -0.1, -0.1, 0.2, 0.4, -0.1, 0.1, 0.2, 0.45, -0.2, -0.1, 0.1, 0, -0.3, -0.25, 0.2]) * np.pi
    for i in range(N):
        if t['name'][i]=='willka_yaku':
            label = 'Willka\nYaku'
            va = 'top'
        else:
            label = t['label'][i]
            va = 'center'
        plt.text(t['rperi'][i,0] * (1.03 + f*np.cos(phase[i])), t['rapo'][i,0] * (1 + f*np.sin(phase[i])), '{:s}'.format(label), fontsize='xx-small', va=va) #, bbox=dict(fc='w', ec='none', alpha=0.3))
    
    # legend
    logm = np.array([3.5, 4, 4.5])
    Nm = len(logm)
    for i in range(Nm):
        sz = 0.15 * logm[i]**5
        plt.scatter(16, 200, c='w', ec='k', linewidths=0.5, s=sz, label='log M$_0$/M$_\odot$ = {:.1f}'.format(logm[i]))
    
    plt.legend(loc=4, frameon=False, fontsize='x-small', handlelength=0.6, scatteryoffsets=[0.75])
    
    
    # add globular clusters
    plt.sca(ax[1])
    gc_sequoia = np.array(['NGC 5466', 'NGC 7006', 'IC 4499'])
    ind_sequoia = np.array([True if x in gc_sequoia else False for x in tgc['name']])
    
    gc_ge = np.array(['NGC 288', 'NGC 362', 'NGC 1261', 'NGC 1851', 'NGC 1904', 'NGC 2298', 'NGC 2808', 'NGC 4147', 'NGC 4833', 'NGC 5286', 'NGC 5897', 'NGC 6205', 'NGC 6235', 'NGC 6284', 'NGC 6341', 'NGC 6779', 'NGC 6864', 'NGC 7089', 'NGC 7099', 'NGC 7492'])
    ind_ge = np.array([True if x in gc_ge else False for x in tgc['name']])
    
    gc_helmi = np.array(['NGC 4590', 'NGC 5024', 'NGC 5053', 'NGC 5272', 'NGC 6981'])
    ind_helmi = np.array([True if x in gc_helmi else False for x in tgc['name']])
    
    gc_sgr = np.array(['NGC 2419', 'NGC 5824', 'NGC 6715', 'Pal 12', 'Ter 7', 'Ter 8', 'Arp 2', 'Whiting 1'])
    ind_sgr = np.array([True if x in gc_sgr else False for x in tgc['name']])

    gc_kraken = np.array(['NGC 5946', 'NGC 5986', 'NGC 6093', 'NGC 6121', 'NGC 6144', 'NGC 6254', 'NGC 6273', 'NGC 6287', 'NGC 6541', 'NGC 6544', 'NGC 6681', 'NGC 6712', 'NGC 6809'])
    ind_kraken = np.array([True if x in gc_kraken else False for x in tgc['name']])
    
    ind_smooth = ~ind_sequoia & ~ind_ge & ~ind_helmi & ~ind_sgr & ~ind_kraken
    
    #plt.plot(tgc['rperi'][ind_sequoia], tgc['rapo'][ind_sequoia], 'k^', ms=8)
    #plt.plot(tgc['rperi'][ind_ge], tgc['rapo'][ind_ge], 'ko', ms=8)
    #plt.plot(tgc['rperi'][ind_helmi], tgc['rapo'][ind_helmi], 'kp', ms=8)
    #plt.plot(tgc['rperi'][ind_sgr], tgc['rapo'][ind_sgr], 'kD', ms=8)
    #plt.plot(tgc['rperi'][ind_kraken], tgc['rapo'][ind_kraken], 'ks', ms=8)
    #plt.scatter(tgc['rperi'][ind_smooth], tgc['rapo'][ind_smooth], c=tgc['lz'][ind_smooth], vmin=-2, vmax=2, zorder=10, cmap='PuOr', s=20, edgecolors='k', linewidths=0.5, label='')
    
    labels = ['Sequoia', 'Gaia - Enceladus', 'Helmi', 'Sagittarius', 'Kraken', 'Unassociated']
    indices = [ind_sequoia, ind_ge, ind_helmi, ind_sgr, ind_kraken, ind_smooth]
    markers = ['^', 'H', 'p', 'D', 's', 'o']
    sizes = [60, 80, 80, 40, 60, 20]

    for i in range(6):
        plt.scatter(tgc['rperi'][indices[i]], tgc['rapo'][indices[i]], c=tgc['lz'][indices[i]], vmin=-2, vmax=2, zorder=10, cmap='PuOr', s=sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label='')
    
        # legend entries
        #if i!=4:
        plt.scatter(0, 0, facecolors='w', s=sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label=labels[i])
    
    # label globular clusters
    if label_all:
        for i in range(len(tgc)):
            if (tgc['rperi'][i]>2) & (tgc['rapo'][i]>6) & (tgc['rapo'][i]<100):
                plt.text(tgc['rperi'][i], tgc['rapo'][i], '{:s}'.format(tgc['name'][i]), fontsize=8)
    else:
        gc_tolabel = ['NGC 288', 'NGC 5024', 'NGC 5824', 'NGC 3201', 'NGC 5272', 'NGC 4590']
        for i in range(len(tgc)):
            if tgc['name'][i] in gc_tolabel:
                print(tgc['name'][i], tb['RAJ2000'][i], tb['DEJ2000'][i])
                plt.text(tgc['rperi'][i]*1.01, tgc['rapo'][i]*1.01, '{:s}'.format(tgc['name'][i]), fontsize='xx-small', zorder=15)

    plt.legend(loc=4, frameon=False, fontsize='x-small', handlelength=0.6, title='Globular clusters', title_fontsize='x-small')
    #print(np.isin(tgc['name'], ind_sequoia))
    #print(tgc['name']==ind_sequoia[0])
    
    if label_all:
        plt.savefig('../plots/streams_peri_apo_labeled.png')
    else:
        plt.savefig('../plots/streams_peri_apo.png')
        plt.savefig('../paper/figures/streams_peri_apo.pdf', bbox_incehs='tight')


def rapo_ecc():
    """"""
    t = Table.read('../data/overall_summary.fits')
    
    plt.close()
    fig, ax = plt.subplots(1,1,figsize=(10,10))
    
    s = 0.05*t['logM0']**6
    im = plt.scatter(t['rapo'][:,0], t['ecc'][:,0], c=t['Lz'], vmin=-2, vmax=2, zorder=1, cmap='PuOr', s=s, edgecolors='k', linewidths=0.5)
    #plt.errorbar(t['ecc'][:,0], t['rapo'][:,0], xerr=(t['ecc'][:,1], t['ecc'][:,2]), yerr=(t['rapo'][:,1], t['rapo'][:,2]), fmt='none', zorder=0, color='k', lw=0.5)
    
    #f = 0.05
    #N = len(t)
    #np.random.seed(12846)
    #phase = (np.random.rand(N) + -1) * 0.5*np.pi
    #for i in range(N):
        #plt.text(t['ecc'][i,0] * (1 + f*np.cos(phase[i])), t['rapo'][i,0] * (1 + f*np.sin(phase[i])), '${:s}$'.format(t['label'][i]), fontsize='x-small') #, bbox=dict(fc='w', ec='none', alpha=0.3))
    
    #plt.ylim(1, 100)
    #plt.gca().set_aspect('equal')
    #plt.gca().set_xscale('log')
    #plt.gca().set_yscale('log')
    plt.ylabel('Eccentricity')
    plt.xlabel('Apocenter [kpc]')
    plt.tight_layout()
    
    pos = ax.get_position()
    print(pos)
    cax = plt.axes([pos.x1+0.02, pos.y0, 0.03, pos.y1-pos.y0])
    plt.colorbar(im, cax=cax, label='z angular momentum [kpc$^2$ Myr$^{-1}$]')

def rapo_lz():
    """"""
    t = Table.read('../data/overall_summary.fits')
    
    plt.close()
    plt.figure()
    
    plt.plot(t['Lz'], t['rapo'][:,0], 'ko')
    
    plt.tight_layout()

def lz():
    """"""
    t = Table.read('../data/overall_summary.fits')
    print(np.median(t['Lz']))
    
    plt.close()
    plt.figure()
    
    plt.hist(t['Lz'], bins=20)
    plt.axvline(np.median(t['Lz']), color='tab:red')
    
    plt.tight_layout()


def orbital_phase():
    """Print current galactocentric radius, as well as the peri and apocenter"""
    
    t = Table.read('../data/overall_summary.fits')
    #t.pprint()
    
    for e, name in enumerate(t['name'][:]):
        stream = Stream(name)
        sampler = pickle.load(open('../data/fits/mcmc_{:s}.pkl'.format(stream.savename), 'rb'))
        chain = sampler.chain[:,256:,:]
        flatchain = np.reshape(chain,(-1,5))
        pbest = np.median(flatchain, axis=0)
        
        p0 = [x*y.unit for x, y in zip(pbest, stream.p0)]
        dec, dist, pmra, pmdec, vr = p0
        
        c = coord.SkyCoord(ra=stream.ra0*u.deg, dec=dec, distance=dist, pm_ra_cosdec=pmra, pm_dec=pmdec, radial_velocity=vr, frame='icrs')
        w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)

        orbit = stream.ham.integrate_orbit(w0, dt=stream.dt, n_steps=stream.nstep)
        #model = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=stream.gc_frame)
        
        rcurrent = np.median(orbit.spherical.pos.distance)
        rperi = t['rperi'][e][0]*u.kpc
        rapo = t['rapo'][e][0]*u.kpc
        f = (rcurrent - rperi) / (rapo - rperi)
        print('{:s} {:.2f}'.format(name, f))
        #print(, , t['rapo'][e][0])
        #print(np.median(orbit.cartesian.pos.norm()))



# follow up

def elz(labeled=0):
    """"""
    t = Table.read('../data/overall_summary.fits')
    tgc = Table.read('../data/gc_orbits.fits')
    
    print(len(t))
    
    # add globular clusters
    gc_sequoia = np.array(['NGC 5466', 'NGC 7006', 'IC 4499'])
    ind_sequoia = np.array([True if x in gc_sequoia else False for x in tgc['name']])
    
    gc_ge = np.array(['NGC 288', 'NGC 362', 'NGC 1261', 'NGC 1851', 'NGC 1904', 'NGC 2298', 'NGC 2808', 'NGC 4147', 'NGC 4833', 'NGC 5286', 'NGC 5897', 'NGC 6205', 'NGC 6235', 'NGC 6284', 'NGC 6341', 'NGC 6779', 'NGC 6864', 'NGC 7089', 'NGC 7099', 'NGC 7492'])
    ind_ge = np.array([True if x in gc_ge else False for x in tgc['name']])
    
    gc_helmi = np.array(['NGC 4590', 'NGC 5024', 'NGC 5053', 'NGC 5272', 'NGC 6981'])
    ind_helmi = np.array([True if x in gc_helmi else False for x in tgc['name']])
    
    gc_sgr = np.array(['NGC 2419', 'NGC 5824', 'NGC 6715', 'Pal 12', 'Ter 7', 'Ter 8', 'Arp 2', 'Whiting 1'])
    ind_sgr = np.array([True if x in gc_sgr else False for x in tgc['name']])

    gc_kraken = np.array(['NGC 5946', 'NGC 5986', 'NGC 6093', 'NGC 6121', 'NGC 6144', 'NGC 6254', 'NGC 6273', 'NGC 6287', 'NGC 6541', 'NGC 6544', 'NGC 6681', 'NGC 6712', 'NGC 6809'])
    ind_kraken = np.array([True if x in gc_kraken else False for x in tgc['name']])
    
    ind_smooth = ~ind_sequoia & ~ind_ge & ~ind_helmi & ~ind_sgr & ~ind_kraken
    
    plt.close()
    fig, ax = plt.subplots(1,1,figsize=(12,12))
    
    labels = ['Sequoia', 'Gaia - Enceladus', 'Helmi', 'Sagittarius', 'Kraken', 'Unassociated']
    indices = [ind_sequoia, ind_ge, ind_helmi, ind_sgr, ind_kraken, ind_smooth]
    markers = ['^', 'H', 'p', 'D', 's', 'o']
    sizes = [60, 80, 80, 40, 60, 20]
    
    uscale = (1*u.kpc**2*u.Myr**-2).to(u.km**2*u.s**-2).value * 1e-6
    #print(uscale)
    uscale = 1.

    for i in range(6):
        plt.scatter(tgc['lz'][indices[i]], tgc['etot'][indices[i]]*uscale, c=0*tgc['lz'][indices[i]], vmin=-2, vmax=2, zorder=10, cmap='PuOr', s=1.5*sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label='')
    
        # legend entries
        plt.scatter(0, 0, facecolors='w', s=sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label=labels[i])
    
    # label globular clusters
    if labeled:
        for i in range(len(tgc)):
            if (tgc['lz'][i]>-4) & (tgc['lz'][i]<4) & (tgc['etot'][i]>-0.17) & (tgc['etot'][i]<0):
                plt.text(tgc['lz'][i], tgc['etot'][i]*uscale, '{:s}'.format(tgc['name'][i]), fontsize='xx-small')

    plt.legend(loc=4, frameon=False, fontsize='small', handlelength=0.6, title='Globular clusters', title_fontsize='small')

    # streams
    s = 0.15*t['logM0']**5
    plt.scatter(t['lz'], t['etot'], c=0*t['lz']+0.5, vmin=-2, vmax=2, zorder=1, cmap='PuOr', s=s, edgecolors='k', linewidths=0.5, label='')
    
    if labeled:
        for i in range(len(t)):
            if t['name'][i]=='willka_yaku':
                label = 'Willka\nYaku'
                va = 'top'
            else:
                label = t['label'][i]
                va = 'center'
            plt.text(t['lz'][i], t['etot'][i], '{:s}'.format(label), va=va, fontsize='x-small')
    
    plt.xlim(-4,4)
    plt.ylim(-0.265, -0.01)
    
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$E_{tot}$ [kpc$^2$ Myr$^{-2}$]')
    
    plt.tight_layout()
    if labeled:
        plt.savefig('../plots/elz_gc_streams_labeled.png')
    else:
        plt.savefig('../plots/elz_gc_streams.png')

def elz_field():
    """"""
    t = Table.read('../data/overall_summary.fits')
    tgc = Table.read('../data/gc_orbits.fits')
    
    # add globular clusters
    gc_sequoia = np.array(['NGC 5466', 'NGC 7006', 'IC 4499'])
    ind_sequoia = np.array([True if x in gc_sequoia else False for x in tgc['name']])
    
    gc_ge = np.array(['NGC 288', 'NGC 362', 'NGC 1261', 'NGC 1851', 'NGC 1904', 'NGC 2298', 'NGC 2808', 'NGC 4147', 'NGC 4833', 'NGC 5286', 'NGC 5897', 'NGC 6205', 'NGC 6235', 'NGC 6284', 'NGC 6341', 'NGC 6779', 'NGC 6864', 'NGC 7089', 'NGC 7099', 'NGC 7492'])
    ind_ge = np.array([True if x in gc_ge else False for x in tgc['name']])
    
    gc_helmi = np.array(['NGC 4590', 'NGC 5024', 'NGC 5053', 'NGC 5272', 'NGC 6981'])
    ind_helmi = np.array([True if x in gc_helmi else False for x in tgc['name']])
    
    gc_sgr = np.array(['NGC 2419', 'NGC 5824', 'NGC 6715', 'Pal 12', 'Ter 7', 'Ter 8', 'Arp 2', 'Whiting 1'])
    ind_sgr = np.array([True if x in gc_sgr else False for x in tgc['name']])

    gc_kraken = np.array(['NGC 5946', 'NGC 5986', 'NGC 6093', 'NGC 6121', 'NGC 6144', 'NGC 6254', 'NGC 6273', 'NGC 6287', 'NGC 6541', 'NGC 6544', 'NGC 6681', 'NGC 6712', 'NGC 6809'])
    ind_kraken = np.array([True if x in gc_kraken else False for x in tgc['name']])
    
    ind_smooth = ~ind_sequoia & ~ind_ge & ~ind_helmi & ~ind_sgr & ~ind_kraken
    
    plt.close()
    fig, ax = plt.subplots(1,2,figsize=(12,6), sharex=True, sharey=True)
    
    labels = ['Sequoia', 'Gaia - Enceladus', 'Helmi', 'Sagittarius', 'Kraken', 'Unassociated']
    indices = [ind_sequoia, ind_ge, ind_helmi, ind_sgr, ind_kraken, ind_smooth]
    markers = ['^', 'H', 'p', 'D', 's', 'o']
    sizes = [60, 80, 80, 40, 60, 20]

    plt.sca(ax[0])
    for i in range(6):
        plt.scatter(tgc['lz'][indices[i]], tgc['etot'][indices[i]], c=tgc['lz'][indices[i]], vmin=-2, vmax=2, zorder=10, cmap='PuOr', s=1.5*sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label='')
    
        # legend entries
        #if i!=4:
        plt.scatter(0, 0, facecolors='w', s=sizes[i], marker=markers[i], edgecolors='k', linewidths=0.5, label=labels[i])
    
    # label globular clusters
    for i in range(len(tgc)):
        if (tgc['lz'][i]>2) & (tgc['lz'][i]<4):
            plt.text(tgc['lz'][i], tgc['etot'][i], '{:s}'.format(tgc['name'][i]), fontsize='x-small')

    plt.legend(loc=4, frameon=False, fontsize='x-small', handlelength=0.6, title='Globular clusters', title_fontsize='x-small')

    # streams
    s = 0.15*t['logM0']**5
    plt.scatter(t['Lz'], t['Etot'], c=t['Lz'], vmin=-2, vmax=2, zorder=1, cmap='PuOr', s=s, edgecolors='k', linewidths=0.5, label='')
    
    for i in range(len(t)):
        if t['name'][i]=='willka_yaku':
            label = 'Willka\nYaku'
            va = 'top'
        else:
            label = t['label'][i]
            va = 'center'
        plt.text(t['Lz'][i], t['Etot'][i], '{:s}'.format(label), va=va, fontsize='x-small')
    
    plt.xlim(-4,4)
    plt.ylim(-0.265, -0.01)
    plt.ylim(-0.165, -0.01)
    
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    plt.ylabel('$E_{tot}$ [kpc$^2$ Myr$^{-2}$]')
    
    # field stars
    plt.sca(ax[1])
    
    th = Table.read('../data/h3_giants.fits')
    print(len(th))
    ind_mpoor = (th['FeH']<-1.) #& (th['SNR']>10)
    th = th[ind_mpoor]
    plt.scatter(th['Lz'], th['E_tot_pot1'], c=th['FeH'], vmin=-3, cmap='magma', s=10)
    
    plt.xlabel('$L_z$ [kpc$^2$ Myr$^{-1}$]')
    
    cbaxes = inset_axes(ax[1], width='100%', height='100%', bbox_to_anchor=(1.02,0,0.04,1), bbox_transform=ax[1].transAxes, borderpad=0) 
    cb = plt.colorbar(cax=cbaxes, orientation='vertical', label='[Fe/H]', ticklocation='right')
    cb.ax.tick_params(labelsize='small')
    cb.ax.set_ylabel(cb.ax.get_ylabel(), fontsize='small')
    
    plt.tight_layout()
    plt.savefig('../plots/elz_gc_streams_h3.png')

    

def mdf():
    """"""
    th = Table.read('../data/h3_giants.fits')
    #ind_mpoor = (th['FeH']<-1.5) #& (th['SNR']>10)
    #th = th[ind_mpoor]
    
    plt.close()
    plt.figure()
    
    feh_bins = np.linspace(-3, -1.5, 20)
    ind_retrograde = (th['Lz']>1.7) & (th['E_tot_pot1']>-0.09)
    ind_prograde = (th['Lz']<-1.7) & (th['E_tot_pot1']>-0.09)
    ind_midretrograde = (th['Lz']>1.) & (th['Lz']<1.7) & (th['E_tot_pot1']>-0.11)
    
    plt.hist(th['FeH'], bins=feh_bins, histtype='step', density=True, color='0.2')
    plt.hist(th['FeH'][ind_retrograde], bins=feh_bins, histtype='step', density=True)
    plt.hist(th['FeH'][ind_prograde], bins=feh_bins, histtype='step', density=True)
    #plt.hist(th['FeH'][ind_midretrograde], bins=feh_bins, histtype='step', density=True)
    #plt.hist(th['FeH'][ind_retrograde], bins=feh_bins, histtype='step')
    #plt.hist(th['FeH'][ind_midretrograde], bins=feh_bins, histtype='step')
    
    
    plt.tight_layout()

def h3_giants():
    """"""
    t = Table.read('/home/ana/data/rcat.fits')
    ind_giant = (t['logg']<3.5) & (t['FLAG']==0)
    tg = t[ind_giant]
    
    tg['Lz'] = (tg['Lz']*u.kpc*u.km/u.s).to(u.kpc**2*u.Myr**-1)
    tg['E_tot_pot1'] = (tg['E_tot_pot1']*u.km**2*u.s**-2).to(u.kpc**2*u.Myr**-2)
    
    tg.write('../data/h3_giants.fits', overwrite=True)
