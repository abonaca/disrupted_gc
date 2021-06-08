from orbit_fits import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.style.use('urf')


def get_names():
    """Get names of streams in the sample"""
    
    streams = ['ophiuchus', 'gd1', 'svol', 'leiptr', 'gjoll', 'fjorm', 'fimbulthul', 'ylgr', 'sylgr', 'slidr', 'phlegethon', 'phoenix', 'turranburra', 'indus', 'elqui', 'jhelum', 'atlas', 'aliqa_uma', 'ravi', 'wambelong', 'willka_yaku', 'turbio']
    
    return sorted(streams)

def proposal():
    """"""
    
    text_color = '0.95'
    
    # set up figure geometry
    plt.close()
    fig = plt.figure(figsize=(12,9.4), facecolor='k')
    
    x0, x1 = 0.04, 0.98
    y0, y1, y2, y3 = 0.02, 0.48, 0.51, 0.98
    ddx = 0.02
    dx = (x1 - x0 - ddx*2) / 3
    
    ax0 = plt.axes([x0, y2, x1-x0, y3-y2])
    ax_ = [plt.axes([x0+i*(dx + ddx), y0, dx, y1-y0]) for i in range(3)]
    
    ax = [ax0] + ax_
    
    ###########
    # Top panel
    # plot all streams
    plt.sca(ax[0])
    d = pickle.load(open('../data/galactocentric_streams_orbits.pkl', 'rb'))
    names = np.array(get_names())
    Nstreams = len(names)
    
    for j in range(Nstreams):
        ms = 2
        alpha = 0.15
        ds = d[names[j]]
        cgal = ds['stream']
        
        phi = 180*u.deg
        stream_x = cgal.x * np.cos(phi) + cgal.y * np.sin(phi)
        stream_y = - cgal.x * np.sin(phi) + cgal.y * np.cos(phi)
        
        #ind_plot = stream_y < rho_in + i*dr
        plt.plot(stream_x, cgal.z, 'o', color=ds['color'], ms=ms)
        #if names[j] in ['gd1', 'jhelum', 'ophiuchus']:
            #plt.text(stream_x[0].value, cgal.z[0].value, names[j], color=ds['color'], fontsize='xx-small')
    
    # annotate Galactic center
    plt.plot(0, 0, 'wx')
    plt.text(0,-1, 'Galactic\nCentre', ha='center', va='top', color=text_color, fontsize='x-small')
    
    # annotate the Sun
    xsun = -gc_frame.galcen_distance * np.cos(phi)
    zsun = gc_frame.z_sun.to(u.kpc)
    plt.plot(xsun, zsun, color='w', marker='$\odot$', ms=10, mew=0.4)
    plt.text(xsun.value+1, zsun.value-1, 'Sun', color=text_color, fontsize='x-small')
    
    text_intro = '''Stellar streams orbiting the Milky Way
    are gravitational antennae that have
    picked up curious signals:'''
    plt.text(0.2, 0.87, text_intro, color=text_color, wrap=True, transform=plt.gca().transAxes, va='top', ha='center')
    plt.xlim(-40,60)
    plt.ylim(-30,13)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().set_facecolor('k')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    plt.ylabel('TODAY', color=text_color)
    for spine in plt.gca().spines.values():
        spine.set_edgecolor('0.9')
        spine.set_linewidth(4)

    # show GD-1
    gd1_color = '#13d100'
    g = Table.read('/home/ana/projects/legacy/GD1-DR2/output/gd1_members.fits')

    ax00 = inset_axes(ax0, '100%', '100%', bbox_to_anchor=[0.58, 0.69, 0.38, 0.2], loc=3, bbox_transform=ax0.transAxes)
    plt.sca(ax00)
    plt.scatter(g['phi1'], g['phi2'], s=g['pmem']*1.5, c=g['pmem'], cmap=mpl.cm.binary, vmin=0.5, vmax=1.1, zorder=0, label='', rasterized=True)
    #plt.plot(g['phi1'], g['phi2'], 'k.', ms=1, alpha=0.5)
    plt.text(-40,-2, 'gap', fontsize='x-small', ha='center')
    plt.text(-32, 2, 'spur', fontsize='x-small', ha='center')
    
    plt.xlim(-60, -20)
    plt.ylim(-4,4)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    for spine in plt.gca().spines.values():
        spine.set_edgecolor(gd1_color)
        spine.set_linewidth(2)
    
    plt.text(0, 1.1, 'Price-Whelan & Bonaca (2018)', transform=plt.gca().transAxes, color=gd1_color, fontsize='small')
    plt.text(0, -0.1, 'Were the GD-1 gap and spur created by an impact\nof a compact object, like a dark-matter subhalo?', va='top', transform=plt.gca().transAxes, fontsize='small', color=gd1_color)
    
    plt.sca(ax[0])
    plt.annotate('', xy=(0,1), xycoords=ax00.transAxes, xytext=(13.5, 5.5), textcoords=ax[0].transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=gd1_color, alpha=0.3))
    plt.annotate('', xy=(0,0), xycoords=ax00.transAxes, xytext=(13.5, 5.5), textcoords=ax[0].transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=gd1_color, alpha=0.3))
    
    # show Pal 5
    tp = Table.read('/home/ana/projects/legacy/Pal5sBiggestFan/data/pal5_likely_members.fits')
    pal5_color = '#00c7ff'
    
    ax01 = inset_axes(ax0, '100%', '100%', bbox_to_anchor=[0.62, 0.02, 0.2, 0.43], loc=3, bbox_transform=ax0.transAxes)
    plt.sca(ax01)
    plt.plot(tp['ra'], tp['dec'], 'k.', ms=0.5, alpha=0.1, rasterized=True)
    plt.text(223.5, -5, 'fan', fontsize='x-small')
    plt.text(228.5, 0.5, 'Pal 5', fontsize='x-small')
    
    plt.xlim(242,222)
    plt.ylim(-10,10)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    for spine in plt.gca().spines.values():
        spine.set_edgecolor(pal5_color)
        spine.set_linewidth(2)
    
    plt.text(0, 1.05, 'Bonaca et al. (2020)', transform=plt.gca().transAxes, color=pal5_color, fontsize='small')
    pal5_text = '''Was the leading tail
of the Palomar 5
stream broadened
into a fan by a
large perturbation,
like the rotating
Galactic bar?'''
    plt.text(1.05, 1, pal5_text, va='top', transform=plt.gca().transAxes, fontsize='small', color=pal5_color)
    
    plt.sca(ax[0])
    plt.annotate('', xy=(0,1), xycoords=ax01.transAxes, xytext=(6, 2), textcoords=ax[0].transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=pal5_color, alpha=0.3), zorder=0)
    plt.annotate('', xy=(0,0), xycoords=ax01.transAxes, xytext=(6, 2), textcoords=ax[0].transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=pal5_color, alpha=0.3), zorder=0)
    
    # show Jhelum
    jhelum_color = '#ff9e00'
    tj = Table.read('/home/ana/projects/legacy/jhelum/data/jhelum_likely_members.fits')
    
    ax02 = inset_axes(ax0, '100%', '100%', bbox_to_anchor=[0.02, 0.2, 0.3, 0.3], loc=3, bbox_transform=ax0.transAxes)
    plt.sca(ax02)
    plt.plot(tj['phi1'], tj['phi2'], 'k.', ms=1.5, alpha=0.3, rasterized=True)
    plt.text(10, -3.5, 'wide component', fontsize='x-small', ha='center', bbox=dict(facecolor='w', alpha=0.4, edgecolor='none'))
    plt.text(-2, 2, 'narrow component', fontsize='x-small', bbox=dict(facecolor='w', alpha=0.4, edgecolor='none'))
    
    plt.ylim(-5,5)
    plt.xlim(-4, 24)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    for spine in plt.gca().spines.values():
        spine.set_edgecolor(jhelum_color)
        spine.set_linewidth(2)
    plt.text(0, 1.05, 'Bonaca et al. (2019)', transform=plt.gca().transAxes, color=jhelum_color, fontsize='small')
    plt.text(0, -0.1, "Is Jhelum's wide envelope a signature of\nits nature set in the early universe or has\nit been nurtured by the Milky Way?", va='top', transform=plt.gca().transAxes, fontsize='small', color=jhelum_color)
    
    plt.sca(ax[0])
    plt.annotate('', xy=(1,1), xycoords=ax02.transAxes, xytext=(-1, -10), textcoords=ax[0].transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=jhelum_color, alpha=0.3), zorder=0)
    plt.annotate('', xy=(1,0), xycoords=ax02.transAxes, xytext=(-1, -10), textcoords=ax[0].transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=jhelum_color, alpha=0.3), zorder=0)
    
    #############
    # Bottom left
    tsgr = Table.read('/home/ana/projects/h3/data/SgrTriax_DYN.dat.gz', format='ascii')
    tsgr = tsgr[::10]
    c_sgr = coord.ICRS(ra=tsgr['ra']*u.deg, dec=tsgr['dec']*u.deg, distance=tsgr['dist']*u.kpc)
    sgr_gal = c_sgr.transform_to(coord.Galactocentric)
    
    d = pickle.load(open('../data/galactocentric_streams_orbits_long.pkl', 'rb'))
    ds = d['slidr']
    cgal = ds['stream']
    orbit = ds['orbit']
    
    plt.sca(ax[1])
    plt.text(0.5,0.9, 'What is the origin of halo stars?', ha='center', wrap=True, transform=plt.gca().transAxes, color=text_color)
    plt.text(0.04, 0.87, '$_{^{\\blacktriangleright}}$ I will derive orbits of halo stars and\n associate structures of common origin', va='top', wrap=True, transform=plt.gca().transAxes, color=jhelum_color, fontsize='small')
    
    plt.plot(sgr_gal.x, sgr_gal.z, 'o', color='#ffc15d', ms=1)
    plt.plot(cgal.x, cgal.z, 'o', color='orangered', label='')
    plt.plot(orbit.x, orbit.z, '-', color='orangered', lw=2)
    
    # annotate orientation marks
    plt.plot(0, 0, 'wx')
    plt.plot(-gc_frame.galcen_distance, gc_frame.z_sun.to(u.kpc), color='w', marker='$\odot$', ms=10, mew=0.4)
    
    slidr_text = '''Slidr, a globular cluster stream, has an orbit\nsimilar to Sagittarius, indicating a possible link'''
    plt.text(0.5, 0.02, slidr_text, transform=plt.gca().transAxes, fontsize='x-small', color=text_color, ha='center', va='bottom') #, bbox=dict(facecolor='w', alpha=0.2))
    plt.text(55,-45, 'Sagittarius', fontsize='x-small', color='#ffc15d')
    plt.text(-45,35, 'Slidr', fontsize='x-small', color='orangered')
    
    plt.ylim(-115,135)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().set_facecolor('k')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    plt.ylabel('FUTURE RESEARCH PLAN', color=text_color)
    for spine in plt.gca().spines.values():
        spine.set_edgecolor(jhelum_color)
        spine.set_linewidth(4)

    ###############
    # Bottom center
    dl = pickle.load(open('../data/galactocentric_lmc_elqui.pkl', 'rb'))
    orbit_on = dl['orbits'][:,1]
    orbit_off = dl['orbits0'][:,1]
    ds = d['elqui']
    cgal = ds['stream']
    
    plt.sca(ax[2])
    plt.text(0.5,0.9, 'How did our Galaxy evolve?', ha='center', wrap=True, transform=plt.gca().transAxes, color=text_color)
    plt.text(0.04, 0.87, '$_{^{\\blacktriangleright}}$ I will quantify the dynamical impact\nof the most massive progenitors', va='top', wrap=True, transform=plt.gca().transAxes, color=pal5_color, fontsize='small', bbox=dict(facecolor='k', alpha=0.5))
    
    plt.plot(orbit_on.x, orbit_on.z, '-', color='#3dd4ff', lw=2)
    plt.plot(orbit_off.x, orbit_off.z, '--', color='#0081de', lw=2)
    plt.plot(cgal.x, cgal.z, 'o', color='#3dd4ff', label='')
    
    # annotate orientation marks
    plt.plot(0, 0, 'wx')
    plt.plot(-gc_frame.galcen_distance, gc_frame.z_sun.to(u.kpc), color='w', marker='$\odot$', ms=10, mew=0.4)
    
    elqui_text = '''
Massive satellites like the LMC change orbits of
less massive satellites, like the Elqui stream'''
    plt.text(0.5, 0.02, elqui_text, transform=plt.gca().transAxes, fontsize='x-small', color=text_color, ha='center', va='bottom')
    plt.text(15,-44, 'w/o LMC', fontsize='x-small', color=text_color)
    plt.text(11.6,-20, 'with LMC', fontsize='x-small', color=text_color)
    
    plt.xlim(-25,25)
    plt.ylim(-90,110)
    #plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().set_facecolor('k')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    for spine in plt.gca().spines.values():
        spine.set_edgecolor(pal5_color)
        spine.set_linewidth(4)

    ##############
    # Bottom right
    plt.sca(ax[3])
    plt.text(0.5,0.9, 'What is the nature of dark matter?', ha='center', wrap=True, transform=plt.gca().transAxes, color=text_color, bbox=dict(facecolor='k', alpha=0.5))
    plt.text(0.04, 0.87, "$_{^{\\blacktriangleright}}$ I will measure the granularity of the\nMilky Way's dark matter halo", va='top', transform=plt.gca().transAxes, color=gd1_color, fontsize='small', bbox=dict(facecolor='k', alpha=0.5))
    
    ds = d['fimbulthul']
    cgal = ds['stream']
    orbit = ds['orbit']
    plt.plot(cgal.x, cgal.z, 'o', color='#4fff3d', label='')
    Nf = 3000
    plt.plot(orbit.x[:Nf], orbit.z[:Nf], '-', color='#4fff3d', lw=1)
    
    ds = d['willka_yaku']
    cgal = ds['stream']
    orbit = ds['orbit']
    plt.plot(cgal.x, cgal.z, 'o', color='#0e9800', label='')
    plt.plot(orbit.x, orbit.z, '-', color='#0e9800', lw=1)
    
    # annotate orientation marks
    plt.plot(0, 0, 'wx')
    plt.plot(-gc_frame.galcen_distance, gc_frame.z_sun.to(u.kpc), color='w', marker='$\odot$', ms=10, mew=0.4)
    
    fimbulthul_text = '''
If dark matter is fuzzy, nearby
streams like Fimbulthul will be
heated by quantum turbulence'''
    plt.text(0.98, 0.47, fimbulthul_text, transform=plt.gca().transAxes, fontsize='x-small', color=text_color, ha='right', va='top')
    
    willka_text = '''If dark matter is cold, distant streams like Willka\nYaku will reveal gaps from subhalo impacts'''
    plt.text(0.5, 0.02, willka_text, transform=plt.gca().transAxes, fontsize='x-small', color=text_color, ha='center', va='bottom')
    
    plt.xlim(-30,30)
    #plt.ylim(50,-50)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().set_facecolor('k')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    for spine in plt.gca().spines.values():
        spine.set_edgecolor(gd1_color)
        spine.set_linewidth(4)
    
    # panel labels
    plt.sca(ax[0])
    plt.text(0.0125, 0.945, 'a', fontsize='small', transform=plt.gca().transAxes, color=text_color)
    
    for i, label in enumerate(['b', 'c', 'd']):
        plt.sca(ax[i+1])
        plt.text(0.05, 0.71, label, fontsize='small', transform=plt.gca().transAxes, color=text_color)

    
    plt.savefig('../plots/proposal_summary.png', facecolor=fig.get_facecolor())
    plt.savefig('../plots/proposal_summary.pdf', facecolor=fig.get_facecolor())

