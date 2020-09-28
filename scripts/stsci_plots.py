from orbit_fits import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.style.use('lgray_ucon')

def gallery(step=0):
    """"""
    
    text_color = '0.95'
    
    # set up figure geometry
    plt.close()
    fig = plt.figure(figsize=(16,9), facecolor='k')
    ax = plt.axes([0,0,1,1])
    
    # plot all streams
    plt.sca(ax)
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
    plt.text(0,-1, 'Galactic\nCenter', ha='center', va='top', color=text_color, fontsize='small')
    
    # annotate the Sun
    xsun = -gc_frame.galcen_distance * np.cos(phi)
    zsun = gc_frame.z_sun.to(u.kpc)
    plt.plot(xsun, zsun, color='w', marker='$\odot$', ms=10, mew=0.4)
    plt.text(xsun.value+1, zsun.value-1, 'Sun', color=text_color, fontsize='small')
    
    #text_intro = '''Stellar streams orbiting the Milky Way
    #are gravitational antennae that have
    #picked up curious signals:'''
    #plt.text(0.2, 0.87, text_intro, color=text_color, wrap=True, transform=plt.gca().transAxes, va='top', ha='center')
    plt.xlim(-40,60)
    plt.ylim(-30,13)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().set_facecolor('k')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    #plt.ylabel('TODAY', color=text_color)
    #for spine in plt.gca().spines.values():
        #spine.set_edgecolor('0.9')
        #spine.set_linewidth(4)

    if step>2:
        # show GD-1
        gd1_color = '#13d100'
        g = Table.read('/home/ana/projects/legacy/GD1-DR2/output/gd1_members.fits')

        ax00 = inset_axes(ax, '100%', '100%', bbox_to_anchor=[0.58, 0.73, 0.38, 0.15], loc=3, bbox_transform=ax.transAxes)
        plt.sca(ax00)
        plt.scatter(g['phi1'], g['phi2'], s=g['pmem']*1.5, c=g['pmem'], cmap=mpl.cm.binary, vmin=0.5, vmax=1.1, zorder=0, label='', rasterized=True)
        #plt.plot(g['phi1'], g['phi2'], 'k.', ms=1, alpha=0.5)
        plt.text(-40,-2, 'gap', fontsize='small', ha='center')
        plt.text(-32, 2, 'spur', fontsize='small', ha='center')
        
        plt.xlim(-60, -20)
        plt.ylim(-4,4)
        plt.gca().set_aspect('equal', adjustable='datalim')
        plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
        for spine in plt.gca().spines.values():
            spine.set_edgecolor(gd1_color)
            spine.set_linewidth(2)
        
        plt.text(0, 1.1, 'Price-Whelan & Bonaca (2018)', transform=plt.gca().transAxes, color=gd1_color, fontsize='medium')
        plt.text(0, -0.1, 'Were the GD-1 gap and spur created by an impact of a compact\nobject, like a dark-matter subhalo?', va='top', transform=plt.gca().transAxes, fontsize='small', color=gd1_color)
        
        plt.sca(ax)
        plt.annotate('', xy=(0,1), xycoords=ax00.transAxes, xytext=(13.5, 5.5), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=gd1_color, alpha=0.3))
        plt.annotate('', xy=(0,0), xycoords=ax00.transAxes, xytext=(13.5, 5.5), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=gd1_color, alpha=0.3))
        
    if step>1:
        # show Pal 5
        tp = Table.read('/home/ana/projects/legacy/Pal5sBiggestFan/data/pal5_likely_members.fits')
        pal5_color = '#00c7ff'
        
        ax01 = inset_axes(ax, '100%', '100%', bbox_to_anchor=[0.64, 0.2, 0.3, 0.31], loc=3, bbox_transform=ax.transAxes)
        plt.sca(ax01)
        plt.plot(tp['ra'], tp['dec'], 'k.', ms=0.5, alpha=0.1, rasterized=True)
        plt.text(223.5, -5, 'fan', fontsize='small')
        plt.text(228.5, 0.5, 'Pal 5', fontsize='small')
        
        plt.xlim(242,222)
        plt.ylim(-10,10)
        plt.gca().set_aspect('equal', adjustable='datalim')
        plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
        for spine in plt.gca().spines.values():
            spine.set_edgecolor(pal5_color)
            spine.set_linewidth(2)
        
        plt.text(0, 1.05, 'Bonaca et al. (2020a)', transform=plt.gca().transAxes, color=pal5_color, fontsize='medium')
        pal5_text = '''Was the leading tail of the Palomar 5 stream
broadened into a fan by a large perturbation, like
the rotating Galactic bar?
[see also Starkman et al. (2020), #poster-starkman]'''
        plt.text(0, -0.1, pal5_text, va='top', transform=plt.gca().transAxes, fontsize='small', color=pal5_color)
        
        plt.sca(ax)
        plt.annotate('', xy=(0,1), xycoords=ax01.transAxes, xytext=(6, 2), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=pal5_color, alpha=0.3), zorder=0)
        plt.annotate('', xy=(0,0), xycoords=ax01.transAxes, xytext=(6, 2), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=pal5_color, alpha=0.3), zorder=0)
    
    if step>0:
        # show Jhelum
        jhelum_color = '#ff9e00'
        tj = Table.read('/home/ana/projects/legacy/jhelum/data/jhelum_likely_members.fits')
        
        ax02 = inset_axes(ax, '100%', '100%', bbox_to_anchor=[0.01, 0.3, 0.34, 0.17], loc=3, bbox_transform=ax.transAxes)
        plt.sca(ax02)
        plt.plot(tj['phi1'], tj['phi2'], 'k.', ms=1.5, alpha=0.3, rasterized=True)
        plt.text(10, -3.5, 'wide component', fontsize='small', ha='center', bbox=dict(facecolor='w', alpha=0.4, edgecolor='none'))
        plt.text(-2, 2, 'narrow component', fontsize='small', bbox=dict(facecolor='w', alpha=0.4, edgecolor='none'))
        
        plt.ylim(-5,5)
        plt.xlim(-4, 24)
        plt.gca().set_aspect('equal', adjustable='datalim')
        plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
        for spine in plt.gca().spines.values():
            spine.set_edgecolor(jhelum_color)
            spine.set_linewidth(2)
        plt.text(0, 1.05, 'Bonaca et al. (2019b)', transform=plt.gca().transAxes, color=jhelum_color, fontsize='medium')
        jhelum_text='''Is Jhelum's wide envelope a signature of its nature set
in the early universe or has it been nurtured by
the Milky Way?
[see also Malhan et al. (2020), #poster-malhan]
'''
        plt.text(0, -0.1, jhelum_text, va='top', transform=plt.gca().transAxes, fontsize='small', color=jhelum_color)
        
        plt.sca(ax)
        plt.annotate('', xy=(1,1), xycoords=ax02.transAxes, xytext=(-1, -10), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=jhelum_color, alpha=0.3), zorder=0)
        plt.annotate('', xy=(1,0), xycoords=ax02.transAxes, xytext=(-1, -10), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=jhelum_color, alpha=0.3), zorder=0)
    
    plt.savefig('../plots/gallery_perturbed_{:d}.png'.format(step), dpi=120, facecolor='k')

def conclusion(step=0):
    """"""
    
    text_color = '0.95'
    
    # set up figure geometry
    plt.close()
    fig = plt.figure(figsize=(16,9), facecolor='k')
    ax = plt.axes([0,0,1,1])
    
    # plot all streams
    plt.sca(ax)
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
        if (step>0) & (names[j] not in ['gd1', 'jhelum', 'fimbulthul']):
            plt.text(stream_x[-1].value, cgal.z[-1].value, '?', color=ds['color'], fontsize='large', alpha=0.5)
        #if names[j] in ['gd1', 'jhelum', 'ophiuchus']:
            #plt.text(stream_x[0].value, cgal.z[0].value, names[j], color=ds['color'], fontsize='xx-small')
    
    # annotate Galactic center
    plt.plot(0, 0, 'wx')
    plt.text(0,-1, 'Galactic\nCenter', ha='center', va='top', color=text_color, fontsize='small')
    
    # annotate the Sun
    xsun = -gc_frame.galcen_distance * np.cos(phi)
    zsun = gc_frame.z_sun.to(u.kpc)
    plt.plot(xsun, zsun, color='w', marker='$\odot$', ms=10, mew=0.4)
    plt.text(xsun.value+1, zsun.value-1, 'Sun', color=text_color, fontsize='small')
    
    #text_intro = '''Stellar streams orbiting the Milky Way
    #are gravitational antennae that have
    #picked up curious signals:'''
    #plt.text(0.2, 0.87, text_intro, color=text_color, wrap=True, transform=plt.gca().transAxes, va='top', ha='center')
    plt.xlim(-40,60)
    plt.ylim(-30,13)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().set_facecolor('k')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    #plt.ylabel('TODAY', color=text_color)
    #for spine in plt.gca().spines.values():
        #spine.set_edgecolor('0.9')
        #spine.set_linewidth(4)

    if step>-1:
        # show GD-1
        gd1_color = '#13d100'
        g = Table.read('/home/ana/projects/legacy/GD1-DR2/output/gd1_members.fits')

        ax00 = inset_axes(ax, '100%', '100%', bbox_to_anchor=[0.58, 0.73, 0.38, 0.15], loc=3, bbox_transform=ax.transAxes)
        plt.sca(ax00)
        plt.scatter(g['phi1'], g['phi2'], s=g['pmem']*1.5, c=g['pmem'], cmap=mpl.cm.binary, vmin=0.5, vmax=1.1, zorder=0, label='', rasterized=True)
        #plt.plot(g['phi1'], g['phi2'], 'k.', ms=1, alpha=0.5)
        plt.text(-40,-2, 'gap', fontsize='small', ha='center')
        plt.text(-32, 2, 'spur', fontsize='small', ha='center')
        
        plt.xlim(-60, -20)
        plt.ylim(-4,4)
        plt.gca().set_aspect('equal', adjustable='datalim')
        plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
        for spine in plt.gca().spines.values():
            spine.set_edgecolor(gd1_color)
            spine.set_linewidth(2)
        
        plt.text(0, 1.1, 'Price-Whelan & Bonaca (2018)', transform=plt.gca().transAxes, color=gd1_color, fontsize='medium')
        #plt.text(0, -0.1, 'Were the GD-1 gap and spur created by an impact of a compact\nobject, like a dark-matter subhalo?', va='top', transform=plt.gca().transAxes, fontsize='small', color=gd1_color)
        
        plt.sca(ax)
        plt.annotate('', xy=(0,1), xycoords=ax00.transAxes, xytext=(13.5, 5.5), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=gd1_color, alpha=0.3))
        plt.annotate('', xy=(0,0), xycoords=ax00.transAxes, xytext=(13.5, 5.5), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=gd1_color, alpha=0.3))
        
    if step>-1:
        # show Pal 5
        tp = Table.read('/home/ana/projects/legacy/Pal5sBiggestFan/data/pal5_likely_members.fits')
        pal5_color = '#00c7ff'
        
        ax01 = inset_axes(ax, '100%', '100%', bbox_to_anchor=[0.03, 0.52, 0.3, 0.31], loc=3, bbox_transform=ax.transAxes)
        plt.sca(ax01)
        plt.plot(tp['ra'], tp['dec'], 'k.', ms=0.5, alpha=0.1, rasterized=True)
        plt.text(223.5, -5, 'fan', fontsize='small')
        plt.text(228.5, 0.5, 'Pal 5', fontsize='small')
        
        plt.xlim(242,222)
        plt.ylim(-10,10)
        plt.gca().set_aspect('equal', adjustable='datalim')
        plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
        for spine in plt.gca().spines.values():
            spine.set_edgecolor(pal5_color)
            spine.set_linewidth(2)
        
        plt.text(0, 1.05, 'Bonaca et al. (2020a)', transform=plt.gca().transAxes, color=pal5_color, fontsize='medium')
        pal5_text = '''Was the leading tail of the Palomar 5 stream
broadened into a fan by a large perturbation, like
the rotating Galactic bar?
[see also Starkman et al. (2020), #poster-starkman]'''
        #plt.text(0, -0.1, pal5_text, va='top', transform=plt.gca().transAxes, fontsize='small', color=pal5_color)
        
        plt.sca(ax)
        plt.annotate('', xy=(1,1), xycoords=ax01.transAxes, xytext=(6, 2), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=pal5_color, alpha=0.3), zorder=0)
        plt.annotate('', xy=(1,0), xycoords=ax01.transAxes, xytext=(6, 2), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=pal5_color, alpha=0.3), zorder=0)
    
    if step>-1:
        # show Jhelum
        jhelum_color = '#ff9e00'
        tj = Table.read('/home/ana/projects/legacy/jhelum/data/jhelum_likely_members.fits')
        
        ax02 = inset_axes(ax, '100%', '100%', bbox_to_anchor=[0.03, 0.12, 0.34, 0.17], loc=3, bbox_transform=ax.transAxes)
        plt.sca(ax02)
        plt.plot(tj['phi1'], tj['phi2'], 'k.', ms=1.5, alpha=0.3, rasterized=True)
        plt.text(10, -3.5, 'wide component', fontsize='small', ha='center', bbox=dict(facecolor='w', alpha=0.4, edgecolor='none'))
        plt.text(-2, 2, 'narrow component', fontsize='small', bbox=dict(facecolor='w', alpha=0.4, edgecolor='none'))
        
        plt.ylim(-5,5)
        plt.xlim(-4, 24)
        plt.gca().set_aspect('equal', adjustable='datalim')
        plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
        for spine in plt.gca().spines.values():
            spine.set_edgecolor(jhelum_color)
            spine.set_linewidth(2)
        plt.text(0, 1.05, 'Bonaca et al. (2019b)', transform=plt.gca().transAxes, color=jhelum_color, fontsize='medium')
        jhelum_text='''Is Jhelum's wide envelope a signature of its nature set
in the early universe or has it been nurtured by
the Milky Way?
[see also Malhan et al. (2020), #poster-malhan]
'''
        #plt.text(0, -0.1, jhelum_text, va='top', transform=plt.gca().transAxes, fontsize='small', color=jhelum_color)
        
        plt.sca(ax)
        plt.annotate('', xy=(1,1), xycoords=ax02.transAxes, xytext=(-1, -10), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=jhelum_color, alpha=0.3), zorder=0)
        plt.annotate('', xy=(0,1), xycoords=ax02.transAxes, xytext=(-1, -10), textcoords=ax.transData, arrowprops=dict(arrowstyle='-', lw=1.5, color=jhelum_color, alpha=0.3), zorder=0)
    
    plt.savefig('../plots/gallery_conclusion_{:d}.png'.format(step), dpi=120, facecolor='k')


def animate():
    """"""
    d = pickle.load(open('../data/galactocentric_streams_orbits.pkl', 'rb'))
    names = np.array(get_names())
    r = np.array([d[x]['r'].value for x in names]) * u.kpc
    isort = np.argsort(r)
    
    Nstreams = len(names)
    Nind = 80
    Ntot = Nstreams * Nind
    
    x0_in, x1_in, y0_in, y1_in = -7.25, -3.75, 1, 3
    x0_fin, x1_fin, y0_fin, y1_fin = -88, 88, -60, 40
    
    dx0 = x0_fin - x0_in
    dx1 = x1_fin - x1_in
    dy0 = y0_fin - y0_in
    dy1 = y1_fin - y1_in
    
    rho_in, rho_out = 2*u.kpc, 50*u.kpc
    dr = (rho_out - rho_in)/Ntot
    
    dphi = np.pi*0.5/Ntot * u.radian
    
    text_color = '0.9'
    
    for i in range(Ntot):
        plt.close()
        plt.figure(figsize=(16,9))
        ax = plt.axes([0,0,1,1])
        
        ind_focus = int(i/Nind)
        
        for j in range(Nstreams):
            if ind_focus==j:
                ms = 6
                alpha = 0.5
            else:
                ms = 2
                alpha = 0.15
            
            ds = d[names[isort][j]]
            cgal = ds['stream']
            
            phi = dphi * i
            stream_x = cgal.x * np.cos(phi) + cgal.y * np.sin(phi)
            stream_y = - cgal.x * np.sin(phi) + cgal.y * np.cos(phi)
            
            ind_plot = stream_y < rho_in + i*dr
            plt.plot(stream_x[ind_plot], cgal.z[ind_plot], 'o', color=ds['color'], ms=ms)
            
            orbit = ds['orbit']
            orbit_x = orbit.x * np.cos(phi) + orbit.y * np.sin(phi)
            orbit_y = - orbit.x * np.sin(phi) + orbit.y * np.cos(phi)
            
            plt.plot(orbit_x, orbit.z, '-', color=ds['color'], alpha=alpha, lw=1, zorder=0)
        
        ## annotate Galactic center
        #plt.plot(0, 0, 'wx')
        #plt.text(0,-1, 'Galactic\nCenter', ha='center', va='top', color=text_color, fontsize='small')
        
        ## annotate the Sun
        #xsun = -gc_frame.galcen_distance * np.cos(phi)
        #zsun = gc_frame.z_sun.to(u.kpc)
        #plt.plot(xsun, zsun, color='w', marker='$\odot$', ms=10, mew=0.4)
        #plt.text(xsun.value+1, zsun.value-1, 'Sun', color=text_color, fontsize='small')
        
        
        x0 = x0_in + dx0 * i/Ntot
        x1 = x1_in + dx1 * i/Ntot
        y0 = y0_in + dy0 * i/Ntot
        y1 = y1_in + dy1 * i/Ntot
        
        # show scale
        scale_values = np.array([1, 2, 5, 10, 15])
        delta = 0.1*(x1 - x0)
        ind = np.argmin(np.abs(scale_values - delta))
        delta = scale_values[ind]
        
        xd0 = x0 + 0.05*(x1 - x0)
        xd1 = xd0 + delta
        xd = np.array([xd0, xd1])
        yd = np.ones(2)*(y0 + 0.05*(y1 - y0))
        
        plt.plot(xd, yd, '-', color=text_color, lw=2)
        plt.text(xd[0], y0 + 0.07*(y1 - y0), '{:d} kpc'.format(delta), fontsize='medium', color=text_color, va='bottom')
        
        plt.xlim(x0, x1)
        plt.ylim(y0, y1)
        plt.gca().tick_params(labelbottom=False, labelleft=False)
        plt.gca().set_facecolor('k')
        
        plt.savefig('../plots/stream_sample_wide/stream_sample.{:04d}.png'.format(i), dpi=120, facecolor='k', bbox_inches='tight', pad_inches=0)

