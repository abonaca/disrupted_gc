from orbit_fits import *

def save_galactocentric():
    """"""
    
    names = get_names()
    np.random.seed(2752)
    
    plt.close()
    plt.figure(figsize=(12,12))
    ax = plt.axes([0,0,1,1])
    
    outdict = {}
    
    for name in names[:]:
        stream = Stream(name, ham=ham, save_ext='')

        # best-fit orbit
        res = pickle.load(open('../data/fits/minimization_{:s}.pkl'.format(stream.savename), 'rb'))
        orbit = stream.orbital_properties(pbest=res.x, t=stream.tstream)

        # orbital distance
        orbit_eq = orbit.to_coord_frame(coord.ICRS)
        if stream.fra:
            model_x = orbit_eq.ra.wrap_at(stream.wangle)
            ix = 0
        else:
            model_x = orbit_eq.dec
            ix = -1

        pp = np.polyfit(model_x.value, orbit_eq.distance.value, 3)
        interp = np.poly1d(pp)
        dist = interp(stream.data['dec'][ix].value)*u.kpc

        # galactocentric coordinates
        c = coord.SkyCoord(ra=stream.data['dec'][0], dec=stream.data['dec'][1], distance=dist, frame='icrs')
        cgal = c.transform_to(coord.Galactocentric)
        
        # long orbit
        orbit = stream.orbital_properties(pbest=res.x, t=20*stream.tstream)
        lz = np.median(orbit.angular_momentum()[:,2])
        ecc = orbit.eccentricity()
        
        #color = mpl.cm.YlOrRd(0.8*np.random.rand())
        
        if lz<0:
            icolor = 0.5 - 0.5*ecc
            color = mpl.cm.PuOr(icolor)
        else:
            icolor = 0.5 + 0.5*ecc
            color = mpl.cm.BrBG(icolor)
        
        plt.plot(cgal.x, cgal.z, 'o', color=color, zorder=1, ms=2)
        plt.plot(orbit.pos.x, orbit.pos.z, '-', color=color, alpha=0.15, lw=1, zorder=0)
        
        xoff = -6*u.kpc
        zoff = 2*u.kpc
        outdict[name] = dict(stream=cgal, orbit=orbit.pos, color=color, r=np.median(np.sqrt((cgal.x-xoff)**2 + cgal.y**2 + (cgal.z-zoff)**2)))

    pickle.dump(outdict, open('../data/galactocentric_streams_orbits.pkl', 'wb'))
    
    plt.xlim(-50,50)
    plt.ylim(-60,40)
    plt.gca().set_facecolor('k')
    #plt.gca().set_aspect('equal')
    #plt.tight_layout()

def animate():
    """"""
    d = pickle.load(open('../data/galactocentric_streams_orbits.pkl', 'rb'))
    names = np.array(get_names())
    r = np.array([d[x]['r'].value for x in names]) * u.kpc
    isort = np.argsort(r)
    
    Nstreams = len(names)
    Nind = 80
    Ntot = Nstreams * Nind
    
    x0_in, x1_in, y0_in, y1_in = -7, -4, 0, 3
    x0_fin, x1_fin, y0_fin, y1_fin = -50, 50, -60, 40
    
    dx0 = x0_fin - x0_in
    dx1 = x1_fin - x1_in
    dy0 = y0_fin - y0_in
    dy1 = y1_fin - y1_in
    
    rho_in, rho_out = 2*u.kpc, 50*u.kpc
    dr = (rho_out - rho_in)/Ntot
    
    dphi = np.pi*0.5/Ntot * u.radian
    
    for i in range(Ntot):
        plt.close()
        plt.figure(figsize=(12,12))
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
        
        x0 = x0_in + dx0 * i/Ntot
        x1 = x1_in + dx1 * i/Ntot
        y0 = y0_in + dy0 * i/Ntot
        y1 = y1_in + dy1 * i/Ntot
        
        plt.xlim(x0, x1)
        plt.ylim(y0, y1)
        plt.gca().tick_params(labelbottom=False, labelleft=False)
        plt.gca().set_facecolor('k')
        
        plt.savefig('../plots/stream_sample/stream_sample.{:04d}.png'.format(i), dpi=90, facecolor='k')
    

#################
# proposal figure

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def save_galactocentric_long():
    """"""
    
    names = get_names()
    np.random.seed(2752)
    
    plt.close()
    plt.figure(figsize=(12,12))
    ax = plt.axes([0,0,1,1])
    
    outdict = {}
    
    for name in names[:]:
        stream = Stream(name, ham=ham, save_ext='')

        # best-fit orbit
        res = pickle.load(open('../data/fits/minimization_{:s}.pkl'.format(stream.savename), 'rb'))
        orbit = stream.orbital_properties(pbest=res.x, t=stream.tstream)

        # orbital distance
        orbit_eq = orbit.to_coord_frame(coord.ICRS)
        if stream.fra:
            model_x = orbit_eq.ra.wrap_at(stream.wangle)
            ix = 0
        else:
            model_x = orbit_eq.dec
            ix = -1

        pp = np.polyfit(model_x.value, orbit_eq.distance.value, 3)
        interp = np.poly1d(pp)
        dist = interp(stream.data['dec'][ix].value)*u.kpc

        # galactocentric coordinates
        c = coord.SkyCoord(ra=stream.data['dec'][0], dec=stream.data['dec'][1], distance=dist, frame='icrs')
        cgal = c.transform_to(coord.Galactocentric)
        
        # long orbit
        orbit = stream.orbital_properties(pbest=res.x, t=5*u.Gyr)
        lz = np.median(orbit.angular_momentum()[:,2])
        ecc = orbit.eccentricity()
        
        #color = mpl.cm.YlOrRd(0.8*np.random.rand())
        
        if lz<0:
            icolor = 0.5 - 0.5*ecc
            color = mpl.cm.PuOr(icolor)
        else:
            icolor = 0.5 + 0.5*ecc
            color = mpl.cm.BrBG(icolor)
        
        plt.plot(cgal.x, cgal.z, 'o', color=color, zorder=1, ms=2)
        plt.plot(orbit.pos.x, orbit.pos.z, '-', color=color, alpha=0.15, lw=1, zorder=0)
        
        xoff = -6*u.kpc
        zoff = 2*u.kpc
        outdict[name] = dict(stream=cgal, orbit=orbit.pos, color=color, r=np.median(np.sqrt((cgal.x-xoff)**2 + cgal.y**2 + (cgal.z-zoff)**2)))

    pickle.dump(outdict, open('../data/galactocentric_streams_orbits_long.pkl', 'wb'))
    
    plt.xlim(-50,50)
    plt.ylim(-60,40)
    plt.gca().set_facecolor('k')
    #plt.gca().set_aspect('equal')
    #plt.tight_layout()

def plot_long_orbits():
    """"""
    
    tsgr = Table.read('/home/ana/projects/h3/data/SgrTriax_DYN.dat.gz', format='ascii')
    tsgr = tsgr[::10]
    c_sgr = coord.ICRS(ra=tsgr['ra']*u.deg, dec=tsgr['dec']*u.deg, distance=tsgr['dist']*u.kpc)
    sgr_gal = c_sgr.transform_to(coord.Galactocentric)
    
    d = pickle.load(open('../data/galactocentric_streams_orbits_long.pkl', 'rb'))
    names = np.array(get_names())
    Nstreams = len(names)
    
    for i in range(Nstreams):
        ds = d[names[i]]
        cgal = ds['stream']
        orbit = ds['orbit']
        lw = 0.8

        plt.close()
        fig, ax = plt.subplots(1,2,figsize=(10,5))
        
        plt.sca(ax[0])
        plt.plot(sgr_gal.x, sgr_gal.z, 'ko', ms=1, label='Sgr')
        plt.plot(cgal.x, cgal.z, 'ro', label='')
        plt.plot(orbit.x, orbit.z, 'r-', lw=lw, label=names[i])
        plt.gca().set_aspect('equal', adjustable='datalim')
        plt.xlabel('x [kpc]')
        plt.ylabel('z [kpc]')
        plt.legend(fontsize='small', handlelength=0.5)
        
        plt.sca(ax[1])
        plt.plot(sgr_gal.y, sgr_gal.z, 'ko', ms=1)
        plt.plot(cgal.y, cgal.z, 'ro')
        plt.plot(orbit.y, orbit.z, 'r-', lw=lw)
        plt.gca().set_aspect('equal', adjustable='datalim')
        plt.xlabel('y [kpc]')
        plt.ylabel('z [kpc]')

        plt.tight_layout()
        plt.savefig('../plots/diag/sgr_{:s}.png'.format(names[i]))

from gala.units import galactic
def lmc_effect():
    """"""
    
    c_lmc = coord.ICRS(ra=78.76*u.deg, dec=-69.19*u.deg, distance=(10**(0.2*18.50+1)*u.pc).to(u.kpc), radial_velocity=262.2*u.km/u.s, pm_ra_cosdec=1.91*u.mas/u.yr, pm_dec=0.229*u.mas/u.yr)
    lmc_gal = c_lmc.transform_to(coord.Galactocentric)
    print(lmc_gal)
    
    name = 'elqui'
    stream = Stream(name, ham=ham, save_ext='')
    res = pickle.load(open('../data/fits/minimization_{:s}.pkl'.format(stream.savename), 'rb'))
    dec, d, pmra, pmdec, vr = res.x
    c_stream = coord.ICRS(ra=stream.ra0*u.deg, dec=dec*u.deg, distance=d*u.kpc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=vr*u.km/u.s)
    
    c = coord.concatenate([c_lmc, c_stream])
    w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)
    #print(w0)
    
    masses = np.array([1e10, 0]) * u.Msun
    b = 2.5*u.kpc
    potentials = [gp.PlummerPotential(x, b, units=galactic) for x in masses]

    nbody = gd.DirectNBody(w0, potentials, external_potential=ham.potential, units=galactic, save_all=True)
    orbits = nbody.integrate_orbit(dt=-1*u.Myr, n_steps=5000)
    
    plt.close()
    plt.figure()
    
    plt.plot(orbits.x, orbits.z, '-', lw=1)
    plt.plot(orbits.x[-1,:], orbits.z[-1, :], 'k.')
    
    masses = np.array([0, 0]) * u.Msun
    b = 2.5*u.kpc
    potentials = [gp.PlummerPotential(x, b, units=galactic) for x in masses]
    #print(w0)

    nbody0 = gd.DirectNBody(w0, potentials, external_potential=ham.potential, units=galactic, save_all=True)
    orbits0 = nbody0.integrate_orbit(dt=-1*u.Myr, n_steps=5000)
    
    plt.plot(orbits0.pos.x, orbits0.z, '--', lw=1)
    plt.plot(orbits0.x[-1,:], orbits0.z[-1, :], 'k.')
    
    outdict = dict(orbits=orbits.pos, orbits0=orbits0.pos)
    pickle.dump(outdict, open('../data/galactocentric_lmc_elqui.pkl', 'wb'))
    

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
    plt.scatter(g['phi1'], g['phi2'], s=g['pmem']*1.5, c=g['pmem'], cmap=mpl.cm.binary, vmin=0.5, vmax=1.1, zorder=0, label='')
    #plt.plot(g['phi1'], g['phi2'], 'k.', ms=1, alpha=0.5)
    
    plt.xlim(-60, -20)
    plt.ylim(-4,4)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    for spine in plt.gca().spines.values():
        spine.set_edgecolor(gd1_color)
        spine.set_linewidth(2)
    
    plt.text(0, 1.1, 'Price-Whelan & Bonaca (2018)', transform=plt.gca().transAxes, color=gd1_color, fontsize='small')
    plt.text(0, -0.1, 'Were the GD-1 gap and spur created by an impact\nof a compact object, like a dark-matter subhalo?', va='top', transform=plt.gca().transAxes, fontsize='small', color=gd1_color)
    
    # show Pal 5
    tp = Table.read('/home/ana/projects/legacy/Pal5sBiggestFan/data/pal5_likely_members.fits')
    pal5_color = '#00c7ff'
    
    ax01 = inset_axes(ax0, '100%', '100%', bbox_to_anchor=[0.62, 0.02, 0.2, 0.43], loc=3, bbox_transform=ax0.transAxes)
    plt.sca(ax01)
    plt.plot(tp['ra'], tp['dec'], 'k.', ms=0.5, alpha=0.1, rasterized=True)
    
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
    
    # show Jhelum
    jhelum_color = '#ff9e00'
    tj = Table.read('/home/ana/projects/legacy/jhelum/data/jhelum_likely_members.fits')
    
    ax02 = inset_axes(ax0, '100%', '100%', bbox_to_anchor=[0.02, 0.2, 0.3, 0.3], loc=3, bbox_transform=ax0.transAxes)
    plt.sca(ax02)
    plt.plot(tj['phi1'], tj['phi2'], 'k.', ms=1.5, alpha=0.3)
    
    plt.ylim(-5,5)
    plt.xlim(-4, 24)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    for spine in plt.gca().spines.values():
        spine.set_edgecolor(jhelum_color)
        spine.set_linewidth(2)
    plt.text(0, 1.05, 'Bonaca et al. (2019)', transform=plt.gca().transAxes, color=jhelum_color, fontsize='small')
    plt.text(0, -0.1, "Is Jhelum's wide envelope a signature of\nits nature set in the early universe or has\nit been nurtured by the Milky Way?", va='top', transform=plt.gca().transAxes, fontsize='small', color=jhelum_color)
    
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
    
    plt.ylim(-115,135)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().set_facecolor('k')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    plt.ylabel('THE  URF  PLAN', color=text_color)
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


def plot_gd1():
    """"""
    g = Table.read('/home/ana/projects/legacy/GD1-DR2/output/gd1_members.fits')

    #ax00 = inset_axes(ax0, '100%', '100%', bbox_to_anchor=[0.58, 0.69, 0.38, 0.2], loc=3, bbox_transform=ax0.transAxes)
    #plt.sca(ax00)
    plt.close()
    fig, ax = plt.subplots(1,1,figsize=(10,2), facecolor='k')
    
    plt.scatter(g['phi1'], g['phi2'], s=g['pmem']*1.5, c=g['pmem'], cmap=mpl.cm.binary, vmin=0.5, vmax=1.1, zorder=0, label='', rasterized=True)
    #plt.plot(g['phi1'], g['phi2'], 'k.', ms=1, alpha=0.5)
    plt.text(-40,-3, 'gap', fontsize='small', ha='center')
    plt.text(-32, 3, 'spur', fontsize='small', ha='center')
    
    plt.xlim(-80, 0)
    plt.ylim(-4,4)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    
    plt.savefig('../plots/fancy_gd1.png', facecolor=fig.get_facecolor())

def plot_jhelum():
    """"""
    tj = Table.read('/home/ana/projects/legacy/jhelum/data/jhelum_likely_members.fits')
    
    plt.close()
    fig, ax = plt.subplots(1,1,figsize=(10,4), facecolor='k')
    
    plt.plot(tj['phi1'], tj['phi2'], 'k.', ms=4, alpha=0.4, rasterized=True)
    plt.text(10, -3.5, 'wide component', fontsize='small', ha='center', bbox=dict(facecolor='w', alpha=0.6, edgecolor='none'))
    plt.text(-2, 2, 'narrow component', fontsize='small', bbox=dict(facecolor='w', alpha=0.6, edgecolor='none'))
    
    plt.ylim(-5,5)
    plt.xlim(-4, 24)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)

    plt.savefig('../plots/fancy_jhelum.png', facecolor=fig.get_facecolor())

def plot_pal5():
    """"""
    tp = Table.read('/home/ana/projects/legacy/Pal5sBiggestFan/data/pal5_likely_members.fits')
    pal5_color = '#00c7ff'
    
    plt.close()
    fig, ax = plt.subplots(1,1,figsize=(8,6), facecolor='k')
    
    plt.plot(tp['ra'], tp['dec'], 'k.', ms=0.5, alpha=0.3, rasterized=True)
    plt.text(223.5, -5, 'fan', fontsize='small', bbox=dict(facecolor='w', alpha=0.4, edgecolor='none'))
    plt.text(228.5, 0.5, 'Pal 5', fontsize='small', bbox=dict(facecolor='w', alpha=0.4, edgecolor='none'))
    
    plt.xlim(243,223)
    plt.ylim(-10,10)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.gca().tick_params(labelbottom=False, labelleft=False, top=False, bottom=False, right=False, left=False)
    
    plt.savefig('../plots/fancy_pal5.png', facecolor=fig.get_facecolor())
