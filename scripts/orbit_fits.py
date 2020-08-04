import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.table import Table, QTable, hstack, vstack
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import ascii

import gala.coordinates as gc
import gala.potential as gp
import gala.dynamics as gd
#import galstreams

from scipy.optimize import minimize
from scipy.interpolate import InterpolatedUnivariateSpline

import pickle

ham = gp.Hamiltonian(gp.MilkyWayPotential())
ham_bovy = gp.Hamiltonian(gp.BovyMWPotential2014())
ham_heavy = gp.Hamiltonian(gp.MilkyWayPotential(nucleus=dict(m=0), halo=dict(c=0.95, m=7e11), bulge=dict(m=4e9), disk=dict(m=5.5e10)))

coord.galactocentric_frame_defaults.set('v4.0')
gc_frame = coord.Galactocentric()

class GeneralStream():
    def __init__(self, name, label='', wangle=360*u.deg, ra0=np.nan*u.deg, tstream=50*u.Myr, dt=-0.5*u.Myr, vnorm=1., pmnorm=1., minra=True, dra=0.5*u.deg, ham=ham, gc_frame=gc_frame, save_ext=''):
        self.name = name
        if len(save_ext):
            self.savename = '{:s}_{:s}'.format(self.name, save_ext)
        else:
            self.savename = self.name
        if len(label):
            self.label = label
        else:
            self.label = self.name
        
        self.data = pickle.load(open('../data/streams/data_{:s}.pkl'.format(self.name), 'rb'))
        
        self.wangle = wangle
        if ~np.isfinite(ra0.value):
            self.ra0 = self.get_ra0(minra=min, dra=dra)
        else:
            self.ra0 = ra0
        
        self.dt = dt
        self.tstream = tstream
        self.nstep = int((self.tstream/np.abs(self.dt)).decompose())
        
        self.ham = ham
        self.gc_frame = gc_frame
            
    def get_ra0(self, minra=True, dra=0.5*u.deg):
        """Select min/max RA as the orbital fiducial point"""
        
        if minra:
            ra0 = np.min(data['phi2'][0]) - dra
        else:
            ra0 = np.max(data['phi2'][1]) + dra
        
        return ra0
    
    def rm_dataunits(self):
        """"""
        self.data_nounits = dict()
        for k in self.data.keys():
            self.data_nounits[k] = [x.value for x in self.data[k]]
    
    def orbit_minimize(self, p0=[], save=True):
        """Find best-fitting orbit by maximizing log likelihood"""
        
        if len(p0)==0:
            p0 = self.p0
        
        self.rm_dataunits()
        p0_input = [x_.value for x_ in p0]
        
        res = minimize(lambda *x: -ln_likelihood_icrs(*x), x0=p0_input, args=(self.ra0, self.data_nounits, self.nstep, self.dt, self.wangle, self.ham, self.gc_frame))
        self.pbest = res.x
        
        if save:
            pickle.dump(res, open('../data/fits/minimization_{:s}.pkl'.format(self.savename), 'wb'))
        
        return res
    
    def orbital_properties(self, pbest=[], t=5*u.Gyr, save=True):
        """"""
        if len(pbest)==0:
            pbest = self.pbest
        
        dec, d, pmra, pmdec, vr = pbest
        
        c = coord.ICRS(ra=self.ra0*u.deg, dec=dec*u.deg, distance=d*u.kpc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=vr*u.km/u.s)
        w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)
        
        n_long = int((t/np.abs(self.dt)).decompose())
        long_orbit = self.ham.integrate_orbit(w0, dt=self.dt, n_steps=n_long)
        
        if save:
            name = np.array([self.name])
            rperi = long_orbit.pericenter()
            rperi = np.array([rperi.value]) * rperi.unit
            rapo = long_orbit.apocenter()
            rapo = np.array([rapo.value]) * rapo.unit
            ecc = np.array([long_orbit.eccentricity()])
            vcirc = self.ham.potential.circular_velocity(np.array([long_orbit.apocenter().to(u.kpc).value, 0, 0]))
            vcirc = np.array([vcirc.value]) * vcirc.unit
            
            tout = Table([name, rperi, rapo, ecc, vcirc], names=('name', 'rperi', 'rapo', 'ecc', 'vcirc'))
            tout.write('../data/fits/minimization_orbit_{:s}.fits'.format(self.savename), overwrite=True)
        
        return long_orbit

class Stream(GeneralStream):
    def __init__(self, name, dt=-0.5*u.Myr ,ham=ham, gc_frame=gc_frame, save_ext=''):
        prop = get_properties(name)
        self._prop = prop
        
        self.name = name
        if len(save_ext):
            self.savename = '{:s}_{:s}'.format(self.name, save_ext)
        else:
            self.savename = self.name
        self.label = prop['label']
        self.data = pickle.load(open('../data/streams/data_{:s}.pkl'.format(self.name), 'rb'))
        
        self.wangle = prop['wangle']
        self.ra0 = prop['ra0'].value
        self.p0 = [prop[x] for x in ['dec0', 'd0', 'pmra0', 'pmdec0', 'vr0']]
        
        self.dt = dt
        self.tstream = prop['tstream']
        self.nstep = int((self.tstream/np.abs(self.dt)).decompose())
        
        self.ham = ham
        self.gc_frame = gc_frame

def ln_likelihood_icrs(p, ra_0, data, n_steps, dt, wangle, ham, gc_frame):
    # initial conditions at ra_0
    dec, d, pmra, pmdec, vr = p
    
    if (d<0) | (np.abs(vr)>500) | (dec<-90) | (dec>90):
        return -np.inf
    
    wdeg = wangle.to(u.deg).value
    
    c = coord.ICRS(ra=ra_0*u.deg, dec=dec*u.deg, distance=d*u.kpc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=vr*u.km/u.s)
    w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)
    
    if ham.energy(w0)>0:
        return -np.inf
    
    orbit = ham.integrate_orbit(w0, dt=dt, n_steps=n_steps)
    model_stream = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=gc_frame)
    model_x = model_stream.ra.wrap_at(wangle).degree
    if model_x[-1] < wdeg - 360:
        return -np.inf
    
    model_dec = model_stream.dec.degree
    model_dist = model_stream.distance.to(u.kpc).value
    model_pmra = model_stream.pm_ra_cosdec.to(u.mas/u.yr).value
    model_pmdec = model_stream.pm_dec.to(u.mas/u.yr).value
    model_vr = model_stream.radial_velocity.to(u.km/u.s).value

    ix = np.argsort(model_x)
    model_x = model_x[ix]
    
    # define interpolating functions
    order = 3
    bbox = [wdeg - 360, wdeg]
    
    interp = {}
    interp['dec'] = InterpolatedUnivariateSpline(model_x, model_dec[ix], k=order, bbox=bbox)
    interp['dist'] = InterpolatedUnivariateSpline(model_x, model_dist[ix], k=order, bbox=bbox)
    interp['pmra'] = InterpolatedUnivariateSpline(model_x, model_pmra[ix], k=order, bbox=bbox)
    interp['pmdec'] = InterpolatedUnivariateSpline(model_x, model_pmdec[ix], k=order, bbox=bbox)
    interp['vr'] = InterpolatedUnivariateSpline(model_x, model_vr[ix], k=order, bbox=bbox)
    
    # model smoothing
    isigma = {}
    isigma['dec'] = 0.2 # deg
    isigma['dist'] = 0.1 # kpc
    isigma['pmra'] = 0. # mas/yr
    isigma['pmdec'] = 0. # mas/yr
    isigma['vr'] = 1 # km/s
    
    chi2 = 0
    keys = data.keys()
    for k in keys:
        sigma = np.sqrt(isigma[k]**2 + data[k][2]**2)
        chi2 += np.sum(-(interp[k](data[k][0]) - data[k][1])**2 / sigma**2 - 2*np.log(sigma))
    
    return chi2


##################################
# Save individual stream data sets

def prep_ophiuchus():
    """Store dictionary with Ophiuchus data"""
    
    # read in data from Caldwell et al. (2020)
    tdata = Table.read('../data/streams/docs/temp2_oph_members.tab', format='ascii.commented_header')
    N = len(tdata)
    
    # uncertainties
    w0 = 0.05*u.deg
    d0 = 0.1*u.kpc
    #vr0 = 10*u.km/u.s

    w = np.ones(N) * w0
    derr = np.ones(N) * d0
    #verr = np.ones(N) * vr0
    verr = tdata['CZERR'] * u.km/u.s
    
    # construct the data dictionary
    data = dict()
    data['dec'] = (tdata['RA']*u.deg, tdata['DEC']*u.deg, w)
    data['dist'] = (tdata['RA']*u.deg, tdata['d']*u.kpc, derr)
    data['pmra'] = (tdata['RA']*u.deg, tdata['pmra']*u.mas/u.yr, tdata['pmra_error']*u.mas/u.yr)
    data['pmdec'] = (tdata['RA']*u.deg, tdata['pmdec']*u.mas/u.yr, tdata['pmdec_error']*u.mas/u.yr)
    data['vr'] = (tdata['RA']*u.deg, tdata['VELOCITY']*u.km/u.s, verr)
    
    pickle.dump(data, open('../data/streams/data_ophiuchus.pkl', 'wb'))

def prep_gd1():
    """Store dictionary with GD-1 data"""
    
    # track from PWB 2018
    track = Table.read('../data/streams/docs/dr2_stream_track.fits')
    track['phi1'] *= u.deg

    # Koposov et al. 
    kop_pm = ascii.read("""phi1 pm_phi1 pm_phi2 err
    -55.00 -13.60 -5.70 1.30
    -45.00 -13.10 -3.30 0.70
    -35.00 -12.20 -3.10 1.00
    -25.00 -12.60 -2.70 1.40
    -15.00 -10.80 -2.80 1.00""")

    kop_phi2 = ascii.read("""phi1 phi2 err
    -60.00 -0.64 0.15
    -56.00 -0.89 0.27
    -54.00 -0.45 0.15
    -48.00 -0.08 0.13
    -44.00 0.01 0.14
    -40.00 -0.00 0.09
    -36.00 0.04 0.10
    -34.00 0.06 0.13
    -32.00 0.04 0.06
    -30.00 0.08 0.10
    -28.00 0.03 0.12
    -24.00 0.06 0.05
    -22.00 0.06 0.13
    -18.00 -0.05 0.11
    -12.00 -0.29 0.16
    -2.00 -0.87 0.07""")

    kop_dist = ascii.read("""phi1 dist err
    -55.00 7.20 0.30
    -45.00 7.59 0.40
    -35.00 7.83 0.30
    -25.00 8.69 0.40
    -15.00 8.91 0.40
    0.00 9.86 0.50""")

    kop_vr = ascii.read("""phi1 phi2 vr err
    -45.23 -0.04 28.8 6.9
    -43.17 -0.09 29.3 10.2
    -39.54 -0.07 2.9  8.7
    -39.25 -0.22 -5.2 6.5
    -37.95 0.00 1.1   5.6
    -37.96 -0.00 -11.7 11.2
    -35.49 -0.05 -50.4 5.2
    -35.27 -0.02 -30.9 12.8
    -34.92 -0.15 -35.3 7.5
    -34.74 -0.08 -30.9 9.2
    -33.74 -0.18 -74.3 9.8
    -32.90 -0.15 -71.5 9.6
    -32.25 -0.17 -71.5 9.2
    -29.95 -0.00 -92.7 8.7
    -26.61 -0.11 -114.2 7.3
    -25.45 -0.14 -67.8 7.1
    -24.86 0.01 -111.2 17.8
    -21.21 -0.02 -144.4 10.5
    -14.47 -0.15 -179.0 10.0
    -13.73 -0.28 -191.4 7.5
    -13.02 -0.21 -162.9 9.6
    -12.68 -0.26 -217.2 10.7
    -12.55 -0.23 -172.2 6.6""")
    
    # LAMOST radial velocities
    tl = Table.read('../data/streams/docs/lamost_vr.fits')
    tl.keep_columns(['phi1', 'phi2', 'vr', 'err'])
    
    t = Table.read('../data/streams/docs/members_catalog.fits')
    ind = ((t['field']==1) | (t['field']==3) | (t['field']==7) | (t['field']==8)) & (t['std_Vrad']<5)
    tchelle = t[ind]
    
    tchelle.keep_columns(['phi1', 'phi2', 'Vrad', 'std_Vrad'])
    tchelle.rename_columns(['Vrad', 'std_Vrad'], ['vr', 'err'])
    
    tchelle['phi1'].unit = None
    tchelle['phi2'].unit = None
    tvr = vstack([kop_vr, tl, tchelle])
    
    # convert to equatorial coordinates
    ctrack = gc.GD1(phi1=track['phi1'], phi2=track['phi2'], pm_phi1_cosphi2=track['pm_phi1_cosphi2'], pm_phi2=track['pm_phi2'])
    ctrack_eq = ctrack.transform_to(coord.ICRS)
    
    cvr = gc.GD1(phi1=tvr['phi1'], phi2=tvr['phi2'])
    cvr_eq = cvr.transform_to(coord.ICRS)
    
    interp_track = InterpolatedUnivariateSpline(kop_phi2['phi1'], kop_phi2['phi2'], k=3, bbox=[-180,180])
    kop_dist_phi2 = interp_track(kop_dist['phi1'])
    ckop = gc.GD1(phi1=kop_dist['phi1']*u.deg, phi2=kop_dist_phi2*u.deg)
    ckop_eq = ckop.transform_to(coord.ICRS)

    # construct data dictionary
    data = dict()
    data['dec'] = (ctrack_eq.ra, ctrack_eq.dec, track['w'].quantity)
    data['dist'] = (ckop_eq.ra, kop_dist['dist']*u.kpc, kop_dist['err']*u.kpc)
    data['pmra'] = (ctrack_eq.ra, ctrack_eq.pm_ra_cosdec, track['pm_phi1_cosphi2_error'].quantity)
    data['pmdec'] = (ctrack_eq.ra, ctrack_eq.pm_dec, track['pm_phi2_error'].quantity)
    data['vr'] = (cvr_eq.ra, tvr['vr'].quantity, tvr['err'].quantity)
    
    pickle.dump(data, open('../data/streams/data_gd1.pkl', 'wb'))


def test_oph():
    """"""
    
    oph = Stream('ophiuchus')
    p0 = [-7.3*u.deg, 10*u.kpc, -4*u.mas/u.yr, -4.5*u.mas/u.yr, 270*u.km/u.s]
    res = oph.orbit_minimize(p0=p0, save=True)
    
    print(res.x)

def diag_oph():
    """"""
    
    oph = Stream('ophiuchus')
    res = pickle.load(open('../data/fits/minimization_ophiuchus.pkl', 'rb'))
    
    # find good initial guess
    p0 = [-7.3*u.deg, 10*u.kpc, -4*u.mas/u.yr, -4.5*u.mas/u.yr, 270*u.km/u.s]
    p0_fit = [x*y.unit for x,y in zip(res.x, p0) ]
    dec, dist, pmra, pmdec, vr = p0_fit

    c = coord.SkyCoord(ra=oph.ra0*u.deg, dec=dec, distance=dist, pm_ra_cosdec=pmra, pm_dec=pmdec, radial_velocity=vr, frame='icrs')
    w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)

    orbit = ham.integrate_orbit(w0, dt=oph.dt, n_steps=oph.nstep)
    model = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=gc_frame)
    
    # visualize data and the first attempt at an orbit
    plt.close()
    fig, ax = plt.subplots(5, 1, figsize=(7,11), sharex=True)

    fields = ['dec', 'dist', 'pmra', 'pmdec', 'vr']
    labels = ['Dec [deg]', 'Distance [kpc]', '$\mu_\\alpha$ [mas yr$^{-1}$]', '$\mu_\delta$ [mas yr$^{-1}$]',
            '$V_r$ [km s$^{-1}$]']
    model_fields = [model.dec, model.distance, model.pm_ra_cosdec, model.pm_dec, model.radial_velocity]
    istart, iend = 4, -4

    for i in range(5):
        plt.sca(ax[i])
        
        plt.plot(oph.data[fields[i]][0], oph.data[fields[i]][1], 'k.', label='Data')
        plt.errorbar(oph.data[fields[i]][0].value, oph.data[fields[i]][1].value, yerr=oph.data[fields[i]][2].value,
                    fmt='none', color='k', alpha=0.7, label='')
        
        plt.plot(model.ra[istart:iend], model_fields[i][istart:iend], '-', color='tab:blue', label='Best-fit orbit')
        
        plt.ylabel(labels[i])
        if i==0:
            plt.legend(loc=4, fontsize='small', handlelength=1)

    plt.ylim(210,360)
    plt.minorticks_on()
    plt.gca().invert_xaxis()
    plt.xlabel('R.A. [deg]')

    plt.tight_layout(h_pad=0)
    
def orbit_oph():
    """"""
    oph = Stream('ophiuchus')
    res = pickle.load(open('../data/fits/minimization_ophiuchus.pkl', 'rb'))
    pfit = res.x
    
    oph.orbital_properties(pfit=pfit)
    
    t = Table.read('../data/fits/minimization_orbit_ophiuchus.fits')
    t.pprint()



def initialize():
    """Construct a table with stream parameters and initial orbital guesses"""
    
    names = get_names()
    
    # initialize table
    t0 = get_properties(names[0])
    t0['name'] = names[0]
    
    name = np.array([t0['name']])
    label = np.array([t0['label']])
    wangle = np.array([t0['wangle'].value])*t0['wangle'].unit
    ra0 = np.array([t0['ra0'].value])*t0['ra0'].unit
    dec0 = np.array([t0['dec0'].value])*t0['dec0'].unit
    d0 = np.array([t0['d0'].value])*t0['d0'].unit
    pmra0 = np.array([t0['pmra0'].value])*t0['pmra0'].unit
    pmdec0 = np.array([t0['pmdec0'].value])*t0['pmdec0'].unit
    vr0 = np.array([t0['vr0'].value])*t0['vr0'].unit
    tstream = np.array([t0['tstream'].value])*t0['tstream'].unit
    
    tout = Table([name, label, wangle, ra0, dec0, d0, pmra0, pmdec0, vr0, tstream], names=('name', 'label', 'wangle', 'ra0', 'dec0', 'd0', 'pmra0', 'pmdec0', 'vr0', 'tstream'))
    
    # add subsequent rows
    for i in range(1, len(names)):
        t = get_properties(names[i])
        t['name'] = names[i]
        tout.add_row(t)
    
    tout.pprint()
    tout.write('../data/streams/initialize.fits', overwrite=True)

def get_names():
    """Get names of streams in the sample"""
    
    streams = ['ophiuchus', 'gd1']
    
    return sorted(streams)

def get_properties(name):
    """Return initial positions"""
    
    props = {}

    props['ophiuchus'] = dict(label='Ophiuchus', wangle=360*u.deg, ra0=240.5*u.deg, dec0=-7.3*u.deg, d0=10*u.kpc, pmra0=-4*u.mas/u.yr, pmdec0=-4.5*u.mas/u.yr, vr0=270*u.km/u.s, tstream=13*u.Myr)
    props['gd1'] = dict(label='GD-1', wangle=360*u.deg, ra0=123*u.deg, dec0=-10*u.deg, d0=9*u.kpc, pmra0=-2*u.mas/u.yr, pmdec0=-7*u.mas/u.yr, vr0=300*u.km/u.s, tstream=110*u.Myr)
    
    return props[name]

def test(name, dra=2, best=True):
    """"""
    stream = Stream(name)
    
    if best:
        res = pickle.load(open('../data/fits/minimization_{:s}.pkl'.format(stream.savename), 'rb'))
        p0 = [x*y.unit for x, y in zip(res.x, stream.p0)]
        dec, dist, pmra, pmdec, vr = p0
    else:
        dec, dist, pmra, pmdec, vr = stream.p0
        
    print(ham_heavy.potential.parameters)
    
    c = coord.SkyCoord(ra=stream.ra0*u.deg, dec=dec, distance=dist, pm_ra_cosdec=pmra, pm_dec=pmdec, radial_velocity=vr, frame='icrs')
    w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)

    orbit = stream.ham.integrate_orbit(w0, dt=stream.dt, n_steps=stream.nstep)
    model = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=stream.gc_frame)
    
    # plot data
    plt.close()
    fig, ax = plt.subplots(5, 1, figsize=(7,11), sharex=True)

    fields = ['dec', 'dist', 'pmra', 'pmdec', 'vr']
    labels = ['Dec [deg]', 'Distance [kpc]', '$\mu_\\alpha$ [mas yr$^{-1}$]', '$\mu_\delta$ [mas yr$^{-1}$]',
            '$V_r$ [km s$^{-1}$]']
    model_fields = [model.dec, model.distance, model.pm_ra_cosdec, model.pm_dec, model.radial_velocity]
    istart, iend = 0, -1

    for i in range(5):
        plt.sca(ax[i])
        
        if fields[i] in stream.data.keys():
            plt.plot(stream.data[fields[i]][0], stream.data[fields[i]][1], 'k.', label='Data')
            plt.errorbar(stream.data[fields[i]][0].value, stream.data[fields[i]][1].value, yerr=stream.data[fields[i]][2].value, fmt='none', color='k', alpha=0.7, label='')
            
        plt.plot(model.ra[istart:iend], model_fields[i][istart:iend], '-', color='tab:blue', label='Initial orbit')
        
        plt.ylabel(labels[i])
        if i==0:
            plt.legend(loc=4, fontsize='small', handlelength=1)

    plt.minorticks_on()
    plt.xlim(np.min(stream.data['dec'][0].to(u.deg).value)-dra, np.max(stream.data['dec'][0].to(u.deg).value)+dra)
    plt.xlabel('R.A. [deg]')

    plt.tight_layout(h_pad=0)
    if best:
        plt.savefig('../plots/diag/best_{:s}.png'.format(stream.name))
    
def fit_stream(name):
    """"""
    
    stream = Stream(name, ham=ham, save_ext='')
    res = stream.orbit_minimize(save=True)
    stream.orbital_properties(save=True)
    
    t = Table.read('../data/fits/minimization_orbit_{:s}.fits'.format(name))
    t.pprint()
    
    stream = Stream(name, ham=ham_bovy, save_ext='bovy')
    res = stream.orbit_minimize(save=True)
    stream.orbital_properties(save=True)
    
    t = Table.read('../data/fits/minimization_orbit_{:s}_bovy.fits'.format(name))
    t.pprint()
    
    stream = Stream(name, ham=ham_heavy, save_ext='heavy')
    res = stream.orbit_minimize(save=True)
    stream.orbital_properties(save=True)
    
    t = Table.read('../data/fits/minimization_orbit_{:s}_heavy.fits'.format(name))
    t.pprint()

    