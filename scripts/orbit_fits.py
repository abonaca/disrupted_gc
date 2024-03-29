import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.table import Table, QTable, hstack, vstack
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import ascii
from astropy.coordinates import frame_transform_graph
from astropy.coordinates.matrix_utilities import matrix_transpose

import gala.coordinates as gc
import gala.potential as gp
import gala.dynamics as gd
#import galstreams

from scipy.optimize import minimize
from scipy.interpolate import InterpolatedUnivariateSpline

import pickle
import emcee
from multiprocessing import Pool
import corner


ham = gp.Hamiltonian(gp.MilkyWayPotential())
ham_bovy = gp.Hamiltonian(gp.BovyMWPotential2014())
ham_heavy = gp.Hamiltonian(gp.MilkyWayPotential(nucleus=dict(m=0), halo=dict(c=0.95, m=7e11), bulge=dict(m=4e9), disk=dict(m=5.5e10)))

coord.galactocentric_frame_defaults.set('v4.0')
gc_frame = coord.Galactocentric()


class AAU(coord.BaseCoordinateFrame):
    """
    A Heliocentric spherical coordinate system defined by the track of the ATLAS Aliqa Uma stream

    For more information about this class, see the Astropy documentation
    on coordinate frames in :mod:`~astropy.coordinates`.

    Parameters
    ----------
    representation : :class:`~astropy.coordinates.BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)

    phi1 : angle_like, optional, must be keyword
        The longitude-like angle corresponding to AAU's orbit.
    phi2 : angle_like, optional, must be keyword
        The latitude-like angle corresponding to AAU's orbit.
    distance : :class:`~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.

    pm_phi1_cosphi2 : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in the longitude-like direction corresponding to
        the AAU stream's orbit.
    pm_phi2 : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in the latitude-like direction perpendicular to the
        AAU stream's orbit.
    radial_velocity : :class:`~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.

    """
    default_representation = coord.SphericalRepresentation
    default_differential = coord.SphericalCosLatDifferential

    frame_specific_representation_info = {
        coord.SphericalRepresentation: [
            coord.RepresentationMapping('lon', 'phi1'),
            coord.RepresentationMapping('lat', 'phi2'),
            coord.RepresentationMapping('distance', 'distance')],
        coord.SphericalCosLatDifferential: [
            coord.RepresentationMapping('d_lon_coslat', 'pm_phi1_cosphi2'),
            coord.RepresentationMapping('d_lat', 'pm_phi2'),
            coord.RepresentationMapping('d_distance', 'radial_velocity')],
        coord.SphericalDifferential: [
            coord.RepresentationMapping('d_lon', 'pm_phi1'),
            coord.RepresentationMapping('d_lat', 'pm_phi2'),
            coord.RepresentationMapping('d_distance', 'radial_velocity')]
    }

    frame_specific_representation_info[coord.UnitSphericalRepresentation] = \
        frame_specific_representation_info[coord.SphericalRepresentation]
    frame_specific_representation_info[coord.UnitSphericalCosLatDifferential] = \
        frame_specific_representation_info[coord.SphericalCosLatDifferential]
    frame_specific_representation_info[coord.UnitSphericalDifferential] = \
        frame_specific_representation_info[coord.SphericalDifferential]

@frame_transform_graph.transform(coord.StaticMatrixTransform, coord.ICRS, AAU)
def icrs_to_aau():
    """ Compute the transformation from Galactic spherical to
        heliocentric AAU coordinates.
    """
    rotmat = np.array([[0.83697865, 0.29481904, -0.4610298], [0.51616778, -0.70514011, 0.4861566], [0.18176238, 0.64487142, 0.74236331]])

    return rotmat

@frame_transform_graph.transform(coord.StaticMatrixTransform, AAU, coord.ICRS)
def aau_to_icrs():
    """ Compute the transformation from heliocentric AAU coordinates to
        spherical Galactic.
    """
    return matrix_transpose(icrs_to_aau())

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
        
        self.rm_dataunits()
            
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
        
        #self.rm_dataunits()
        p0_input = [x_.value for x_ in p0]
        
        res = minimize(lambda *x: -ln_likelihood_icrs(*x), x0=p0_input, args=(self.x0, self.data_nounits, self.nstep, self.dt, self.wangle, self.ham, self.gc_frame, self.fra))
        self.pbest = res.x
        
        if save:
            pickle.dump(res, open('../data/fits/minimization_{:s}.pkl'.format(self.savename), 'wb'))
        
        return res
    
    def orbital_properties(self, pbest=[], t=5*u.Gyr, save=True):
        """"""
        if len(pbest)==0:
            pbest = self.pbest
        
        if self.fra:
            ra = self.x0
            dec, d, pmra, pmdec, vr = pbest
        else:
            dec = self.x0
            ra, d, pmra, pmdec, vr = pbest
        
        c = coord.ICRS(ra=ra*u.deg, dec=dec*u.deg, distance=d*u.kpc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=vr*u.km/u.s)
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
        self.fra = prop['fra']
        self.ra0 = prop['ra0'].value
        self.dec0 = prop['dec0'].value
        if self.fra:
            self.x0 = prop['ra0'].value
            self.p0 = [prop[x] for x in ['dec0', 'd0', 'pmra0', 'pmdec0', 'vr0']]
        else:
            self.x0 = prop['dec0'].value
            self.p0 = [prop[x] for x in ['ra0', 'd0', 'pmra0', 'pmdec0', 'vr0']]
        
        self.dt = dt
        self.tstream = prop['tstream']
        self.nstep = int((self.tstream/np.abs(self.dt)).decompose())
        
        self.ham = ham
        self.gc_frame = gc_frame
        
        self.rm_dataunits()

def ln_likelihood_icrs(p, x_0, data, n_steps, dt, wangle, ham, gc_frame, fra):
    wdeg = wangle.to(u.deg).value
    # initial conditions at ra_0
    if fra:
        ra = x_0
        dec, d, pmra, pmdec, vr = p
    else:
        dec = x_0
        ra, d, pmra, pmdec, vr = p
    
    if (d<0) | (np.abs(vr)>500) | (dec<-90) | (dec>90) | (ra<wdeg-360) | (ra>wdeg):
        return -np.inf
    
    
    c = coord.ICRS(ra=ra*u.deg, dec=dec*u.deg, distance=d*u.kpc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=vr*u.km/u.s)
    w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)
    orbit = ham.integrate_orbit(w0, dt=dt, n_steps=n_steps)
    model_stream = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=gc_frame)
    
    model_ra = model_stream.ra.wrap_at(wangle).degree
    if model_ra[-1] < wdeg - 360:
        return -np.inf
    model_dec = model_stream.dec.degree
    
    if fra:
        model_x = model_ra
        model_y = model_dec
        indx = 0
        bbox = [wdeg - 360, wdeg]
    else:
        model_x = model_dec
        model_y = model_ra
        indx = -1
        bbox = [-90, 90]
        
        # switch data order
        data['dec'][1] = data['dec'][0]
    
    model_dist = model_stream.distance.to(u.kpc).value
    model_pmra = model_stream.pm_ra_cosdec.to(u.mas/u.yr).value
    model_pmdec = model_stream.pm_dec.to(u.mas/u.yr).value
    model_vr = model_stream.radial_velocity.to(u.km/u.s).value

    ix = np.argsort(model_x)
    model_x = model_x[ix]
    
    # define interpolating functions
    order = 3
    
    interp = {}
    interp['dec'] = InterpolatedUnivariateSpline(model_x, model_y[ix], k=order, bbox=bbox)
    interp['dist'] = InterpolatedUnivariateSpline(model_x, model_dist[ix], k=order, bbox=bbox)
    interp['pmra'] = InterpolatedUnivariateSpline(model_x, model_pmra[ix], k=order, bbox=bbox)
    interp['pmdec'] = InterpolatedUnivariateSpline(model_x, model_pmdec[ix], k=order, bbox=bbox)
    interp['vr'] = InterpolatedUnivariateSpline(model_x, model_vr[ix], k=order, bbox=bbox)
    
    # model smoothing
    isigma = {}
    isigma['dec'] = 0.01 # deg
    isigma['dist'] = 0.1 # kpc
    isigma['pmra'] = 0. # mas/yr
    isigma['pmdec'] = 0. # mas/yr
    isigma['vr'] = 1 # km/s
    
    #isigma['dec'] = 0. # deg
    #isigma['dist'] = 0. # kpc
    #isigma['pmra'] = 0. # mas/yr
    #isigma['pmdec'] = 0. # mas/yr
    #isigma['vr'] = 0. # km/s
    
    chi2 = 0
    keys = data.keys()
    for k in keys:
        sigma = np.sqrt(isigma[k]**2 + data[k][2]**2)
        chi2 += np.sum(-(interp[k](data[k][indx]) - data[k][1])**2 / sigma**2 - 2*np.log(sigma))
    
    return chi2

def ln_likelihood_icrs_mcmc(p, x_0, data, n_steps, dt, wangle, fra):
    wdeg = wangle.to(u.deg).value
    # initial conditions at ra_0
    if fra:
        ra = x_0
        dec, d, pmra, pmdec, vr = p
    else:
        dec = x_0
        ra, d, pmra, pmdec, vr = p
    
    if (d<0) | (np.abs(vr)>500) | (dec<-90) | (dec>90) | (ra<wdeg-360) | (ra>wdeg):
        return -np.inf
    
    
    c = coord.ICRS(ra=ra*u.deg, dec=dec*u.deg, distance=d*u.kpc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=vr*u.km/u.s)
    ## initial conditions at ra_0
    #dec, d, pmra, pmdec, vr = p
    
    #if (d<0) | (np.abs(vr)>500) | (dec<-90) | (dec>90):
        #return -np.inf
    
    #wdeg = wangle.to(u.deg).value
    
    #c = coord.ICRS(ra=ra_0*u.deg, dec=dec*u.deg, distance=d*u.kpc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=vr*u.km/u.s)
    w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)
    orbit = ham.integrate_orbit(w0, dt=dt, n_steps=n_steps)
    model_stream = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=gc_frame)
    
    model_ra = model_stream.ra.wrap_at(wangle).degree
    if model_ra[-1] < wdeg - 360:
        return -np.inf
    model_dec = model_stream.dec.degree
    
    if fra:
        model_x = model_ra
        model_y = model_dec
        indx = 0
        bbox = [wdeg - 360, wdeg]
    else:
        model_x = model_dec
        model_y = model_ra
        indx = -1
        bbox = [-90, 90]
        
        # switch data order
        data['dec'][1] = data['dec'][0]
    
    model_dist = model_stream.distance.to(u.kpc).value
    model_pmra = model_stream.pm_ra_cosdec.to(u.mas/u.yr).value
    model_pmdec = model_stream.pm_dec.to(u.mas/u.yr).value
    model_vr = model_stream.radial_velocity.to(u.km/u.s).value

    ix = np.argsort(model_x)
    model_x = model_x[ix]
    
    # define interpolating functions
    order = 3
    
    interp = {}
    interp['dec'] = InterpolatedUnivariateSpline(model_x, model_y[ix], k=order, bbox=bbox)
    interp['dist'] = InterpolatedUnivariateSpline(model_x, model_dist[ix], k=order, bbox=bbox)
    interp['pmra'] = InterpolatedUnivariateSpline(model_x, model_pmra[ix], k=order, bbox=bbox)
    interp['pmdec'] = InterpolatedUnivariateSpline(model_x, model_pmdec[ix], k=order, bbox=bbox)
    interp['vr'] = InterpolatedUnivariateSpline(model_x, model_vr[ix], k=order, bbox=bbox)
    
    # model smoothing
    isigma = {}
    isigma['dec'] = 0.01 # deg
    isigma['dist'] = 0.1 # kpc
    isigma['pmra'] = 0. # mas/yr
    isigma['pmdec'] = 0. # mas/yr
    isigma['vr'] = 1 # km/s
    
    #isigma['dec'] = 0. # deg
    #isigma['dist'] = 0. # kpc
    #isigma['pmra'] = 0. # mas/yr
    #isigma['pmdec'] = 0. # mas/yr
    #isigma['vr'] = 0. # km/s
    
    chi2 = 0
    keys = data.keys()
    for k in keys:
        sigma = np.sqrt(isigma[k]**2 + data[k][2]**2)
        chi2 += np.sum(-(interp[k](data[k][indx]) - data[k][1])**2 / sigma**2 - 2*np.log(sigma))
    
    return chi2

def ln_likelihood_icrs_mcmc_bovy(p, ra_0, data, n_steps, dt, wangle, fra):
    # initial conditions at ra_0
    dec, d, pmra, pmdec, vr = p
    
    if (d<0) | (np.abs(vr)>500) | (dec<-90) | (dec>90):
        return -np.inf
    
    wdeg = wangle.to(u.deg).value
    
    c = coord.ICRS(ra=ra_0*u.deg, dec=dec*u.deg, distance=d*u.kpc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=vr*u.km/u.s)
    w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)
    orbit = ham_bovy.integrate_orbit(w0, dt=dt, n_steps=n_steps)
    model_stream = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=gc_frame)
    
    model_ra = model_stream.ra.wrap_at(wangle).degree
    if model_ra[-1] < wdeg - 360:
        return -np.inf
    model_dec = model_stream.dec.degree
    
    if fra:
        model_x = model_ra
        model_y = model_dec
        indx = 0
        bbox = [wdeg - 360, wdeg]
    else:
        model_x = model_dec
        model_y = model_ra
        indx = -1
        bbox = [-90, 90]
        
        # switch data order
        data['dec'][1] = data['dec'][0]
    
    model_dist = model_stream.distance.to(u.kpc).value
    model_pmra = model_stream.pm_ra_cosdec.to(u.mas/u.yr).value
    model_pmdec = model_stream.pm_dec.to(u.mas/u.yr).value
    model_vr = model_stream.radial_velocity.to(u.km/u.s).value

    ix = np.argsort(model_x)
    model_x = model_x[ix]
    
    # define interpolating functions
    order = 3
    
    interp = {}
    interp['dec'] = InterpolatedUnivariateSpline(model_x, model_y[ix], k=order, bbox=bbox)
    interp['dist'] = InterpolatedUnivariateSpline(model_x, model_dist[ix], k=order, bbox=bbox)
    interp['pmra'] = InterpolatedUnivariateSpline(model_x, model_pmra[ix], k=order, bbox=bbox)
    interp['pmdec'] = InterpolatedUnivariateSpline(model_x, model_pmdec[ix], k=order, bbox=bbox)
    interp['vr'] = InterpolatedUnivariateSpline(model_x, model_vr[ix], k=order, bbox=bbox)
    
    # model smoothing
    isigma = {}
    isigma['dec'] = 0.01 # deg
    isigma['dist'] = 0.1 # kpc
    isigma['pmra'] = 0. # mas/yr
    isigma['pmdec'] = 0. # mas/yr
    isigma['vr'] = 1 # km/s
    
    #isigma['dec'] = 0. # deg
    #isigma['dist'] = 0. # kpc
    #isigma['pmra'] = 0. # mas/yr
    #isigma['pmdec'] = 0. # mas/yr
    #isigma['vr'] = 0. # km/s
    
    chi2 = 0
    keys = data.keys()
    for k in keys:
        sigma = np.sqrt(isigma[k]**2 + data[k][2]**2)
        chi2 += np.sum(-(interp[k](data[k][indx]) - data[k][1])**2 / sigma**2 - 2*np.log(sigma))
    
    return chi2

def ln_likelihood_icrs_mcmc_heavy(p, ra_0, data, n_steps, dt, wangle, fra):
    # initial conditions at ra_0
    dec, d, pmra, pmdec, vr = p
    
    if (d<0) | (np.abs(vr)>500) | (dec<-90) | (dec>90):
        return -np.inf
    
    wdeg = wangle.to(u.deg).value
    
    c = coord.ICRS(ra=ra_0*u.deg, dec=dec*u.deg, distance=d*u.kpc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=vr*u.km/u.s)
    w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)
    orbit = ham_heavy.integrate_orbit(w0, dt=dt, n_steps=n_steps)
    model_stream = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=gc_frame)
    
    model_ra = model_stream.ra.wrap_at(wangle).degree
    if model_ra[-1] < wdeg - 360:
        return -np.inf
    model_dec = model_stream.dec.degree
    
    if fra:
        model_x = model_ra
        model_y = model_dec
        indx = 0
        bbox = [wdeg - 360, wdeg]
    else:
        model_x = model_dec
        model_y = model_ra
        indx = -1
        bbox = [-90, 90]
        
        # switch data order
        data['dec'][1] = data['dec'][0]
    
    model_dist = model_stream.distance.to(u.kpc).value
    model_pmra = model_stream.pm_ra_cosdec.to(u.mas/u.yr).value
    model_pmdec = model_stream.pm_dec.to(u.mas/u.yr).value
    model_vr = model_stream.radial_velocity.to(u.km/u.s).value

    ix = np.argsort(model_x)
    model_x = model_x[ix]
    
    # define interpolating functions
    order = 3
    
    interp = {}
    interp['dec'] = InterpolatedUnivariateSpline(model_x, model_y[ix], k=order, bbox=bbox)
    interp['dist'] = InterpolatedUnivariateSpline(model_x, model_dist[ix], k=order, bbox=bbox)
    interp['pmra'] = InterpolatedUnivariateSpline(model_x, model_pmra[ix], k=order, bbox=bbox)
    interp['pmdec'] = InterpolatedUnivariateSpline(model_x, model_pmdec[ix], k=order, bbox=bbox)
    interp['vr'] = InterpolatedUnivariateSpline(model_x, model_vr[ix], k=order, bbox=bbox)
    
    # model smoothing
    isigma = {}
    isigma['dec'] = 0.01 # deg
    isigma['dist'] = 0.1 # kpc
    isigma['pmra'] = 0. # mas/yr
    isigma['pmdec'] = 0. # mas/yr
    isigma['vr'] = 1 # km/s
    
    #isigma['dec'] = 0. # deg
    #isigma['dist'] = 0. # kpc
    #isigma['pmra'] = 0. # mas/yr
    #isigma['pmdec'] = 0. # mas/yr
    #isigma['vr'] = 0. # km/s
    
    chi2 = 0
    keys = data.keys()
    for k in keys:
        sigma = np.sqrt(isigma[k]**2 + data[k][2]**2)
        chi2 += np.sum(-(interp[k](data[k][indx]) - data[k][1])**2 / sigma**2 - 2*np.log(sigma))
    
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
    data['dec'] = (coord.Longitude(tdata['RA']*u.deg), coord.Latitude(tdata['DEC']*u.deg), w)
    data['dist'] = (coord.Longitude(tdata['RA']*u.deg), tdata['d']*u.kpc, derr)
    data['pmra'] = (coord.Longitude(tdata['RA']*u.deg), tdata['pmra']*u.mas/u.yr, tdata['pmra_error']*u.mas/u.yr)
    data['pmdec'] = (coord.Longitude(tdata['RA']*u.deg), tdata['pmdec']*u.mas/u.yr, tdata['pmdec_error']*u.mas/u.yr)
    data['vr'] = (coord.Longitude(tdata['RA']*u.deg), tdata['VELOCITY']*u.km/u.s, verr)
    
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

def prep_svol():
    """"""
    t1 = Table.read('../data/streams/docs/svol_l_b.csv', format='ascii.no_header', delimiter=',')
    t2 = Table.read('../data/streams/docs/svol_pmra_b.csv', format='ascii.no_header', delimiter=',')
    t3 = Table.read('../data/streams/docs/svol_pmdec_b.csv', format='ascii.no_header', delimiter=',')
    
    tc = Table.read('../data/stream_endpoints_5d.fits')
    ind = tc['name']=='Sv\\"{o}l'
    ts = tc[ind]
    
    # typo in Riley table
    ceq_ = coord.SkyCoord(ra=ts['ra'], dec=ts['dec'], distance=ts['d'], frame='icrs')[0]
    cgal_ = ceq_.transform_to(coord.Galactic)
    cgal_end = coord.SkyCoord(l=cgal_.b.degree*u.deg, b=cgal_.l.degree*u.deg, distance=cgal_.distance, frame='galactic')
    ceq_end = cgal_end.transform_to(coord.ICRS)
    
    # uncertainties
    l_err = np.mean(ts['dec_err'])*ts['dec_err'].unit
    d_err = np.mean(ts['d_err'])*ts['d_err'].unit
    pm_err = np.mean(ts['pm_err'])*ts['pm_err'].unit
    #pm_err = 0.5*u.mas/u.yr
    
    # convert to equatorial coordinates
    c = coord.Galactic(l=t1['col2']*u.deg, b=t1['col1']*u.deg)
    c_eq = c.transform_to(coord.ICRS)
    
    np.random.seed(193)
    t1['col1'] += np.random.randn(len(t1))*1e-6
    isort = np.argsort(t1['col1'])
    interp_l = InterpolatedUnivariateSpline(t1['col1'][isort], t1['col2'][isort], k=3, bbox=[-90,90])
    
    l_pmra = interp_l(t2['col1'])
    cpmra = coord.Galactic(l=l_pmra*u.deg, b=t2['col1']*u.deg)
    cpmra_eq = cpmra.transform_to(coord.ICRS)
    
    l_pmdec = interp_l(t3['col1'])
    cpmdec = coord.Galactic(l=l_pmdec*u.deg, b=t3['col1']*u.deg)
    cpmdec_eq = cpmdec.transform_to(coord.ICRS)
    
    data = dict()
    data['dec'] = (c_eq.ra, c_eq.dec, np.ones(len(t1))*l_err)
    data['dist'] = (ceq_end.ra, ceq_end.distance, np.ones(np.size(ceq_end.ra))*d_err)
    data['pmra'] = (cpmra_eq.ra, t2['col2']*u.mas/u.yr, np.ones(len(t2))*pm_err)
    data['pmdec'] = (cpmdec_eq.ra, t3['col2']*u.mas/u.yr, np.ones(len(t3))*pm_err)
    
    #data['dec'] = (ceq_end.ra, ceq_end.dec, np.ones(np.size(ceq_end.ra))*l_err)
    #data['pmra'] = (ceq_end.ra, ts['pmra'].quantity[0], np.ones(np.size(ceq_end.ra))*pm_err)
    #data['pmdec'] = (ceq_end.ra, ts['pmdec'].quantity[0], np.ones(np.size(ceq_end.ra))*pm_err)
    
    pickle.dump(data, open('../data/streams/data_svol.pkl', 'wb'))

def reily_name(name):
    """Maps our simplified names to full Reily names"""
    
    names = dict(svol='Sv\\"{o}l', leiptr='Leiptr', gjoll='Gj\\"{o}ll', fjorm='Fj\\"{o}rm', fimbulthul='Fimbulthul', ylgr='Ylgr', sylgr='Sylgr', slidr='Slidr', phlegethon='Phlegethon', aliqa_uma='Aliqa Uma', atlas='ATLAS', elqui='Elqui', indus='Indus', phoenix='Phoenix', turranburra='Turranburra', jhelum='Jhelum', ravi='Ravi', turbio='Turbio', wambelong='Wambelong', willka_yaku='Willka Yaku')
    
    return names[name]

def prep_ibata_l(name, graph=False):
    """"""
    t1 = Table.read('../data/streams/docs/{:s}_b_l.csv'.format(name), format='ascii.no_header', delimiter=',')
    t2 = Table.read('../data/streams/docs/{:s}_pmra_l.csv'.format(name), format='ascii.no_header', delimiter=',')
    t3 = Table.read('../data/streams/docs/{:s}_pmdec_l.csv'.format(name), format='ascii.no_header', delimiter=',')
    
    tc = Table.read('../data/stream_endpoints_5d.fits')
    ind = tc['name']==reily_name(name)
    ts = tc[ind]
    ceq_end = coord.SkyCoord(ra=ts['ra'], dec=ts['dec'], distance=ts['d'], frame='icrs')[0]

    # mean uncertainties
    l_err = np.nanmean(ts['dec_err'])*ts['dec_err'].unit
    d_err = np.nanmean(ts['d_err'])*ts['d_err'].unit
    pm_err = np.nanmean(ts['pm_err'])*ts['pm_err'].unit
    
    # galactic latitude interpolation
    if name=='fjorm':
        k = 6
    else:
        k = 3
    pp = np.polyfit(t1['col1'], t1['col2'], k)
    interp_b = np.poly1d(pp)
    
    if graph:
        plt.close()
        plt.plot(t1['col1'], t1['col2'], 'ko')
        x = np.linspace(np.min(t1['col1']), np.max(t1['col1']), 100)
        y = interp_b(x)
        plt.plot(x, y, 'r-')
    
    # convert to equatorial coordinates
    c = coord.Galactic(l=t1['col1']*u.deg, b=t1['col2']*u.deg)
    c_eq = c.transform_to(coord.ICRS)
    
    b_pmra = interp_b(t2['col1'])
    cpmra = coord.Galactic(l=t2['col1']*u.deg, b=b_pmra*u.deg)
    cpmra_eq = cpmra.transform_to(coord.ICRS)
    
    b_pmdec = interp_b(t3['col1'])
    cpmdec = coord.Galactic(l=t3['col1']*u.deg, b=b_pmdec*u.deg)
    cpmdec_eq = cpmdec.transform_to(coord.ICRS)
    
    data = dict()
    data['dec'] = [c_eq.ra, c_eq.dec, np.ones(len(t1))*l_err, c_eq.dec]
    data['dist'] = [ceq_end.ra, ceq_end.distance, np.ones(np.size(ceq_end.ra))*d_err, ceq_end.dec]
    data['pmra'] = [cpmra_eq.ra, t2['col2']*u.mas/u.yr, np.ones(len(t2))*pm_err, cpmra_eq.dec]
    data['pmdec'] = [cpmdec_eq.ra, t3['col2']*u.mas/u.yr, np.ones(len(t3))*pm_err, cpmdec_eq.dec]
    
    pickle.dump(data, open('../data/streams/data_{:s}.pkl'.format(name), 'wb'))

def get_gaia_shipp():
    """Download Gaia DR2 data for likely members of DES streams"""
    t = Table.read('../data/streams/docs/shipp_members.txt', format='ascii.commented_header', delimiter=',')
    t.pprint()
    
    ids = '(' + ' '.join(list('{:d},'.format(x) for x in t['Gaia_Source_ID']))[:-1] + ')'
    q_base ='''SELECT * FROM gaiadr2.gaia_source WHERE source_id IN {:s}'''.format(ids)
    print(q_base)

def prep_shipp_members(name):
    """"""
    ti = Table.read('../data/streams/docs/shipp_members.txt', format='ascii.commented_header', delimiter=',')
    isort = np.argsort(ti['Gaia_Source_ID'])
    ti = ti[isort]
    t = Table.read('../data/streams/docs/shipp_members_gdr2.gz')
    
    if name=='jhelum':
        ind = (ti['Stream']=='Jhelum-a') | (ti['Stream']=='Jhelum-b')
    else:
        ind = ti['Stream']==reily_name(name)
    t = t[ind]
    
    props = get_properties(name)
    wangle = props['wangle']
    
    tc = Table.read('../data/stream_endpoints_5d.fits')
    ind = tc['name']==reily_name(name)
    ts = tc[ind]
    ceq_end = coord.SkyCoord(ra=ts['ra'], dec=ts['dec'], distance=ts['d'], frame='icrs')[0]

    # mean uncertainties
    l_err = np.nanmean(ts['dec_err'])*ts['dec_err'].unit
    d_err = np.nanmean(ts['d_err'])*ts['d_err'].unit
    pm_err = np.nanmean(ts['pm_err'])*ts['pm_err'].unit
    
    data = dict()
    data['dec'] = (coord.Longitude(t['ra'].quantity).wrap_at(wangle), t['dec'].quantity, np.ones(len(t))*l_err)
    data['dist'] = (ceq_end.ra.wrap_at(wangle), ceq_end.distance, np.ones(np.size(ceq_end.ra))*d_err)
    data['pmra'] = (coord.Longitude(t['ra'].quantity).wrap_at(wangle), t['pmra'].quantity, t['pmra_error'].quantity)
    data['pmdec'] = (coord.Longitude(t['ra'].quantity).wrap_at(wangle), t['pmdec'].quantity, t['pmdec_error'].quantity)
    
    pickle.dump(data, open('../data/streams/data_{:s}.pkl'.format(name), 'wb'))

def prep_aau(name):
    """"""
    ti = Table.read('../data/streams/docs/shipp_members.txt', format='ascii.commented_header', delimiter=',')
    isort = np.argsort(ti['Gaia_Source_ID'])
    ti = ti[isort]
    t = Table.read('../data/streams/docs/shipp_members_gdr2.gz')
    
    ind = ti['Stream']==reily_name(name)
    t = t[ind]
    
    props = get_properties(name)
    wangle = props['wangle']
    
    # endpoints
    tc = Table.read('../data/stream_endpoints_5d.fits')
    ind = tc['name']==reily_name(name)
    ts = tc[ind]
    ceq_end = coord.SkyCoord(ra=ts['ra'], dec=ts['dec'], distance=ts['d'], frame='icrs')[0]

    # mean uncertainties
    l_err = np.nanmean(ts['dec_err'])*ts['dec_err'].unit
    d_err = np.nanmean(ts['d_err'])*ts['d_err'].unit
    pm_err = np.nanmean(ts['pm_err'])*ts['pm_err'].unit
    vr_err = 4.8*u.km/u.s # velocity dispersion in Li et al.
    
    # radial velocities from Li et al. (2020)
    ceq = coord.SkyCoord(ra=t['ra'], dec=t['dec'], frame='icrs')
    caau = ceq.transform_to(AAU)
    phi10 = caau.phi1.wrap_at(wangle).deg/10
    vgsr = (-131.33 + 0.07*phi10 + 5.68*phi10**2)*u.km/u.s
    vr = gc.vgsr_to_vhel(ceq, vgsr)
    
    data = dict()
    data['dec'] = (coord.Longitude(t['ra'].quantity).wrap_at(wangle), coord.Latitude(t['dec'].quantity), np.ones(len(t))*l_err)
    data['dist'] = (ceq_end.ra.wrap_at(wangle), ceq_end.distance, np.ones(np.size(ceq_end.ra))*d_err)
    data['pmra'] = (coord.Longitude(t['ra'].quantity).wrap_at(wangle), t['pmra'].quantity, t['pmra_error'].quantity)
    data['pmdec'] = (coord.Longitude(t['ra'].quantity).wrap_at(wangle), t['pmdec'].quantity, t['pmdec_error'].quantity)
    data['vr'] = (coord.Longitude(t['ra'].quantity).wrap_at(wangle), vr, np.ones(len(t))*vr_err)
    
    pickle.dump(data, open('../data/streams/data_{:s}.pkl'.format(name), 'wb'))

def prep_shipp_avg(name, N=20):
    """"""
    
    # endpoints
    tc = Table.read('../data/stream_endpoints_5d.fits')
    ind = tc['name']==reily_name(name)
    t = tc[ind]
    
    props = get_properties(name)
    wangle = props['wangle']
    
    # find stream reference frame
    ceq_end = coord.SkyCoord(ra=t['ra'], dec=t['dec'], distance=t['d'], frame='icrs')[0]
    stream_pole = gc.pole_from_endpoints(ceq_end[0], ceq_end[1])
    stream_frame = gc.GreatCircleICRSFrame(pole=stream_pole)
    cend = ceq_end.transform_to(stream_frame)
    
    # distribute additional points along the great circle
    phi1 = np.linspace(np.min(cend.phi1), np.max(cend.phi1), N)
    phi2 = np.zeros(N)*u.deg
    call = coord.SkyCoord(phi1=phi1, phi2=phi2, frame=stream_frame)
    ceq_all = call.transform_to(coord.ICRS)
    
    dec_err = np.ones(N) * np.median(t['dec_err']) * t['dec_err'].unit
    d = np.ones(N) * np.median(t['d']) * t['d'].unit
    d_err = np.ones(N) * np.median(t['d_err']) * t['d_err'].unit
    pmra = np.ones(N) * np.median(t['pmra']) * t['pmra'].unit
    pmdec = np.ones(N) * np.median(t['pmdec']) * t['pmdec'].unit
    pm_err = np.ones(N) * np.median(t['pm_err']) * t['pm_err'].unit
    
    data = dict()
    data['dec'] = [ceq_all.ra.wrap_at(wangle), ceq_all.dec, dec_err, ceq_all.dec]
    data['dist'] = [ceq_all.ra.wrap_at(wangle), d, d_err, ceq_all.dec]
    data['pmra'] = [ceq_all.ra.wrap_at(wangle), pmra, pm_err, ceq_all.dec]
    data['pmdec'] = [ceq_all.ra.wrap_at(wangle), pmdec, pm_err, ceq_all.dec]
    
    pickle.dump(data, open('../data/streams/data_{:s}.pkl'.format(name), 'wb'))

def prep_s300():
    """Store dictionary with 300km/s Stream data"""
    
    # read in data (Gaia x Fu)
    ta, ts = pickle.load(open('../data/streams/docs/300S_members.pickle','rb'))
    
    ra = coord.Longitude(np.concatenate([ta['R.A.'], ts['R.A.']]) * u.deg)
    dec = coord.Latitude(np.concatenate([ta['Decl.'], ts['Decl.']]) * u.deg)
    N = np.size(ra)
    
    w0 = 1*u.deg
    w = np.ones(N) * w0
    
    d0 = 15*u.kpc
    derr0 = 3*u.kpc
    dist = np.ones(N) * d0
    dist_err = np.ones(N) * derr0
    
    pmra = np.concatenate([ta['pmra'], ts['pmra']]) * u.mas/u.yr
    pmra_err = np.concatenate([ta['pmra_error'], ts['pmra_error']]) * u.mas/u.yr
    
    pmdec = np.concatenate([ta['pmdec'], ts['pmdec']]) * u.mas/u.yr
    pmdec_err = np.concatenate([ta['pmdec_error'], ts['pmdec_error']]) * u.mas/u.yr
    
    vr = np.concatenate([ta['V_helio'], ts['V_helio']]) * u.km/u.s
    vr_err = np.concatenate([ta['sig_vhelio'], ts['sig_vhelio']]) * u.km/u.s
    
    # construct the data dictionary
    data = dict()
    data['dec'] = (ra, dec, w, dec)
    data['dist'] = (ra, dist, dist_err, dec)
    data['pmra'] = (ra, pmra, pmra_err, dec)
    data['pmdec'] = (ra, pmdec, pmdec_err, dec)
    data['vr'] = (ra, vr, vr_err, dec)
    
    pickle.dump(data, open('../data/streams/data_s300.pkl', 'wb'))


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


######################
# Find best-fit orbits

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
    
    streams = ['turbio', 'phoenix', 'elqui', 'willka_yaku', 'turranburra', 'slidr', 'svol', 'sylgr', 'ylgr']
    streams = ['ophiuchus', 'gd1', 'svol', 'leiptr', 'gjoll', 'fjorm', 'fimbulthul', 'ylgr', 'sylgr', 'slidr', 'phlegethon', 'phoenix', 'turranburra', 'indus', 'elqui', 'jhelum', 'atlas', 'aliqa_uma', 'ravi', 'wambelong', 'willka_yaku', 'turbio', 'triangulum', 's300']
    
    return sorted(streams)

def get_properties(name):
    """Return initial positions"""
    
    props = {}

    props['ophiuchus'] = dict(label='Ophiuchus', wangle=360*u.deg, ra0=240.5*u.deg, dec0=-7.3*u.deg, d0=10*u.kpc, pmra0=-4*u.mas/u.yr, pmdec0=-4.5*u.mas/u.yr, vr0=270*u.km/u.s, tstream=13*u.Myr, fra=True, provenance=[7,7,7])
    
    props['gd1'] = dict(label='GD-1', wangle=360*u.deg, ra0=123*u.deg, dec0=-10*u.deg, d0=9*u.kpc, pmra0=-2*u.mas/u.yr, pmdec0=-7*u.mas/u.yr, vr0=300*u.km/u.s, tstream=110*u.Myr, fra=True, provenance=[5,5,6])
    
    props['svol'] = dict(label='Sv\\"{o}l', wangle=360*u.deg, ra0=250*u.deg, dec0=25*u.deg, d0=8*u.kpc, pmra0=3.5*u.mas/u.yr, pmdec0=-6*u.mas/u.yr, vr0=-150*u.km/u.s, tstream=30*u.Myr, fra=True, provenance=[4,2,np.nan])
    
    props['leiptr'] = dict(label='Leiptr', wangle=360*u.deg, ra0=98*u.deg, dec0=-35*u.deg, d0=8*u.kpc, pmra0=10*u.mas/u.yr, pmdec0=-8*u.mas/u.yr, vr0=250*u.km/u.s, tstream=30*u.Myr, fra=True, provenance=[4,2,np.nan])

    props['gjoll'] = dict(label='Gj\\"{o}ll', wangle=360*u.deg, ra0=90*u.deg, dec0=-21*u.deg, d0=3.5*u.kpc, pmra0=24*u.mas/u.yr, pmdec0=-22*u.mas/u.yr, vr0=150*u.km/u.s, tstream=13*u.Myr, fra=True, provenance=[4,2,np.nan])
    
    props['fjorm'] = dict(label='Fj\\"{o}rm', wangle=360*u.deg, ra0=260*u.deg, dec0=70*u.deg, d0=5*u.kpc, pmra0=6*u.mas/u.yr, pmdec0=3*u.mas/u.yr, vr0=-100*u.km/u.s, tstream=30*u.Myr, fra=True, provenance=[4,2,np.nan])
    
    props['fimbulthul'] = dict(label='Fimbulthul', wangle=360*u.deg, ra0=198*u.deg, dec0=-32*u.deg, d0=4*u.kpc, pmra0=-9*u.mas/u.yr, pmdec0=-9*u.mas/u.yr, vr0=250*u.km/u.s, tstream=20*u.Myr, fra=True, provenance=[4,2,np.nan])
    props['fimbulthul'] = dict(label='Fimbulthul', wangle=360*u.deg, ra0=198*u.deg, dec0=-31*u.deg, d0=5*u.kpc, pmra0=-9*u.mas/u.yr, pmdec0=-7.5*u.mas/u.yr, vr0=400*u.km/u.s, tstream=8*u.Myr, fra=True, provenance=[4,2,np.nan])
    
    props['ylgr'] = dict(label='Ylgr', wangle=360*u.deg, ra0=183*u.deg, dec0=-38*u.deg, d0=9*u.kpc, pmra0=-0.5*u.mas/u.yr, pmdec0=-5*u.mas/u.yr, vr0=320*u.km/u.s, tstream=30*u.Myr, fra=True, provenance=[4,2,np.nan])

    props['sylgr'] = dict(label='Sylgr', wangle=360*u.deg, ra0=164*u.deg, dec0=-13*u.deg, d0=4*u.kpc, pmra0=-25*u.mas/u.yr, pmdec0=-22*u.mas/u.yr, vr0=-200*u.km/u.s, tstream=15*u.Myr, fra=True, provenance=[4,2,np.nan])
    
    props['slidr'] = dict(label='Slidr', wangle=360*u.deg, ra0=148*u.deg, dec0=17*u.deg, d0=3.5*u.kpc, pmra0=-28*u.mas/u.yr, pmdec0=-10*u.mas/u.yr, vr0=-50*u.km/u.s, tstream=20*u.Myr, fra=True, provenance=[4,2,np.nan])

    props['phlegethon'] = dict(label='Phlegethon', wangle=360*u.deg, ra0=299*u.deg, dec0=-61*u.deg, d0=3.5*u.kpc, pmra0=-12*u.mas/u.yr, pmdec0=-25*u.mas/u.yr, vr0=250*u.km/u.s, tstream=60*u.Myr, fra=False, provenance=[8,2,np.nan])
    
    props['phoenix'] = dict(label='Phoenix', wangle=360*u.deg, ra0=28*u.deg, dec0=-42*u.deg, d0=16*u.kpc, pmra0=2.6*u.mas/u.yr, pmdec0=-0.2*u.mas/u.yr, vr0=50*u.km/u.s, tstream=30*u.Myr, fra=True, provenance=[1,2,np.nan])
    
    props['turranburra'] = dict(label='Turranburra', wangle=360*u.deg, ra0=56*u.deg, dec0=-16*u.deg, d0=21*u.kpc, pmra0=0.6*u.mas/u.yr, pmdec0=-1.1*u.mas/u.yr, vr0=-50*u.km/u.s, tstream=140*u.Myr, fra=True, provenance=[1,2,np.nan])
    
    props['indus'] = dict(label='Indus', wangle=180*u.deg, ra0=4*u.deg, dec0=-66*u.deg, d0=14*u.kpc, pmra0=5.*u.mas/u.yr, pmdec0=-3*u.mas/u.yr, vr0=20*u.km/u.s, tstream=40*u.Myr, fra=True, provenance=[1,2,np.nan])

    props['elqui'] = dict(label='Elqui', wangle=360*u.deg, ra0=7*u.deg, dec0=-33*u.deg, d0=40*u.kpc, pmra0=0.1*u.mas/u.yr, pmdec0=-0.5*u.mas/u.yr, vr0=-150*u.km/u.s, tstream=150*u.Myr, fra=True, provenance=[1,2,np.nan])
    
    props['jhelum'] = dict(label='Jhelum', wangle=180*u.deg, ra0=5.5*u.deg, dec0=-52*u.deg, d0=10*u.kpc, pmra0=8*u.mas/u.yr, pmdec0=-3*u.mas/u.yr, vr0=-50*u.km/u.s, tstream=30*u.Myr, fra=True, provenance=[1,2,np.nan])

    props['atlas'] = dict(label='ATLAS', wangle=180*u.deg, ra0=9*u.deg, dec0=-20*u.deg, d0=18*u.kpc, pmra0=-0.5*u.mas/u.yr, pmdec0=-1*u.mas/u.yr, vr0=-150*u.km/u.s, tstream=60*u.Myr, fra=True, provenance=[1,2,3])
    
    props['aliqa_uma'] = dict(label='Aliqa Uma', wangle=180*u.deg, ra0=31*u.deg, dec0=-32*u.deg, d0=26*u.kpc, pmra0=0.25*u.mas/u.yr, pmdec0=-0.7*u.mas/u.yr, vr0=-60*u.km/u.s, tstream=40*u.Myr, fra=True, provenance=[1,2,3])
    
    props['ravi'] = dict(label='Ravi', wangle=360*u.deg, ra0=344.1*u.deg, dec0=-59*u.deg, d0=25*u.kpc, pmra0=0.9*u.mas/u.yr, pmdec0=-2.5*u.mas/u.yr, vr0=100*u.km/u.s, tstream=80*u.Myr, fra=False, provenance=[1,2,np.nan])
    
    props['turbio'] = dict(label='Turbio', wangle=360*u.deg, ra0=28.6*u.deg, dec0=-39*u.deg, d0=16*u.kpc, pmra0=2.*u.mas/u.yr, pmdec0=2*u.mas/u.yr, vr0=100*u.km/u.s, tstream=40*u.Myr, fra=False, provenance=[1,2,np.nan])
    
    props['wambelong'] = dict(label='Wambelong', wangle=360*u.deg, ra0=92.5*u.deg, dec0=-46*u.deg, d0=19*u.kpc, pmra0=1.8*u.mas/u.yr, pmdec0=-1.5*u.mas/u.yr, vr0=200*u.km/u.s, tstream=80*u.Myr, fra=True, provenance=[1,2,np.nan])
    props['wambelong'] = dict(label='Wambelong', wangle=360*u.deg, ra0=92.5*u.deg, dec0=-46*u.deg, d0=17*u.kpc, pmra0=2.5*u.mas/u.yr, pmdec0=-2*u.mas/u.yr, vr0=150*u.km/u.s, tstream=80*u.Myr, fra=True, provenance=[1,2,np.nan])
    #props['wambelong'] = dict(label='Wambelong', wangle=360*u.deg, ra0=78.5*u.deg, dec0=-33*u.deg, d0=17*u.kpc, pmra0=2*u.mas/u.yr, pmdec0=-2.5*u.mas/u.yr, vr0=-50*u.km/u.s, tstream=100*u.Myr, fra=True, provenance=[1,2,np.nan])
    
    props['willka_yaku'] = dict(label='Willka Yaku', wangle=360*u.deg, ra0=41*u.deg, dec0=-51*u.deg, d0=35*u.kpc, pmra0=1.1*u.mas/u.yr, pmdec0=0.1*u.mas/u.yr, vr0=100*u.km/u.s, tstream=100*u.Myr, fra=False, provenance=[1,2,np.nan])

    props['triangulum'] = dict(label='Triangulum', wangle=360*u.deg, ra0=21.2*u.deg, dec0=35*u.deg, d0=28*u.kpc, pmra0=0.8*u.mas/u.yr, pmdec0=0.3*u.mas/u.yr, vr0=-68*u.km/u.s, tstream=70*u.Myr, fra=False)

    props['distant'] = dict(label='New', wangle=360*u.deg, ra0=251.6*u.deg, dec0=2*u.deg, d0=87*u.kpc, pmra0=-0.6*u.mas/u.yr, pmdec0=-0.5*u.mas/u.yr, vr0=50*u.km/u.s, tstream=70*u.Myr, fra=True)
    
    props['s300'] = dict(label='S300', wangle=360*u.deg, ra0=152*u.deg, dec0=16*u.deg, d0=18*u.kpc, pmra0=-2.5*u.mas/u.yr, pmdec0=-2.75*u.mas/u.yr, vr0=300*u.km/u.s, tstream=30*u.Myr, fra=True)
    #props['distant'] = dict(label='New', wangle=360*u.deg, ra0=253*u.deg, dec0=-3*u.deg, d0=70*u.kpc, pmra0=-0.5*u.mas/u.yr, pmdec0=-0.5*u.mas/u.yr, vr0=-200*u.km/u.s, tstream=300*u.Myr, fra=False)

    return props[name]

def test(name, dra=2, best=True):
    """"""
    stream = Stream(name)
    
    if best:
        res = pickle.load(open('../data/fits/minimization_{:s}.pkl'.format(stream.savename), 'rb'))
        p0 = [x*y.unit for x, y in zip(res.x, stream.p0)]
        y, dist, pmra, pmdec, vr = p0
        print(p0)
        fit_label = 'Best-fit'
    else:
        y, dist, pmra, pmdec, vr = stream.p0
        fit_label = 'Initial'
        
    if stream.fra:
        ra = stream.x0*u.deg
        dec = y
    else:
        ra = y
        dec = stream.x0*u.deg
    
    c = coord.SkyCoord(ra=ra, dec=dec, distance=dist, pm_ra_cosdec=pmra, pm_dec=pmdec, radial_velocity=vr, frame='icrs')
    w0 = gd.PhaseSpacePosition(c.transform_to(gc_frame).cartesian)

    orbit = stream.ham.integrate_orbit(w0, dt=stream.dt, n_steps=stream.nstep)
    model = orbit.to_coord_frame(coord.ICRS, galactocentric_frame=stream.gc_frame)
    
    lo = stream.ham.integrate_orbit(w0, dt=stream.dt, n_steps=5000)
    print(np.nanmedian(lo.energy()), np.nanmedian(lo.angular_momentum()[2]))
    
    # determine orientation
    if stream.fra:
        model_x = model.ra.wrap_at(stream.wangle)
        model_y = model.dec
        ix = 0
        iy = 1
        xlabel = 'R.A. [deg]'
        ylabel = 'Dec [deg]'
    else:
        model_x = model.dec
        model_y = model.ra.wrap_at(stream.wangle)
        ix = -1
        #iy = 0
        tmp = stream.data['dec'][1]
        stream.data['dec'][1] = stream.data['dec'][0]
        stream.data['dec'][0] = tmp
        xlabel = 'Dec [deg]'
        ylabel = 'R.A. [deg]'
    
    # plot data
    plt.close()
    fig, ax = plt.subplots(5, 1, figsize=(7,11), sharex=True)

    fields = ['dec', 'dist', 'pmra', 'pmdec', 'vr']
    labels = [ylabel, 'Distance [kpc]', '$\mu_\\alpha$ [mas yr$^{-1}$]', '$\mu_\delta$ [mas yr$^{-1}$]',
            '$V_r$ [km s$^{-1}$]']
    model_fields = [model_y, model.distance, model.pm_ra_cosdec, model.pm_dec, model.radial_velocity]
    istart, iend = 0, -1

    for i in range(5):
        plt.sca(ax[i])
        
        if fields[i] in stream.data.keys():
            plt.plot(stream.data[fields[i]][ix], stream.data[fields[i]][1], 'k.', label='Data')
            plt.errorbar(stream.data[fields[i]][ix].value, stream.data[fields[i]][1].value, yerr=stream.data[fields[i]][2].value, fmt='none', color='k', alpha=0.7, label='')
            
        plt.plot(model_x[istart:iend], model_fields[i][istart:iend], '-', color='tab:blue', label='{:s} orbit'.format(fit_label))
        
        plt.ylabel(labels[i])
        if i==0:
            plt.legend(fontsize='small', handlelength=1)

    plt.minorticks_on()
    plt.xlim(np.min(stream.data['dec'][0].to(u.deg).value)-dra, np.max(stream.data['dec'][0].to(u.deg).value)+dra)
    plt.xlabel(xlabel)

    plt.tight_layout(h_pad=0)
    if best:
        plt.savefig('../plots/diag/best_{:s}.png'.format(stream.name))
    
def fit_stream(name, full=False):
    """"""
    
    stream = Stream(name, ham=ham, save_ext='')
    res = stream.orbit_minimize(save=True)
    stream.orbital_properties(save=True)
    
    t = Table.read('../data/fits/minimization_orbit_{:s}.fits'.format(name))
    t.pprint()
    
    if full:
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

def get_ham(pot):
    """Return the hamiltonian and stream extension for a given potential"""
    
    if pot=='bovy':
        ham_ = ham_bovy
        ext_ = 'bovy'
        lhood_ = ln_likelihood_icrs_mcmc_bovy
    elif pot=='heavy':
        ham_ = ham_heavy
        ext_ = 'heavy'
        lhood_ = ln_likelihood_icrs_mcmc_heavy
    else:
        ham_ = ham
        ext_ = ''
        lhood_ = ln_likelihood_icrs_mcmc
    
    return (ham_, ext_, lhood_)

def mcmc_stream(name, seed=249, nwalkers=64, nsteps=512, nth=3, pot='fid'):
    """"""
    
    ham_, ext_, lhood_ = get_ham(pot)
    stream = Stream(name, ham=ham_, save_ext=ext_)
    res = pickle.load(open('../data/fits/minimization_{:s}.pkl'.format(stream.savename), 'rb'))
    p0s = res.x
    
    pool = Pool(nth)
    np.random.seed(seed)
    p0 = emcee.utils.sample_ball(p0s, [1e-3, 1e-3, 1e-3, 1e-3, 1e-3], nwalkers)
    p0[:,1] = np.abs(p0[:,1])
    
    sampler = emcee.EnsembleSampler(nwalkers, p0.shape[1], log_prob_fn=lhood_, pool=pool, args=(stream.x0, stream.data_nounits, stream.nstep, stream.dt, stream.wangle, stream.fra))
    _ = sampler.run_mcmc(p0, nsteps)
    
    pickle.dump(sampler, open('../data/fits/mcmc_{:s}.pkl'.format(stream.savename), 'wb'))
    print(np.median(sampler.flatchain, axis=0))
    
    pool.close()

def plot_chains(sampler, name):
    """Plot chain"""
    names = [r'$\phi_2$', r'd', r'$\mu_{\phi_1}$', r'$\mu_{\phi_2}$', r'$V_r$']
    plt.close()
    fig, ax = plt.subplots(sampler.ndim, figsize=(10,10), sharex=True)

    for k in range(sampler.ndim):
        for walker in sampler.chain[..., k]:
            ax[k].plot(walker, marker='', drawstyle='steps-mid', alpha=0.2)
        ax[k].set_ylabel(names[k])

    plt.sca(ax[sampler.ndim-1])
    plt.xlabel('Step')

    plt.tight_layout(h_pad=0)
    plt.savefig('../plots/diag/chains_{:s}.png'.format(name))

def plot_corner(flatchain, name, bins=25):
    """Corner plot"""
    names = [r'$\phi_2$', r'd', r'$\mu_{\phi_1}$', r'$\mu_{\phi_2}$', r'$V_r$']
    
    plt.close()
    corner.corner(flatchain, bins=bins, labels=names, show_titles=True, title_fmt='.2f', title_kwargs=dict(fontsize='small'))

    plt.tight_layout(h_pad=0.1,w_pad=0.1)
    plt.savefig('../plots/diag/corner_{:s}.png'.format(name))

def plot_models(flatchain, stream, nplot=100, dra=2):
    """"""
    
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
    
    # plot data
    plt.close()
    fig, ax = plt.subplots(5, 1, figsize=(7,11), sharex=True)

    fields = ['dec', 'dist', 'pmra', 'pmdec', 'vr']
    labels = [ylabel, 'Distance [kpc]', '$\mu_\\alpha$ [mas yr$^{-1}$]', '$\mu_\delta$ [mas yr$^{-1}$]', '$V_r$ [km s$^{-1}$]']
    istart, iend = 0, -1

    for i in range(5):
        plt.sca(ax[i])
        
        if fields[i] in stream.data.keys():
            plt.plot(stream.data[fields[i]][ix], stream.data[fields[i]][1], 'k.', label='Data')
            plt.errorbar(stream.data[fields[i]][ix].value, stream.data[fields[i]][1].value, yerr=stream.data[fields[i]][2].value, fmt='none', color='k', alpha=0.7, label='')
            
        plt.ylabel(labels[i])
    
    plt.minorticks_on()
    plt.xlim(np.min(stream.data['dec'][0].to(u.deg).value)-dra, np.max(stream.data['dec'][0].to(u.deg).value)+dra)
    plt.xlabel(xlabel)

    # plot models
    for j in range(nplot):
        p0 = [x*y.unit for x, y in zip(flatchain[j], stream.p0)]
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

            plt.plot(model_x[istart:iend], model_fields[i][istart:iend], '-', color='tab:blue', label='Sample orbit', zorder=0, lw=0.5, alpha=0.3)

            if (i==0) & (j==0):
                plt.legend(fontsize='small', handlelength=1)

    plt.tight_layout(h_pad=0)
    plt.savefig('../plots/diag/stream_models_{:s}.png'.format(stream.savename))

def save_orbits(flatchain, stream):
    """Save orbits for a sample of points in the pdf"""
    
    nsample = np.shape(flatchain)[0]
    tout = Table(names=('ecc', 'rperi', 'rapo', 'vcirc', 'lx', 'ly', 'lz', 'etot'))
    
    # calculate orbits for steps in the chain
    for j in range(nsample):
        orbit = stream.orbital_properties(pbest=flatchain[j])
        
        l = orbit.angular_momentum()
        lx = np.nanmedian(l[0])
        ly = np.nanmedian(l[1])
        lz = np.nanmedian(l[2])
        etot = np.nanmedian(orbit.energy())
        rperi = orbit.pericenter()
        rapo = orbit.apocenter()
        ecc = orbit.eccentricity()
        vcirc = ham.potential.circular_velocity(np.array([rapo.to(u.kpc).value, 0, 0]))[0]
        
        trow = dict(rperi=rperi, rapo=rapo, ecc=ecc, vcirc=vcirc, lx=lx, ly=ly, lz=lz, etot=etot)
        tout.add_row(trow)
    
    # add units
    for k in trow.keys():
        tout[k].unit = trow[k].unit
    
    tout.write('../data/output/orbit_props_{:s}.fits'.format(stream.savename), overwrite=True)

def save_samples(flatchain, stream):
    """Save the flatchain subsample"""
    
    N = np.shape(flatchain)[0]
    
    if stream.fra:
        ra = np.ones(N) * stream.x0 * u.deg
        dec = flatchain[:,0] * u.deg
    else:
        ra = flatchain[:,0] * u.deg
        dec = np.ones(N) * stream.x0 * u.deg
    dist = flatchain[:,1] * u.kpc
    pmra = flatchain[:,2] * u.mas/u.yr
    pmdec = flatchain[:,3] * u.mas/u.yr
    vr = flatchain[:,4] * u.km/u.s
    
    tout = Table([ra, dec, dist, pmra, pmdec, vr], names=('ra', 'dec', 'dist', 'pmra', 'pmdec', 'vr'))
    tout.write('../data/output/samples_{:s}.fits'.format(stream.savename), overwrite=True)

def check_orbit_props(name):
    """"""
    t = Table.read('../data/output/orbit_props_{:}.fits'.format(name))
    #t.pprint()
    
    for k in t.colnames:
        print(k, '{:.3f} {:.3f} {:.3f}'.format(*np.nanpercentile(t[k], [16,50,84])))

def diagnose_mcmc(name, stage=0, pot='fid'):
    """"""
    
    ham_, ext_, lhood_ = get_ham(pot)
    stream = Stream(name, ham=ham_, save_ext=ext_)
    sampler = pickle.load(open('../data/fits/mcmc_{:s}.pkl'.format(stream.savename), 'rb'))
    
    if stage==0:
        plot_chains(sampler, stream.savename)
    
    chain = sampler.chain[:,256:,:]
    flatchain = np.reshape(chain,(-1,5))
    
    if stage==1:
        plot_corner(flatchain, stream.savename)
        
    np.random.seed(391)
    flatchain_short = np.random.permutation(flatchain)[:1000,:]
    
    if stage==2:
        plot_models(flatchain_short, stream, nplot=50)
    
    if stage==3:
        save_orbits(flatchain_short, stream)
    
    if stage==4:
        save_samples(flatchain_short, stream)


#####################
# Orbital systematics

def collate_fits(save_ext=''):
    """"""
    
    names = get_names()
    
    if len(save_ext):
        save_ext = '_' + save_ext
    
    tout = Table()
    
    for stream in names:
        t = Table.read('../data/fits/minimization_orbit_{:s}{:s}.fits'.format(stream, save_ext))
        
        tout = vstack([tout, t])
    
    tout.pprint()
    tout.write('../data/minimization_orbit{:s}.fits'.format(save_ext), overwrite=True)

def collate_mcmc_fits():
    """"""
    names = get_names()
    
    for name in names[:]:
        tfid = Table.read('../data/output/orbit_props_{:s}.fits'.format(name))
        tfid['potential'] = 'fiducial'
        tbovy = Table.read('../data/output/orbit_props_{:s}_bovy.fits'.format(name))
        tbovy['potential'] = 'bovy'
        theavy = Table.read('../data/output/orbit_props_{:s}_heavy.fits'.format(name))
        theavy['potential'] = 'heavy'
        
        tout = vstack([tfid, tbovy, theavy])
        tout.write('../data/output/orbit_props_{:s}_combined.fits'.format(name), overwrite=True)

def potential_comparison():
    """"""
    
    t = Table.read('../data/minimization_orbit.fits')
    t_heavy = Table.read('../data/minimization_orbit_heavy.fits')
    t_bovy = Table.read('../data/minimization_orbit_bovy.fits')
    
    props = ['rperi', 'rapo', 'ecc']
    labels = ['$r_{peri}$', '$r_{apo}$', 'eccentricity']
    units = ['[kpc]', '[kpc]', '']
    potentials = ['fiducial', 'heavy', 'light']
    
    colors = ['k', 'tab:blue', 'tab:orange']
    tables = [t, t_heavy, t_bovy]
    Nbin = 10
    bins = [np.linspace(0,30,Nbin), np.linspace(0,80,Nbin), np.linspace(0,1,Nbin)]
    
    plt.close()
    fig, ax = plt.subplots(3, 3, figsize=(17,8.5), gridspec_kw=dict(height_ratios=[1,3,1]), sharex='col')
    
    for e in range(3):
        plt.sca(ax[0][e])
        for i in range(3):
            plt.hist(tables[i][props[e]], bins=bins[e], histtype='step', lw=2, color=colors[i], zorder=3-i, density=False, label=potentials[i])
        plt.ylabel('Number')
        if e==0:
            plt.legend(loc=1, fontsize='small', frameon=False)
        
        plt.sca(ax[1][e])
        x = np.linspace(np.min(t[props[e]]), np.max(t[props[e]]), 100)
        plt.plot(x, x, 'k-', lw=1, alpha=0.5, label='')
        
        plt.plot(t[props[e]], t_heavy[props[e]], 'o', label='heavy')
        plt.plot(t[props[e]], t_bovy[props[e]], 'o', label='light')
        
        plt.ylabel('Alternative {:s} {:s}'.format(labels[e], units[e]))
        
        plt.sca(ax[2][e])
        plt.axhline(0, color='k', lw=1, alpha=0.5, label='')
        
        f_heavy = 1 - t_heavy[props[e]]/t[props[e]]
        hm, hsig = np.median(f_heavy), np.std(f_heavy)
        plt.plot(t[props[e]], f_heavy, 'o', label='{:.2f}, {:.2f}'.format(hm, hsig))
        
        f_light = 1 - t_bovy[props[e]]/t[props[e]]
        lm, lsig = np.median(f_light), np.std(f_light)
        plt.plot(t[props[e]], f_light, 'o', label='{:.2f}, {:.2f}'.format(lm, lsig))
        
        plt.legend(fontsize='small', frameon=False)
        plt.xlabel('Fiducial {:s} {:s}'.format(labels[e], units[e]))
        plt.ylabel('1 - alt / fid')
    
    plt.tight_layout(h_pad=0)
    plt.savefig('../plots/potential_comparison.png')


#######
# Paper

def ham_params():
    """"""
    print(ham.potential.parameters)
