from orbit_fits import *

def ref_frame():
    """Print properties of the reference frame"""
    print(gc_frame)
    
def potentials():
    """Print properties of the gravitational potentials used"""
    
    pot = [ham, ham_bovy, ham_heavy]
    name = ['fiducial', 'bovy', 'heavy']
    
    #pos = np.array([[0, 0, 25], [0,0,200]]).T * u.kpc
    pos = np.array([[25, 0, 0], [200,0,0]]).T * u.kpc
    mass = np.zeros((3,2)) * u.Msun
    
    for e, p in enumerate(pot):
        print(name[e])
        
        # potential parameters
        keys = p.potential.parameters.keys()
        for k in keys:
            print(k, p.potential.parameters[k])
        
        # enclosed mass
        mass[e] = p.potential.mass_enclosed(pos)
        print(mass[e])
    
    print(mass[0])
    print(mass[1]/mass[0])
    print(mass[2]/mass[0])

def plot_enclosed_mass():
    """Plot the ratio of enclosed mass for the adopted potentials"""
    
    pot = [ham, ham_bovy, ham_heavy]
    name = ['Fiducial', 'MWPotential2014', 'Price-Whelan & Bonaca (2018)']
    colors = ['k', 'tab:blue', 'tab:red']
    
    pos = np.zeros((3,100)) * u.kpc
    pos[0] = np.logspace(np.log10(20./100.), np.log10(20*10.), pos.shape[1]) * u.kpc
    
    mass = np.zeros((3,100))
    
    for e, p in enumerate(pot):
        mass[e] = p.potential.mass_enclosed(pos)
    
    plt.close()
    plt.figure(figsize=(8,6))
    
    for i in range(3):
        plt.plot(pos[0], mass[i]/mass[0], '-', color=colors[i], label=name[i])
    
    plt.axvline(25, color='k', ls=':')
    plt.legend(fontsize='small', loc=0)
    plt.ylim(0.7, 1.3)
    plt.xlim(0,200)
    
    plt.xlabel('r [kpc]')
    plt.ylabel('M(<r) / M$_{fid}$(<r)')
    
    plt.tight_layout()
    plt.savefig('../plots/response/enclosed_mass_potentials.png')
