import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d, interp2d
from scipy.special import kn as besselk

from scipy.optimize import minimize

from mpl_toolkits.mplot3d import Axes3D, axes3d, art3d
from matplotlib.patches import PathPatch
from matplotlib.text import TextPath
from matplotlib.transforms import TransformedBbox, Affine2D

def get_graphene_vectors(strain, carbon_carbon_distance=1.42, poisson_ratio=0.165):
    δ = strain;
    η = poisson_ratio;
    a0 = carbon_carbon_distance;
    R_δ = δ * np.array([[-η,0],[0,1]]);

    Am_δ0 = (a0/2)*np.array([np.sqrt(3),3]); # isotropic
    Am = (a0/2)*np.array([np.sqrt(3)*(1 - δ*η), 3*(1 + δ)]); # with strain

    An_δ0 = (a0/2)*np.array([-np.sqrt(3),3]); # isotropic
    An = (a0/2)*np.array([-np.sqrt(3)*(1 - δ*η),3*(1 + δ)]); # with strain

    Am_δ0_60 = np.array([a0*np.sqrt(3),0]); # Am_δ0 rotated 60 degrees to sit on x-axis
    Am_60 = Am_δ0_60 + (R_δ@Am_δ0_60); # rotated with strain
    An_60_δ0 = (a0/2)*np.array([np.sqrt(3),3]); # An_δ0 rotated 60
    An_60 = An_60_δ0 + (R_δ@An_60_δ0); # rotated with strain


    # basis vectors
    b1 = a0*np.array([0.0, 1 + δ]);
    b2 = a0*np.array([0.0, 2*(1 + δ)]);

    # reciprocal lattice vectors

    gm = (2*np.pi/3/a0)*np.array([3/(1 - δ*η), np.sqrt(3)/(1 + δ)]);
    gn = np.array([-gm[0],gm[1]]);

    gm_60 = (2*np.pi/3/a0)*np.array([np.sqrt(3)/(1 - δ*η), -1/(1 + δ)]);
    gn_60 = (4*np.pi/3/a0)*np.array([0, 1/(1 + δ)]);
    return Am_60, An_60, b1, b2, gm_60, gn_60

def get_graphene_c_one_third_vectors(strain, carbon_carbon_distance=1.42, poisson_ratio=0.165):
    Am, An, b1, b2, gm, gn = get_graphene_vectors(strain, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)
    Am_c_one_third = 3 * Am;
    An_c_one_third = np.array([3,1]) * An;
    return Am_c_one_third, An_c_one_third

def get_graphene_carbon_atoms(strain, box_dims, carbon_carbon_distance=1.42, poisson_ratio=0.165):
    if (type(box_dims) is float) or (type(box_dims) is int):
        box_dims = [-box_dims/2,box_dims/2,-box_dims/2,box_dims/2]
    elif (len(box_dims) == 2) and (np.size(box_dims) == 2):
        box_dims = [-box_dims[0]/2,box_dims[0]/2,-box_dims[1]/2,box_dims[1]/2]
    elif (len(box_dims) != 4) or (np.size(box_dims) != 4):
        raise ValueError("box_dims must be float, 1D array-like with 2 elements, or 1D array-like with 4 elements")

    Am, An, b1, b2, gm, gn = get_graphene_vectors(strain, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)

    box_dims = np.array(box_dims).reshape((2,2))
    box_x = box_dims[0,:]
    box_y = box_dims[1,:]

    max_j = max(int(np.round((box_y[1] - b1[1])/An[1])),int(np.round((box_y[1] - b2[1])/An[1])))
    min_j = min(int(np.round((box_y[0] - b1[1])/An[1])),int(np.round((box_y[0] - b2[1])/An[1])))

    max_i = max(int(np.round((box_x[1] - min_j*An[0] - b1[0])/Am[0])),int(np.round((box_x[1] - min_j*An[0] - b2[0])/Am[0])))
    min_i = min(int(np.round((box_x[0] - max_j*An[0] - b1[0])/Am[0])),int(np.round((box_x[0] - max_j*An[0] - b2[0])/Am[0])))

    grid = np.mgrid[min_i:max_i:(1 + max_i - min_i)*1j,min_j:max_j:(1 + max_j - min_j)*1j]

    a_site = grid[0,:,:,None]*Am[:,None].T + grid[1,:,:,None]*An[:,None].T + b1
    b_site = grid[0,:,:,None]*Am[:,None].T + grid[1,:,:,None]*An[:,None].T + b2
    a_site_idx = (a_site[:,:,0] <= box_x[1] + 0.001) & (a_site[:,:,0] >= box_x[0] - 0.001) & (a_site[:,:,1] <= box_y[1] + 0.001) & (a_site[:,:,1] >= box_y[0] - 0.001)
    b_site_idx = (b_site[:,:,0] <= box_x[1] + 0.001) & (b_site[:,:,0] >= box_x[0] - 0.001) & (b_site[:,:,1] <= box_y[1] + 0.001) & (b_site[:,:,1] >= box_y[0] - 0.001)

    _a_site = np.array([a_site[:,:,0][a_site_idx], a_site[:,:,1][a_site_idx]]).T
    _b_site = np.array([b_site[:,:,0][b_site_idx], b_site[:,:,1][b_site_idx]]).T
    return _a_site, _b_site

def get_graphene_c_one_third_phase(strain, box_dims, carbon_carbon_distance=1.42, poisson_ratio=0.165):
    if (type(box_dims) is float) or (type(box_dims) is int):
        box_dims = [-box_dims/2,box_dims/2,-box_dims/2,box_dims/2]
    elif (len(box_dims) == 2) and (np.size(box_dims) == 2):
        box_dims = [-box_dims[0]/2,box_dims[0]/2,-box_dims[1]/2,box_dims[1]/2]
    elif (len(box_dims) != 4) or (np.size(box_dims) != 4):
        raise ValueError("box_dims must be float, 1D array-like with 2 elements, or 1D array-like with 4 elements")

    Am, An = get_graphene_c_one_third_vectors(strain, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)

    box_dims = np.array(box_dims).reshape((2,2))
    box_x = box_dims[0,:]
    box_y = box_dims[1,:]

    max_j = max(int(np.round((box_y[1])/An[1])),int(np.round((box_y[1])/An[1])))
    min_j = min(int(np.round((box_y[0])/An[1])),int(np.round((box_y[0])/An[1])))

    max_i = max(int(np.round((box_x[1] - min_j*An[0])/Am[0])),int(np.round((box_x[1] - min_j*An[0])/Am[0])))
    min_i = min(int(np.round((box_x[0] - max_j*An[0])/Am[0])),int(np.round((box_x[0] - max_j*An[0])/Am[0])))

    grid = np.mgrid[min_i:max_i:(1 + max_i - min_i)*1j,min_j:max_j:(1 + max_j - min_j)*1j]

    a_site = grid[0,:,:,None]*Am[:,None].T + grid[1,:,:,None]*An[:,None].T
    a_site_idx = (a_site[:,:,0] <= box_x[1]) & (a_site[:,:,0] >= box_x[0]) & (a_site[:,:,1] <= box_y[1]) & (a_site[:,:,1] >= box_y[0])

    _a_site = np.array([a_site[:,:,0][a_site_idx], a_site[:,:,1][a_site_idx]]).T
    return _a_site

def get_graphene_carbon_atom_NN_distace(strain, carbon_carbon_distance=1.42, poisson_ratio=0.165):
    δ = strain;
    η = poisson_ratio;
    a0 = carbon_carbon_distance;
    armchair_NN_distance = a0*(1 + δ)
    zigzag_NN_distance = np.sqrt(np.sum(((a0/2)*np.array([(1 - δ*η)*np.sqrt(3),(1 + δ)]))**2))
    return armchair_NN_distance, zigzag_NN_distance

def plot_graphene_lattice(strain, box_dims, fig_ax=None, carbon_carbon_distance=1.42, poisson_ratio=0.165):
    if fig_ax:
        fig, ax = fig_ax
    else:
        fig, ax = plt.subplots()
    Am, An, b1, b2, gm, gn = get_graphene_vectors(strain, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)
    a_site, b_site = get_graphene_carbon_atoms(strain,box_dims, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)

    if (type(box_dims) is float) or (type(box_dims) is int):
        box_dims = [-box_dims/2,box_dims/2,-box_dims/2,box_dims/2]
    elif (len(box_dims) == 2) and (np.size(box_dims) == 2):
        box_dims = [-box_dims[0]/2,box_dims[0]/2,-box_dims[1]/2,box_dims[1]/2]
    elif (len(box_dims) != 4) or (np.size(box_dims) != 4):
        raise ValueError("box_dims must be float, 1D array-like with 2 elements, or 1D array-like with 4 elements")

    box_dims = np.array(box_dims).reshape((2,2))
    box_x = box_dims[0,:]
    box_y = box_dims[1,:]

    ax.plot(a_site[:,0],a_site[:,1],marker="o",ms=4,linestyle="None",color="#555555", alpha=0.6)
    ax.plot(b_site[:,0],b_site[:,1],marker="o",ms=4,linestyle="None",color="#555555", alpha=0.6)

#    rect = Rectangle( (box_x[0],box_y[0]), box_x[1] - box_x[0], box_y[1] - box_y[0], linestyle = "dashed", facecolor = "None", edgecolor="k", clip_on=False)
#    ax.add_patch(rect)

#    ax.set_xlim([box_x[0] - 0.1*abs(box_x[0]),box_x[1] + 0.1*abs(box_x[1])])
#    ax.set_ylim([box_y[0] - 0.1*abs(box_y[0]),box_y[1] + 0.1*abs(box_y[1])])
    ax.set_xlim([box_x[0] ,box_x[1] ])
    ax.set_ylim([box_y[0] ,box_y[1] ])

    ax.set_aspect("equal")
    ax.set_xlabel(r"$x\ \mathrm{[\AA]}$")
    ax.set_ylabel(r"$y\ \mathrm{[\AA]}$")
    return fig, ax

def plot_graphene_lattice_with_c_one_third(strain, box_dims, fig_ax=None, carbon_carbon_distance=1.42, poisson_ratio=0.165):
    if fig_ax:
        fig, ax = fig_ax
    else:
        fig, ax = plt.subplots()
    c_one_third_lattice_points = get_graphene_c_one_third_phase(strain,box_dims, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)

    fig, ax = plot_graphene_lattice(strain, box_dims, fig_ax=(fig,ax), carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)

    ax.plot(c_one_third_lattice_points[:,0],c_one_third_lattice_points[:,1],marker="o",ms=8,linestyle="None",color="C2")
    return fig, ax

def c_one_third_commensurate_command(m,n,strain, carbon_carbon_distance=1.42, poisson_ratio=0.165):
    """Prints QMC simulation cell parameters for `2*m*n` C1/3 adsorption sites."""
    N = int(2*m*n)
    Am_c_on_third, An_c_on_third = get_graphene_c_one_third_vectors(strain, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)
    Lx = Am_c_on_third[0]*m
    Ly = An_c_on_third[1]*2*n
    print("-N {} --Lx {} --Ly {}".format(N,Lx,Ly))

def c_one_third_commensurate_command_plot(m,n,strain,fig_ax=None, carbon_carbon_distance=1.42, poisson_ratio=0.165):
    """Plots graphene lattice with `2*m*n` C1/3 adsorption sites."""
    Am_c_on_third, An_c_on_third = get_graphene_c_one_third_vectors(strain, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)
    box_dims = [Am_c_on_third[0]*m + 0.00000001,An_c_on_third[1]*2*n + 0.00000001]
    fig, ax = plot_graphene_lattice_with_c_one_third(strain, box_dims, fig_ax=fig_ax, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)
    return fig, ax

def roughly_square(n,strain, carbon_carbon_distance=1.42, poisson_ratio=0.165):
    """Prints QMC parameters for `(2n)^2` C1/3 adsorption sites with a roughly square simulation cell (for isotropic graphene, considerably less square with strain)."""
    c_one_third_commensurate_command(n,2*n,strain,carbon_carbon_distance=carbon_carbon_distance,poisson_ratio=poisson_ratio)

def roughly_square_plot(n,strain,fig_ax=None, carbon_carbon_distance=1.42, poisson_ratio=0.165):
    """Plots graphene lattice with `(2n)^2` C1/3 adsorption sites for a roughly square simulation cell (for isotropic graphene, considerably less square with strain)."""
    fig, ax = c_one_third_commensurate_command_plot(n,2*n,strain,fig_ax=fig_ax,carbon_carbon_distance=carbon_carbon_distance,poisson_ratio=poisson_ratio)
    return fig, ax

def get_LJ_parameters(strain,conventional=False, **kwargs):
    """Interpolate Lennard-Jones parameters based on parameters generated by Nichols et al. (https://doi.org/10.1103/PhysRevB.93.205412) for helium interacting with uniaxially strained graphene. Passes `**kwargs` to `scipy.interpolate.interp1d`."""
    if conventional:
        return 2.74, 16.2463

    optimized_lj_parameters = np.array([
        [0.00, 2.6426682611845180, 16.961130160262375],
        [0.05, 2.6691003148269594, 16.757965573159975],
        [0.10, 2.6999208415947040, 16.529721451450662],
        [0.15, 2.7381160674611458, 16.256015311769104],
        [0.16, 2.7485482036865276, 16.169791193974070],
        [0.17, 2.7585554088442380, 16.086298168462033],
        [0.18, 2.7688858608199607, 16.000654315229383],
        [0.19, 2.7800711002084983, 15.910775851211197],
        [0.20, 2.7906774893502124, 15.824759255848775],
        [0.21, 2.8040035563982340, 15.719998352422472],
        [0.22, 2.8171361080910140, 15.617586197043124],
        [0.23, 2.8311584765852196, 15.509815911364557],
        [0.24, 2.8456701010275367, 15.399913806990275],
        [0.25, 2.8580790323978964, 15.303052988043815],
        [0.26, 2.8784702596739513, 15.156337176229170],
        [0.27, 2.8964649920317695, 15.025373911567815],
        [0.28, 2.9165976207500730, 14.882223175940368],
        [0.29, 2.9377464777914035, 14.732693536111157],
        [0.34, 3.0653537015148230, 13.962563962971998]
    ])

    get_σ = interp1d(optimized_lj_parameters[:,0],optimized_lj_parameters[:,1],**kwargs)
    get_ε = interp1d(optimized_lj_parameters[:,0],optimized_lj_parameters[:,2],**kwargs)
    return np.array([get_σ(strain), get_ε(strain)])

# Vz
def Vz_64(z,sigma):
    return (sigma**4*(2*sigma**6 - 5*z**6))/(5*z**10)

def gradVz_x_64(z,sigma):
    return 0.0

def gradVz_y_64(z,sigma):
    return 0.0

def gradVz_z_64(z,sigma):
    return (4*sigma**4*(-sigma**6 + z**6))/z**11

def grad2Vz_64(z,sigma):
    return (4*sigma**4*(11*sigma**6 - 5*z^6))/z**12

# Vg
def Vg_64(x, y, z, sigma, g, gi, gj, bi, bj):
    return (g**2*sigma**4*(-480*z**3*besselk(2, g*z) + g**3*sigma**6*besselk(5, g*z))*np.cos(gi*(bi + x) + gj*(bj + y)))/(960*z**5)

def gradVg_x_64(x, y, z, sigma, g, gi, gj, bi, bj):
    return -(g**2*gi*sigma**4*(-480*z**3*besselk(2, g*z) + g**3*sigma**6*besselk(5, g*z))*sin(gi*(bi + x) + gj*(bj + y)))/(960*z**5)

def gradVg_y_64(x, y, z, sigma, g, gi, gj, bi, bj):
    return -(g**2*gj*sigma**4*(-480*z**3*besselk(2, g*z) + g**3*sigma**6*besselk(5, g*z))*sin(gi*(bi + x) + gj*(bj + y)))/(960*z**5)

def gradVg_z_64(x, y, z, sigma, g, gi, gj, bi, bj):
    return -(g**3*sigma**4*(10*(g**2*sigma**6*z - 48*z**5)*besselk(3, g*z) + g*sigma**6*(80 + g**2*z**2)*besselk(4, g*z))*cos(gi*(bi + x) + gj*(bj + y)))/(960*z**7)

def grad2Vg_64(x, y, z, sigma, g, gi, gj, bi, bj):
    return (g**2*sigma**4*(-120*g*z**6*(g*z*besselk(0, g*z) + 8*besselk(1, g*z)) + z*(19*g**4*sigma**6*z**2 + 480*z**4*(-6 + (gi**2 + gj**2)*z**2) - 8*g**2*(45*z**6 + sigma**6*(-110 + (gi**2 + gj**2)*z**2)))*besselk(2, g*z) + g*(-1680*z**6 + sigma**6*(5280 + z**2*(-48*(gi**2 + gj**2) + g**4*z**2 + g**2*(224 - (gi**2 + gj**2)*z**2))))*besselk(3, g*z))*cos(gi*(bi + x) + gj*(bj + y)))/(960*z**9)

#V
def V_64(strain, sigma, epsilon, x, y, z, carbon_carbon_distance=1.42, poisson_ratio=0.165, k_max=10000,potential="V"):
    if potential not in ["V", "gradVx", "gradVy", "gradVz", "grad2V"]:
        raise ValueError("potential option must be V, gradVx, gradVy, gradVz, or grad2V")
    if potential == "V":
        Vz_func = Vz_64
        Vg_func = Vg_64
    elif potential == "gradVx":
        Vz_func = gradVz_x_64
        Vg_func = gradVg_x_64
    elif potential == "gradVy":
        Vz_func = gradVz_y_64
        Vg_func = gradVg_y_64
    elif potential == "gradVz":
        Vz_func = gradVz_z_64
        Vg_func = gradVg_z_64
    elif potential == "grad2V":
        Vz_func = grad2Vz_64
        Vg_func = grad2Vg_64
    Am, An, b1, b2, gm, gn = get_graphene_vectors(strain, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)

    area_lattice = np.cross(Am,An)
    gi_array, gj_array = np.mgrid[-200:200:401j,-200:200:401j]
    g_magnitude_array = np.sqrt(((gi_array[:,:,None]*gm + gj_array[:,:,None]*gn)**2).sum(axis=-1))
    g_magnitude_array = g_magnitude_array.ravel()
    gi_array = gi_array.ravel()
    gj_array = gj_array.ravel()
    g_mag_sort_idx = np.argsort(g_magnitude_array)
    gi_array = gi_array[g_mag_sort_idx]
    gj_array = gj_array[g_mag_sort_idx]
    g_magnitude_array = g_magnitude_array[g_mag_sort_idx]

    flag_1 = False
    flag_2 = False
    flag_3 = False

    pf = 2*np.pi*epsilon*(sigma**2)/area_lattice
    _V = 0.0
    for b in [b1, b2]:
        _V += Vz_func(z,sigma)

    _V_old = 0.0

    for k in range(g_magnitude_array.shape[0] - 1):
        k = k + 1
        gi = gi_array[k]
        gj = gj_array[k]
        g = (gi * gm) + (gj * gn)
        g_magnitude = g_magnitude_array[k]
        for b in [b1, b2]:
            _V += Vg_func(x, y, z, sigma, g_magnitude, g[0], g[1], b[0], b[1])

        if (_V == _V_old) and (g_magnitude != g_magnitude_array[k+1]) and (not flag_1):
            _V_old = _V;
            flag_1 = True;
        elif (_V == _V_old) and (g_magnitude != g_magnitude_array[k+1]) and (not flag_2):
            _V_old = _V;
            flag_2 = True;
        elif (_V == _V_old) and (g_magnitude != g_magnitude_array[k+1]) and (not flag_3):
            _V_old = _V;
            flag_3 = True;
        elif (_V == _V_old) and (g_magnitude != g_magnitude_array[k+1]) and flag_1 and flag_2 and flag_3:
            _V_old = _V;
            break
        elif (k >= k_max) and (g_magnitude != g_magnitude_array[k+1]):
            #println("WARNING: k_max reached")
            break
        else:
            _V_old = _V;
    return pf * _V

def generate_V1D(x,y,z,strain=0.00,sigma=None,epsilon=None,conventional=False,carbon_carbon_distance=1.42,poisson_ratio=0.165,k_max=10000,potential="V"):
    
    if np.size(x) != 1:
        raise ValueError("x must be float or int")
    if np.size(y) != 1:
        raise ValueError("y must be float or int")

    _x = float(x)
    _y = float(y)
    V = np.zeros_like(z)

    _sigma, _epsilon = get_LJ_parameters(strain,conventional=conventional)
    if sigma is None:
        sigma=_sigma
    if epsilon is None:
        epsilon = _epsilon

    for i, _z in enumerate(z):
        V[i] = V_64(strain, sigma, epsilon, _x, _y, _z, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio, k_max=k_max, potential=potential)
    
    return V

def plot_V1D(x,y,z,strains=0.00,V_array=None,sigmas=None,epsilons=None,conventional=False,carbon_carbon_distance=1.42,poisson_ratio=0.165,k_max=10000,potential="V",mplstylefile="default",dpi=None,xlim=None,ylim=None,linestyle = "-",label = r"$\delta={:.2f}$",fig_ax=None,xlabel=r"$z\ \mathrm{[\AA]}$",ylabel=r"$V\ \mathrm{[K]}$"):
    if (type(strains) is float) or (type(strains) is int):
        strains = np.array([strains])
    strains = np.array(strains)

    if sigmas is not None:
        if (type(sigmas) is float) or (type(sigmas) is int):
            sigmas = np.ones_like(strains)*sigmas
        elif np.shape(sigmas) != np.shape(strains):
            raise ValueError("sigmas must be array like with same shape as strains")
        sigmas = np.array(sigmas)

    if epsilons is not None:
        if (type(epsilons) is float) or (type(epsilons) is int):
            epsilons = np.ones_like(strains)*epsilons
        elif np.shape(epsilons) != np.shape(strains):
            raise ValueError("epsilons must be array like with same shape as strains")
        epsilons = np.array(epsilons)

    if V_array is None:
        V_array = []
        for i,strain in enumerate(strains):
            if sigmas is not None:
                sigma = sigmas[i]
            else:
                sigma = None
            if epsilons is not None:
                epsilon = epsilons[i]
            else:
                epsilon = None
            V = generate_V1D(x,y,z,strain=strain,sigma=sigma,epsilon=epsilon,conventional=conventional,carbon_carbon_distance=carbon_carbon_distance,poisson_ratio=poisson_ratio,k_max=k_max,potential=potential)
            V_array.append(V)
        V_array = np.array(V_array)
    else:
        if V_array.shape[0] != strains.shape[0]:
            raise ValueError("length of V_array must be same as length of strains")
        if V_array.shape[1] != np.shape(z)[0]:
            raise ValueError("V_array.shape[1] must be same as length of z")
    with plt.style.context(mplstylefile):
        if fig_ax:
            fig, ax = fig_ax
        else:
            fig,ax = plt.subplots(dpi=dpi)

        for i, strain in enumerate(strains):
            V = V_array[i,:]
            ax.plot(z,V,label = label.format(strain),linestyle=linestyle)

        if ylim is None:
            ax.set_ylim(V_array.min()*1.1,50)
        else:
            ax.set_ylim(ylim)

        if xlim is not None:
            ax.set_xlim(xlim)
        ax.legend()
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    return fig, ax, V_array

def get_V2D(fn,box_dims,resolution = 1001,**kwargs):
    if (type(box_dims) is float) or (type(box_dims) is int):
        box_dims = [-box_dims/2,box_dims/2,-box_dims/2,box_dims/2]
    elif (len(box_dims) == 2) and (np.size(box_dims) == 2):
        box_dims = [-box_dims[0]/2,box_dims[0]/2,-box_dims[1]/2,box_dims[1]/2]
    elif (len(box_dims) != 4) or (np.size(box_dims) != 4):
        raise ValueError("box_dims must be float, 1D array-like with 2 elements, or 1D array-like with 4 elements")

    extent = np.array(box_dims)
    box_dims = extent.reshape((2,2))
    box_x = box_dims[0,:]
    box_y = box_dims[1,:]

    data = np.load(fn)
    cell_length_a,cell_length_b,cell_angle_gamma = data["uc_info"]
    side = np.array([cell_length_a, cell_length_b])
    uc_x = data["uc_x"]
    uc_y = data["uc_y"]
    V = data["potential"]
    A = data["A_xy_to_uc"]

    uc_to_V = interp2d(uc_x,uc_y,V,**kwargs)

    big_xy_x, big_xy_y = np.mgrid[box_x[0]:box_x[1]:resolution*1j,box_y[0]:box_y[1]:resolution*1j]
    big_uc_x, big_uc_y = cartesian_to_uc(big_xy_x, big_xy_y, A[0,0], A[0,1], A[1,0], A[1,1])

    big_uc_x, big_uc_y = put_in_BC(big_uc_x,big_uc_y,side)

    big_V = np.zeros_like(big_uc_x)
    for i in range(big_uc_x.shape[0]):
        for j in range(big_uc_x.shape[1]):
            _x = big_uc_x[i,j]
            _y = big_uc_y[i,j]
            big_V[i,j] = uc_to_V(_x,_y)

    return big_V, big_xy_x, big_xy_y, extent

def generate_V2D(x,y,z,strain,sigma,epsilon,carbon_carbon_distance=1.42,poisson_ratio=0.165,k_max=10000,potential="V"):
    if np.ndim(x) != 2:
        raise ValueError("x and y must be 2D array like")
    if np.shape(x) != np.shape(y):
        raise ValueError("x and y must have same shape.")
    if np.shape(z):
        raise ValueError("z must be int or float")
    x = np.array(x)
    y = np.array(y)
    V = np.zeros_like(x)
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            _x = x[i,j]
            _y = y[i,j]
            V[i,j] = V_64(strain, sigma, epsilon, _x, _y, z, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio, k_max=k_max,potential=potential)
    return V

def generate_V2D_uc(z, strain=0.00, resolution=101, fn_prefix="./helium-graphene-2D-unit-cell", carbon_carbon_distance=1.42, poisson_ratio=0.165, potential="V", k_max=10000, conventional=False, sigma=None, epsilon=None, **kwargs):
    Am, An, b1, b2, gm, gn = get_graphene_vectors(strain, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)
    cell_length_a = calculate_magnitude(Am)
    cell_length_b = calculate_magnitude(An)
    cell_angle_gamma= calculate_angle(Am,An);
    cell_angle_gamma_degrees = cell_angle_gamma*180/np.pi;
    area_lattice = cell_length_a*cell_length_b*np.sin(cell_angle_gamma);

    # Set up transfer matrices to transfer from unit cell coordinates to Cartesian coordinates and vice versa
    B = np.array([
        [1, np.cos(cell_angle_gamma)],
        [0, np.sin(cell_angle_gamma)]
    ]); #transfer to Cartesian

    A = np.array([
        [1, -1/np.tan(cell_angle_gamma)],
        [0, 1/np.sin(cell_angle_gamma)]
    ]); # transfer to unit cell coords

    uc_x, uc_y = np.mgrid[0:cell_length_a:resolution*1j,0:cell_length_b:resolution*1j]
    xy_x, xy_y = uc_to_cartesian(uc_x,uc_y,B[0,0],B[0,1],B[1,0],B[1,1])

    _sigma, _epsilon = get_LJ_parameters(strain,conventional=conventional, **kwargs)
    if sigma is None:
        sigma=_sigma
    if epsilon is None:
        epsilon = _epsilon

    V = generate_V2D(xy_x,xy_y,z,strain,sigma,epsilon,carbon_carbon_distance=carbon_carbon_distance,poisson_ratio=poisson_ratio,k_max=k_max,potential=potential)
    fn = "{}_{}_strain_{:.5f}_z_{:.5f}_res_{}.npz".format(fn_prefix,potential,strain,z,resolution)
    print("Saving 2D lookup table to {}".format(fn))
    np.savez(fn,
             potential=V,
             uc_x = uc_x,
             uc_y = uc_y,
             A_xy_to_uc = A,
             B_uc_to_xy = B,
             LJ_parameters = np.array([strain,sigma,epsilon]),
             uc_info = np.array([cell_length_a,cell_length_b,cell_angle_gamma]),
             lattice_vector_parameters = np.array([strain,carbon_carbon_distance,poisson_ratio])
            )
    data = np.load(fn)
    return data

def plot_V2D(fn,box_dims,V2D_data=None,resolution=1001,graphene_lattice=True,plot_filename=None,fig_ax=None,mplstylefile="default",dpi=None,imshow_kwargs = {},**kwargs):
    if V2D_data is None:
        big_V, big_xy_x, big_xy_y, extent = get_V2D(fn,box_dims,resolution=resolution,**kwargs)
    else:
        big_V, big_xy_x, big_xy_y, extent = V2D_data
        
    V2D_data = (big_V, big_xy_x, big_xy_y, extent)
    with plt.style.context(mplstylefile):
        if fig_ax:
            fig, ax = fig_ax
        else:
            fig,ax = plt.subplots(dpi=dpi)
     
        _imshow_kwargs = {"origin":"lower", "cmap":"viridis", "interpolation":"none", "rasterized":False, "extent":extent}
        for key in imshow_kwargs.keys():
            _imshow_kwargs[key] = imshow_kwargs[key]
        im = ax.imshow(big_V.T,**_imshow_kwargs)
        cb = fig.colorbar(im)

        if graphene_lattice:
            data = np.load(fn)
            strain, carbon_carbon_distance, poisson_ratio = data["lattice_vector_parameters"]
            armchair_NN_distance, zigzag_NN_distance = get_graphene_carbon_atom_NN_distace(strain, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)
            a_site_atoms, b_site_atoms = get_graphene_carbon_atoms(strain, extent, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)
            graphene_atoms = np.vstack((a_site_atoms,b_site_atoms))

            for _i in range(graphene_atoms.shape[0]):
                atom1 = graphene_atoms[_i,:]
                for _j in range(_i + 1,graphene_atoms.shape[0]):
                    atom2 = graphene_atoms[_j,:]
                    draw_line = False
                    atom_distance = np.sqrt(np.sum((atom2 - atom1)**2))
                    if (((atom_distance <= 1.01*armchair_NN_distance) and (atom_distance >= 0.99*armchair_NN_distance)) or
                        ((atom_distance <= 1.01*zigzag_NN_distance) and (atom_distance >= 0.99*zigzag_NN_distance))
                       ):
                        draw_line = True
                    if draw_line:
                        ax.plot([atom1[0],atom2[0]],[atom1[1],atom2[1]], color = "k",alpha = 0.5)
            ax.plot(graphene_atoms[:,0],graphene_atoms[:,1], linestyle = "None", color = "k", marker="o",alpha = 1.0)

        ax.set_xlabel(r"$x\ \mathrm{[\AA]}$")
        ax.set_ylabel(r"$y\ \mathrm{[\AA]}$")
        cb.set_label(r"$V\ \mathrm{[K]}$")
        if plot_filename:
            fig.savefig(plot_filename,dpi=dpi)

    return fig, ax, V2D_data

def calculate_angle_with_magnitude(x,y,magnitude_x,magnitude_y):
    return np.arccos(np.dot(x,y)/magnitude_x/magnitude_y)

def calculate_angle(x,y):
    return calculate_angle_with_magnitude(x,y,calculate_magnitude(x),calculate_magnitude(y))

def calculate_magnitude(vec):
    return np.sqrt(np.sum(vec**2))

def put_in_BC(x,y,side):
    x -= side[0] * np.floor(x/side[0])
    y -= side[1] * np.floor(y/side[1])
    return x, y

def cartesian_to_uc(x, y, A00, A01, A10, A11):
    return A00*x + A01*y, A10*x + A11*y

def uc_to_cartesian(uc_x, uc_y, B00, B01, B10,B11):
    return cartesian_to_uc(uc_x, uc_y, B00, B01, B10, B11)

def get_z_min(x=0.0,y=0.0,z0=3.0,strain=0.00,sigma=None,epsilon=None,conventional=False,carbon_carbon_distance=1.42,poisson_ratio=0.165,k_max=10000,potential="V"):
    #minimize_me = lambda z: generate_V1D(x,y,z,strain=strain,sigma=sigma,epsilon=epsilon,conventional=conventional,carbon_carbon_distance=carbon_carbon_distance,poisson_ratio=poisson_ratio,k_max=k_max,potential=potential)
    result = minimize(lambda z: generate_V1D(x,y,z,strain=strain,sigma=sigma,epsilon=epsilon,conventional=conventional,carbon_carbon_distance=carbon_carbon_distance,poisson_ratio=poisson_ratio,k_max=k_max,potential=potential),z0)
    return result.x[0], result.fun

def get_z_mins(strains=0.0,**kwargs):
    if (type(strains) is float) or (type(strains) is int):
        strains = np.array([strains])
    strains = np.array(strains)
    z_mins = np.zeros_like(strains)
    V_mins = np.zeros_like(strains)
    for i, strain in enumerate(strains):
        z_mins[i], V_mins[i] = get_z_min(strain=strain,**kwargs)
    return z_mins, V_mins

# 3D plotting tools
def lims(mplotlims):
    scale = 1.021
    offset = (mplotlims[1] - mplotlims[0])*scale
    return mplotlims[1] - offset, mplotlims[0] + offset

def text_on_2D(ax, s, xy = (0.5,0.95), xlim=None, ylim=None, size=None, angle=0.0, usetex=False, center=True, aspect_ratio=1.0, **kwargs):
    #FIXME expand and add to text_on_the_wall
    #xlim = ax.get_xlim()
    #ylim = ax.get_ylim()
    #aspect_ratio = (ylim[1] - ylim[0])/(ylim[1] - ylim[0])
    #size = 0.8
    #center = True
    #angle = 0.0
    #usetex=False

    #s = r"$\delta=0.00$"
    

    if xlim == None:
        xlim = ax.get_xlim()
    if ylim == None:
        ylim = ax.get_ylim()
    x = xlim[0] + (xlim[1] - xlim[0])*xy[0]
    y = ylim[0] + (ylim[1] - ylim[0])*xy[1]

    xy1 = (x, y)

    text_path = TextPath((0, 0), s, size=size, usetex=usetex)

    #Scale by aspect ratio
    scale_path = Affine2D().scale(1.0, sy=aspect_ratio)
    text_path = scale_path.transform_path(text_path)

    #Get bbox to center text
    bbox = text_path.get_extents()
    bbox_points = bbox.get_points()
    _b = bbox_points.sum(axis=0)/2

    if center == True:
        trans = Affine2D().rotate(angle).translate(xy1[0] - _b[0], xy1[1] - _b[1])
    else:
        trans = Affine2D().rotate(angle).translate(xy1[0], xy1[1] - _b[1])

    tp = trans.transform_path(text_path)
    p1 = PathPatch(tp, **kwargs)
    ax.add_patch(p1)
    
def text_on_the_wall(ax, xyz, s, zdir="z", size=None, angle=0, usetex=False, center=True, aspect_ratio=1.0, **kwargs):
    x, y, z = xyz
    if zdir == "y":
        xy1, z1 = (x, z), y
    elif zdir == "x":
        xy1, z1 = (y, z), x
    else:
        xy1, z1 = (x, y), z

    text_path = TextPath((0, 0), s, size=size, usetex=usetex)
    
    #Scale by aspect ratio
    scale_path = Affine2D().scale(1.0, sy=aspect_ratio)
    text_path = scale_path.transform_path(text_path)

    #Get bbox to center text
    bbox = text_path.get_extents()
    bbox_points = bbox.get_points()
    _b = bbox_points.sum(axis=0)/2
    
    if center == True:
        trans = Affine2D().rotate(angle).translate(xy1[0] - _b[0], xy1[1] - _b[1])
    else:
        trans = Affine2D().rotate(angle).translate(xy1[0], xy1[1] - _b[1])
  
    tp = trans.transform_path(text_path)
    p1 = PathPatch(tp, **kwargs)
    ax.add_patch(p1)

    art3d.pathpatch_2d_to_3d(p1, z=z1, zdir=zdir)
    
def plot_V3D(fn,box_dims,V2D_data=None,resolution=1001,graphene_lattice=True,plot_filename=None,fig_ax=None,mplstylefile="default",dpi=None,add_pane_edge=True,pane_edge_color="#c8c8c8",pane_face_color="#f0f0f5",horizontal_shift=.5,vertical_shift=-0.12,labelpad=[10.0,10.0,17.0],ztickpad=10.0,figsize=None,xlim=None,ylim=None,zlim=None,surf_kwargs={},**kwargs):
    if V2D_data is None:
        big_V, big_xy_x, big_xy_y, extent = get_V2D(fn,box_dims,resolution=resolution,**kwargs)
    else:
        big_V, big_xy_x, big_xy_y, extent = V2D_data

    V2D_data = (big_V, big_xy_x, big_xy_y, extent)

    data = np.load(fn)
    strain, carbon_carbon_distance, poisson_ratio = data["lattice_vector_parameters"]

    #https://stackoverflow.com/questions/47334760/3d-figures-from-matplotlib-visibility-of-pane-edge
    #https://stackoverflow.com/questions/37324837/tweaking-axis-labels-and-names-orientation-for-3d-plots-in-matplotlib

    with plt.style.context(mplstylefile):
        if fig_ax:
            fig, ax = fig_ax
        else:
            fig = plt.figure(figsize=figsize,dpi=dpi)
            ax = fig.add_subplot(111, projection='3d', proj_type = 'ortho')
        az,el = ax.azim, ax.elev

        #ax.dist = 10
        #ax.view_init(azim=az, elev= el)

        # Remove grid lines
        ax.grid(False)

        # Remove tick labels
        #ax.set_xticklabels([])
        ##ax.set_yticklabels([])
        #ax.set_zticklabels([])

        # No ticks
        #ax.set_xticks([]) 
        #ax.set_yticks([]) 
        #ax.set_zticks([])
        
        _surf_kwargs = {"cmap":"viridis", "linewidth":0, "antialiased":False, "alpha":0.5}
        for key in surf_kwargs.keys():
            _surf_kwargs[key] = surf_kwargs[key]
        surf = ax.plot_surface(big_xy_x, big_xy_y, big_V, **_surf_kwargs)
        
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)    
        if zlim:
            ax.set_zlim(zlim)
            
        ax.set_xlabel(r"$x\ \mathrm{[\AA]}$", labelpad=labelpad[0])
        ax.set_ylabel(r"$y\ \mathrm{[\AA]}$", labelpad=labelpad[1])
        ax.set_zlabel(r"$V\ \mathrm{[K]}$", labelpad=labelpad[2])

        ax.tick_params(axis='z', which='major', pad=ztickpad)

        xlims2, ylims2, zlims2 = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()

        if graphene_lattice:
            armchair_NN_distance, zigzag_NN_distance = get_graphene_carbon_atom_NN_distace(strain, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)
            a_site_atoms, b_site_atoms = get_graphene_carbon_atoms(strain, extent, carbon_carbon_distance=carbon_carbon_distance, poisson_ratio=poisson_ratio)
            graphene_atoms = np.vstack((a_site_atoms,b_site_atoms))

            for _i in range(graphene_atoms.shape[0]):
                atom1 = graphene_atoms[_i,:]
                for _j in range(_i + 1,graphene_atoms.shape[0]):
                    atom2 = graphene_atoms[_j,:]
                    draw_line = False
                    atom_distance = np.sqrt(np.sum((atom2 - atom1)**2))
                    if (((atom_distance <= 1.01*armchair_NN_distance) and (atom_distance >= 0.99*armchair_NN_distance)) or
                        ((atom_distance <= 1.01*zigzag_NN_distance) and (atom_distance >= 0.99*zigzag_NN_distance))
                       ):
                        draw_line = True
                    if draw_line:
                        ax.plot([atom1[0],atom2[0]],[atom1[1],atom2[1]], zlims2[0], color = "k",alpha = 0.5)
            ax.plot(graphene_atoms[:,0],graphene_atoms[:,1], zlims2[0], linestyle = "None", color = "k", marker="o",alpha = 1.0)

        if add_pane_edge:
            xlims, ylims, zlims = lims(ax.get_xlim()), lims(ax.get_ylim()), lims(ax.get_zlim())
            i = np.array([xlims[0], ylims[0], zlims[0]])
            f = np.array([xlims[0], ylims[0], zlims[1]])
            p = art3d.Poly3DCollection(np.array([[i, f]]))

            p.set_color(pane_edge_color)
            ax.add_collection3d(p)

        #Set pane edge color
        ax.xaxis.pane.set_edgecolor(pane_edge_color)
        ax.yaxis.pane.set_edgecolor(pane_edge_color)
        ax.zaxis.pane.set_edgecolor(pane_edge_color)

        #Set pane face color
        ax.xaxis.pane.set_facecolor(pane_face_color)
        ax.yaxis.pane.set_facecolor(pane_face_color)
        ax.zaxis.pane.set_facecolor(pane_face_color)

        ax.xaxis.pane.set_alpha(1)
        ax.yaxis.pane.set_alpha(1)
        ax.zaxis.pane.set_alpha(1)
        ax.xaxis.pane.fill = True
        ax.yaxis.pane.fill = True
        ax.zaxis.pane.fill = True

        # Get aspect ratio to fix text on wall
        aspect_ratio_xz = (zlims2[1] - zlims2[0])/(xlims2[1] - xlims2[0])
        aspect_ratio_yz = (zlims2[1] - zlims2[0])/(ylims2[1] - ylims2[0])
        size = 0.5
        text_on_the_wall(ax, ((xlims2[1] + xlims2[0])/2, ylims2[1], zlims2[0] + (zlims2[1] - zlims2[0])*0.95), r"$\mathrm{Adsorption\ Potential\ [K]}$", size=size, zdir="y", usetex=plt.rcParams["text.usetex"], aspect_ratio=aspect_ratio_xz, ec="none", fc="k")
        text_on_the_wall(ax, (xlims2[0], (ylims2[1] + ylims2[0])/2, zlims2[0] + (zlims2[1] - zlims2[0])*0.95), r"$\mathrm{{Strain\ \delta = {:.2f}}}$".format(strain), size=size, zdir="x", usetex=plt.rcParams["text.usetex"], aspect_ratio=aspect_ratio_yz, ec="none", fc="k")

        # Set ticks inward
        ax.xaxis._axinfo['tick']['outward_factor'] = 0.3
        ax.xaxis._axinfo['tick']['inward_factor'] = 0.0
        ax.yaxis._axinfo['tick']['outward_factor'] = 0.3
        ax.yaxis._axinfo['tick']['inward_factor'] = 0.0
        ax.zaxis._axinfo['tick']['outward_factor'] = 0.3
        ax.zaxis._axinfo['tick']['inward_factor'] = 0.0
        fig.tight_layout()
        fig.subplots_adjust(top=1, bottom=0, left=0, right=1, wspace=0, hspace=0)

        if plot_filename:
            # Shift figure to cut out whitespace
            figsize = fig.get_size_inches()
            bbox = fig.bbox_inches.from_bounds(horizontal_shift, vertical_shift, figsize[0], figsize[1])
            fig.savefig(plot_filename,dpi=dpi,bbox_inches=bbox)

        return fig, ax, V2D_data

