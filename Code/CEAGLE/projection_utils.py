
import numpy as np
import yt 
yt.set_log_level(40)


def fibonacci_sphere(samples: int = 10) -> np.ndarray:
    """
    Generate approximately uniformly spaced directions on a unit sphere
    using the Fibonacci lattice method.
    """
    indices = np.arange(0, samples, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/samples)
    theta = np.pi * (1 + 5**0.5) * indices  # golden angle

    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)
    return np.stack((x, y, z), axis=1)



def compute_bounding_box(data_all, length_unit='Mpc', margin=0.05, cubic=False, verbose=True):
    """
    Compute a bounding box encompassing all particle types in a flat yt-style dictionary.

    Parameters
    ----------
    data_all : dict
        Dictionary with keys like ("star","particle_position_x").
    length_unit : str, optional
        Unit for yt.load_particles. Default: 'Mpc'.
    mass_unit : str, optional
        Unit for yt.load_particles. Default: 'Msun'.
    margin : float, optional
        Fractional padding to extend the bounding box beyond min/max values (default 5%).
    cubic : bool, optional
        If True, make the bounding box cubic and centered on the system.
    verbose : bool, optional
        If True, print information about the box.

    Returns
    -------
    bbox : np.ndarray
        The 3×2 bounding box array.
    width : float
        The full width (largest side length) of the box.
    """

    # collect all positions across particle types
    x_vals, y_vals, z_vals = [], [], []

    for (ptype, field) in data_all.keys():
        
        if field == "particle_position_x":
            x_vals.append(data_all[(ptype, field)])
        elif field == "particle_position_y":
            y_vals.append(data_all[(ptype, field)])
        elif field == "particle_position_z":
            z_vals.append(data_all[(ptype, field)])

    if not x_vals or not y_vals or not z_vals:
        raise ValueError("No particle position fields found in data_all.")

    # concatenate all particle positions
    x = np.concatenate(x_vals)
    y = np.concatenate(y_vals)
    z = np.concatenate(z_vals)

    # global min/max for each axis
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    zmin, zmax = np.min(z), np.max(z)

    print("data shape: ",x.shape,xmin, xmax, ymin, ymax, zmin, zmax)

   

   # find farthest extent from center (0,0,0)
    width = max(np.max(np.abs(x)), np.max(np.abs(y)), np.max(np.abs(z)))

   
    # add margin
    width *= (1 + margin)

    print(f"The width is : +- {width}")
     
    bbox = np.array([
        [-width, width],
        [-width, width],
        [-width, width]
    ])
    
    
    if verbose:
        print(f"✅ Bounding box (in {length_unit}):\n{bbox}")
        print(f"→ Box width: {width:.4f} {length_unit}")
        
    return bbox, width


# Function to make fits file containing the projected positions of the clusters
# called from savefits function
def make_fits(comp, proj, ds, width, res, kind):
    """
    Create a particle projection for a specific component (star, dm, gas)
    and export it to FITS format.
    """
    # field_name = f"{kind}_{comp}_density_{proj}"
    
    # ParticleProjectionPlot for a single type
    prjpd_plot = yt.ParticleProjectionPlot(
        ds,
        proj,
        (comp, "particle_mass"),       # component-specific
        weight_field=(comp, "particle_mass"),
        width=(2 * width, "Mpc")
    )

    prjpd_plot.set_buff_size((res, res))

    # Access the FRB for this specific component
    frb = prjpd_plot.frb[(comp, "particle_mass")]

    return frb.v

# helper to build datasets rotated to view_dir
def make_ds_for_view(data_dict=None, view_dir='x', bbox=None):
        # rotate particle positions so view_dir aligns with x-axis
        posx = data_dict[("star","particle_position_x")]
        posy = data_dict[("star","particle_position_y")]
        posz = data_dict[("star","particle_position_z")]
        positions = np.vstack([posx, posy, posz]).T
        rot = rotate_vectors(positions, view_dir)
        data_rot = data_dict.copy()
        data_rot[("star","particle_position_x")] = rot[:,0]
        data_rot[("star","particle_position_y")] = rot[:,1]
        data_rot[("star","particle_position_z")] = rot[:,2]
        # also rotate dm and gas if present
        for comp in ['dm','gas']:
            keyx = (comp,'particle_position_x')
            if keyx in data_dict:
                px = data_dict[keyx]
                py = data_dict[(comp,'particle_position_y')]
                pz = data_dict[(comp,'particle_position_z')]
                ppos = np.vstack([px,py,pz]).T
                prot = rotate_vectors(ppos, view_dir)
                data_rot[(comp,'particle_position_x')] = prot[:,0]
                data_rot[(comp,'particle_position_y')] = prot[:,1]
                data_rot[(comp,'particle_position_z')] = prot[:,2]
        
        # print number of data_rot particles

        if bbox is None:
            bbox,width = compute_bounding_box(data_rot, verbose=False)
            ds = yt.load_particles(data_rot, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
            return ds,width, bbox
        
        else:
            ds = yt.load_particles(data_rot, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
            return ds
            
# Utility: rotate a set of 3D particle positions so that the new viewing
# axis (the line-of-sight) aligns with the +x axis used by yt projections.
def rotate_vectors(positions, view_dir):
    """
    Rotate positions so that view_dir (3,) becomes the x-axis.
    positions: (N,3) array
    view_dir: iterable length 3
    Returns rotated positions (N,3).
    """
    v = np.array(view_dir, dtype=float)
    if np.allclose(v, 0):
        raise ValueError("view_dir must be non-zero")
    v = v / np.linalg.norm(v)

    # target x-axis
    t = np.array([1.0, 0.0, 0.0])
    # if v is already the x-axis, no rotation needed
    if np.allclose(v, t):
        return positions.copy()

    # rotation axis = cross(v, t)
    axis = np.cross(v, t)
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-12:
        # v is parallel or anti-parallel to t
        # if anti-parallel, rotate 180 degrees about y-axis
        if np.dot(v, t) < 0:
            R = np.array([[-1,0,0],[0,1,0],[0,0,-1]])
            return positions.dot(R.T)
        else:
            return positions.copy()

    axis = axis / axis_norm
    angle = np.arccos(np.clip(np.dot(v, t), -1.0, 1.0))

    # Rodrigues' rotation formula
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K.dot(K))

    return positions.dot(R.T)
