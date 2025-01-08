

def sersic_profile(r, I0, R0, n):
    """ Sérsic profile function. """
    return I0 * np.exp(-2 * (np.abs(r / R0) ** (1/np.abs(n)) - 1))

def fit_sersic_profile(r, I):
    """ Fit a Sérsic profile to the data and return the parameters. """
    # Initial guess for parameters
    initial_guess = [np.max(I), np.median(r[r > 0]), 2.0]  # Start with a slightly larger n
    
    # Set bounds to ensure R0 and n stay positive
    bounds = (0, [np.inf, np.inf, 10])  # Example: allow n to vary up to 10
    
    try:
        # Fit the Sérsic profile to the data
        popt, _ = curve_fit(sersic_profile, r, I, p0=initial_guess, bounds=bounds, maxfev=10000)
        I0, R0, n = popt
        # print(f"Fitted Sérsic index (n): {n:.2f}")
        return I0, R0, n
    except RuntimeError as e:
        # print(f"Optimal parameters not found: {e}")
        return None, None, None
    

def project_3d_to_2d(x, y, z):
    """ Project 3D particle data to 2D by summing along z-axis. """
    r = np.sqrt(x**2 + y**2)
    hist, edges = np.histogram(r, bins=100, range=(0, np.max(r)))
    bin_centers = (edges[:-1] + edges[1:]) / 2
    return bin_centers, hist




def get_morphology(self):
       
        
        if self.star_mass[:].shape[0]>0 and self.galID !='ICL':
            # Calculate the center of mass (COM) for the chosen particle type
            com = np.average(self.star_pos, weights=self.star_mass, axis=0)

            # Calculate the distance from COM for each particle
            pos_gal = self.star_pos[:] - com
            
            

            x = pos_gal[:, 0]
            y = pos_gal[:, 1]
            z = pos_gal[:, 2]

            r, intensity = project_3d_to_2d(x, y, z)
            

            # Fit the Sérsic profile
            I0, R0, n_ser = fit_sersic_profile(r, intensity)

            # Print the Sérsic index
            return n_ser
