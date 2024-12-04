import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# Set a random seed for reproducibility
seed_value = 42
np.random.seed(seed_value)

# Simulation parameters
n_photons = 500000         # Increased number of photons for better statistics
voxel_size = 0.1          # Size of each voxel in mm
weight_threshold = 1e-9   # Threshold for photon weight
survival_chance = 0.1     # Survival probability in Russian Roulette

# Emitter properties
emitter_width = 0.2      # Reduced width in mm
emitter_height = 0.2     # Reduced height in mm
emitter_position = (-3.22, 0.0)  # 3.22 mm to the left of the center

# Detector properties
detector_width = 1.29       # Reduced width in mm
detector_height = 2.69      # Reduced height in mm
detector_position = (0.0, 0.0)  # Center of the volume

# Fixed grid size (dimensions of the tissue volume)
grid_size_x = 100  # Number of voxels in X direction
grid_size_y = 100  # Number of voxels in Y direction
grid_size_z = 162  # Number of voxels in Z direction

grid_size = (grid_size_x, grid_size_y, grid_size_z)

# Define six layers with their optical properties and thicknesses
layers = [
    {'mu_a': 0.12, 'mu_s': 5.3,  'g': 0.715, 'thickness': 1.6, 'n': 1.36},  # Skin
    {'mu_a': 0.118, 'mu_s': 26.8, 'g': 0.85, 'thickness': 0.3, 'n': 1.37},  # SCF
    {'mu_a': 0.13, 'mu_s': 1.18, 'g': 0.732, 'thickness': 1.15, 'n': 1.37}, # Platysma Muscle
    {'mu_a': 0.118, 'mu_s': 26.8, 'g': 0.85, 'thickness': 0.3, 'n': 1.37},  # DCF
    {'mu_a': 0.13, 'mu_s': 1.18, 'g': 0.732, 'thickness': 5.8, 'n': 1.37},  # SCM
    {'mu_a': 0.118, 'mu_s': 26.8, 'g': 0.85, 'thickness': 7, 'n': 1.37},   # PF
]

# Refractive indices of surrounding media
n_air = 1.0
n_bottom = 1.0

# Compute total tissue depth and layer interfaces
def compute_total_tissue_depth(layers):
    return sum(layer['thickness'] for layer in layers)

total_tissue_depth = compute_total_tissue_depth(layers)
print(f"Total tissue depth: {total_tissue_depth} mm")

def compute_layer_interfaces(layers):
    interfaces = []
    current_depth = 0.0
    for layer in layers[:-1]:
        current_depth += layer['thickness']
        interfaces.append(current_depth)
    return interfaces

layer_interfaces = compute_layer_interfaces(layers)

# Initialize counters and grids
interface_counts = [0 for _ in layer_interfaces]
photons_exited_top = 0
photons_exited_bottom = 0
photons_exited_lateral = 0
photons_detected = 0  # Photons detected by the detector

heatmap_resolution = (grid_size[0], grid_size[1])
interface_heatmaps = [np.zeros(heatmap_resolution) for _ in layer_interfaces]
tissue_grid = np.zeros(grid_size)

def launch_photon_rectangular():
    x = emitter_position[0] + (np.random.random() - 0.5) * emitter_width
    y = emitter_position[1] + (np.random.random() - 0.5) * emitter_height
    z = 0.0  # At the bottom layer
    return x, y, z

def hop_step(mu_t):
    RND = np.random.random()
    step_size = -np.log(RND) / mu_t
    return step_size

def drop_weight(weight, mu_a, mu_t, step_size):
    absorbed = weight * (mu_a / mu_t) * (1 - np.exp(-mu_t * step_size))
    weight -= absorbed
    return weight, absorbed

def spin_photon(ux, uy, uz, g):
    RND = np.random.random()
    if g == 0:
        costheta = 2 * RND - 1
    else:
        temp = (1 - g**2) / (1 - g + 2 * g * RND)
        costheta = (1 + g**2 - temp**2) / (2 * g)
    sintheta = np.sqrt(1 - costheta**2)
    phi = 2 * np.pi * np.random.random()
    cosphi = np.cos(phi)
    sinphi = np.sin(phi)

    # Rotate the photon direction
    if abs(uz) > 0.99999:
        ux_new = sintheta * cosphi
        uy_new = sintheta * sinphi
        uz_new = costheta * np.sign(uz)
    else:
        denom = np.sqrt(1 - uz**2) + 1e-12  # Added epsilon to prevent division by zero
        ux_new = (sintheta * (ux * uz * cosphi - uy * sinphi)) / denom + ux * costheta
        uy_new = (sintheta * (uy * uz * cosphi + ux * sinphi)) / denom + uy * costheta
        uz_new = -sintheta * cosphi * denom + uz * costheta

    # Normalize the new direction cosines
    norm = np.sqrt(ux_new**2 + uy_new**2 + uz_new**2) + 1e-12  # Added epsilon to prevent division by zero
    ux_new /= norm
    uy_new /= norm
    uz_new /= norm

    return ux_new, uy_new, uz_new

def update_energy_grid(grid, x, y, z, absorbed):
    ix = int(np.floor((x + (grid.shape[0] * voxel_size) / 2) / voxel_size))
    iy = int(np.floor((y + (grid.shape[1] * voxel_size) / 2) / voxel_size))
    iz = int(np.floor(z / voxel_size))

    # Boundary checks with clipping
    ix = np.clip(ix, 0, grid.shape[0]-1)
    iy = np.clip(iy, 0, grid.shape[1]-1)
    iz = np.clip(iz, 0, grid.shape[2]-1)

    grid[ix, iy, iz] += absorbed

def get_tissue_properties(z):
    current_depth = 0.0
    for layer in layers:
        if current_depth <= z < current_depth + layer['thickness']:
            return layer
        current_depth += layer['thickness']
    return layers[-1]

def check_layer_crossing(z_prev, z_new):
    crossed_indices = []
    for idx, interface_z in enumerate(layer_interfaces):
        # Check if the photon crosses the interface
        if (z_prev - interface_z) * (z_new - interface_z) <= 0 and z_prev != z_new:
            crossed_indices.append(idx)
    return crossed_indices

def handle_interface_crossing(ux, uy, uz, n1, n2):
    cos_theta_i = abs(uz)
    sin_theta_i = np.sqrt(1 - cos_theta_i**2)
    sin_theta_t = (n1 / n2) * sin_theta_i

    if sin_theta_t > 1.0:
        return ux, uy, -uz, True  # Total internal reflection
    else:
        cos_theta_t = np.sqrt(1 - sin_theta_t**2)
        Rs = ((n1 * cos_theta_i - n2 * cos_theta_t) / (n1 * cos_theta_i + n2 * cos_theta_t)) ** 2
        Rp = ((n1 * cos_theta_t - n2 * cos_theta_i) / (n1 * cos_theta_t + n2 * cos_theta_i)) ** 2
        R = (Rs + Rp) / 2
        if np.random.random() < R:
            return ux, uy, -uz, True  # Photon is reflected
        else:
            ux_t = ux * (n1 / n2)
            uy_t = uy * (n1 / n2)
            uz_t = np.sign(uz) * cos_theta_t
            norm = np.sqrt(ux_t**2 + uy_t**2 + uz_t**2) + 1e-12  # Added epsilon
            ux_t /= norm
            uy_t /= norm
            uz_t /= norm
            return ux_t, uy_t, uz_t, False

def monte_carlo_simulation():
    global photons_exited_top, photons_exited_bottom, photons_exited_lateral, photons_detected

    for _ in tqdm(range(n_photons), desc="Simulating Photons", unit="photons"):
        x, y, z = launch_photon_rectangular()
        ux, uy, uz = 0.0, 0.0, 1.0  # Moving upwards
        weight = 1.0

        while True:
            props = get_tissue_properties(z)
            mu_a = props['mu_a']
            mu_s = props['mu_s']
            mu_t = mu_a + mu_s
            g = props['g']
            n_current = props['n']

            step_size = hop_step(mu_t)
            x_prev, y_prev, z_prev = x, y, z
            x += step_size * ux
            y += step_size * uy
            z += step_size * uz

            # Check for layer crossing
            crossed_layers = check_layer_crossing(z_prev, z)
            for idx in crossed_layers:
                interface_counts[idx] += 1
                # Update heatmap
                ix = int(np.floor((x + (heatmap_resolution[0] * voxel_size) / 2) / voxel_size))
                iy = int(np.floor((y + (heatmap_resolution[1] * voxel_size) / 2) / voxel_size))
                if 0 <= ix < heatmap_resolution[0] and 0 <= iy < heatmap_resolution[1]:
                    interface_heatmaps[idx][ix, iy] += 1

                # Handle interface crossing
                if z_prev < z:  # Moving upwards
                    layer_below = get_tissue_properties(layer_interfaces[idx] - 1e-6)
                    layer_above = get_tissue_properties(layer_interfaces[idx] + 1e-6)
                else:  # Moving downwards
                    layer_below = get_tissue_properties(layer_interfaces[idx] + 1e-6)
                    layer_above = get_tissue_properties(layer_interfaces[idx] - 1e-6)
                n1 = layer_below['n']
                n2 = layer_above['n']
                ux, uy, uz, is_reflected = handle_interface_crossing(ux, uy, uz, n1, n2)
                if is_reflected:
                    z = z_prev  # Reflect back to previous position
                    break
                else:
                    z = z  # Continue with updated direction

            # Boundary checks
            if z >= total_tissue_depth:
                # Photon exits the tissue at the top
                ux, uy, uz, is_reflected = handle_interface_crossing(ux, uy, uz, n_current, n_air)
                if is_reflected:
                    z = total_tissue_depth - 1e-6  # Reflect back into tissue
                else:
                    photons_exited_top += 1
                    break
            elif z < 0:
                # Photon exits the tissue at the bottom
                ux, uy, uz, is_reflected = handle_interface_crossing(ux, uy, uz, n_current, n_bottom)
                if is_reflected:
                    z = 1e-6  # Reflect back into tissue
                else:
                    photons_exited_bottom += 1
                    # Check if photon hits the detector
                    if (detector_position[0] - detector_width / 2 <= x <= detector_position[0] + detector_width / 2 and
                        detector_position[1] - detector_height / 2 <= y <= detector_position[1] + detector_height / 2):
                        photons_detected += 1
                    break
            elif (abs(x) > (grid_size[0] * voxel_size) / 2 or
                  abs(y) > (grid_size[1] * voxel_size) / 2):
                photons_exited_lateral += 1
                break

            # Update weight due to absorption
            weight, absorbed = drop_weight(weight, mu_a, mu_t, step_size)
            update_energy_grid(tissue_grid, x, y, z, absorbed)

            # Photon scattering
            ux, uy, uz = spin_photon(ux, uy, uz, g)

            # Russian Roulette
            if weight < weight_threshold:
                if np.random.random() < survival_chance:
                    weight /= survival_chance
                else:
                    break

def visualize_depth_sum(grid):
    depth_sum = np.sum(grid, axis=2)
    extent = [
        -grid_size[0] * voxel_size / 2,
        grid_size[0] * voxel_size / 2,
        -grid_size[1] * voxel_size / 2,
        grid_size[1] * voxel_size / 2,
    ]
    plt.figure(figsize=(8, 6))
    plt.imshow(
        depth_sum.T,
        cmap='hot',
        extent=extent,
        origin='lower',
        interpolation='nearest',
        aspect='auto',
        norm=plt.Normalize(vmin=0, vmax=np.max(depth_sum))
    )
    plt.colorbar(label='Total Absorbed Energy')
    plt.title('Summed Absorption Over Depth')
    plt.xlabel('X Position (mm)')
    plt.ylabel('Y Position (mm)')

    # Plot emitter and detector positions
    plt.scatter(emitter_position[0], emitter_position[1], color='blue', marker='o', label='Emitter')
    plt.scatter(detector_position[0], detector_position[1], color='green', marker='s', label='Detector')

    plt.legend()
    plt.tight_layout()
    plt.show()

def visualize_interface_heatmaps(interface_heatmaps, voxel_size, interface_depths):
    for idx, heatmap in enumerate(interface_heatmaps):
        extent = [
            -grid_size[0] * voxel_size / 2,
            grid_size[0] * voxel_size / 2,
            -grid_size[1] * voxel_size / 2,
            grid_size[1] * voxel_size / 2,
        ]
        plt.figure(figsize=(8, 6))
        plt.imshow(
            heatmap.T,
            cmap='hot',
            extent=extent,
            origin='lower',
            interpolation='nearest',
            aspect='auto',
            norm=plt.Normalize(vmin=0, vmax=np.max(heatmap))
        )
        plt.colorbar(label='Photon Crossing Counts')
        plt.title(f'Photon Crossings at Interface {idx+1} (Depth {interface_depths[idx]:.2f} mm)')
        plt.xlabel('X Position (mm)')
        plt.ylabel('Y Position (mm)')

        # Optionally, mark emitter and detector positions
        plt.scatter(emitter_position[0], emitter_position[1], color='blue', marker='o', label='Emitter')
        plt.scatter(detector_position[0], detector_position[1], color='green', marker='s', label='Detector')
        plt.legend()

        plt.tight_layout()
        plt.show()

# Run the simulation
monte_carlo_simulation()

# Display results
print(f"\nSimulation Results:")
print(f"Photons detected by the detector: {photons_detected}")
print(f"Photons exited at the top surface: {photons_exited_top}")
print(f"Photons exited at the bottom: {photons_exited_bottom}")
print(f"Photons exited laterally: {photons_exited_lateral}")

# Visualize the results
visualize_depth_sum(tissue_grid)
visualize_interface_heatmaps(interface_heatmaps, voxel_size, layer_interfaces)

# Display the number of photons crossing each interface
print("\nPhotons Crossing Each Interface:")
for idx, count in enumerate(interface_counts):
    interface_depth = layer_interfaces[idx]
    print(f"Interface at depth {interface_depth:.2f} mm: {count} photons crossed")
