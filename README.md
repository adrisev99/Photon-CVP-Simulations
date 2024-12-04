# Photon-CVP-Simulations

This repository contains Python code developed for simulating photon propagation in infant neck tissue using Monte Carlo methods. The work is part of the study detailed in the thesis **"Photon Propagation in Infant Neck Tissue: A Monte Carlo and Ray-Tracing Simulation Study"**, which explores the optimization of non-invasive central venous pressure (CVP) monitoring devices.

The code is designed to model light-tissue interactions in a six-layer anatomical model of infant neck tissue. The study employs these simulations to analyze photon scattering, absorption, and detection to enhance the performance of optical biomonitoring devices, particularly focusing on signal optimization for capturing photon interactions with the internal jugular vein (IJV).

---

## Project Background

Non-invasive CVP monitoring is a critical aspect of managing pediatric and neonatal hemodynamic conditions, where traditional invasive techniques pose significant risks. Optical methods, such as photoplethysmography (PPG), offer a safe alternative, but the complex interactions of light with multi-layered biological tissues challenge signal detection. 

This project specifically investigates photon behavior in a 10 mm × 10 mm × 16.15 mm tissue model using Monte Carlo simulations. The simulations incorporate anatomical and optical properties of skin, fascia, muscle, and other tissues to evaluate photon penetration, interaction, and capture efficiency at a wavelength of 630 nm. 

The insights derived from these simulations aim to inform the design of medical devices that optimize photon detection from deep tissue structures like the IJV, critical for accurate CVP monitoring.

---

## Features of the Code

1. **Monte Carlo Simulation Framework**:
   - Tracks individual photons as they propagate, scatter, and are absorbed within a 3D voxelated tissue model.
   - Includes photon boundary interactions, interface crossings, and detection at the tissue boundaries.

2. **Detailed Anatomical Modeling**:
   - Represents the infant neck tissue as a six-layer model: skin, superficial cervical fascia (SCF), platysma muscle, deep cervical fascia (DCF), sternocleidomastoid muscle (SCM), and prevertebral fascia (PF).
   - Each layer is assigned specific optical properties (absorption, scattering, anisotropy, and refractive index).

3. **Photon Interaction Metrics**:
   - Calculates the number of photons crossing layer interfaces, absorbed energy distribution, and detection efficiency at the bottom detector.
   - Tracks photons exiting at the tissue's top, bottom, and lateral boundaries.

4. **Emitter-Detector Configuration**:
   - Models a rectangular emitter-detector setup with customizable geometry and positions.
   - Evaluates the effect of emitter-detector spacing on signal capture efficiency and signal-to-noise ratio (SNR).

5. **Energy Deposition and Heatmap Generation**:
   - Produces 2D and 3D heatmaps to visualize energy absorbed by tissue layers.
   - Provides detailed spatial data on photon crossing counts at each layer interface.

---

## Code Highlights

### Monte Carlo Methodology
The Monte Carlo algorithm simulates photon-tissue interactions, including:
- Randomized photon step size and scattering direction based on optical properties.
- Interaction probabilities computed using tissue-specific absorption and scattering coefficients.
- Boundary interactions handled via Fresnel equations for reflection and refraction.

### Voxelated Tissue Model
The 3D tissue model is discretized into fine-resolution voxels (0.1 mm). This approach:
- Enables high spatial accuracy in photon path tracking.
- Supports detailed analysis of photon absorption and scattering within specific regions of tissue.

### Emitter and Detector Configuration
- Emitter: Positioned at 3.22 mm offset from the center with a rectangular emitting area.
- Detector: Located centrally on the bottom boundary with a customizable rectangular area.
- Both components are adjustable, allowing studies on their spatial influence on photon detection.

### Simulation Results
- Photon absorption and scattering are quantified within each tissue layer.
- Photon detection rates and interaction efficiencies with the IJV are calculated.
- Visualization tools highlight energy deposition and layer interface crossings.

---

## Applications

The simulations and resulting insights from this code are directly applicable to:
- Designing and optimizing optical biomonitoring devices for neonatal and pediatric CVP monitoring.
- Understanding light transport in multi-layered biological tissues for medical imaging.
- Developing non-invasive diagnostic tools with improved photon penetration and capture efficiency.

---

## Future Extensions

The current code is a robust foundation for further enhancements:
- Incorporating additional anatomical features, such as dynamic blood vessel modeling.
- Expanding to simulate a wider range of wavelengths for improved tissue penetration.
- Implementing advanced statistical techniques for noise reduction in simulated detection.

This repository provides a comprehensive framework for studying photon-tissue interactions, supporting the development of innovative non-invasive monitoring technologies.
