# Host Galaxy Association

## Overview

Identifying the host galaxy of a transient event is a fundamental step in
transient astronomy. The host provides redshift (distance), environment,
and progenitor constraints. prost uses a two-stage approach:

1. **Geometric ranking** via directional light radius (DLR)
2. **Bayesian scoring** with offset, redshift, and magnitude likelihoods

## Directional Light Radius (DLR)

The DLR measures the angular separation between a transient and a galaxy
center, normalized by the galaxy's effective radius *in the direction of
the transient*. This accounts for galaxy morphology — a transient along
the major axis of an edge-on galaxy may be physically closer than one
at the same angular separation along the minor axis.

### Geometry

Given a galaxy at position \((\alpha_g, \delta_g)\) with elliptical
profile characterized by semi-major axis \(a\), semi-minor axis \(b\),
and position angle \(\theta_{PA}\):

1. **Tangent-plane offsets** (arcsec):

    \[
    \Delta x = (\alpha_t - \alpha_g) \cos(\delta_g) \times 3600
    \]
    \[
    \Delta y = (\delta_t - \delta_g) \times 3600
    \]

2. **Rotate into galaxy frame**:

    \[
    x_{\text{maj}} = \Delta x \cos\theta_{PA} + \Delta y \sin\theta_{PA}
    \]
    \[
    y_{\text{min}} = -\Delta x \sin\theta_{PA} + \Delta y \cos\theta_{PA}
    \]

3. **Directional radius** — the galaxy's effective radius at the angle
   toward the transient:

    \[
    r(\phi) = \frac{ab}{\sqrt{(b\cos\phi)^2 + (a\sin\phi)^2}}
    \]

    where \(\phi = \arctan(y_{\text{min}} / x_{\text{maj}})\).

4. **Fractional offset**:

    \[
    d_{FR} = \frac{\text{separation}}{r(\phi)}
    \]

A fractional offset \(d_{FR} < 1\) means the transient is within the
galaxy's effective radius along that direction.

## Bayesian Posterior

### Offset Likelihood

The fractional offset distribution for true host galaxies follows a
Gamma distribution with shape parameter \(a = 0.75\)
(Gagliano et al. 2021):

\[
\mathcal{L}(d_{FR}) = \frac{d_{FR}^{-0.25} \, e^{-d_{FR}}}{\Gamma(0.75)}
\]

where \(\Gamma(0.75) \approx 1.2254\). This distribution peaks near
zero and has a heavy tail, reflecting that most transients occur close
to their host center but some are found at large offsets.

### Redshift Likelihood

When both transient and galaxy have measured redshifts:

\[
\mathcal{L}(z) = \exp\left(-\frac{(z_g - z_t)^2}{2(\sigma_{z_g}^2 + \sigma_{z_t}^2)}\right)
\]

This Gaussian likelihood strongly downweights galaxies at discrepant
redshifts, breaking positional degeneracies.

### Prior

A uniform prior on fractional offset:

\[
\pi(d_{FR}) = \frac{1}{d_{FR,\max}} \quad \text{for } d_{FR} \in [0, d_{FR,\max}]
\]

### Null Hypotheses

The posterior includes probability mass for three alternative scenarios:

| Hypothesis | Default | Description |
|-----------|---------|-------------|
| \(P_{\text{outside}}\) | 0.01 | True host beyond search radius |
| \(P_{\text{unobserved}}\) | 0.01 | True host too faint for catalog |
| \(P_{\text{hostless}}\) | 0.005 | Transient is genuinely hostless |

### Normalization

The posterior for candidate \(i\) is:

\[
P_i = \frac{\mathcal{L}_i \, \pi_i}{\sum_j \mathcal{L}_j \, \pi_j + P_{\text{outside}} + P_{\text{unobserved}} + P_{\text{hostless}}}
\]

ensuring \(\sum_i P_i + P_{\text{none}} = 1\).

## Tractor Shape Parameters

The Legacy Survey Tractor catalog parameterizes galaxy shapes as:

- `shape_r` — effective radius (arcsec)
- `shape_e1`, `shape_e2` — ellipticity components

Conversion to ellipse parameters:

\[
e = \sqrt{e_1^2 + e_2^2}, \quad
q = \frac{1-e}{1+e}, \quad
\theta_{PA} = \frac{1}{2} \arctan\left(\frac{e_2}{e_1}\right)
\]

\[
a = \texttt{shape\_r}, \quad b = \max(a \cdot q, \, b_{\min})
\]

## Source Extraction

When catalog morphology is unavailable, prost can extract sources
directly from image cutouts using a SExtractor-inspired pipeline:

1. **Background estimation** — The image is divided into mesh cells
   (default 64 px). Sigma-clipped statistics are computed per cell,
   then bilinearly interpolated to full resolution.

2. **Thresholding** — Pixels exceeding `thresh_sigma * rms` above
   the local background are flagged.

3. **Connected-component labeling** — Flagged pixels are grouped
   using 8-connectivity (union-find algorithm).

4. **Shape measurement** — Flux-weighted second moments yield the
   centroid, semi-axes, and position angle for each source.

## Sersic Profile Fitting

For galaxies where accurate morphology is required (or catalog shapes
are unreliable), prost fits a 2D Sersic profile:

\[
I(r) = I_e \exp\left(-b_n \left[\left(\frac{r}{r_e}\right)^{1/n} - 1\right]\right)
\]

where \(b_n \approx 1.9992n - 0.3271\) (Ciotti & Bertin 1999), and
\(r\) is the elliptical radius:

\[
r = \sqrt{x_{\text{rot}}^2 + (y_{\text{rot}} / q)^2}
\]

The fit uses Levenberg-Marquardt optimization with 7 free parameters:
centroid (\(x_c, y_c\)), surface brightness (\(I_e\)), effective
radius (\(r_e\)), Sersic index (\(n\)), axis ratio (\(q\)), and
position angle (\(\theta\)). The Jacobian is computed numerically.

The fitted parameters map to `GalaxyCandidate` shape via:

- `a_arcsec` = \(r_e / q\) (semi-major of the half-light ellipse)
- `b_arcsec` = \(r_e\) (semi-minor)
- `pa_deg` = \(\theta\)

## References

1. Gagliano, A. et al. (2021). "PROST: A Probabilistic Method for
   Associating Transients to Their Host Galaxies."
   [*ApJ*, 908, 170](https://doi.org/10.3847/1538-4357/abd02b).

2. Gagliano, A. et al. (2021). "GHOST: Using Host Galaxy Properties
   to Inform the Classification of Supernovae."
   [*ApJL*, 911, L14](https://doi.org/10.3847/2041-8213/abf1e4).

3. Sullivan, M. et al. (2006). "Rates and Properties of Type Ia Supernovae
   as a Function of Mass and Star Formation in Their Host Galaxies."
   [*ApJ*, 648, 868](https://doi.org/10.1086/506137).

4. Ciotti, L. & Bertin, G. (1999). "Analytical properties of the
   R^{1/m} law." [*A&A*, 352, 447](https://ui.adsabs.harvard.edu/abs/1999A%26A...352..447C).
