Simulate Ronchigrams on the fly using WebAssembly (Fast) or Javascript (Slow). See it in action at http://ronchigram.com
WebAssembly mode uses KissFFT (https://github.com/mborgerding/kissfft) and JS mode uses mrquincle's JS FFT (https://gist.github.com/mrquincle/b11fff96209c9d1396b0) and math.js (https://github.com/josdejong/mathjs)

## Ronchigram Simulation Algorithm

A Ronchigram is an in-focus diffraction pattern formed by a STEM (Scanning Transmission Electron Microscope) when a broad, quasi-parallel electron beam passes through a thin amorphous sample. It encodes the aberrations of the electron optical system and is widely used for aberration diagnosis and correction.

### Step 1 — Electron Wavelength

The relativistically corrected de Broglie wavelength of the electron beam is:

$$\lambda = \frac{h}{\sqrt{2\,m_0\,eV\!\left(1 + \dfrac{eV}{2\,m_0 c^2}\right)}}$$

which simplifies (in practical units) to:

$$\lambda \;[\text{m}] = \frac{12.3986}{\sqrt{(2\times511+E)\times E}} \times 10^{-10}$$

where $E$ is the beam energy in keV, $m_0 c^2 = 511$ keV is the electron rest energy, $h$ is Planck's constant, and $e$ is the elementary charge. This is implemented in `calculateLambda()`.

### Step 2 — Polar Mesh and Objective Aperture

A 2-D grid is created in angular (reciprocal) space. Each pixel $(i,\,j)$ is mapped to Cartesian angular coordinates:

$$x_i = \frac{i - N/2}{N/2}\,\alpha_{\max}, \qquad y_j = \frac{j - N/2}{N/2}\,\alpha_{\max}$$

and then converted to polar form:

$$\alpha_{ij} = \sqrt{x_i^2 + y_j^2}, \qquad \Phi_{ij} = \mathrm{atan2}(y_j,\, x_i)$$

where $\alpha$ is the semi-angle (mrad, same units as $\alpha_{\max}$) and $\Phi$ is the azimuthal angle. The objective aperture mask is:

$$\mathrm{oapp}(\alpha) = \begin{cases} 1 & \alpha \le \alpha_{\mathrm{ap}} \\ 0 & \alpha > \alpha_{\mathrm{ap}} \end{cases}$$

This is implemented in `polarMeshnOapp()`.

### Step 3 — Aberration Phase Function χ₀

The geometrical aberration function in the Krivanek notation is:

$$\chi_0(\alpha,\Phi) = \frac{2\pi}{\lambda} \sum_{k} \frac{C_{n_k,m_k}}{n_k+1}\;\alpha^{n_k+1} \cos\!\bigl[m_k\,(\Phi - \Phi_{n_k,m_k})\bigr]$$

where $C_{n,m}$ (meters) is the aberration coefficient magnitude, $\Phi_{n,m}$ (radians) is its orientation, and the indices $(n,\,m)$ label the order and symmetry of each term.

The 14 aberration terms supported (up to 5th order) are:

| k  | n | m | Name                      | Krivanek symbol |
|----|---|---|---------------------------|-----------------|
| 0  | 1 | 0 | Defocus                   | C₁₀             |
| 1  | 1 | 2 | 2-fold astigmatism        | C₁₂             |
| 2  | 2 | 1 | Axial coma                | C₂₁             |
| 3  | 2 | 3 | 3-fold astigmatism        | C₂₃             |
| 4  | 3 | 0 | 3rd-order spherical       | C₃₀             |
| 5  | 3 | 2 | 3rd-order axial star      | C₃₂             |
| 6  | 3 | 4 | 4-fold astigmatism        | C₃₄             |
| 7  | 4 | 1 | 4th-order axial coma      | C₄₁             |
| 8  | 4 | 3 | 3-lobe aberration         | C₄₃             |
| 9  | 4 | 5 | 5-fold astigmatism        | C₄₅             |
| 10 | 5 | 0 | 5th-order spherical       | C₅₀             |
| 11 | 5 | 2 | 5th-order axial star      | C₅₂             |
| 12 | 5 | 4 | 5th-order rosette         | C₅₄             |
| 13 | 5 | 6 | 6-fold astigmatism        | C₅₆             |

This is implemented in `calculateChi0()`.

### Step 4 — Complex Pupil Function χ

The complex pupil function (exit-wave in aperture space) is:

$$\chi(\alpha,\Phi) = e^{-i\,\chi_0(\alpha,\Phi)}$$

The modulus is identically 1; only the phase encodes the lens aberrations. This is implemented in `calculateChi()`.

### Step 5 — Sample and Transmission Function

An amorphous sample is modelled as a random noisy grating. The sample transmission function under the phase-object approximation is:

$$t(\mathbf{r}) = \exp\!\left(-i\,\frac{\pi}{4}\,\sigma\,V(\mathbf{r})\right)$$

where $V(\mathbf{r})$ is the random projected potential and $\sigma$ is the relativistic interaction parameter:

$$\sigma = \frac{2\pi}{\lambda E} \cdot \frac{m_0 c^2 + eV}{2\,m_0 c^2 + eV}$$

In the code $\sigma$ is normalised to its value at 300 kV so that the default coefficients are on a human-readable scale. This is implemented in `generateSample()`, `calculateInteractionParam()`, and `generateTransmissionFn()`.

### Step 6 — Ronchigram (Diffraction Intensity)

The Ronchigram intensity is computed as:

$$I(\mathbf{q}) = \Bigl|\mathcal{F}\bigl[t(\mathbf{r})\cdot\mathcal{F}[\chi(\boldsymbol{\alpha})]\bigr](\mathbf{q})\Bigr|^{2} \times \mathrm{oapp}(\mathbf{q})$$

That is:
1. $\mathcal{F}[\chi(\boldsymbol{\alpha})]$ — Fourier-transform the pupil function into real space.
2. $t(\mathbf{r})\cdot(\ldots)$ — Multiply by the sample transmission function.
3. $\mathcal{F}[\cdots](\mathbf{q})$ — Fourier-transform back to reciprocal space.
4. $|\cdots|^{2}\times\mathrm{oapp}$ — Take the squared modulus and apply the objective aperture mask.

This is implemented in `calcDiffract()`.

### Step 7 — π/4 Phase Limit Map

Pixels where the wave-front error exceeds the Rayleigh quarter-wave criterion are masked to zero:

$$\mathrm{mask}(\alpha,\Phi) = \begin{cases} 255 & |\chi_0(\alpha,\Phi)| \le \dfrac{\pi}{4} \\[4pt] 0 & |\chi_0(\alpha,\Phi)| > \dfrac{\pi}{4} \end{cases}$$

The radius $r_{\max}$ of the largest circle that fits entirely inside the white (valid) region is then:

$$r_{\max} = \min\bigl\{\alpha \;:\; |\chi_0(\alpha,\Phi)| > \tfrac{\pi}{4}\bigr\}$$

This is implemented in `maskChi0()`.

### Step 8 — Electron Probe

The real-space probe intensity formed by the aberrated lens is:

$$P(\mathbf{r}) = \Bigl|\mathcal{F}\bigl[\chi(\boldsymbol{\alpha})\cdot A(\boldsymbol{\alpha})\bigr](\mathbf{r})\Bigr|^{2}$$

where $A(\boldsymbol{\alpha})$ is the circular aperture function. The result is fftshift-ed (so the probe centre appears at the image centre) and normalised to 0–255. This is implemented in `probeGeneration()`.

### Step 9 — Strehl Ratio

The Strehl ratio $S$ measures how close the aberrated probe is to the diffraction limit (ideal unaberrated lens):

$$S = \left(\frac{\max\!\left|\mathcal{F}[\chi\cdot A_s]\right|}{\max\!\left|\mathcal{F}[A_s]\right|}\right)^{2}$$

where $A_s$ is a circular aperture of radius $r_s$. A value $S = 1$ corresponds to a perfect lens; $S \ge 0.8$ is the Maréchal criterion for diffraction-limited performance.

The code also performs a bisection search to find the largest aperture semi-angle $r_s$ that still satisfies $S \ge 0.8$:

$$r_{0.8} = \arg\max_{r_s}\, r_s \quad \text{subject to} \quad S(r_s) \ge 0.8$$

This is implemented in `singleStrehl()` and `pointEightStrehlSearch()`.

### Data Flow Summary

```
keV  ──► λ
                    ┌─ chi0 = Σ C_{n,m} terms ──► χ = exp(-i χ0) ─┬─► Ronchigram = |FFT(t × FFT(χ))|² × oapp
α_max, obj_ap_r ──►│ polar mesh (α, Φ), oapp                       ├─► π/4 map  (|χ0| threshold)
aberration coefs ──┘                                                ├─► Electron probe = |FFT(χ × oapp)|²
                    sample ──► t(r) = exp(-i σ V)                  └─► Strehl ratio
```

### Source Files

| File | Description |
|------|-------------|
| `emsrc/calculateRonch.cpp` | Core C++ simulation compiled to WebAssembly via Emscripten |
| `emsrc/kiss_fft*.{c,h}` | KissFFT library used for 2-D FFTs in the C++ path |
| `ronch_calc.js` | Pure-JavaScript fallback implementation |
| `index.js` | Emscripten-generated WebAssembly loader |
| `emsrc/makefile` | Build instructions (`make all` with Emscripten) |
