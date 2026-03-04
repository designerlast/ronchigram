Simulate Ronchigrams on the fly using WebAssembly (Fast) or Javascript (Slow). See it in action at http://ronchigram.com
WebAssembly mode uses KissFFT (https://github.com/mborgerding/kissfft) and JS mode uses mrquincle's JS FFT (https://gist.github.com/mrquincle/b11fff96209c9d1396b0) and math.js (https://github.com/josdejong/mathjs)

## Ronchigram Simulation Algorithm

A Ronchigram is an in-focus diffraction pattern formed by a STEM (Scanning Transmission Electron Microscope) when a broad, quasi-parallel electron beam passes through a thin amorphous sample. It encodes the aberrations of the electron optical system and is widely used for aberration diagnosis and correction.

### Step 1 — Electron Wavelength

The relativistically corrected de Broglie wavelength of the electron beam is:

```
λ = 12.3986 / sqrt((2 × 511 + E) × E) × 10⁻¹⁰  [m]
```

where `E` is the beam energy in keV. This is implemented in `calculateLambda()`.

### Step 2 — Polar Mesh and Objective Aperture

A 2-D grid is created in angular (reciprocal) space, parameterised by semi-angle `α` (radians) and azimuthal angle `Φ`. Each pixel `(i, j)` maps to:

```
α = sqrt(x² + y²),   Φ = atan2(y, x)
```

where `x` and `y` run from −α_max to +α_max. A circular objective aperture mask `oapp` is set to 1 inside radius `obj_ap_r` and 0 outside. This is implemented in `polarMeshnOapp()`.

### Step 3 — Aberration Phase Function χ₀

The geometrical aberration function in the Krivanek notation is:

```
χ₀(α, Φ) = (2π/λ) × Σ_{k} C_{n,m} × α^(n+1) × cos(m × (Φ − Φ_{n,m})) / (n+1)
```

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

```
χ(α, Φ) = exp(−i × χ₀(α, Φ))
```

This is implemented in `calculateChi()`.

### Step 5 — Sample and Transmission Function

An amorphous sample is modelled as a random noisy grating. The sample transmission function is:

```
t(r) = exp(−i × (π/4) × σ × V(r))
```

where `V(r)` is the random projected potential and `σ` is the interaction parameter (relativistically corrected, normalised to 300 kV). This is implemented in `generateSample()`, `calculateInteractionParam()`, and `generateTransmissionFn()`.

### Step 6 — Ronchigram (Diffraction Intensity)

The Ronchigram intensity is computed as:

```
I(q) = |FFT( t(r) × FFT( χ(α) ) )|² × oapp(q)
```

That is:
1. Fourier-transform the pupil function χ into real space.
2. Multiply by the sample transmission function t(r).
3. Fourier-transform back to reciprocal space.
4. Take the squared modulus and apply the objective aperture mask.

This is implemented in `calcDiffract()`.

### Step 7 — π/4 Phase Limit Map

Pixels where `|χ₀| > π/4` are masked to zero (shown in white), revealing the isochronous zone where the wave front aberration stays within the Rayleigh quarter-wave criterion. The radius of the largest circle that fits inside this zone is `r_max` (in mrad). This is implemented in `maskChi0()`.

### Step 8 — Electron Probe

The real-space probe intensity is:

```
probe(r) = |FFT( χ(α) × aperture(α) )|²
```

The probe is shifted (fftshift) and normalised to 0–255. This is implemented in `probeGeneration()`.

### Step 9 — Strehl Ratio

The Strehl ratio measures how close the aberrated probe is to the diffraction limit:

```
S = ( |FFT(χ × A_s)|_max / |FFT(A_s)|_max )²
```

where `A_s` is a circular aperture of radius `r_strehl`. The code also searches for the aperture semi-angle that maximises the Strehl ratio to 0.8 (the Maréchal criterion). This is implemented in `singleStrehl()` and `pointEightStrehlSearch()`.

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
