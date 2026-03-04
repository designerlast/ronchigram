#include <math.h>
#include <stdint.h>
#include <iostream>
#include <stdlib.h>
#include <complex>
#include <time.h>
#include "kiss_fftnd.h"

using namespace std;
extern "C" {
int sub2ind(int x, int y, int z, int dimX, int dimY, int dimZ) {
    return dimX * dimY * z + dimX * y + x;
}

float* getMinMax(float* base, int dimX, int dimY) {
    float* minMax = new float[2];
    minMax[0] = base[0]; //min
    minMax[1] = base[0]; //max

    for (int i = 0; i < dimX * dimY; i++) {
        if (base[i] < minMax[0]) {
            minMax[0] = base[i];
        } else if (base[i] > minMax[1]) {
            minMax[1] = base[i];
        }
    }
    return minMax;
}

float* normalize(float* base, float scale, int dimX, int dimY) {
    float* minMax = getMinMax(base, dimX, dimY);
    for (int i = 0; i < dimX * dimY; i++) {
        base[i] = (base[i] - minMax[0]) / minMax[1] * scale;
    }
    return base;
}
// Relativistically corrected de Broglie wavelength for electrons.
// Formula: λ = h / sqrt(2 m0 eV (1 + eV / (2 m0 c²)))
// Simplified to: λ [m] = 12.3986 / sqrt((2*511 + keV)*keV) * 1e-10
float calculateLambda(float keV) {
    float lambda = 12.3986 / sqrt((2 * 511 + keV) * keV) * 1e-10;
    return lambda;
}

// Build a 2-D polar mesh in angular (reciprocal) space.
// Each pixel (i,j) maps to semi-angle r = sqrt(x²+y²) and azimuth p = atan2(y,x),
// where x and y run from -r_max to +r_max (in mrad).
// The objective aperture mask oapp is 1 inside obj_ap_r and 0 outside.
int polarMeshnOapp(float* rr, float* pp, float* oapp, float r_max, float obj_ap_r, int numPx) {
    float center = numPx / 2;
    int idx;
    float xval;
    float yval;
    for (int j = 0; j < numPx; j++) {
        for (int i = 0; i < numPx; i++) {
            idx = i + numPx * j;
            xval = (i - center) / center * r_max;
            yval = (j - center) / center * r_max;
            rr[idx] = sqrt(pow(xval, 2) + pow(yval, 2));
            pp[idx] = atan2(yval, xval);
            if (rr[idx] > obj_ap_r) {
                oapp[idx] = 0;
            } else {
                oapp[idx] = 1;
            }
        }
    }
    return 0;
}

// Generate a random noisy grating as a simple amorphous-sample model.
// Values are uniform random [0,1).  A coarser sub-sampled grid is tiled
// at 'scaleFactor' to avoid aliasing at high resolution.
float* noisyGrating(int dimX, int dimY) {
    srand(time(NULL));
    float* vals = new float[dimX * dimY];
    for (int i = 0; i < dimX * dimY; i++) {
        vals[i] = rand() / float(RAND_MAX);
    }
    return vals;
}

float* generateSample(int dimX, int dimY, int scaleFactor) {
    float* subsample = noisyGrating(dimX / scaleFactor, dimX / scaleFactor);
    float* supersample = new float[dimX * dimY];
    for (int j = 0; j < dimY; j++) {
        for (int i = 0; i < dimX; i++) {
            int idx = sub2ind(i, j, 0, dimX, dimY, 1);
            int idxsub = sub2ind(i / scaleFactor, j / scaleFactor, 0, dimX / scaleFactor, dimY / scaleFactor, 1);
            supersample[idx] = subsample[idxsub];
        }
    }
    return supersample;
}

// Build the sample transmission function t(r) = exp(-i * π/4 * σ * V(r))
// where V(r) is the projected potential (random sample values) and σ is
// the interaction parameter (relativistically corrected, normalised to 300 kV).
complex<float>* generateTransmissionFn(float* sample, int dimX, int dimY, float interactionParam) {
    complex<float>* trans = new complex<float> [dimX * dimY];
    complex<float> imag(0.0, 1.0);
    for (int i = 0; i < dimX * dimY; i++) {
        complex<float> real_part(M_PI / 4 * sample[i] * interactionParam, 0.0);
        trans[i] = exp(-imag * real_part);
    }
    return trans;
}

float* complexToReal(complex<float>* orig, int dimX, int dimY) {
    float* real = new float[dimX * dimY];
    for (int i = 0; i < dimX * dimY; i++) {
        real[i] = abs(orig[i]);
    }
    return real;
}

complex<float>* realToComplex(float* orig, int dimX, int dimY) {
    complex<float>* comp = new complex<float> [dimX * dimY];
    for (int i = 0; i < dimX * dimY; i++) {
        comp[i].real(orig[i]);
        comp[i].imag(0);
    }

    return comp;
}

kiss_fft_cpx* complexToKiss(complex<float>* orig, int dimX, int dimY) {
    kiss_fft_cpx* kisscpx = new kiss_fft_cpx[dimX * dimY];
    kiss_fft_cpx holder;
    for (int i = 0; i < dimX * dimY; i++) {
        holder.r = real(orig[i]);
        holder.i = imag(orig[i]);
        kisscpx[i] = holder;
    }
    return kisscpx;
}

complex<float>* kissToComplex(kiss_fft_cpx* orig, int dimX, int dimY) {
    complex<float>* comp = new complex<float> [dimX * dimY];
    for (int i = 0; i < dimX * dimY; i++) {
        //comp[i] = (orig[i].r,orig[i].i);
        comp[i].real(orig[i].r);
        comp[i].imag(orig[i].i);
    }
    return comp;
}

float* packageOutput(float* base1, float* base2, float* base3, float* scalars, int dimX, int dimY, int nScal) {
    int sz = dimX * dimY;
    float* imageStack = new float[sz * 2 + nScal];
    for (size_t i = 0; i < sz; i++) {
        imageStack[i] = base1[i]; // copy the allocated memory
        imageStack[sz + i] = base2[i];
        imageStack[2 * sz + i] = base3[i];
    }
    for (size_t i = 0; i < nScal; i++) {
        imageStack[3 * sz + i] = scalars[i];
    }
    return imageStack;
}

// Compute the geometrical aberration phase χ₀(α, Φ) in the Krivanek notation:
//   χ₀ = (2π/λ) × Σ_k  C_{n,m} × α^(n+1) × cos(m×(Φ−Φ_{n,m})) / (n+1)
// magptr[k] = aberration magnitude C_{n,m} (metres)
// angleptr[k] = aberration orientation Φ_{n,m} (radians)
// 14 terms are supported, covering aberrations up to 5th order.
float* calculateChi0(float* magptr, float* angleptr, float* alrr, float* alpp, int numPx, int numAbs, float keV) {
    float* chi0 = new float[numPx * numPx];
    int n[14] = {1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5};
    int m[14] = {0, 2, 1, 3, 0, 2, 4, 1, 3, 5, 0, 2, 4, 6};
    float lambda = calculateLambda(keV);

    for (int i = 0; i < numPx*numPx; i++) {
        chi0[i] = 0;
        for (int k = 0; k < numAbs; k++) {
            chi0[i] = chi0[i] + 2 * M_PI / lambda * magptr[k] * pow(alrr[i], n[k] + 1) * cos(m[k] * (alpp[i] - angleptr[k])) / (n[k] + 1);
        }
    }
    return chi0;
}

// Convert the real-valued aberration phase χ₀ into the complex pupil function
//   χ(α, Φ) = exp(−i × χ₀(α, Φ))
// This represents the phase-only transfer function of the aberrated lens.
complex<float>* calculateChi(float* chi0, int numPx) {
    complex<float>* chi = new complex<float> [numPx * numPx];
    complex<float> imag(0.0, 1.0);
    for (int i = 0; i < numPx * numPx; i++) {
        chi[i] = exp(-imag * chi0[i]);
    }
    return chi;
}

// Build the π/4 phase-limit mask.
// Pixels where |χ₀| > threshold (default π/4) are set to 0 (outside the
// isochronous zone); others are set to 255.
// Also returns the radius of the largest circle fitting inside the zone (r_max, mrad × 1000).
float maskChi0(float* chi0, float* alrr, int numPx, float threshold) {
    float rmax = 1e5;
    for (int i = 0; i < numPx * numPx; i++) {
        chi0[i] = 255 * ((abs(chi0[i]) > threshold) ? 0 : 1);
        float cv = (1 - chi0[i]) * alrr[i];
        if (cv > 0 && cv < rmax) {
            rmax = cv;
        }
    }
    return rmax * 1000;
}

float* fftShift(float* original, int numPx) {
    // only for square matrices...
    float* shifted = new float[numPx * numPx];

    for (int j = 0; j < numPx / 2; j++) {
        for (int i = 0; i < numPx / 2; i++) {
            int idx_q1 = sub2ind(i, j, 0, numPx, numPx, 1);
            int idx_q3 = sub2ind(i + numPx / 2, j + numPx / 2, 0, numPx, numPx, 1);
            int idx_q2 = sub2ind(i + numPx / 2, j, 0, numPx, numPx, 1);
            int idx_q4 = sub2ind(i, j + numPx / 2, 0, numPx, numPx, 1);
            shifted[idx_q1] = original[idx_q3];
            shifted[idx_q3] = original[idx_q1];
            shifted[idx_q4] = original[idx_q2];
            shifted[idx_q2] = original[idx_q4];
        }
    }
    return shifted;
}

complex<float>* cmplxFFT(complex<float>* comp, int dimX, int dimY) {
    int isInverseFFT = 0;
    int ndims = 2;
    int dims[2];
    dims[0] = dimX;
    dims[1] = dimY;

    kiss_fft_cpx* cxin;
    kiss_fft_cpx* cxout;

    kiss_fftnd_cfg cfg = kiss_fftnd_alloc(dims, ndims, isInverseFFT, 0, 0);
    // preallocation
    cxin = complexToKiss(comp, dimX, dimY);
    cxout = complexToKiss(comp, dimX, dimY);
    kiss_fftnd(cfg, cxin, cxout);
    complex<float>* fftResult = kissToComplex(cxout, dimX, dimY);

    return fftResult;
}

// Compute the Ronchigram diffraction intensity:
//   I(q) = |FFT( t(r) × FFT( χ(α) ) )|² × oapp(q)
// Steps:
//  1. FFT(χ) → real-space probe (before aperture)
//  2. Multiply by sample transmission function t(r)
//  3. FFT back to reciprocal space
//  4. |·|² and apply objective aperture mask
float* calcDiffract(complex<float>* chi, complex<float>* trans, float* oapp, int numPx) {
    // want: abs(fft(trans*fft(chi)))

    //in place b/c chi won't be reused
    chi = cmplxFFT(chi, numPx, numPx);
    for (int i = 0; i < numPx * numPx; i++) {
        chi[i] = chi[i] * trans[i];
    }
    chi = cmplxFFT(chi, numPx, numPx);
    float* diffInt = new float[numPx * numPx];
    diffInt = complexToReal(chi, numPx, numPx);
    for (int i = 0; i < numPx * numPx; i++) {
        diffInt[i] = diffInt[i] * oapp[i] * diffInt[i];
    }
    return diffInt;
}

// Relativistically corrected interaction parameter σ, normalised to the
// value at 300 kV, used as the multiplicative factor in the transmission
// function phase: σ = (2π / λE) × (m₀c² + eV) / (2m₀c² + eV).
float calculateInteractionParam(float keV) {
    float keV_300 = 300;
    float c = 3e8;
    float mass_e = 9.11e-31;
    float charge_e = 1.602e-19;
    float lambda = calculateLambda(keV);
    float lambda_300 = calculateLambda(keV_300);
    float param = 2 * M_PI / (lambda * keV / charge_e * 1000) * (mass_e * c * c + keV * 1000) / (2 * mass_e * c * c + keV * 1000);
    float param_300 = 2 * M_PI / (lambda_300 * keV_300 / charge_e * 1000) * (mass_e * c * c + keV_300 * 1000) / (2 * mass_e * c * c + keV_300 * 1000);
    return param / param_300;
}

// Compute the Strehl ratio for a given aperture radius r_strehl (mrad × 1000).
// S = ( max|FFT(χ × A_s)| / max|FFT(A_s)| )²
// where A_s is a circular aperture of radius r_strehl.
// A Strehl ratio of 1 indicates a perfect (diffraction-limited) lens.
float singleStrehl(float rmax, complex<float>* chi, float al_max, int numPx){
    float strr[numPx * numPx];
    float stpp[numPx * numPx];
    float sapp[numPx * numPx];
    float  strehl_radius = rmax/1000;
    polarMeshnOapp(strr, stpp, sapp, al_max, strehl_radius, numPx);
    float* strehl_aperture_real = sapp;
    complex<float>* strehl_aperture = realToComplex(strehl_aperture_real, numPx, numPx);
    complex<float>* strehl_inner = new complex<float> [numPx * numPx];
     for(int i = 0; i < numPx * numPx; i++){
         strehl_inner[i] = chi[i] * strehl_aperture[i];
     };
    complex<float>* strehl_inner_fft = cmplxFFT(strehl_inner, numPx, numPx);
    float* strehl_inner_fft_real = complexToReal(strehl_inner_fft, numPx, numPx);
    float strehl_inner_max = getMinMax(strehl_inner_fft_real, numPx, numPx)[1];
    complex<float>*strehl_bottom = cmplxFFT(strehl_aperture, numPx, numPx);
    float* strehl_bottom_real = complexToReal(strehl_bottom, numPx, numPx);
    float strehl_bottom_max = getMinMax(strehl_bottom_real, numPx, numPx)[1];
    float strehl = strehl_inner_max/strehl_bottom_max;
    strehl = pow(strehl,2);
    return strehl;
}

// Generate the real-space electron probe intensity:
//   probe(r) = |FFT( χ(α) × aperture(α) )|²
// The result is fftshift-ed and normalised to 0–255.
float* probeGeneration(complex<float>* chi, int numPx, float* oapp){
    complex<float>* obj_aperture = realToComplex(oapp, numPx, numPx);
    complex<float>* probe = new complex<float> [numPx * numPx];
     for(int i = 0; i < numPx * numPx; i++){
         probe[i] = chi[i] * obj_aperture[i];
     };
    complex<float>* probe_fft = cmplxFFT(probe, numPx, numPx);
    float* probe_abs = complexToReal(probe_fft, numPx, numPx);
    float* probe_shift = fftShift(probe_abs, numPx);
    float* probe_out = normalize(probe_shift, 255, numPx, numPx);
    return probe_out;
}

// Binary search for the aperture semi-angle that achieves a Strehl ratio of
// exactly 0.8 (the Maréchal criterion for diffraction-limited imaging).
// Starts from the π/4 radius (pi_radius) and searches within [1×, 1.8×] that range.
float* pointEightStrehlSearch(float pi_radius, complex<float>* chi, float al_max, int numPx){

    float upper_bound = 1.8 * pi_radius/1000;
    float lower_bound = pi_radius/1000;
    float upper_strehl = 1;
    float lower_strehl = 0.5;
    float tolerance = 0.01;
    float middle;
    float ans;
    float count = 0;
    while((upper_strehl-lower_strehl)>=tolerance){
        middle = (upper_bound+lower_bound)/2;
        ans = singleStrehl(1000*middle, chi, al_max, numPx);
        if(ans > 0.8){
            lower_bound = middle;
            upper_strehl = ans;
        }
        if(ans < 0.8){
            upper_bound = middle;
            lower_strehl = ans;
        }
        count = count + 1;
        if(count == 10){
            break;
        }
    }
    middle = (upper_bound+lower_bound)/2;
    float r_strehl = singleStrehl(1000*middle, chi, al_max, numPx);
    float* returns = new float[3];
    returns[0] = middle;
    returns[1] = r_strehl;
    returns[2] = count;
    return returns;
}


// Main entry point called from JavaScript / WebAssembly.
// buffer layout:
//   [0]      numPx          – image size (pixels, square)
//   [1]      al_max         – display semi-angle (mrad)
//   [2]      obj_ap_r       – objective aperture semi-angle (mrad)
//   [3]      scalefactor    – sample coarsening factor
//   [4]      keV            – beam energy (keV)
//   [5]      calcStrehl     – 1 = compute Strehl search, 0 = skip
//   [6..19]  aberration magnitudes C_{n,m} (metres, 14 terms)
//   [20..33] aberration orientations Φ_{n,m} (radians, 14 terms)
//
// Returns a flat float array:
//   [0 .. numPx²-1]             Ronchigram image (0–255)
//   [numPx² .. 2·numPx²-1]      π/4 phase-limit map (0 or 255)
//   [2·numPx² .. 3·numPx²-1]    Electron probe image (0–255)
//   [3·numPx²]                  r_max (mrad × 1000) – π/4 aperture radius
//   [3·numPx²+1]                Strehl ratio at r_max
//   [3·numPx²+2]                Strehl ratio at optimal aperture (−1 if skipped)
//   [3·numPx²+3]                Optimal aperture semi-angle (mrad × 1000)
//   [3·numPx²+4]                Number of bisection iterations used
float* calcRonch(float* buffer, int bufSize) {
    int numPx = static_cast < int > (buffer[0]);
    float al_max = buffer[1]; //mrad
    float obj_ap_r = buffer[2]; //mrad
    int scalefactor = buffer[3];
    float keV = buffer[4];
    float calcStrehl = buffer[5];
    float* outputScalars = new float[5];
    float alrr[numPx * numPx];
    float alpp[numPx * numPx];
    float oapp[numPx * numPx];
    polarMeshnOapp(alrr, alpp, oapp, al_max, obj_ap_r, numPx);

    float* sample = generateSample(numPx, numPx, scalefactor);

    complex<float>* trans = generateTransmissionFn(sample, numPx, numPx, calculateInteractionParam(keV));

    float* chi0 = calculateChi0( & buffer[6], & buffer[20], alrr, alpp, numPx, 14, keV);

    complex<float>* chi = calculateChi(chi0, numPx);
    // Calculate r_max and return,  and turn chi0 into pi/4map Normalized
    outputScalars[0] = maskChi0(chi0, alrr, numPx, M_PI / 4);
    
    float strehl = singleStrehl(outputScalars[0], chi, al_max, numPx);
    outputScalars[1] = strehl;

    if(calcStrehl == true){
    float* strehl_search_results = pointEightStrehlSearch(outputScalars[0], chi, al_max, numPx);

    outputScalars[2] = strehl_search_results[1]; //strehl ratio achieved
    outputScalars[3] = strehl_search_results[0]*1000; //strehl radius
    outputScalars[4] = strehl_search_results[2]; //number of iterations needed to find strehl
    }
    else{
        outputScalars[2] = -1; //strehl ratio achieved
        outputScalars[3] = 0; //strehl radius
        outputScalars[4] = 0; //count
    }

    float strr[numPx * numPx];
    float stpp[numPx * numPx];
    float sapp[numPx * numPx];
    float  strehl_radius = outputScalars[3]/1000;


   if(calcStrehl == true){
    polarMeshnOapp(strr, stpp, sapp, al_max, strehl_radius, numPx);
    }
    else{
    polarMeshnOapp(strr, stpp, sapp, al_max, outputScalars[0]/1000, numPx);
    
    }
    
    float* probe = probeGeneration(chi, numPx, sapp);

    // Normalize to 0-255 for output
    float* ronch = normalize(calcDiffract(chi, trans, oapp, numPx), 255, numPx, numPx);
    // Package results and return
    auto arrayPtr = packageOutput(ronch, chi0, probe, outputScalars, numPx, numPx, 5);
    return arrayPtr;
}

}