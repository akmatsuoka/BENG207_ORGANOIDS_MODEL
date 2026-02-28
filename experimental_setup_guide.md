# Experimental Setup: Acoustic Tweezer Organoid Deformation Imaging

## Equipment Inventory

| Component         | Model                                       | Do we have it |
| ----------------- | ------------------------------------------- | ------------- |
| Microscope        | Olympus BX51WI (upright, water immersion)   | BSB 4041      |
| Fluorescence lamp | Olympus U-RFL-T (Hg arc, ~8200 hrs)         | BSB 4041      |
| Transmitted light | Olympus TH4-100 (halogen)                   | BSB 4041      |
| Shutter           | Uniblitz electronic + VCM-D1 driver         | BSB 4041      |
| Stage platform    | Siskiyou MXMS-125                           | BSB 4041      |
| Camera            | acA1920-155um (IMX174, global shutter,USB3) | NO            |
| C-mount adapter   | 0.63× (Basler, 1/1.2"U-TV0.63XC or generic  | NO            |

| Optical table     | Breadboard on vibration-isolation legs              | TONS BSB 4041      

## Additional Equipment Needed

| Component                       | Purpose                                 | Suggested model / spec                              |
| ------------------------------- | --------------------------------------- | --------------------------------------------------- |
| Function generator              | RF signal source, 20 MHz                | We have this                                        |
| RF amplifier                    | Drive transducer, 50 dB gain, DC–50 MHz | We have this                                        |
| Delay / pulse generator         | Sync US pulses with camera frames       | We have this                                        |
| Transducer                      | 20 MHz unfocused or focused piston      | We have this                                        |
| Couplant                        | Acoustic coupling to glass              | US gel (Aquasonic 100) will steal it from Perlman   |
| Microfluidic chip / chamber     | Hold organoids + create acoustic cavity | Custom glass–glass or glass–Si bonded (Alexi knows) |
| BNC cables + impedance matching | Connect RF chain                        | 50 Ω throughout, everywhere in my lab               |

---

## System Architecture

```
                    TRIGGER BUS (TTL)
                         │
        ┌────────────────┼──────────────────┐
        │                │                  │
        ▼                ▼                  ▼
 ┌──────────┐    ┌──────────────┐    ┌───────────┐
 │  Delay/  │    │   Function   │    │  Basler   │
 │  Pulse   │───▶│  Generator   │    │  Camera   │
 │  Gen     │    │  (20 MHz CW) │    │ acA1920   │
 └──────────┘    └──────┬───────┘    └─────┬─────┘
                        │                  │
                        ▼                  │ USB3
                 ┌──────────────┐          │
                 │ RF Amplifier │          │
                 │ (50 dB gain) │          │
                 └──────┬───────┘          │
                        │ coax             │
                        ▼                  │
                 ┌──────────────┐          │
                 │  Transducer  │          │
                 │  (20 MHz)    │          │
                 └──────┬───────┘          │
                   US gel│couplant         │
                        ▼                  │
              ┌─────────────────┐          │
              │  Glass slide    │          │
              │  ─ ─ ─ ─ ─ ─ ─  │          │
              │  Fluid channel  │◀── Organoid (200 µm)
              │  ─ ─ ─ ─ ─ ─ ─  │          │
              │  Coverslip      │          │
              └─────────────────┘          │
                        │                  │
                        ▼                  ▼
              ┌─────────────────────────────────┐
              │  Olympus BX51WI                 │
              │  Objective (10× or 20×)         │
              │  → Basler camera (USB3)         │
              │  Siskiyou MXMS-125 stage        │
              └─────────────────────────────────┘
```

---

## Timing Diagram

The key challenge: synchronize acoustic pulses with camera frames so every deformation event is captured.

```
Frame period = 1/fps
For 155 fps (full frame): T_frame = 6.45 ms
For 300 fps (ROI crop):   T_frame = 3.33 ms

TRIGGER SEQUENCE (one deformation cycle):

         T_on = 500 ms              T_off = 2000 ms
    ◀─────────────────────▶   ◀──────────────────────────▶

    ┌─────────────────────┐
US  │█████████████████████│
    └─────────────────────┘───────────────────────────────

    ┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐
CAM ││││││││││││││││││││││││││││││││││││││││││││││││││││
    └┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘
    |← pre-trigger →|← US ON →|← recovery imaging  →|
       ~50 frames     ~77 frames     ~310 frames

                 6.45 ms per frame at 155 fps
                 ──▶│◀──

SHUTTER (Uniblitz): OPEN throughout acquisition
     (or sync to camera if photobleaching is a concern, we need to do
 thsi first beads of course)
```

### Timing configuration

| Parameter          | Value                                      | Notes                                              |
| ------------------ | ------------------------------------------ | -------------------------------------------------- |
| Camera fps         | 155 (full) or 300+ (ROI)                   | Global shutter eliminates motion blur              |
| Camera exposure    | 2–4 ms                                     | Shorter than frame period; brightfield OK at ~2 ms |
| US pulse ON        | 500 ms (primary), also test 100 ms and 2 s | Captures both τ₁ and τ₂                            |
| US pulse OFF       | 2–5 s (short), 30 s (long recovery)        | Must exceed τ₂ = 10 s for full recovery test       |
| Pre-trigger frames | 50 frames (~320 ms)                        | Baseline before deformation                        |
| Number of pulses   | 5–10 per acquisition                       | Track residual strain accumulation                 |
| Total acquisition  | 25–55 s per pulse train                    | Depends on protocol                                |

---

## Basler acA1920-155um Camera Settings

### Pylon configuration

```
Pixel Format:        Mono8 (fastest) or Mono12 (if dynamic range needed)
Sensor readout:      Global shutter (default for IMX174)
Trigger mode:        On (external trigger from delay generator)
  or                 Free-run at fixed fps (simpler, timestamp each frame)
Trigger source:      Line1 (hardware trigger via GPIO)
Exposure time:       2000–4000 µs
Gain:                0 dB (brightfield) or auto (fluorescence)
```

### ROI for higher frame rate

The IMX174 can exceed 155 fps with a reduced region of interest. For a 200 µm organoid:

| Objective           | Organoid on sensor | Recommended ROI | Approx fps |
| ------------------- | ------------------ | --------------- | ---------- |
| 10× (0.65 µm/px)    | ~308 px diameter   | 512 × 512       | ~350 fps   |
| 20× (0.325 µm/px)   | ~615 px diameter   | 800 × 800       | ~250 fps   |
| 10× with 1.6× relay | ~492 px diameter   | 640 × 640       | ~300 fps   |

For measuring ~12 µm diameter changes at 20×: that's ~37 pixels of deformation — well within sub-pixel tracking accuracy.

### Data rate and storage

| Config                     | Data rate | 60 s acquisition |
| -------------------------- | --------- | ---------------- |
| 1920×1200 @ 155 fps, Mono8 | ~358 MB/s | ~21 GB           |
| 800×800 @ 250 fps, Mono8   | ~160 MB/s | ~9.6 GB          |
| 512×512 @ 350 fps, Mono8   | ~92 MB/s  | ~5.5 GB          |

Use SSD (NVMe preferred). USB3 bandwidth is ~350 MB/s max.

---

## Acoustic Channel Design (from Echo Chamber Model)

Based on the model results at 20 MHz:

| Parameter          | Recommended                    | Why                                          |
| ------------------ | ------------------------------ | -------------------------------------------- |
| Channel length L   | 2.5 mm                         | Δf = 300 kHz spacing (robust to drift)       |
| Channel depth      | 100–200 µm                     | Fits 200 µm organoid; supports standing wave |
| End condition      | Glass–glass bonded (rigid)     |                                              |
| Channel width      | 500 µm–1 mm                    | Allows organoid entry + rotation             |
| Couplant thickness | < 7.5 µm (λ/10 at 20 MHz)      | Minimizes coupling sensitivity               |
| Couplant type      | US gel (Aquasonic 100)         | Or UV epoxy for stability                    |
| Substrate          | Borosilicate glass, 1 mm       | Matches model parameters                     |
| Cover              | Borosilicate coverslip, 150 µm | Matches model; fits BX51WI WD                |

### Standing wave nodes in channel

At 20 MHz in water: λ/2 = 37.1 µm. In a 200 µm deep channel, you get ~5 pressure nodes across the depth. The organoid (200 µm) spans the full channel depth, so it experiences the integrated radiation pressure — consistent with the "effective stress" approach in the organoid model.

---

## Experimental Protocols

### Protocol A: Fast timescale identification (τ₁)--

```
Purpose:   Capture fast viscoelastic response (cortex deformation)
US pulse:  100 ms ON, 2 s OFF, × 5 pulses
Camera:    300+ fps (use ROI), ~15 s total acquisition
Objective: 20× (water immersion) for maximum resolution
Analysis:  Fit exponential rise/fall to get τ₁ and A₁ = σ₀/E₁
Expected:  τ₁ ~ 0.1 s → ~30 frames to capture rise
```

### Protocol B: Slow timescale identification (τ₂)

```
Purpose:   Capture slow relaxation (ECM, cell rearrangement)
US pulse:  2 s ON, 30 s OFF, × 3 pulses
Camera:    50 fps (full frame OK), ~100 s total acquisition
Objective: 10× (wider field, track any drift)
Analysis:  Fit long recovery tail to get τ₂ and A₂ = σ₀/E₂
Expected:  τ₂ ~ 10 s → need >30 s recovery window
```

### Protocol C: Residual strain accumulation

```
Purpose:   Demonstrate two-timescale behavior
US pulse:  500 ms ON, 2 s OFF, × 10 pulses
Camera:    155 fps (full frame), ~25 s total
Objective: 20× for precision
Analysis:  Track rising baseline between pulses
Expected:  Single KV → no baseline rise
            Modified KV → baseline increases by ~1–2% over 10 pulses
```

### Protocol D: Stress calibration (drag balance)

```
Purpose:   Calibrate σ₀ independently
Method:    Apply flow at known velocity U while US is on
           At equilibrium: F_ac = 6πμRU
           σ₀ = F_ac / (πR²)
Camera:    155 fps, track organoid displacement
Flow:      Syringe pump, 0.1–10 µL/min through channel
Analysis:  Find U where organoid just escapes trap
```

---

## Signal Chain: Wiring Diagram

```
DELAY/PULSE GENERATOR (BNC 575 or SRS DG535)
├── Ch A: Master trigger (start of acquisition)
│         → Camera trigger input (GPIO Line1)
│         → Also starts recording software
│
├── Ch B: US gate pulse (delayed from A)
│         Duration: T_on (100 ms / 500 ms / 2 s)
│         Repeat: N_pulses at T_on + T_off period
│         → Function generator GATE/TRIG input
│
├── Ch C: (optional) Uniblitz shutter sync
│         → VCM-D1 trigger input (BNC on back)
│
└── Ch D: (optional) End-of-sequence marker
          → Camera GPIO Line2 (metadata stamp)

FUNCTION GENERATOR (20 MHz CW)
├── Output: gated by Ch B from delay gen
│         Amplitude: ~100 mVpp into 50 Ω
│         → RF amplifier input
│
└── Sync out → oscilloscope (monitoring)

RF AMPLIFIER
├── Input: from function generator
├── Output: 50 Ω coax → transducer
│          Power: 1–5 W typical
└── Monitor: directional coupler → oscilloscope

TRANSDUCER → [US gel] → glass slide → channel → coverslip
                                                    │
                                                    ▼
                                        Olympus BX51WI objective
                                                    │
                                              Basler camera
                                                    │
                                              USB3 → PC
                                         (Pylon Viewer or
                                          custom Python/MATLAB or whatever)
```

---

## Image Analysis Pipeline

### Step 1: Frame-by-frame contour extraction

```
Tool:    Fiji/ImageJ + TrackMate-Cellpose (same as in the paper)
         or Python + OpenCV contour detection
Input:   Image stack (TIFF or raw from Pylon)
Output:  Major axis a(t), minor axis b(t) per frame
```

### Step 2: Strain computation

```python
# From contour axes
a0 = np.mean(a[pre_trigger_frames])   # baseline major axis
b0 = np.mean(b[pre_trigger_frames])   # baseline minor axis

epsilon_a = (a - a0) / a0              # axial strain (compression)
epsilon_b = (b - b0) / b0              # transverse strain (extension)

# Volume conservation check
V_proxy = a * b**2                     # prolate spheroid proxy
V_ratio = V_proxy / V_proxy[0]        # should stay ~1.0
```

### Step 3: Fit KV model parameters

```python
from scipy.optimize import curve_fit

def modified_kv_pulse(t, A1, tau1, A2, tau2, T_on):
    """Two-branch KV response to rectangular pulse."""
    eps = np.zeros_like(t)
    on = t < T_on
    off = ~on
    eps[on]  = A1*(1-np.exp(-t[on]/tau1)) + A2*(1-np.exp(-t[on]/tau2))
    e1T = A1*(1-np.exp(-T_on/tau1))
    e2T = A2*(1-np.exp(-T_on/tau2))
    eps[off] = e1T*np.exp(-(t[off]-T_on)/tau1) + e2T*np.exp(-(t[off]-T_on)/tau2)
    return eps

# Fit
popt, pcov = curve_fit(lambda t, A1, t1, A2, t2: modified_kv_pulse(t, A1, t1, A2, t2, T_on=0.5),
                       t_data, epsilon_data,
                       p0=[0.03, 0.1, 0.03, 10],
                       bounds=([0,0.01,0,0.5], [0.5,2,0.5,100]))
A1, tau1, A2, tau2 = popt
```

### Step 4: Convert to material properties (after σ₀ calibration)

```python
E1 = sigma0 / A1      # Pa
E2 = sigma0 / A2      # Pa
eta1 = E1 * tau1       # Pa·s
eta2 = E2 * tau2       # Pa·s
```

---

## Checklist Before Experiment

- [ ] Transducer impedance matched (check with network analyzer if available)
- [ ] US gel applied thin (<10 µm) — use minimal amount, press firmly
- [ ] Channel filled, no air bubbles (critical for standing wave)-ask Spencer (He is the king of getting rid of bubbles)
- [ ] Organoid loaded, settled to channel bottom or suspended at node
- [ ] Camera ROI set, exposure verified, frames saving to SSD
- [ ] Delay generator programmed: pre-trigger → US gate → N pulses
- [ ] Oscilloscope monitoring: RF amplitude at transducer
- [ ] Baseline images acquired (50+ frames, US off)
- [ ] Test pulse at low power to verify acoustophoresis before full experiment
- [ ] Temperature noted (affects c_f and viscosity)

---

## Expected Results (from Model Predictions)

| Measurement                           | Expected value      | Model source             |
| ------------------------------------- | ------------------- | ------------------------ |
| Peak strain (500 ms pulse, 500 kPa)   | ~6%                 | Organoid KV model        |
| Diameter change                       | ~12 µm              | Organoid KV model        |
| Fast time constant τ₁                 | 0.05–0.3 s          | Cortex mechanics         |
| Slow time constant τ₂                 | 5–30 s              | ECM/rearrangement        |
| Residual strain after 10 pulses       | 1–3%                | Modified KV accumulation |
| Minimum detectable ΔD (20× obj)       | ~0.3 µm (sub-pixel) | Camera resolution        |
| Minimum acoustic pressure for 1 µm ΔD | ~200 kPa            | Stress sweep (fig 7)     |
