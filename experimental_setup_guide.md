# Experimental Setup: Acoustic Tweezer Organoid Deformation Imaging

This guide covers two configurations:

- **Part A — Single Transducer (original):** 20 MHz single-source standing wave. Organoid deformation under uniaxial radiation stress. Useful for initial validation and parameter fitting.
- **Part B — Dual PZT Config B:** Two PZT-5A transducers immersed in water on opposing sides of the channel, counter-propagating at 741 kHz. Perfect standing wave (infinite contrast), biaxial radiation stress, phase-steerable nodes. This is the target configuration for the biohybrid cochlear implant prototype.

---

## Equipment Inventory

### Part A — Single Transducer

| Component         | Model                                               | Location / Status                     |
| ----------------- | --------------------------------------------------- | ------------------------------------- |
| Microscope        | Olympus BX51WI (upright, water immersion)           | BSB 4041                              |
| Fluorescence lamp | Olympus U-RFL-T (Hg arc, ~3.3 hrs)                  | BSB 4041                              |
| Transmitted light | Olympus TH4-100 (halogen)                           | BSB 4041                              |
| Shutter           | Uniblitz electronic + VCM-D1 driver                 | BSB 4041                              |
| Stage platform    | Siskiyou MXMS-125                                   | BSB 4041                              |
| Camera            | Basler acA1920-155um (IMX174, global shutter, USB3) | NO — order from Edmund Optics #89-985 |
| C-mount adapter   | U-TV1X-2 (1×) already installed under U-CMAD3       | BSB 4041 (no new adapter needed)      |
| Optical table     | Breadboard on vibration-isolation legs              | BSB 4041                              |

### Part B — Dual PZT Config B (additional / different items)

| Component                | Model / Spec                                                     | Location / Status                   |
| ------------------------ | ---------------------------------------------------------------- | ----------------------------------- |
| Transducers (×2)         | PZT-5A, thickness-mode resonant near 741 kHz                     | Need to order                       |
| Function generator ×2    | Two-channel or two units — 741 kHz CW, independent phase control | We have this (verify channel count) |
| RF amplifier ×2          | One per transducer, 50 Ω, sufficient power for 500 kPa           | We have this                        |
| Phase controller         | 0–360° independent phase offset per channel                      | Verify on function generator        |
| Water bath / tank        | PZTs immersed in water — no couplant gel needed                  | Custom fabrication                  |
| PDMS side walls          | ≤3 mm thick (thinned from current 10 mm)                         | Fabricate in lab                    |
| PDMS ceiling             | 1–2 mm thick (thinned from current 9 mm for objective WD)        | Fabricate in lab                    |
| Borosilicate glass slide | 1 mm, underneath channel — provides TIR confinement              | BSB 4041                            |
| All other imaging        | Same as Part A (BX51WI, Basler, shutter, stage)                  | BSB 4041                            |

---

## Additional Equipment Needed

### Part A

| Component                       | Purpose                                 | Suggested model / spec                              |
| ------------------------------- | --------------------------------------- | --------------------------------------------------- |
| Function generator              | RF signal source, 20 MHz                | We have this                                        |
| RF amplifier                    | Drive transducer, 50 dB gain, DC–50 MHz | We have this                                        |
| Delay / pulse generator         | Sync US pulses with camera frames       | We have this                                        |
| Transducer                      | 20 MHz unfocused or focused piston      | We have this                                        |
| Couplant                        | Acoustic coupling to glass              | US gel (Aquasonic 100) — borrow from Perlman lab    |
| Microfluidic chip / chamber     | Hold organoids + create acoustic cavity | Custom glass–glass or glass–Si bonded (Alexi knows) |
| BNC cables + impedance matching | Connect RF chain                        | 50 Ω throughout                                     |

### Part B (additional)

| Component                                           | Purpose                                                     | Notes                                                                                                                                                                                                                                                                                                                                                                               |
| --------------------------------------------------- | ----------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Second RF amplifier                                 | Drive right PZT independently                               | Match gain to left channel for amplitude balance                                                                                                                                                                                                                                                                                                                                    |
| Phase splitter / controller                         | Set φ_L − φ_R from 0 to π                                   | Some dual-channel function generators support this directly                                                                                                                                                                                                                                                                                                                         |
| Water-immersion tank                                | PZTs sit in water on both sides of PDMS walls               | Custom acrylic or 3D-printed holder                                                                                                                                                                                                                                                                                                                                                 |
| PDMS mold (revised)                                 | Side walls ≤3 mm; ceiling 1–2 mm                            | Rework existing mold; thin walls improve T_power to 77–92%                                                                                                                                                                                                                                                                                                                          |
| Amplitude attenuator (×1)                           | Balance p_L = p_R to within ~5% for symmetric nodes         | 3–6 dB fixed pad on one channel if gains differ                                                                                                                                                                                                                                                                                                                                     |
| Hydrophone (required long-term, optional near-term) | Measure actual p_ac in channel; calibrate σ₀ for Protocol D | **Near-term:** use Basler camera + beads as substitution detector (free). **Long-term:** Precision Acoustics HP1 or HP2 (1–2 mm PVDF, NPL-calibrated 30 kHz–40 MHz, ~£3,000–4,500 complete system). Standard Onda HNA/HNP/HNR probes start calibration at 1 MHz — above 741 kHz; request extended low-frequency calibration if ordering from Onda (~$2,000–3,500 + AH-2010 preamp). |

---

## System Architecture

### Part A — Single Transducer

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

### Part B — Dual PZT Config B

```
                         TRIGGER BUS (TTL)
                                │
        ┌───────────────────────┼────────────────────────┐
        │                       │                        │
        ▼                       ▼                        ▼
 ┌──────────┐         ┌──────────────────┐        ┌───────────┐
 │  Delay/  │         │  Dual-Channel    │        │  Basler   │
 │  Pulse   │────────▶│  Function Gen    │        │  Camera   │
 │  Gen     │         │  741 kHz CW      │        │ acA1920   │
 └──────────┘         │  Ch1 ──── Ch2    │        └─────┬─────┘
                      │  (phase offset φ)│              │
                      └──────┬──────┬───┘              │ USB3
                             │      │                   │
                             ▼      ▼                   │
                      ┌──────┐    ┌──────┐              │
                      │ AMP  │    │ AMP  │              │
                      │ Left │    │Right │              │
                      └──┬───┘    └──┬───┘              │
                         │ coax      │ coax             │
                         ▼           ▼                  │
              ┌──────────────────────────────┐          │
              │  PZT-5A     Water    PZT-5A  │          │
              │  (left) ←─ 5mm ─→  (right)  │          │
              │         ┌────────┐           │          │
              │ PDMS    │ Water  │  PDMS     │          │
              │ wall    │channel │  wall     │          │
              │ ≤3mm    │ 20mm   │  ≤3mm     │          │
              │         │~20     │           │          │
              │  p_L →  │organoids← p_R     │          │
              │         │1mm apart           │          │
              │         └────────┘           │          │
              │    Borosilicate glass (1mm)  │          │
              │    TIR θ_c = 15.3° ✓         │          │
              │    PDMS ceiling: 1–2mm       │          │
              └──────────────────────────────┘          │
                              │                         │
                              ▼                         ▼
              ┌─────────────────────────────────────────┐
              │  Olympus BX51WI (upright, images from above)
              │  Objective: UMPLFLN10XW or UMPLFLN20XW  │
              │  WD = 3.5 mm — fits through 1–2 mm PDMS │
              │  → Basler camera (USB3)                  │
              │  Siskiyou MXMS-125 stage                 │
              └─────────────────────────────────────────┘
```

**Key differences from Part A:**

- No couplant gel — PZTs are immersed directly in water (eliminates BENG207_1 failure mode)
- Two independent RF channels with phase control
- PDMS side walls must be ≤3 mm (currently 10 mm — must thin)
- PDMS ceiling must be 1–2 mm (currently 9 mm — must thin for BX51WI WD = 3.5 mm)
- Glass slide underneath is essential — provides TIR acoustic confinement (θ_c = 15.3°)
- Perfect standing wave contrast (infinite) regardless of PDMS wall reflectivity

---

## Timing Diagram

### Part A — Single Transducer

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
     (or sync to camera if photobleaching is a concern —
      test with beads first before organoids)
```

### Part B — Dual PZT Config B

Same camera timing as Part A. The additional consideration is the two-channel RF trigger — both amplifiers must be gated simultaneously so p_L and p_R are on/off together.

```
TRIGGER SEQUENCE (Part B, one deformation cycle):

         T_on = 500 ms              T_off = 2000 ms
    ◀─────────────────────▶   ◀──────────────────────────▶

    ┌─────────────────────┐
L   │████ PZT Left ███████│                               (741 kHz CW)
    └─────────────────────┘───────────────────────────────

    ┌─────────────────────┐
R   │████ PZT Right ██████│                               (741 kHz CW, offset φ)
    └─────────────────────┘───────────────────────────────

    ┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐┌┐
CAM ││││││││││││││││││││││││││││││││││││││││││││││││││││
    └┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘└┘
    |← pre-trigger →|← US ON →|← recovery imaging  →|

Phase sweep protocol (additional, Part B only):
    Step φ from 0° to 180° in N steps between acquisitions.
    At each φ: acquire one full pulse train, extract ε_peak(φ).
    Compare to model prediction: σ₀(φ) = σ₀_max · cos²(φ/2)
```

### Timing configuration

| Parameter          | Part A value                               | Part B value                             | Notes                                             |
| ------------------ | ------------------------------------------ | ---------------------------------------- | ------------------------------------------------- |
| Camera fps         | 155 (full) or 300+ (ROI)                   | Same                                     | Global shutter eliminates motion blur             |
| Camera exposure    | 2–4 ms                                     | Same                                     | Brightfield OK at ~2 ms                           |
| US pulse ON        | 500 ms (primary), also test 100 ms and 2 s | Same                                     | Captures both τ₁ and τ₂                           |
| US pulse OFF       | 2–5 s (short), 30 s (long recovery)        | Same                                     | Must exceed τ₂ = 10 s for full recovery test      |
| Pre-trigger frames | 50 frames (~320 ms)                        | Same                                     | Baseline before deformation                       |
| Number of pulses   | 5–10 per acquisition                       | Same                                     | Track residual strain accumulation                |
| Total acquisition  | 25–55 s per pulse train                    | Same                                     | Depends on protocol                               |
| Phase steps        | N/A                                        | 0°, 30°, 60°, 90°, 120°, 150°, 180°      | One full acquisition per phase step               |
| RF channels        | 1                                          | 2 (synchronized gate, independent phase) | Both channels triggered from same delay generator |

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

For measuring diameter changes at 20×: ~37 pixels per 12 µm (Part A) or ~60 pixels per 19 µm (Part B biaxial) — well within sub-pixel tracking accuracy.

### Data rate and storage

| Config                     | Data rate | 60 s acquisition |
| -------------------------- | --------- | ---------------- |
| 1920×1200 @ 155 fps, Mono8 | ~358 MB/s | ~21 GB           |
| 800×800 @ 250 fps, Mono8   | ~160 MB/s | ~9.6 GB          |
| 512×512 @ 350 fps, Mono8   | ~92 MB/s  | ~5.5 GB          |

Use SSD (NVMe preferred). USB3 bandwidth is ~350 MB/s max.

---

## Acoustic Channel Design

### Part A — Single Transducer (20 MHz)

Based on model results at 20 MHz:

| Parameter          | Recommended                    | Why                                          |
| ------------------ | ------------------------------ | -------------------------------------------- |
| Channel length L   | 2.5 mm                         | Δf = 300 kHz spacing (robust to drift)       |
| Channel depth      | 100–200 µm                     | Fits 200 µm organoid; supports standing wave |
| End condition      | Glass–glass bonded (rigid)     | Maximizes reflection coefficient             |
| Channel width      | 500 µm–1 mm                    | Allows organoid entry + rotation             |
| Couplant thickness | < 7.5 µm (λ/10 at 20 MHz)      | Minimizes coupling sensitivity               |
| Couplant type      | US gel (Aquasonic 100)         | Or UV epoxy for stability                    |
| Substrate          | Borosilicate glass, 1 mm       | Matches model parameters                     |
| Cover              | Borosilicate coverslip, 150 µm | Matches model; fits BX51WI WD                |

At 20 MHz in water: λ/2 = 37.1 µm. In a 200 µm deep channel, ~5 pressure nodes span the depth. The organoid (200 µm diameter) spans the full channel depth, experiencing integrated radiation pressure — consistent with the effective stress approach in the organoid model.

### Part B — Dual PZT Config B (741 kHz)

Based on model results at 741 kHz (BENG207_2 Config B):

| Parameter          | Recommended                                          | Why                                                                 |
| ------------------ | ---------------------------------------------------- | ------------------------------------------------------------------- |
| Channel length L   | 20 mm                                                | Fits ~20 organoids at λ/2 = 1 mm spacing                            |
| Channel depth      | 1 mm                                                 | Current design; depth >> organoid diameter (200 µm)                 |
| PDMS side walls    | ≤3 mm (thin from current 10 mm)                      | T_power = 77% at 3 mm vs 43% at 10 mm (see BENG207_2 Fig 12)        |
| PDMS ceiling       | 1–2 mm (thin from current 9 mm)                      | BX51WI objective WD = 3.5 mm — must image through ceiling           |
| Bottom substrate   | Borosilicate glass, 1 mm                             | θ_c = 15.3° → TIR for lateral wave (essential, not just structural) |
| PZT position       | Immersed in water, 5 mm water gap                    | No couplant gel needed; no aging/degradation                        |
| PZT type           | PZT-5A, thickness-mode resonant ~741 kHz             | LiNbO₃ underpowered at this frequency                               |
| Node spacing       | λ/2 = 1.001 mm                                       | Matches CI electrode spacing target                                 |
| Number of nodes    | ~20 (20 mm / 1 mm)                                   | One organoid per node                                               |
| Transducer spacing | 20 mm + 2×(3 mm PDMS) + 2×(5 mm water) = 36 mm total | Center-to-center across channel                                     |

**PDMS transmission at 741 kHz (5 dB/cm/MHz attenuation model):**

| Wall thickness  | Loss per wall | T_power | Force retained | Recommendation        |
| --------------- | ------------- | ------- | -------------- | --------------------- |
| 10 mm (current) | 3.7 dB        | 0.43    | 43%            | Too thick — must thin |
| 5 mm            | 1.9 dB        | 0.65    | 65%            | Acceptable minimum    |
| 3 mm            | 1.1 dB        | 0.77    | 77%            | **Recommended**       |
| 1 mm            | 0.4 dB        | 0.92    | 92%            | Ideal if moldable     |

**Note on PDMS attenuation measurement:** The table above uses the model value of 5 dB/cm/MHz (Carugo et al. 2012), but actual attenuation varies 3–12 dB/cm/MHz depending on PDMS base:crosslinker ratio, cure temperature, and age. Measure after fabrication using the substitution method before committing to a wall thickness:

1. Drive one PZT into water only (no PDMS wall) at a fixed voltage. Image bead node density or standing wave contrast with the Basler camera.
2. Insert the PDMS wall. Repeat at the same drive voltage.
3. The ratio of node densities or contrast values gives T_power directly — no hydrophone needed.

If measured T_power is lower than expected (e.g., T_power < 0.6 at 3 mm wall), attenuation is higher than 5 dB/cm/MHz. Options: thin the wall further, or remake PDMS at 20:1 base:crosslinker ratio instead of 10:1 for lower attenuation. A calibrated hydrophone (see Instrumentation Note in the README) gives absolute pressure and is needed for Protocol D drag-balance stress calibration, but is not required for the T_power characterization.

---

## Experimental Protocols

### Part A Protocols

#### Protocol A: Fast timescale identification (τ₁)

```
Purpose:   Capture fast viscoelastic response (cortex deformation)
Config:    Part A (single transducer, 20 MHz)
US pulse:  100 ms ON, 2 s OFF, × 5 pulses
Camera:    300+ fps (use ROI), ~15 s total acquisition
Objective: 20× water immersion (UMPLFLN20XW) for maximum resolution
Analysis:  Fit exponential rise/fall to get τ₁ and A₁ = σ₀/E₁
Expected:  τ₁ ~ 0.1 s → ~30 frames to capture rise at 300 fps
```

#### Protocol B: Slow timescale identification (τ₂)

```
Purpose:   Capture slow relaxation (ECM, cell rearrangement)
Config:    Part A (single transducer, 20 MHz)
US pulse:  2 s ON, 30 s OFF, × 3 pulses
Camera:    50 fps (full frame OK), ~100 s total acquisition
Objective: 10× (wider field, track any drift)
Analysis:  Fit long recovery tail to get τ₂ and A₂ = σ₀/E₂
Expected:  τ₂ ~ 10 s → need >30 s recovery window
```

#### Protocol C: Residual strain accumulation

```
Purpose:   Demonstrate two-timescale ratcheting behavior
Config:    Part A (single transducer, 20 MHz)
US pulse:  500 ms ON, 2 s OFF, × 10 pulses
Camera:    155 fps (full frame), ~25 s total
Objective: 20× for precision
Analysis:  Track rising baseline between pulses
Expected:  Single KV → no baseline rise
            Modified KV → baseline increases by ~1–2% over 10 pulses
```

#### Protocol D: Stress calibration (drag balance)

```
Purpose:   Calibrate σ₀ independently of model
Config:    Part A (single transducer, 20 MHz)
Method:    Apply flow at known velocity U while US is on.
           At equilibrium: F_ac = 6πμRU
           σ₀ = F_ac / (πR²)
Camera:    155 fps, track organoid displacement
Flow:      Syringe pump, 0.1–10 µL/min through channel
Analysis:  Find U where organoid just escapes trap
```

### Part B Protocols

#### Protocol E: Dual-PZT standing wave verification (beads first)

```
Purpose:   Verify counter-propagating standing wave before using organoids.
           Always validate with polystyrene beads (Phi ~ 0.22) before
           switching to organoids (Phi ~ 0.059 — 4× harder to trap).
Config:    Part B, φ = 0° (symmetric), 741 kHz
Sample:    10 µm polystyrene beads in water
US pulse:  CW (continuous), low power (start at 10% drive)
Camera:    155 fps, full frame
Objective: 10× (wide field to see all ~20 nodes simultaneously)
Expected:  Beads migrate to pressure nodes at 1 mm spacing within 1–5 s.
           All nodes should fill simultaneously and symmetrically.
           If nodes are uneven → check amplitude balance (p_L ≠ p_R).
           If no nodes form → check frequency, impedance matching, phase.
```

#### Protocol F: Biaxial deformation measurement

```
Purpose:   Measure biaxial KV deformation in dual-PZT Config B.
           Compare D_x (compression) and D_y (extension) independently.
Config:    Part B, φ = 0°, p_L = p_R = 500 kPa each
US pulse:  500 ms ON, 2 s OFF, × 10 pulses
Camera:    250 fps (ROI 800×800), ~25 s total
Objective: 20× water immersion — image from above through PDMS ceiling (≤2 mm)
Analysis:  Extract major axis a(t) [D_x, along acoustic axis] and
           minor axis b(t) [D_y, transverse axis] independently.
           Verify D_y extension ~ 2× larger than Part A (biaxial vs uniaxial).
Expected:  σ₀ = 26.8 Pa, peak ε ~ 9.7% (500 ms pulse)
           D_x compressed by ~19.4 µm
           D_y extended by ~19.4 µm (biaxial disk geometry)
           Contrast with Part A: D_y extended by only ~9.7 µm (uniaxial)
```

#### Protocol G: Phase sweep — node calibration

```
Purpose:   Verify node position control and calibrate σ₀(φ).
           Use organoid deformation as a readout of acoustic field strength.
Config:    Part B, step φ from 0° to 180° in 30° increments (7 steps)
US pulse:  500 ms ON at each phase step (single pulse per step)
Camera:    155 fps, track ε_peak at each φ
Objective: 20× for single-organoid tracking
Analysis:  Plot ε_peak vs φ. Fit to:
           ε_peak(φ) = ε_max · cos²(φ/2)
           ε_max should match model: 9.7% at φ=0°, 0% at φ=180°.
Expected:  φ = 0°:   ε_peak ~ 9.7%,  node shift = 0
           φ = 90°:  ε_peak ~ 4.85%, node shift = 0.25 mm
           φ = 180°: ε_peak ~ 0%,    node shift = 0.5 mm (λ/4)
           This is the primary calibration experiment for CI electrode alignment.
```

#### Protocol H: Amplitude imbalance sensitivity

```
Purpose:   Characterize effect of p_L ≠ p_R (asymmetric drive).
           From BENG207_2 Fig 11: ratio < 0.8 degrades contrast below ~10.
Config:    Part B, fix p_R = 500 kPa, vary p_L / p_R from 0.5 to 1.0
US pulse:  500 ms ON, 2 s OFF, × 5 pulses per amplitude ratio
Camera:    155 fps
Analysis:  Track node uniformity (all 20 organoids should deform equally).
           Uneven deformation → amplitude imbalance → re-balance drive.
Expected:  Ratio > 0.9: all nodes uniform (good)
           Ratio < 0.8: some nodes weaker → organoid loss at those nodes
```

---

## Signal Chain: Wiring Diagram

### Part A — Single Transducer

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

TRANSDUCER → [US gel couplant] → glass slide → fluid channel
                                                    │
                                                    ▼
                                        Olympus BX51WI objective
                                                    │
                                              Basler camera
                                                    │
                                              USB3 → PC
                                         (Pylon Viewer or
                                          custom Python/MATLAB)
```

### Part B — Dual PZT Config B

```
DELAY/PULSE GENERATOR (BNC 575 or SRS DG535)
├── Ch A: Master trigger
│         → Camera trigger input (GPIO Line1)
│
├── Ch B: US gate pulse
│         Duration: T_on, repeated N_pulses
│         → Function generator Ch1 GATE (left PZT)
│         → Function generator Ch2 GATE (right PZT)
│         (Both channels gated simultaneously from same TTL)
│
├── Ch C: (optional) Uniblitz shutter sync
│
└── Ch D: (optional) Phase step trigger
          → External phase controller input (if separate hardware)

DUAL-CHANNEL FUNCTION GENERATOR (741 kHz)
├── Ch1 output: 741 kHz CW, phase = 0°
│              Amplitude: ~100 mVpp → RF Amp Left
│
├── Ch2 output: 741 kHz CW, phase = φ (0° to 180°, set manually per step)
│              Amplitude: ~100 mVpp → RF Amp Right
│
└── Sync out → oscilloscope (monitor both channels, verify phase offset)

RF AMPLIFIER — LEFT
├── Input:  Function generator Ch1
├── Output: 50 Ω coax → PZT-5A Left (immersed in water, left side)
└── Monitor: directional coupler → oscilloscope Ch1

RF AMPLIFIER — RIGHT
├── Input:  Function generator Ch2
├── Output: 50 Ω coax → PZT-5A Right (immersed in water, right side)
└── Monitor: directional coupler → oscilloscope Ch2

PZT Left  → [5 mm water gap] → PDMS wall (≤3 mm) → fluid channel
PZT Right → [5 mm water gap] → PDMS wall (≤3 mm) → fluid channel (opposite side)

                     Counter-propagating standing wave
                     ← p_L        p_R →
                     ~20 organoids at nodes (1 mm spacing)

                     Borosilicate glass slide (1 mm) below channel
                     TIR θ_c = 15.3° → perfect vertical confinement

                                    │
                                    ▼
                        Olympus BX51WI (upright)
                        Images from ABOVE through PDMS ceiling (1–2 mm)
                        UMPLFLN20XW objective (WD = 3.5 mm)
                                    │
                              Basler camera
                                    │
                              USB3 → PC
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
a0 = np.mean(a[pre_trigger_frames])   # baseline major axis (acoustic axis = D_x)
b0 = np.mean(b[pre_trigger_frames])   # baseline minor axis (transverse = D_y)

epsilon_a = (a - a0) / a0              # axial strain (compression, negative)
epsilon_b = (b - b0) / b0              # transverse strain (extension, positive)

# Volume conservation check
# Part A (uniaxial):  V_proxy = a * b^2   → should stay ~1.0
# Part B (biaxial):   V_proxy = a * b^2   → same check (b = D_y = D_z by symmetry)
V_proxy = a * b**2
V_ratio = V_proxy / V_proxy[0]        # should stay ~1.0 if incompressible

# Part B geometry check: transverse strain should be ~ 2x Part A value
# For same sigma0: Part A gives epsilon_b ~ epsilon_a/2
#                  Part B gives epsilon_b ~ epsilon_a  (biaxial disk)
# If epsilon_b < epsilon_a/2 in Part B: possible asymmetry or wall reflection
```

### Step 3: Fit KV model parameters

```python
from scipy.optimize import curve_fit

def modified_kv_pulse(t, A1, tau1, A2, tau2, T_on):
    """Two-branch KV response to rectangular pulse."""
    eps = np.zeros_like(t)
    on  = t < T_on
    off = ~on
    eps[on]  = A1*(1 - np.exp(-t[on]/tau1)) + A2*(1 - np.exp(-t[on]/tau2))
    e1T = A1*(1 - np.exp(-T_on/tau1))
    e2T = A2*(1 - np.exp(-T_on/tau2))
    eps[off] = e1T*np.exp(-(t[off]-T_on)/tau1) + e2T*np.exp(-(t[off]-T_on)/tau2)
    return eps

# Fit — same function for both Part A and Part B
# (KV constitutive equation is configuration-independent)
popt, pcov = curve_fit(
    lambda t, A1, t1, A2, t2: modified_kv_pulse(t, A1, t1, A2, t2, T_on=0.5),
    t_data, epsilon_data,
    p0=[0.03, 0.1, 0.03, 10],
    bounds=([0, 0.01, 0, 0.5], [0.5, 2, 0.5, 100])
)
A1, tau1, A2, tau2 = popt
```

### Step 4: Convert to material properties (after σ₀ calibration)

```python
# After sigma0 is known (from Protocol D drag balance, or from model):
E1   = sigma0 / A1    # Pa   — fast branch stiffness
E2   = sigma0 / A2    # Pa   — slow branch stiffness
eta1 = E1 * tau1       # Pa·s — fast branch viscosity
eta2 = E2 * tau2       # Pa·s — slow branch viscosity

# Part B note: use sigma0_dual = Phi * (p_L + p_R)^2 / (rho_f * c_f^2)
# with Phi = 0.059, p_L = p_R = 500 kPa → sigma0_dual = 26.83 Pa
# This is 57% larger than Part A estimate (17.05 Pa), so E1, E2 will be
# correspondingly larger when computed from the same measured A1, A2.
```

### Step 5: Phase sweep analysis (Part B only)

```python
# Fit sigma0(phi) curve from Protocol G data
phi_deg   = np.array([0, 30, 60, 90, 120, 150, 180])
eps_peak  = np.array([...])  # measured peak strain at each phase

phi_rad   = phi_deg * np.pi / 180
sigma0_max_fit = sigma0_dual  # from model, or fit as free parameter

eps_model = (sigma0_max_fit / E_eff) * (1 - np.exp(-T_on / tau_eff)) \
            * np.cos(phi_rad / 2)**2

# Overlay measured vs model — deviation indicates:
#   Systematic offset → organoid not initially at node center
#   Asymmetric curve  → p_L ≠ p_R (amplitude imbalance)
#   Zero crossing ≠ 180° → verify actual wavelength / frequency
```

---

## Checklist Before Experiment

### Part A Checklist

- [ ] Transducer impedance matched (check with network analyzer if available)
- [ ] US gel applied thin (<10 µm) — use minimal amount, press firmly
- [ ] Channel filled, no air bubbles (critical for standing wave) — ask Spencer
- [ ] Organoid loaded, settled to channel bottom or suspended at node
- [ ] Camera ROI set, exposure verified, frames saving to SSD
- [ ] Delay generator programmed: pre-trigger → US gate → N pulses
- [ ] Oscilloscope monitoring: RF amplitude at transducer
- [ ] Baseline images acquired (50+ frames, US off)
- [ ] Test pulse at low power to verify acoustophoresis before full experiment
- [ ] Temperature noted (affects c_f and viscosity)
- [ ] Hg lamp hours checked — replace at >300 hrs (currently ~3.3 hrs, well within rated life)

### Part B Checklist (additional)

- [ ] PDMS side walls thinned to ≤3 mm — measure with calipers before mounting
- [ ] PDMS ceiling thinned to 1–2 mm — verify BX51WI objective clears with WD = 3.5 mm
- [ ] Both PZTs immersed in water, no air gaps between PZT face and water
- [ ] Left and right RF amplifier gains matched — measure p_L = p_R within ±5%
- [ ] Phase offset verified on oscilloscope: Ch1 vs Ch2 of function generator
- [ ] Standing wave confirmed with polystyrene beads BEFORE switching to organoids (Protocol E)
- [ ] All ~20 nodes fill uniformly with beads — if uneven, re-balance amplitudes
- [ ] Borosilicate glass slide in place underneath channel (TIR confinement)
- [ ] Imaging path confirmed: upright BX51WI objective through PDMS ceiling
- [ ] Phase sweep protocol programmed (φ = 0° → 180° in steps)
- [ ] Channel filled, no air bubbles (same as Part A)
- [ ] Baseline images acquired (50+ frames, both PZTs off)

---

## Expected Results (from Model Predictions)

### Part A — Single Transducer (σ₀ = 17.1 Pa at 500 kPa, uniaxial)

| Measurement                         | Expected value      | Model source             |
| ----------------------------------- | ------------------- | ------------------------ |
| Peak strain (500 ms pulse, 500 kPa) | ~6.2%               | Organoid KV model, Fig 1 |
| D_x diameter change (compression)   | ~12.4 µm            | Organoid KV model, Fig 6 |
| D_y diameter change (extension)     | ~6.2 µm             | Uniaxial: ε/2            |
| Fast time constant τ₁               | 0.05–0.3 s          | Cortex mechanics         |
| Slow time constant τ₂               | 5–30 s              | ECM/rearrangement        |
| Residual strain after 10 pulses     | 1–3%                | Modified KV accumulation |
| Minimum detectable ΔD (20× obj)     | ~0.3 µm (sub-pixel) | Camera resolution        |
| Minimum p_ac for 1 µm ΔD (uniaxial) | ~200 kPa            | Stress sweep, Fig 7      |

### Part B — Dual PZT Config B (σ₀ = 26.8 Pa at 500 kPa per PZT, biaxial)

| Measurement                              | Expected value         | Model source                  |
| ---------------------------------------- | ---------------------- | ----------------------------- |
| Peak strain (500 ms pulse, φ=0°)         | ~9.7%                  | Organoid KV model, Fig 9a     |
| D_x diameter change (compression)        | ~19.4 µm               | Biaxial: D₀·ε                 |
| D_y diameter change (extension)          | ~19.4 µm               | Biaxial: D₀·ε (disk)          |
| Stress increase vs Part A                | ×1.573 (57% more)      | Gorkov + superposition, Fig 8 |
| Node shift at φ=90°                      | 0.25 mm                | Phase sweep, Fig 10           |
| Node shift at φ=180°                     | 0.50 mm (λ/4)          | Phase sweep, Fig 10           |
| ε_peak at φ=90°                          | ~4.85% (half max)      | cos²(φ/2) model, Fig 10       |
| ε_peak at φ=180°                         | ~0% (organoid at null) | Node off organoid             |
| Residual strain after 10 pulses (Part B) | 1.5–5% (amplified)     | Higher σ₀ → more ratcheting   |
| Minimum p_ac for 1 µm ΔD (biaxial)       | ~115 kPa per PZT       | Stress sweep, Fig 8           |

---

## Notes and Known Issues

- **Hg lamp (U-RFL-T):** Currently at ~3.3 hours on a 300-hour rated bulb — well within safe operating life. Log hours after each session. Replace with USH-103OL ($270) at or before 300 hours, or upgrade to LED (Thorlabs SOLIS-470C, ~$1,000 including driver and adapter) when the bulb is due.
- **PDMS ceiling (Part B):** Must be thinned to 1–2 mm before Part B imaging is possible. Current 9 mm ceiling exceeds the BX51WI water immersion objective working distance (3.5 mm). Plan fabrication before scheduling imaging sessions.
- **PDMS attenuation (Part B):** Model assumes 5 dB/cm/MHz (Carugo et al. 2012). Actual value varies 3–12 dB/cm/MHz with cure ratio, temperature, and age — must be measured after fabricating thinned walls. Use the substitution method first (Basler camera + beads, no hydrophone needed): compare standing wave contrast with and without the PDMS wall at matched drive voltage to get T_power. If T_power < 0.6 at 3 mm wall, attenuation is above 8 dB/cm/MHz — thin further to 1–2 mm or remake at 20:1 base:crosslinker ratio. For absolute pressure calibration (Protocol D), a hydrophone with sub-1 MHz coverage is required: **Precision Acoustics HP1 or HP2** (1–2 mm PVDF probe, NPL-calibrated 30 kHz–40 MHz, ~£3,000–4,500 for complete system including preamplifier and DC coupler). Standard Onda HNA/HNP/HNR probes are not suitable at 741 kHz without a custom extended calibration — confirm with Onda before ordering.
- **C-mount adapter (Part A):** The U-TV1X-2 (1× relay) is already installed under U-CMAD3 on the BX51WI. The Basler acA1920-155um screws directly onto the U-CMAD3 — no additional adapter needed. This saves $250–500 (no need to buy U-TV0.63XC).
- **Amplitude balance (Part B):** Amplitude ratio p_L/p_R < 0.8 degrades standing wave contrast below ~10 (from infinite). Check balance by imaging bead distribution — all nodes should be equally populated.
- **Air bubbles:** Any air bubble in the channel will destroy the standing wave. Degas water/media before loading. Consult Spencer for bubble removal protocol.
