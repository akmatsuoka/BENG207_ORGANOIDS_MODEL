<pre>
Title: BENG207_Organoid_KV_Model
Author: BENG207 students (Instructor: Akihiro J. Matsuoka, M.D., D.M.Sc., Ph.D., FACS, Co-instructor: James Friend, Ph.D.)
Course: BENG207 Winter/Spring Semester (2026)
Status: Draft
Type: Analytical model
License: Public domain
Further discussion: TBD
</pre>

Instruction to BENG207 students: Please review all of the variables/parameters.

---

## Table of Contents

* [Abstract](#abstract)
* [Background and Significance](#background-and-significance)
* [Deliverable](#deliverable)
* [BENG207 Organoid KV Model](#beng207-organoid-kv-model-analytical-model)
* [Parameters](#parameters)
* [Output Figures: Captions](#output-figures-captions)
* [References](#references)

---

## Abstract

TBW

## Background and Significance

### Why Model Organoid Deformation?

The BENG207_1 transfer-matrix model describes how acoustic energy enters the fluid channel through the vertical stack. The BENG207_2 rigid cavity resonator model describes how that energy is organized into a lateral standing wave with pressure nodes at ~1 mm spacing, where iPSC organoids form. However, neither model addresses what happens to the organoid *itself* once it is trapped at a pressure node.

An organoid sitting at a pressure node experiences a time-varying acoustic radiation stress. Unlike a rigid bead, the organoid is a soft, viscoelastic biological structure — it deforms under this stress. The deformation has practical consequences for the biohybrid cochlear implant:

* **Too much deformation** may disrupt cell-cell junctions, collapse internal lumens, or trigger mechanotransduction pathways that alter differentiation — potentially converting neural progenitors into unintended cell types.
* **Too little deformation** may indicate the acoustic field is not coupling efficiently to the organoid, or that the organoid has stiffened pathologically (e.g., fibrosis).
* **Residual (non-recovered) deformation** after repeated acoustic pulses may indicate permanent structural damage or plastic remodeling of the extracellular matrix (ECM).

This model provides the quantitative framework to predict organoid deformation as a function of acoustic parameters (pressure amplitude, pulse duration, duty cycle) and organoid mechanical properties (elastic modulus, viscosity, relaxation times). It answers the design question: **for a given acoustic protocol, how much does the organoid deform, and does it fully recover?**

### The Kelvin-Voigt Viscoelastic Framework

Biological tissues are viscoelastic — they exhibit both elastic (spring-like, instantaneous) and viscous (dashpot-like, time-dependent) behavior. The simplest model that captures both is the Kelvin-Voigt (KV) element: a spring (modulus E) in parallel with a dashpot (viscosity η). Under a step stress σ₀, the KV element produces a creep response:

$$
\varepsilon(t) = \frac{\sigma_0}{E}\left(1 - e^{-t/\tau}\right), \quad \tau = \frac{\eta}{E}
$$

where τ is the **relaxation time** — the characteristic timescale for the organoid to approach its equilibrium deformation. When the stress is removed, the strain recovers exponentially with the same time constant:

$$
\varepsilon(t) = \varepsilon_T \, e^{-(t-T)/\tau}
$$

where ε_T is the strain at the moment the stress was removed (t = T). This complete recovery is a defining feature of the KV model — it predicts **no permanent deformation**, only delayed elastic response.

### Single KV vs. Modified (Two-Branch) KV

A single KV element has one relaxation time τ. Real organoids, however, exhibit at least two distinct deformation timescales:

1. **Fast response (τ₁ ~ 0.1 s):** This reflects the deformation of individual cells — cortical actin, cell membrane bending, and cytoplasmic flow. When you poke a single cell with an AFM tip, it deforms within ~100 ms.
2. **Slow response (τ₂ ~ 10 s):** This reflects the rearrangement of cells relative to each other and the remodeling of the extracellular matrix. The organoid as a multicellular structure takes much longer to reach mechanical equilibrium than a single cell.

The modified KV model places two KV elements in series (technically, the stress is shared and strains add):

$$
\varepsilon(t) = \frac{\sigma_0}{E_1}\left(1 - e^{-t/\tau_1}\right) + \frac{\sigma_0}{E_2}\left(1 - e^{-t/\tau_2}\right)
$$

The fast branch (E₁ = 500 Pa, τ₁ = 0.1 s) captures the initial rapid deformation. The slow branch (E₂ = 30 Pa, τ₂ = 10 s) captures the continued creep that occurs over seconds to tens of seconds. The slow branch is softer (lower E₂) because intercellular rearrangement offers less resistance than intracellular deformation.

This two-timescale behavior has practical consequences for acoustic pulse design:

* **Short pulses (T_on < τ₁):** Only the fast branch responds. The organoid barely deforms.
* **Moderate pulses (τ₁ < T_on < τ₂):** The fast branch is fully loaded; the slow branch is partially loaded. Most of the deformation comes from the fast branch.
* **Long pulses (T_on > τ₂):** Both branches reach equilibrium. Maximum deformation occurs.
