# 3BAtwin Microwave Tomography Hackathon

**Challenge:** implement a MATLAB microwave tomographic imaging algorithm and validate it against **real measured S-parameter data**.

**Development:** 25–26 November 2026 · **Final evaluation and pitches:** 27 November 2026  
**Submission deadline:** **Friday, 27 November 2026, 08:59 CET (Europe/Rome)** — strict code freeze.  
**Submission repository:** **[REPOSITORY LINK — TO BE INSERTED]**

## Challenge and teams

Develop an algorithm that turns measured multistatic scattering-parameter data into a **complex-valued 2D reconstruction** of an unknown target. Teams may use available reference datasets, documented a priori information and organizer-approved additional measurements during development. In the final evaluation, every team receives the **same unseen measurement input** and runs on the **same evaluation workstation**.

There will be **3–4 teams of 3–4 people**. Organizers will assemble mixed-institution teams, aim for gender balance, and assign a team lead from a participating institution, distributing institutional leads as evenly as possible. Mentor assistance is permitted; teams must understand and defend all submitted work.

## Required files

```text
team-repository/
├── main.m                     # entry-point script; runs the complete pipeline
├── reconstruct_tomography.m   # required reconstruction function
├── README.md                  # setup, options, method, toolbox requirements
└── ...                        # helpers, configuration, permitted model files
```

The input is a **Touchstone S-parameter filename** (for example `.s2p` or `.sNp`, as specified by organizers). The required function returns a finite **complex 2D MATLAB numeric matrix**:

```matlab
function image2d = reconstruct_tomography(touchstoneFile, opts)
% touchstoneFile: path to a Touchstone scattering-parameter file
% opts: optional struct with documented, permitted settings
% image2d: finite complex Ny-by-Nx reconstructed image
%          with organizer-defined physical meaning, units and grid
end
```

`main.m` must load the provided file, call the function, and display/save the reconstructed image without interactive code edits. Document any extra parameters in `opts` (such as geometry, frequency selection, background data or regularization). Do not assume extra private test inputs are available. The organizers will publish exact grid, units, target quantity, Touchstone port ordering, frequency conventions and MATLAB compatibility requirements. A usable graphical interface is **optional** and can receive bonus credit.

## Submission rules

1. Submit via the designated repository before **08:59 CET on 27 November 2026**.
2. The repository commit at the deadline is frozen; **no changes after the cutoff** will be accepted.
3. Include `main.m`, `reconstruct_tomography.m`, helper dependencies and a working README.
4. Declare MATLAB version/toolboxes and third-party software, pretrained components and AI-assisted development.
5. Organizers may offer equal-length compatibility test slots on the evaluation computer with **public test data only**.
6. Final tests use identical hidden measurements and identical CPU/GPU, memory, MATLAB, timeout and external-compute policies.

## Scoring (proposed: 100 base + up to 5 bonus)

| Category | Criterion | Max points |
|---|---|---:|
| Algorithm | Image accuracy via NRMSE | 40 |
| Algorithm | Execution time on common PC | 10 |
| Jury | Physical/inverse-problem understanding | 10 |
| Jury | Implementation and code defense | 10 |
| Jury | Interpretation and limitations | 5 |
| Jury | Integration and reproducibility | 5 |
| Jury | Pitch clarity and narrative | 8 |
| Jury | Oral delivery and questions | 7 |
| Jury | Slides and visual communication | 5 |
| **Base total** | | **100** |
| **Optional bonus** | Usable GUI / clear imaging workflow | **+5** |

Three jurors, **one representing each 3BAtwin partner university**, independently score jury criteria; their average determines the jury contribution. Jurors will interview the teams about how their algorithm works and how their code was produced. AI-tool use is not automatically disqualifying, but teams must disclose it and demonstrate genuine understanding.

### Reconstruction metric

The cited paper specifies **normalized** RMSE, calculated separately for the real and imaginary components of the complex reconstructed contrast:

```text
NRMSE(q) = sqrt( mean_j[ (A_j(q) - B_j(q))^2 ] ) / Δ(q)
Δ(q)     = min( max(A(q))-min(A(q)), max(B(q))-min(B(q)) )
q        = real or imaginary component
```

Here `A` is the reconstructed map, and `B` is the independently established reference map. **Lower NRMSE is better.** The paper additionally defines `NRMSE_N = 1 - NRMSE`. The organizers will specify how scores from the two components and multiple blind cases are aggregated, how undefined/near-zero ranges are handled, and how raw metric values map to the 40 available points. Runtime-to-points thresholds will likewise be announced before judging. No unannounced normalization or target-specific tuning is allowed.

**Reference:** A. R. Masaquiza-Caiza, D. O. Rodríguez-Duarte, M. Gugliermino, V. Mariano, and F. Vipiana, “On Anatomical Prior Kernel Distortion for Tomographic Microwave Imaging Stroke Monitoring,” *IEEE Open Journal of Antennas and Propagation*, 2026. https://doi.org/10.1109/OJAP.2026.3724057

## Jury interview and pitch

The jury assesses understanding of the electromagnetic forward model, inversion algorithm, numerical assumptions, code, results and limitations. Questions may cover the exact purpose of a function, justification of regularization, or changes required for a new target. The pitch should communicate the scientific objective, main idea, experimental results, novelty, limitations and future value through clear speaking and well-designed slides.

## Decisions to be confirmed by organizers

- [ ] Insert official repository link.
- [ ] Publish exact Touchstone input contract, output complex quantity, dimensions and grid.
- [ ] Publish image-scoring thresholds, two-component aggregation, zero-range policy and timeout.
- [ ] Confirm team rosters/leads, all three jurors, and time available for interviews and pitches.
- [ ] Confirm bonus GUI policy and official final total convention.
