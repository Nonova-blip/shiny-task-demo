## Prototype versions

- `signal_task_v1_stage2_first.Rmd`  
  Earlier feasibility demo. Generated the fuller Stage 2 structure first, then selected a Stage 1 subset. Useful as design history, but may make Stage 1 too diagnostic because the subset is filtered for slope preservation.

- `signal_task_v2_stage1_first_deming_constrained.R`  
  Current design candidate. Generates Stage 1 first, then adds independent Stage 2 evidence from the same latent line. Uses Deming-regression checks to keep realized slopes close enough to the intended target while avoiding Stage 1 subset-selection artifacts.
