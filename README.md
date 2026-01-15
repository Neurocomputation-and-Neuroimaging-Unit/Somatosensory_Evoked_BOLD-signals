# Somatosensory_Evoked_BOLD-signals
This is the code used for acquiring and analyzing the dataset of the experiment with the same name.

The folder-structure follows this order:
- Acquisition 
  - Test Stimuli
  - Experimental Script
- Analysis
  - A: data conversion
    - Dicom to BIDS (nifti format) conversion
    - Requirements: Dicom to BIDS toolbox
  - B: data preprocessing, including:
    - Realignment
    - Detrending
    - Reordering 
    - Co-registration
    - Normaization
    - Smoothing
    - Requirements: SPM12, hMRI toolbox
  - C: statistical analysis, including:
    - First-level FIR-analysis
    - First-level GLM
    - Second-level one-way ANOVA
    - Computation of 4D permutations
    - Requirements: SPM12, hMRI toolbox
  - D: evaluation of results, including:
    - 4D Permutation statistics
    - ROI-analysis tool
    - Requirements: SPM12, hMRI toolbox
    
