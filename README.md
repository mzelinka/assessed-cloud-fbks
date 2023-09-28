# assessed-cloud-fbks

[![DOI](https://zenodo.org/badge/353136800.svg)](https://zenodo.org/badge/latestdoi/353136800)

## Description
A Jupyter notebook is provided that performs the analysis of [Zelinka et al. (2022)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JD035198). It computes GCM cloud feedback components and compares them to the expert-assessed values from [Sherwood et al. (2020)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019RG000678). Figures are generated in the notebook and also saved to the [figures directory](https://github.com/mzelinka/assessed-cloud-fbks/tree/main/figures).

## Packages Needed
----------
- [xcdat](https://xcdat.readthedocs.io/en/stable/)
- xarray
- numpy
- scipy
- matplotlib<br>
...all of which can be installed via conda:
```
conda create -n <ENV_NAME> -c conda-forge xcdat xesmf xarray numpy matplotlib 
conda activate <ENV_NAME>
```

## References
- Zelinka et al. (2022): [Evaluating climate models’ cloud feedbacks against expert judgement](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JD035198), <em>J. Geophys. Res.</em>, 127, e2021JD035198, doi:10.1029/2021JD035198.

- Sherwood et al. (2020): [A combined assessment of Earth’s climate sensitivity](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019RG000678), <em>Rev. Geophys.</em>, 58, e2019RG000678, doi:10.1029/2019RG000678.
