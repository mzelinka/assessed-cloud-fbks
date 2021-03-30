# assessed-cloud-fbks

This code compares GCM cloud feedback components to expert-assessed feedbacks assessed by [Sherwood et al. (2020)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019RG000678). To use, follow these steps:

1. Install CDAT via conda following [these instructions](https://github.com/CDAT/cdat/wiki/install#installing-latest-cdat---821)

2. Activate this environment:
```
conda activate cdat
```

3. In main.py, update the dictionary so it points to your model's amip and amip-p4K files.

4. Run the code:
```
python main.py
```

5. Inspect the generated figures and tables in the /figures/ directory.
