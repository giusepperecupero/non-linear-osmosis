# Fast and stable schemes for non-linear osmosis filtering

**Author:** Giuseppe Antonio Recupero

**Other Authors:** Luca Calatroni, Serena Morigi, Simone Parisotto

**Version 1.0**

**Date**: 19/06/22

This software implements the algorithms described in the paper [Fast and stable schemes for non-linear osmosis filtering](FILL), also available on [arXiv](https://arxiv.org/abs/2203.15570). Please cite as:

FILL

## Task 1: shadow removal

Given an image with a shadow/light-spot (left) and its relative boundary (middle), the result (right) is an image where the shadow/spot-light has been removed.


![16_soft](https://user-images.githubusercontent.com/103272764/174473968-3c88e96d-b5bb-4fb9-a483-996f82a7b4e1.png)
![16](https://user-images.githubusercontent.com/103272764/174473977-8a9cf7f0-8cd2-407c-99b7-1c912c90bbe3.png)
![16_soft_edge2_p1 99_gSI_tau1000_T100000_eps1e-07](https://user-images.githubusercontent.com/103272764/174473959-0f7525e6-7345-4666-8193-fa7c8daa8787.png)


## Task 2: compact data representation

Given an image (left), the algorithm computes the mask of its edges (middle). Then, using only the information of the drift term stored on the mask, it reconstruct an approximation of the original image (right).



```bash
pip install foobar
```

```python
# returns 'words'
foobar.pluralize('word')
```

## Related code

The provided algorithm implements also the linear osmosis described in [A Fully Discrete Theory for Linear Osmosis Filtering](https://link.springer.com/chapter/10.1007/978-3-642-38267-3_31), by simply switching on the "flag_linear" variable in **main_osmosis_removal** and **main_osmosis_cdr**.

The code related to the comparison model [Anisotropic osmosis filtering for shadow removal in images](https://iopscience.iop.org/article/10.1088/1361-6420/ab08d2/meta) is available at [GitHub](https://github.com/simoneparisotto/anisotropic-osmosis-filter).

## License
[FILL](https://choosealicense.com/licenses/mit/)
