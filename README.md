# Fast and stable schemes for non-linear osmosis filtering

**Author:** Giuseppe Antonio Recupero

**Other Authors:** Luca Calatroni, Serena Morigi, Simone Parisotto

**Version 1.0**

**Date**: 19/06/22

This software implements the algorithms described in the paper [Fast and stable schemes for non-linear osmosis filtering](FILL), also available on [arXiv](https://arxiv.org/abs/2203.15570). Please cite as:

FILL

## Task 1: shadow removal

Given an image with a shadow/light-spot (left) and its relative boundary (middle), the result (right) is an image where the shadow/spot-light has been removed.

![17_soft](https://user-images.githubusercontent.com/103272764/174474017-d12c5095-c2a1-4f40-86a6-7650ff793633.png)
![17](https://user-images.githubusercontent.com/103272764/174474021-c24a49cb-a89a-40ad-91e8-469b10c91e29.png)
![17_soft_edge4_p1 99_gSI_tau1000_T10000_eps1e-07](https://user-images.githubusercontent.com/103272764/174474010-32a4a1c7-604f-429e-9c81-ea864e92cccd.png)

<img src="https://user-images.githubusercontent.com/103272764/174474053-9b472e8f-b83a-4ae5-a074-b8e57c9dfe3f.png" width="100" height="100">
<img src="https://user-images.githubusercontent.com/103272764/174474058-1e9a7047-dc22-4065-ae97-fd1727092fa5.png" width="100" height="100">
<img src="https://user-images.githubusercontent.com/103272764/174474046-31b72fcd-b3bd-4e09-8f3f-ecfd7f135584.png" width="100" height="100">

## Task 2: compact data representation

Given an image (left), the algorithm computes the mask of its edges (middle). Then, using only the information of the drift term stored on the mask, it reconstruct an approximation of the original image (right).

![lisa2](https://user-images.githubusercontent.com/103272764/174474256-c7b2838e-8b4a-4a56-a88c-8f868ce467d8.png)
![lisa2edge3](https://user-images.githubusercontent.com/103272764/174474250-cbb7dc99-2b3f-4b19-b673-7c4b4f38357e.png)
![lisa2_g3_tau1000T74000eps1](https://user-images.githubusercontent.com/103272764/174474237-1f0ba043-6da8-426c-94bf-18005148542b.png)

![fry2000](https://user-images.githubusercontent.com/103272764/174474331-2c6834e4-853f-4e0a-9fac-074603a28540.png)
![fry2000 edge0 008_1](https://user-images.githubusercontent.com/103272764/174474322-3dbca8be-50ae-4bbd-ad15-338cf241fd07.png)
![fry2000 _g1_tau10T140eps0 0001](https://user-images.githubusercontent.com/103272764/174474307-528d9851-882a-4ab4-9b9a-8300d0d0b8a0.png)

## Related code

The provided algorithm implements also the linear osmosis described in [A Fully Discrete Theory for Linear Osmosis Filtering](https://link.springer.com/chapter/10.1007/978-3-642-38267-3_31), by simply switching on the "flag_linear" variable in **main_osmosis_removal** and **main_osmosis_cdr**.

The code related to the comparison model [Anisotropic osmosis filtering for shadow removal in images](https://iopscience.iop.org/article/10.1088/1361-6420/ab08d2/meta) is available at [GitHub](https://github.com/simoneparisotto/anisotropic-osmosis-filter).

## License
[FILL](https://choosealicense.com/licenses/mit/)
