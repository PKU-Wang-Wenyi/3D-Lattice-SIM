# 3D-Lattice-SIM

This is the code repository for the article **"Gentle live-cell super-resolution volumetric imaging with three-dimensional lattice structured illumination microscopy"**.

Currently, the package contains:
- `.m` source code
- `.exe` files containing GUI of the 3D-Lattice-SIM reconstruction code
- We have also developed a reconstruction scheme tailored for orthogonal lattice 3D-SIM data acquired with the ZEISS Elyra 7 system.
## 💿️  System Requirements

The analytical reconstruction software is implemented with MATLAB 2025a and tested in Windows 10 & 11 environments.

## Runtime Requirements
To run `.exe` files, **MATLAB Runtime R2025a (25.1)** is necessary.

**Download link:** [https://ww2.mathworks.cn/products/compiler/matlab-runtime.html](https://ww2.mathworks.cn/products/compiler/matlab-runtime.html)
## 💻 Hardware Requirements
A computing device equipped with an NVIDIA GPU is essential to run this software. Our tests were conducted on a system with an Intel Core i9-14900K CPU, 128 GB of RAM, and an NVIDIA GeForce RTX 5090 D GPU with 32 GB of VRAM. To ensure a fair comparison, other software such as Open3DSIM and FO3DSIM were also executed and benchmarked on the same computer.
For optimal runtime performance, it is strongly recommended to store the data on a solid-state drive (SSD) to accelerate I/O speed. Otherwise, data reading and writing can become a major time bottleneck, sometimes even exceeding the reconstruction time itself.
## 🌐 Datasets
The example 3D-Lattice-SIM imaging data are available at the figshare repository: https://doi.org/10.6084/m9.figshare.30675041.  

For more comprehensive data, including the dual-color long-term imaging datasets, please refer to Zenodo.[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17774000.svg)](https://doi.org/10.5281/zenodo.17774000)
