# Dataset and Code for HighDensity_ToolpathGene4CFRTP

<p align="center">
  <img src="teaser.jpg" alt="Teaser Image">
</p>

This repository accompanies the paper: *Tianyu Zhang, Tao Liu, Neelotpal Dutta, Yongxue Chen, Renbo Su, Zhizhou Zhang, Weiming Wang, and Charlie C.L. Wang*, **"Toolpath generation for high density spatial fiber printing guided by principal stresses,"** Composites Part B: Engineering_, vol. 295, 112154 (17 pages), April 2025. [[paper]](https://www.sciencedirect.com/science/article/pii/S1359836825000447?via%3Dihub) [[Video@YouTube]](https://www.youtube.com/watch?v=ylBgGtqyhDE)

---

## Abstract

While multi-axis 3D printing can align continuous fibers along principal stresses in continuous fiber-reinforced thermoplastic (CFRTP) composites to enhance mechanical strength, existing methods have difficulty generating toolpaths with high fiber coverage. This is mainly due to the orientation consistency constraints imposed by vector-field-based methods and the turbulent stress fields around stress concentration regions.  
This work addresses these challenges by introducing a 2-RoSy representation for computing the direction field, which is then converted into a periodic scalar field to generate partial iso-curves for fiber toolpaths with nearly equal hatching distance. To further improve fiber coverage in stress-concentrated regions, such as around holes, we extend the quaternion-based curved slicing method by incorporating winding compatibility considerations.  
Our proposed method achieves toolpath coverage between 87.5% and 90.6% for continuous fibers with 1.1 mm width. Models fabricated using these toolpaths show up to 84.6% improvement in failure load and a 54.4% increase in stiffness compared to multi-axis 3D printing with sparser fiber layouts.

---

## Repository Content
- DataSet of models (Curved layers \& CCF strip toolpath )
- Toolpath generation code

---

## Usage

The installation steps are the same as those for [[S3-Slicer]](https://github.com/zhangty019/S3_DeformFDM).

- Step 0: Click Button "0. Layer Offset" -> Conduct the offset of the curved layers;
- Step 1: Click Button "1. Read Data" -> Read the curved layers;
- Step 2: Click Button "2. Direction Field Optimization" -> Optimize the 2RoSy direction field;
- Step 3: Click Button "3. Strip Pattern Generation" -> Generate strip pattern;
- Step 4: Click Button "4. Strip Toolpath Generation" -> Extract toolpath from strip;
- Step 5: Click Button "5. Post-processing for Toolpath" -> Improve the continuity and smoothness of the toolpath;
- Step 6: Change the value of "Show Edge" to check the printing simulation of the CCF;

---

## Citation

If you find our work useful, please cite our paper:

```bibtex
@article{zhang2025toolpath,
  title={Toolpath generation for high density spatial fiber printing guided by principal stresses},
  author={Zhang, Tianyu and Liu, Tao and Dutta, Neelotpal and Chen, Yongxue and Su, Renbo and Zhang, Zhizhou and Wang, Weiming and Wang, Charlie C.L.},
  journal={Composites Part B: Engineering},
  volume={295},
  pages={112154},
  year={2025},
  publisher={Elsevier}
}

