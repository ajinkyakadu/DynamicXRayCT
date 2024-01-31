<p align="center">
  <p align="center">
    <h1 align="center">DynamicXRayCT: Tomographic Reconstruction of Discrete Dynamic Objects </h1>
  </p>
  <p align="center" style="font-size:16px">
    <a target="_blank" href="https://ajinkyakadu.github.io/"><strong>Ajinkya Kadu</strong></a>
    ·
    <a target="_blank" href="https://felixlucka.github.io/"><strong>Felix Lucka</strong></a>
    ·
    <a target="_blank" href="https://scholar.google.com/citations?user=YYH0BjEAAAAJ&hl=en"><strong>Joost Batenburg</strong></a>
  </p>
  <h2 align="center">IEEE TCI</h2>
  <div align="center"></div> 

  <p align="center">
    <a href="https://astra-toolbox.com/"><img alt="ASTRA-TOOLBOX" src="https://img.shields.io/badge/astra%20toolbox-8A2BE2"></a>
    <a href="https://www.mathworks.com/products/matlab.html"><img alt="MATLAB" src="https://img.shields.io/badge/MATLAB-FFFF00"></a>
    <br>
    <a href='https://arxiv.org/pdf/2311.05269.pdf'>
      <img src='https://img.shields.io/badge/Paper-PDF-green?style=for-the-badge&logo=arXiv&logoColor=green' alt='Paper PDF'>
    </a>
  </p>
<p align="center">



Welcome to **DynamicXRayCT**! We introduce a novel method that integrates level-set method and compressed sensing for the reconstruction of temporal images in dynamic tomographic imaging.

## 🚀 Dynamic Tomographic Reconstruction

In dynamic tomographic imaging, capturing discrete objects with smooth boundaries that vary over time has always been a challenging task, especially with limited measurements at each time point. **DynamicXRayCT** addresses this challenge by advanced algorithm. Our approach combines spatial and temporal information of dynamic objects, using the level-set method for image segmentation and a sinusoidal basis for motion representation. This results in a computationally efficient and easily optimizable variational framework called Dynamic Shape Sensing (DSS). Our method stands out by enabling the reconstruction of high-quality 2D or 3D image sequences with just a single projection per frame.

## 🌈 Key Features

- **Temporal Image Reconstruction**: Achieve high-resolution imaging with limited measurements.
- **Level-Set Image Segmentation**: Segmentation for smoothly evolving objects.
- **Sinusoidal Motion Capture**: Captures the dynamics of objects in motion.
- **Optimized Variational Framework**: Offers computational efficiency and ease of optimization.
- **Xray Dataset Compatibility**: Demonstrate better performance on both synthetic and real X-ray tomography datasets.

## 🛠 Getting Started

Please refer to the following for setting up and running the code:

- [Installation Instructions](#installation-instructions)
- [Usage Guidelines](#usage-guidelines)

Repo Structure:

- `📁 /src`: the source codes.
- `📁 /examples`: CT applications and examples.
- `📁 /docs`: documentation and resources.

## 📚 Citation  
```bibtex
@article{Kadu:DynamicTomo:2023,
  title={Single-shot Tomography of Discrete Dynamic Objects},
  author={Kadu, Ajinkya and Lucka, Felix and Batenburg, Kees Joost},
  journal={arXiv preprint arXiv:2311.05269},
  year={2023}
}
```

## ©️ License
The code and models are available for use without many restrictions. 
See the [LICENSE](LICENSE) file for details. 


## 💡 Contributions and Feedback
Please contact [Ajinkya Kadu](https://ajinkyakadu.github.io/) for any questions. We welcome contributions and feedback from the community. We invite you to:

- 🐛 Report Issues
- 🌟 Suggest Enhancements
- 🤝 Contribute to the Code

Check our [Contribution Guidelines](#contribution-guidelines) for more.
