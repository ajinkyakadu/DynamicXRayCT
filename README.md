# DynamicXRayCTVision


Welcome to DynamicXRayCTVision! This repository houses the MATLAB codebase for our project titled "Single-shot Tomography of Discrete Dynamic Objects". Our work introduces a novel method for the reconstruction of high-resolution temporal images in dynamic tomographic imaging. See our work [here](https://arxiv.org/abs/2311.05269)

## Dynamic Tomographic Reconstruction

In the realm of dynamic tomographic imaging, capturing discrete objects with smooth boundaries that vary over time has always been a challenging task, especially with limited measurements at each time point. DynamicXRayCTVision addresses this challenge head-on. Our approach synergistically combines spatial and temporal information of dynamic objects, leveraging the level-set method for image segmentation and a sinusoidal basis for motion representation. This fusion results in a computationally efficient and easily optimizable variational framework. Our method stands out by enabling the reconstruction of high-quality 2D or 3D image sequences with just a single projection per frame.

## Key Features

- **Temporal Image Reconstruction**: Enables high-resolution imaging with limited time-point measurements.
- **Level-Set Method for Image Segmentation**: Accurate and efficient segmentation for objects with smooth, evolving boundaries.
- **Sinusoidal Motion Representation**: Captures the dynamics of objects in motion effectively.
- **Variational Framework Optimization**: Offers computational efficiency and ease of optimization.
- **Applicability to Xray Datasets**: Demonstrated superior performance on both synthetic and real X-ray tomography datasets.

## Getting Started

Please refer to the [installation instructions] and [usage guidelines] for setting up and running the code. The repository is structured as follows:

- `/src`: Source code of the proposed method.
- `/examples`: Sample scripts demonstrating the usage of the method on different datasets.
- `/docs`: Documentation and additional resources.

## Contributions and Feedback

We welcome contributions and feedback from the community. Please feel free to raise issues, suggest enhancements, or contribute to the codebase. See our [contribution guidelines] for more information.
