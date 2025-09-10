### 📐 Space Intersection – Photogrammetric 3D Point Computation

This repository implements a basic space intersection workflow using photogrammetric principles. It is designed to compute 3D ground point coordinates from a stereo pair of airborne images given known interior and exterior orientation parameters.

### 🧭 Purpose

The goal is to determine the optimized 3D location of a ground control point (GCP) using its known image coordinates (in both images) and corresponding camera parameters. This is a foundational step in aerial triangulation and 3D reconstruction from stereo imagery.

### 📌 Features

Implements space intersection from scratch, without third-party photogrammetry libraries.

Uses collinearity equations to solve for the object-space coordinates.

Accepts:

Image coordinates (x, y) in each stereo pair

Camera interior orientation parameters (IOP)

Exterior orientation parameters (EOP) from aerial calibration or GNSS/IMU

### 🛠️ Usage
```console
git clone https://github.com/your-username/space-intersection.git 
# Change directory
cd space-intersection 
# Run the intersection script
python space_intersection.py
```
📁 Make sure your input files (IOP, EOP, image coordinates) are correctly formatted.

### 🧮 Method

This code uses the collinearity condition equations, linearized and solved via least squares, to compute the ground coordinates of the point observed in both images.
```python
x = x0 - f * (r11*(X-Xs) + r12*(Y-Ys) + r13*(Z-Zs)) / (r31*(X-Xs) + r32*(Y-Ys) + r33*(Z-Zs))
y = y0 - f * (r21*(X-Xs) + r22*(Y-Ys) + r23*(Z-Zs)) / (r31*(X-Xs) + r32*(Y-Ys) + r33*(Z-Zs))
```
Where:

```(X, Y, Z)``` are object-space coordinates to compute

```(x, y)``` are image coordinates

```f``` is the focal length

```x0, y0``` are the principal point offsets

```R``` is the rotation matrix derived from omega, phi, kappa

```(Xs, Ys, Zs)``` are camera positions

### 🧑‍💻 Author

Developed by Abbas Salehi, as part of academic research on photogrammetry and 3D reconstruction.
