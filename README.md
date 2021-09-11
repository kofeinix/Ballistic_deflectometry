[![LinkedIn][linkedin-shield]][linkedin-url]


<!-- ABOUT THE PROJECT -->
## About The Project

In order to study generated magnetic fields qualitatively and quantitatively various techniques can be
employed. One of them is proton radiography (deflectometry). It allows to measure both electric and magnetic fields and even track their evolution (on a ps scale).
Principle scheme of this diagnostic is shown below:
![scheme](https://user-images.githubusercontent.com/90211042/132966051-479b1cd8-1b48-48dc-8ff9-eb96ad0842e1.png)

The laser irradiates thin foil and thus creates a directed beam of protons. This beam is well-collimated and the proton energy distribution is wide.
The beam passes through the target (field) region and the protons are deflected in some way.
They are recorded by the radiochromic films (RCF) stack. RCF changes color when exposed to ionizing radiation. 
The stack is composed in a way that allows to register protons of different energies, known for each layer.
As proton energy has a direct connection with it's speed and thus with time-of-flight, it means that RCF stack layers contain both temporal and spatial evolution of fields.
The mesh can be placed in a way of a proton beam to produce some structure on the RCF and study it's warping.

To process the experimental data, a simulations must be implemented. By comparing the experimental and synthetic image, the electric and magnetic field parameters may be defined.
This project is the ballistic simulation code with variable parameters. It was used for writing the following article - see [Ref.1](https://iopscience.iop.org/article/10.1088/1742-6596/1686/1/012004).

## Contents of the Project
The project contains synthetic_image_generator.py, which is a scipt for generation of images, My_Colormaps.py, which is a custom colormap for plotting, folders with input (fields and coordinates, explained later) and output data.

## Input data
[![Comsol][Comsol-shield]][Comsol-url]
[![Wolfram Mathematica][Wolfram-shield]][Wolfram-url]

The project is build to work with a strictly defined format of input data. The formats are taken from the software, that was used to produce the field files.
* For the electric field, [Comsol](https://www.comsol.ru/) software (v.3.5a) was used. The output is three files Ex, Ey and Ez, each having x, y, z, Ei columns.
* The magnetic field was calculated using [Wolfram Mathematica](https://www.wolfram.com/mathematica/), add-on [Radia](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia) and a custom script.
The output is a singe file with Bx, By, Bz columns. Also, a separate file with coordinates x, y, z should be used.

## Libraries
The most important library beside standard numpy is a library [numba](http://numba.pydata.org/). It is used to significantly accelerate calculation by using @jit decorator and prange instead of range. This will be discussed further.
## Program initialization
At first, all the input and output file paths should be checked. (2 files for B field and 3 - for E field). The the code checks for infs and NaNs, which may be present and tries to get rid of them.
In creates a single big array for B field and a singe big array for E field. In the code everywhere _meters_ are used.
Important step - to write minimum and maximum coordinates that are used for field creation. Also, the step is important. In the presented example, the step is 10 microns, and this leads to 120 x 100 x 150 point grid.
However, even such not detailed grid still a heavy-weight file.
One of the important parameters is "U", which is a potential of the target, used for electric field calculation in COMSOL. The electric field depends linearly on this potential, thus the same multiplication factor is used in the code.

In the code, the particle parameters are present. They are: mass "m", charge "q", energy in MeV "T".  The time step between each calculation (of speed and acceletion) "dt" can be varied, but a smaller value will lead to increased calculatio time.

## Calculation principle
The proton trajectory is calculated using an equation on motion calculation with defined time-step. The particles fly individually from some source (coordinates defined in the beginning of "calculation" function).
The mesh is set at some distance from the target (defined in "ttomesh" variable). 1500/inch is used as an exapmle in the code, and all particles, that touch it's wires, are deleted.
The RCF (detector) is placed at some distance from (0,0,0) and is defined at "ttorcf". 

When flying through the field area, field in corresponding coordinates is extracted from field arrays. _There is no searching for the coordinates in the array, the necessary row is defined using known x, y and z limits, their sorting and step value_.
This helps to accelerate the code, as any search will increase calculation time drastically.
The trilinear interpolation of field value in the point between grid lines is implemented. The results look cleaner with interpolation turned on, however it takes much more time than taking a closest value. 
The @njit(parallel=True) uses the numba library and accelerates code significantly. In the demo code itself ("calculation" function) the "prange" is used instead of range. 
__The main difference is that prange uses all the CPU cores at 100%! Thus the PC is hard to use while prange is active! Simple range is slower, but allows to use PC as usual__
Without @njit only range can be used, however the calculation takes incomparably more time (x100 or even more). Thus __using [numba](http://numba.pydata.org/) is necessary__.
Note that using prange will mess the progress check print statement up, but this print is useful for a simle range.

Bznach contains desired field values (in the demo case - Bz) in the point (0, 0, 0). Eznach contains desired potential (U) values to be calculated.
The results are saved as a 5000x5000 image file. 

In the end, the any field component can be plotted in either XY or XZ plane at any cross-section or in many cross-sections - PlotField function is used for that. An example of usage for many cross-sections is commented in the end, the example images are presented in the output folder. 


## Notes
_The work is still in progress and some improvements and corrections may be soon done._

## Contact

Iurii Kochetkov -  iu.kochetkov@gmail.com

Project Link: [https://github.com/kofeinix/Ballistic_deflectometry](https://github.com/kofeinix/Ballistic_deflectometry)


<!-- MARKDOWN LINKS & IMAGES -->

[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/iu-kochetkov/
[Comsol-shield]: https://cdn.comsol.com/company/logo/comsol-logo-130x20.png
[Comsol-url]: https://www.comsol.ru/
[Wolfram-shield]: https://www.wolfram.com/common/images/gl-logo-spikey.ru.png
[Wolfram-url]: https://www.wolfram.com/mathematica/
