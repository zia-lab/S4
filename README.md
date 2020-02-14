# S4

S4 is a wonderful tool coded at Stanford some years ago to do RCWA in layered structures. It originally provided both `lua` and `python` interfaces but here I've trimmed down to just the Python extension. It can be finnicki to build, and even more so if the build is made without admin privileges. Here I've forked phoebe's go at it, and added details to build it on Brown's University computer cluster.

The originally provided examples were all in lua, but here I'll work on making their pythonic counterparts.

## Installation Instructions

(Works on the cluster at Brown which as of Feb 2020 is running on Red Hat Enterprise 7.3. Install uses both modules provided by the cluster and packages installed through **conda**. Resulting install has some caveats as explained below.)

```
module load anaconda/3-5.2.0 suiteparse/4.5.4 blas/3.7.0 lapack/3.7.0 fftw/3.3.6
conda create --name python35 python=3.5
source activate python35
conda install -c conda-forge openblas lapack fftw suitesparse
pip install numpy

git clone https://github.com/zia-lab/S4
cd S4
make
make boost
make S4_pyext
```

## If all goes well try this

```
module load anaconda/3-5.2.0
python
```
```
import S4
import numpy as np
import matplotlib.pyplot as plt

S = S4.New(Lattice = ((.24,0),(0,.24)), NumBasis=100)

S.SetMaterial(Name='silicon', Epsilon =12+0.01j)
S.SetMaterial(Name='vacuum', Epsilon=1+0.00j)
S.SetMaterial(Name='b', Epsilon = 12+0.00j)

S.AddLayer(Name = 'incidentMedium', Thickness=0, Material='vacuum');
S.AddLayer(Name = 'pillar', Thickness=0.120, Material='vacuum')

S.SetRegionRectangle('pillar', 'silicon', (0,0), 0, (0.021, 0.021))

S.AddLayer(Name = 'base', Thickness = 10e-3, Material = 'silicon')
S.AddLayer(Name = 'airBellow', Thickness = 0, Material = 'vacuum')

S.SetExcitationPlanewave(IncidenceAngles = (0, 0),
                        pAmplitude = 1 + 0j,
                        sAmplitude = 0 + 0j)

freq = np.linspace(2.5, 1.11111, 10)

reflections = []
for freqI in freq:
    S.SetFrequency(freqI)
    [f_front,b_front] = S.GetPowerFlux('incidentMedium', 0)
    [f_back, b_back] = S.GetPowerFlux('airBellow',0)
    reflection = -b_front.real/f_front.real
    transmission = f_back.real/f_front.real
    reflections.append(reflection)

plt.figure()
plt.plot(freq, reflections)
plt.show()


```



Here's the original description of S4:

-------------------------------------

S4: Stanford Stratified Structure Solver (http://fan.group.stanford.edu/S4/)

A program for computing electromagnetic fields in periodic, layered
structures, developed by Victor Liu (victorliu@alumni.stanford.edu) of the
Fan group in the Stanford Electrical Engineering Department.

See the S4 manual, in doc/index.html, for a complete
description of the package and its user interface, as well as
installation instructions, the license and copyright, contact
addresses, and other important information.

---------------------------------------
