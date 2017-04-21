<!--
# The `ALaDyn` PIC Code
-->

#### P. Londrillo¹, A. Sgattoniº, S. Sinigardi², A. Marocchinoª, F. Mira^, D. Terzani³, D. Palla¤, J. Babaei˜, T. Rovelli², G. Turchetti²  
¹INAF, Osservatorio Astronomico Bologna  
²Dipartimento di Fisica e Astronomia, Università di Bologna, INFN Sezione di Bologna, Via Irnerio 46, I-40126 Bologna (BO), Italy  
³Università Federico II di Napoli  
¤INO/CNR Pisa  
^Università La Sapienza Roma 1  
ªINFN/LNF  
ºLESIA, Observatoire de Paris, CNRS, UPMC, Université Paris Diderot, Meudon, France, LULI-UPMC Univ. Paris 06: Sorbonne Université, CNRS, École Polytechnique, CEA: Univ. Paris-Saclay, France, Istituto Nazionale di Ottica, CNR/INO, U.O.S. Adriano Gozzini, 56127 Pisa  
˜Department of Physics, Faculty of Basic Sciences, University of Mazandaran, Babolsar, Iran

## Overview
`ALaDyn` is a Particle in Cell (PIC) code designed to investigate three main physical regimes:  
- Laser-plasma interaction in under-dense gas targets for electron acceleration (LWFA).  
- Beam-plasma interaction in under-dense gas targets for electron acceleration (PWFA).  
- Laser-plasma interaction in over-dense solid targets for proton(ions) acceleration and related phenomenologies.  

Recent developments and applications are reported in the attached references.

## PIC numerical model
A PIC method is based on a hybrid formal setting, whereby plasma particles are represented
on a Lagrangian framework whereas self-consistent fields are represented on the Eulerian framework
given by Maxwell equations.  
As most other PIC codes, `ALaDyn` discretize particle and field dynamical equations 
by centred finite differences on a staggered  space and time grid (Yee's module) using one step 
second order leap-frog integrator.  
To connect Lagrangian particles to Eulerian fields collocated on the spatial grid, finite order
B-splines are used. B-splines are local polynomials with compact support, allowing to 
represent delta-like point particles on a grid for charge deposition and, by converse,
to assign field grid data to a point particle.  
Energy preserving PIC schemes do not satisfy local charge conservation and the related 
Poisson equation. By converse, using one of the many numerical recipes to enforce the
continuity equations, energy conservation is heavily damaged.  
`ALaDyn` code implements both charge or energy preserving schemes, letting the user
to choose, depending on the problem at hand.  
Besides the standard leap-frog integrator,
`ALaDyn` also implements a fourth order in space and time Runge-Kutta integrator. This scheme
requires larger computational resources, of course, but can be of help to 
improve on accuracy and reduce dispersive effects of wave propagation.  
The code implements also reduced models based on:  
- Envelope (two-scale) approximation of the laser fields and of the particle dynamics;  
- A cold fluid approximation of the wake fields.  

Reduced models are well tested only for (quasi) linear regimes. Serious problems
arise for non linear dynamics.

In all the above configurations field induced ionization (tunnelling) are also implemented
and can be activated if requested by the user.
All ionization models are based on ADK scheme plus barrier suppression (BSI) for higher
*Z* ions. For solid targets, impact (collisional) ionization is under development.

## Implementation
`ALaDyn` is almost completely written in `Fortran 90`, but a couple of utility modules are written in `C++`.
`C` can also be easily used to extend code functionalities.
`Fortran 90` is the most popular computational language in PIC codes, probably because of 
the higher efficiency in handling multidimensional arrays on a grid.
Finite difference integration allows to exploit efficient parallelism by domain decomposition 
using MPI technique to distribute the computational work among CPU units. 

`ALaDyn` has been successfully ported to many HPC architectures, both in Italy at CINECA and in Europe through
PRACE Partnerships. From the CINECA IBM-SP6 system in 2011, to the test system at CINECA based on IBM/BGP in 2012, 
then CINECA FERMI in 2014 and MARCONI in 2016, `ALaDyn` run on a multitude of HPC architectures, always extracting 
top range performances.

A new version has been recently released open source on the web, with a GPLv3 license.
In part it has been rewritten from scratch and it is the basis for future development.
Sources can be found, together with other codes, in our organization GitHub page at
[github.com/ALaDyn](https://github.com/ALaDyn)



### About ALaDyn
- C. Benedetti et al., 
ALaDyn: A High-Accuracy PIC Code for the Maxwell-Vlasov Equations, 
IEEE Transactions on Plasma Science, 36.4 (2008): 1790-1798 
[Paper](http://dx.doi.org/10.1109/TPS.2008.927143)
- P. Londrillo et al., 
Charge preserving high order PIC schemes, 
Nucl. Instr. and Meth. A, 620.1 (2010): 28-35 
[Paper](http://www.sciencedirect.com/science/article/pii/S0168900210001233)
- P. Londrillo et al., 
Numerical investigation of beam-driven PWFA in quasi-nonlinear regime, 
Nucl. Instrum. Meth. A, 740 (2014) 
[Paper](http://www.sciencedirect.com/science/article/pii/S0168900213013740)


### Papers
- C. Benedetti, P. Londrillo, T.V. Liseykina, A. Macchi, A. Sgattoni, G. Turchetti, 
Ion acceleration by petawatt class laser pulses and pellet compression in a fast ignition scenario, 
Nucl. Instr. and Meth. A, Volume 606, Issues 1-2, 11 July 2009, Pages 89-93, ISSN 0168-9002
[Paper](http://www.sciencedirect.com/science/article/pii/S0168900209005531)

- C. Benedetti, P. Londrillo, V. Petrillo, L. Serafini, A. Sgattoni, P. Tomassini, G. Turchetti, 
PIC simulations of the production of high-quality electron beams via laser-plasma interaction, 
Nucl. Instr. and Meth. A, Volume 608, Issue 1, Supplement, 1 September 2009, Pages S94-S98, ISSN 0168-9002
[Paper](http://www.sciencedirect.com/science/article/pii/S0168900209009784)

- A. Sgattoni, C. Benedetti, P. Londrillo, G. Turchetti, 
Simulation of the laser-plasma acceleration for the PLASMONX project with the PIC code ALaDyn, 
Radiation Effects and Defects in Solids Vol. 165, Iss. 6-10, 2010
[Paper](http://www.tandfonline.com/doi/abs/10.1080/10420151003732072)

- G. Turchetti, A. Sgattoni, C. Benedetti, P. Londrillo, L. Di Lucchio, 
Comparison of scaling laws with PIC simulations for proton acceleration with long wavelength pulses, 
Nucl. Instr. and Meth. A, Volume 620, Issue 1, 1 August 2010, Pages 51-55, ISSN 0168-90022
[Paper](http://www.sciencedirect.com/science/article/pii/S0168900210001270)

- A. Sgattoni et al., 
Laser ion acceleration using a solid target coupled with a low-density layer, 
Physical Review E 85.3 (2012) 036405
[Paper](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.85.036405)

- G. Turchetti, S. Sinigardi, P. Londrillo, F. Rossi, M. Sumini, D. Giove and C. De Martinis, 
The LILIA experiment: Energy selection and post-acceleration of laser generated protons, 
AIP Conf. Proc. 1507, 820 (2012)
[Paper](http://scitation.aip.org/content/aip/proceeding/aipcp/10.1063/1.4773804)

- L. A. Gizzi et al., 
Acceleration with self-injection for an all-optical radiation source at LNF, 
Nucl. Instr. and Meth. B, 309 (2013): 202-209
[Paper](http://www.sciencedirect.com/science/article/pii/S0168583X13003017)

- L. Gizzi et al., 
Laser-Plasma Acceleration and Radiation Sources for Applications, 
2013 Conference on Lasers and Electro-Optics Pacific Rim, Optical Society of America
[Paper](http://www.opticsinfobase.org/abstract.cfm?uri=CLEOPR-2013-TuD3_1)

- L. Gizzi et al., 
Laser-Plasma Acceleration with FLAME and ILIL Ultraintense Lasers, 
Appl. Sci. 2013, 3(3), 559-5800
[Paper](http://www.mdpi.com/2076-3417/3/3/559)

- A. Sgattoni et al., 
Laser plasma proton acceleration experiments using foam-covered and grating targets, 
Proc. SPIE 8779, Laser Acceleration of Electrons, Protons, and Ions II; and 
Medical Applications of Laser-Generated Beams of Particles II; and Harnessing 
Relativistic Plasma Waves III, 87790L (May 7, 2013)
[Paper](http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=1686153)

- S. Agosteo et al., 
The LILIA (Light ions laser induced acceleration) experiment at LNF, 
Lasers and Electro-Optics Europe (CLEO EUROPE/IQEC), pp.1, (2013)
[Paper](http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6801173&isnumber=6800590)

- S. Sinigardi, G. Turchetti, P. Londrillo, F. Rossi, D. Giove, C. De Martinis, M. Sumini, 
Transport and energy selection of laser generated protons for postacceleration with a compact linac, 
Phys. Rev. ST Accel. Beams 16, 031301 (2013)
[Paper](http://journals.aps.org/prstab/abstract/10.1103/PhysRevSTAB.16.031301)

- A. Macchi, A. Sgattoni, S. Sinigardi, M. Borghesi, M. Passoni, 
Advanced strategies for ion acceleration using high-power lasers, 
Plasma Phys. Control. Fusion 55 124020 (2013)
[Paper](http://iopscience.iop.org/0741-3335/55/12/124020/)

- A. Bazzani, M. Giovannozzi, P. Londrillo, S. Sinigardi, G. Turchetti, 
Case studies in space charge and plasma acceleration of charged beams, 
CERN-ACC-2014-0073 (2014), Accepted for publication in Comptes Rendus de Mècanique Acadèmie des Sciences de Paris
[Paper](http://cds.cern.ch/record/1712519/files/CERN-ACC-2014-0073.pdf)

- S. Agosteo et al., 
The LILIA (laser induced light ions acceleration) experiment at LNF, 
Nucl. Instr. and Meth. B, 331, 15 (2014)
[Paper](http://www.sciencedirect.com/science/article/pii/S0168583X14001207)

- S. Sinigardi, G. Turchetti, F. Rossi, P. Londrillo, D. Giove, C. De Martinis, P. R. Bolton, 
High quality proton beams from hybrid integrated laser-driven ion acceleration systems, 
Nucl. Instr. and Meth. A, 740, 99 (2014)
[Paper](http://www.sciencedirect.com/science/article/pii/S0168900213014873)

- A. Sgattoni, S. Sinigardi, L. Fedeli, F. Pegoraro, A. Macchi, 
Laser-Driven Rayleigh-Taylor Instability: Plasmonics Effects and Three-Dimensional Structures, 
Phys. Rev. E, 91, 013106 (2015)
[Paper](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.013106)

- A. Sgattoni, S. Sinigardi, A. Macchi, 
High Energy Gain in Three-Dimensional Simulations of Light Sail Acceleration, 
Appl. Phys. Lett. 105, 084105 (2014)
[Paper](http://aip.scitation.org/doi/10.1063/1.4894092)

- L. A. Gizzi et al.,
Role of laser contrast and foil thickness in target normal sheath acceleration,
Nucl. Instr. and Meth. A, 829, 144 (2015)
[Paper](http://www.sciencedirect.com/science/article/pii/S0168900216000528)

- S. Sinigardi, 
Numerical simulations of recent proton acceleration experiments with sub-100 TW laser systems,
Nucl. Instr. and Meth. A, 829, 167 (2015)
[Paper](http://www.sciencedirect.com/science/article/pii/S0168900216301620)

