## ESEPP: an event generator for the elastic scattering of charged leptons on protons ##

ESEPP is a multipurpose event generator developed for the Monte Carlo simulation of unpolarized elastic scattering of charged leptons (electrons, positrons, muons, and antimuons) on protons. The generator takes into account the lowest-order QED radiative corrections to the Rosenbluth cross section including first-order bremsstrahlung without using the soft-photon or ultrarelativistic approximations. A detailed theoretical description of the generator is given in [arXiv:1401.2959](https://arxiv.org/abs/1401.2959) or [J. Phys. G: Nucl. Part. Phys. 41, 115001](https://doi.org/10.1088/0954-3899/41/11/115001) (users of ESEPP are kindly requested to cite this paper).

The event generator is written in the C++ programming language using some [ROOT](http://root.cern.ch) classes. In particular, the [TLorentzVector](http://root.cern.ch/root/html/TLorentzVector.html) class allows us to simplify working with the four-momenta of particles, the [TRandom3](http://root.cern.ch/root/html/TRandom3.html) class is used for pseudorandom number generation, and the [TFoam](http://root.cern.ch/root/html/TFoam.html) class is the basis for event generation and multidimensional phase space integration. The source code of ESEPP is freely available under the GNU GPL license and can be used as a basis for future developments.

To compile the event generator, one should run the **make** command to build an executable file. A computer with ROOT installed is required (including the [MathMore](http://root.cern.ch/drupal/content/mathmore-library) library, the availability of which can be checked with the command **root-config --has-mathmore**). We use ESEPP under the GNU/Linux operating system, but it should also work on MS Windows and OS X platforms.

The source distribution of ESEPP includes the following files:

* DOCUMENTATION
  * [README.md](README.md) (this markdown file)
  * [input_output.pdf](input_output.pdf) (description of the input parameters and the output files format)
  * [calc_lepton.pdf](calc_lepton.pdf) (Mathematica/FeynCalc calculation of the lepton bremsstrahlung term)
  * [calc_proton.pdf](calc_proton.pdf) (Mathematica/FeynCalc calculation of the proton bremsstrahlung term)
  * [calc_interference.pdf](calc_interference.pdf) (Mathematica/FeynCalc calculation of the interference bremsstrahlung term)
  * [LICENSE](LICENSE) (the GNU General Public License)
* SOURCE CODE
  * [Makefile](Makefile) (makefile to build an executable program)
  * [Makefile.arch](Makefile.arch) (architecture-dependent part of the makefile)
  * [esepp.cxx](esepp.cxx) (main source file)
  * [esepp.h](esepp.h) (main header file)
  * [const.h](const.h) (header file containing mathematical and physical constants)
  * [dialog.h](dialog.h) (header file for an initial dialogue)
  * [lepton.h](lepton.h) (header file containing the lepton bremsstrahlung term)
  * [proton.h](proton.h) (header file containing the proton bremsstrahlung term)
  * [interference.h](interference.h) (header file containing the interference bremsstrahlung term)
  * [vpol.dat](vpol.dat) (data file for the vacuum polarization correction, extracted from [here](http://cmd.inp.nsk.su/~ignatov/vpl/vpol_all_bare_sum_v1.dat))
