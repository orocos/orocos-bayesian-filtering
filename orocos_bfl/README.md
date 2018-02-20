# BFL: BAYESIAN FILTERING LIBRARY
The Bayesian Filtering Library (BFL) provides an application
independent framework for inference in Dynamic Bayesian Networks,
i.e., recursive information processing and estimation algorithms based
on Bayes' rule, such as (Extended) Kalman Filters, Particle Filters 
(or Sequential Monte Carlo methods), etc. The Bayesian Filtering
Library (BFL) provides an application independent framework for
inference in Dynamic Bayesian Networks, i.e., recursive information
processing and estimation algorithms based on Bayes' rule, such as
(Extended) Kalman Filters, Particle Filters (or Sequential Monte Carlo
methods), etc. These algorithms can, for example, be run on top of the
Realtime Services, or be used for estimation in Kinematics & Dynamics
applications.

This library encoporates ideas from several available software
libraries:
* Scene (Andrew Davison).  See <http://www.robots.ox.ac.uk/~ajd/Scene/>
* Bayes++ (from ACFR). See <http://www.acfr.usyd.edu.au/> 
* The CES programming library (Sebastian Thrun).  See
  <http://www-2.cs.cmu.edu/afs/cs.cmu.edu/user/thrun/public_html/papers/thrun.ces-tr.html> 

It's most important features are:
* Released under the GNU LGPL licence
* Wrapper around common matrix and RNG libraries, so you can use your own
  favourite matrix/rng library.  Wrappers exist for
  * boost <http://www.boost.org/> Matrix and RNG (default)
  * eigen <http://eigen.tuxfamily.org/> Matrix Library
  * LTIlib <http://ltilib.sourceforge.net/doc/homepage/index.shtml> matrix/RNG library of: a library with algorithms and data structures frequently used in image processing and computer vision.
  * NEWMAT <http://www.robertnz.net/nm_intro.htm> Matrix Library
  * Scythe <http://scythe.berkeley.edu> RNG library
* "Bayesian unifying Design".  This allows to incorporate any Bayesian
  filtering algorithm! Currently the following filter schemes are
  implemented:
  * Standard KF, EKF, IEKF and Non-minimal State KF (See
  <http://people.mech.kuleuven.ac.be/~tlefebvr/publicaties/BayesStat.ps.gz> 
  * Standard Particle filter (arbitrary proposal), BootstrapFilter
  (Proposal = System Model PDF), Auxiliary Particle filter, Extended
  Kalman Particle Filter. 

Tinne De Laet Contributed a tutorial which still needs to be converted
into markdown syntax (WIP)... 
It discusses how to construct your first filter in bfl.

Wim Meeussen and Tinne De Laet contributed a installation guide which
also needs to be converted into markdown...



      










