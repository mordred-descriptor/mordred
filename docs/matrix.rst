Matrix aggregating methods
==========================

SpAbs
-----
.. math::
    {\rm SpAbs} = \sum_{i = 1}^N \left| \lambda_i \right|

SpMax
-----
.. math::
    {\rm SpMax} = \max_{i = 1}^N \lambda_i

SpDiam
------
.. math::
    {\rm SpDiam} = {\rm SpMax} - {\rm SpMin}

SpAD
----
.. math::
    {\rm SpAD} = \sum_{i = 1}^N \left| \lambda_i - \bar{\lambda} \right|

SpMAD
-----
.. math::
    {\rm SpMAD} = \frac{\rm SpAD}{A}

EE
--
.. math::
    {\rm EE} = \log(\sum_{i = 1}^N \exp(\lambda_i))

SM1
---
.. math::
    {\rm SM1} = \sum_{i = 1}^N \lambda_i

VE1
---
.. math::
    {\rm VE1} = \sum_{i = 1}^N \left| \ell_i \right|

VE2
---
.. math::
    {\rm VE2} = \frac{\rm VE1}{A}

VE3
---
.. math::
    {\rm VE3} = \frac{A}{10} \log({\rm VE1})

VR1
---
.. math::
    {\rm VR1} = \sum_{(i, j) \in {\rm binds}} \left( \ell_i \cdot \ell_j \right)^{-1/2}

VR2
---
.. math::
    {\rm VR2} = \frac{\rm VR1}{A}

VR3
---
.. math::
    {\rm VR3} = \frac{A}{10} \log({\rm VR1})
