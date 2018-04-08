.. title:: Publications

.. meta::
   :description: MIRTK Publications. Academic articles/papers about MIRTK tools.
   :keywords: MIRTK Publication, MIRTK Validation, MIRTK FFD, MIRTK Applications.

.. role:: red
.. role:: blue
.. role:: underline

.. |br| raw:: html

   <br />

============
Publications
============


Registration
------------

**Classic free-form deformation**

The following papers introduce the classic B-spline free-form deformation (FFD) algorithm.

.. [Schnabel01] `A generic framework for non-rigid registration based on non-uniform multi-level free-form deformations <http://link.springer.com/chapter/10.1007%2F3-540-45468-3_69>`__. |br|
              J. A. Schnabel, D. Rueckert, M. Quist, J. M. Blackall, A. D. Castellano Smith, T. Hartkens, G. P. Penney,
              W. A. Hall, H. Liu, C. L. Truwit, F.A. Gerritsen, D. L. G. Hill, and D.J. Hawkes.
              In 4th Int. Conf. on Medical Image Computing and Computer-Assisted Intervention (MICCAI '01),
              LNCS 2208:573-581, Utrecht, NL, October 2001.

.. [Rueckert99] `Non-rigid registration using free-form deformations: Application to breast MR images <http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=796284>`__. |br|
           D. Rueckert, L. I. Sonoda, C. Hayes, D. L. G. Hill, M. O. Leach, and D. J. Hawkes.
           IEEE Transactions on Medical Imaging, 18(8):712-721, 1999.

.. [Denton99] `Comparison and evaluation of rigid and non-rigid registration of breast MR images <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.46.5218&rep=rep1&type=pdf>`__. |br|
            E. R. E. Denton, L. I. Sonoda, D. Rueckert, S. C. Rankin, C. Hayes, M. Leach, D. L. G. Hill, and D. J. Hawkes.
            Journal of Computer Assisted Tomography, 23:800-805, 1999.


**Sparse free-from deformation**

The use of the sparsity constraint for the spatio-temporal registration of cardiac image sequences
was first published in the following conference article.

.. [Shi12] `Registration Using Sparse Free-Form Deformations <http://link.springer.com/chapter/10.1007%2F978-3-642-33418-4_81>`__. |br|
              W. Shi, X. Zhuang, L. Pizarro, W. Bai, H. Wang, K.-P. Tung, P. Edwards, and D. Rueckert.
              In 15th Int. Conf. on Medical Image Computing and Computer-Assisted Intervention (MICCAI '12),
              LNCS 7511:659-666, 2012.


**Statistical free-from deformation**

The following publications learn a statistical free-form deformation model from a training dataset
to restrict the deformation on new images to the learned plausible deformations.

.. [Rueckert03] `Automatic Construction of 3D Statistical Deformation Models Using Non-rigid Registration <http://link.springer.com/chapter/10.1007%2F3-540-45468-3_10>`__. |br|
           D. Rueckert, A. F. Frangi, J. A. Schnabel.
           In 4th Int. Conf. on Medical Image Computing and Computer-Assisted Intervention (MICCAI '01),
           LNCS 2208:77-84, Utrecht, NL, October 2001.

.. [Pszczolkowski12] `Nonrigid free-form registration using landmark-based statistical deformation models. <http://pubs.doc.ic.ac.uk/SDM-nonrigid-registration/SDM-nonrigid-registration.pdf>`__. |br|
            S. Pszczolkowski, L. Pizarro, R. Guerrero, D. Rueckert.
            Proc. SPIE 8314, Medical Imaging 2012: Image Processing, 831418, February 2012.


**Symmetric/Diffeomorphic free-form deformation**

The details and application of the symmetric diffeomorphic registration which parameterizes the
non-rigid transformation by a stationary velocity field with B-spline interpolation can be found
in the following workshop paper.

.. [Schuh14] `Construction of a 4D Brain Atlas and Growth Model using Diffeomorphic Registration <http://andreasschuh.com/wp-content/uploads/2015/09/miccai2014-stia.pdf>`__. |br|
            A. Schuh, M. Murgasova, A. Makropoulos, C. Ledig, S. J. Counsell, J. V. Hajnal, P. Aljabar D. Rueckert.
            In 3rd MICCAI Workshop on Spatiotemporal Image Analysis for Longitudinal and Time-Series Image Data (STIA),
            LNCS 8682, pp. 27-37, Boston, MA, 2014.

.. [Schuh18] `Unbiased construction of a temporally consistent morphological atlas of neonatal brain development <https://www.biorxiv.org/content/early/2018/01/28/251512>`__. |br|
            A. Schuh, A. Makropoulos, E. C. Robinson, L. Cordero-Grande, E. Hughes, J. Hutter, A. N. Price, M. Murgasova, R. Pedro A. G. Teixeira, N. Tusor, J. K. Steinweg, S. Victor, M. A. Rutherford, J. V. Hajnal, A. D. Edwards, D. Rueckert.
            bioRxiv, 2018.


Draw-EM
-------

When you use the Developing brain Region Annotation With Expectation-Maximization (Draw-EM) package, please cite the following article:

.. [Makropoulos14] `Automatic whole brain MRI segmentation of the developing neonatal brain <http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6810848>`__. |br|
           A. Makropoulos, I. S. Gousias, C. Ledig, P. Aljabar, A. Serag, J. V. Hajnal, A. D. Edwards, S. J. Counsell, D. Rueckert.
           IEEE Transactions on Medical Imaging, vol. 33 (9), pp. 1818-1831, 2014.


Deformable
----------

**Neonatal cortex reconstruction**

The reconstruction of the neonatal cortex from T2-weighted (and T1-weighted) brain MR images is detailed in:

.. [Schuh17] `A deformable model for the reconstruction of the neonatal cortex <https://doi.org/10.1109/ISBI.2017.7950639>`__. |br|
            A. Schuh, A. Makropoulos, R. Wright, E. C. Robinson, N. Tusor, J. Steinweg, E. Hughes, L. Cordero-Grande, A. N. Price, J. Hutter, J. V. Hajnal, D. Rueckert.
            In 14th IEEE International Symposium on Biomedical Imaging (ISBI), April 2017.

.. [Makropoulos18] `The developing human connectome project: A minimal processing pipeline for neonatal cortical surface reconstruction <https://doi.org/10.1101/125526>`__. |br|
            A. Makropoulos, E. C Robinson, A. Schuh, R. Wright, S. Fitzgibbon, J. Bozek, S. J. Counsell, J. Steinweg, J. Passerat-Palmbach, G. Lenz, F. Mortari, T. Tenev, E. P. Duff, M. Bastiani, L. Cordero-Grande, E. Hughes, N. Tusor, J.-D. Tournier, J. Hutter, A. N. Price, M. Murgasova, C. Kelly, M. A. Rutherford, S. M. Smith, A. D. Edwards, J. V. Hajnal, M. Jenkinson, D. Rueckert.
            NeuroImage, vol. 173, pp. 88-112, 2018.


Point Set
---------

**Spectral surface matching**

The spectral surface matching implementations are based on the MATLAB code provided by
`Herve Lombaert <http://step.polymtl.ca/~rv101/projects.php>`__.

.. [Lombaert13] `Diffeomorphic Spectral Matching of Cortical Surfaces. <http://link.springer.com/chapter/10.1007%2F978-3-642-38868-2_32#page-1>`__
                H. Lombaert, J. Sporring, and K. Siddiqi." << endl;
                In 23rd Int. Conf. on Image Processing in Medical Imaging (IPMI '13)
                LNCS 7917:376-389, Asilomar, CA, USA, June 28-July 3, 2013.
