===
FAQ
===

Q1. No reads captured from bam file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A: If there is no read captured from the bam file, there can be multiple 
reasons:

* your reads don't have UMI tag (please ``--UMItag None``) or the UMI tag is not
  ``UR`` (please specify)

