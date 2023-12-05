..
   News
   ====

 
Before 2023-11-26
    Since v1.2.3, ``UB``, instead of ``UR``, is used as default UMI tag when 
    barcodes are given.

    The ``Too many open files`` issue has been fixed (since v1.2.0). 
    The issue is commonly caused by exceeding the `RLIMIT_NOFILE`_ resource 
    limit (i.e. the max number of files allowed to be opened by system for 
    single process), which is typically 1024. 
    Specifically, in the case of ``M`` input files and ``N`` threads, 
    *cellsnp-lite* would open in total about ``M*N`` files.
    So the issue would more likely happen when large M or N is given. 
    In order to fix it, cellsnp-lite would firstly try to increase the limit to 
    the max possible value (which is typically 4096) and then use a fail-retry 
    strategy to auto detect the most suitable number of threads (which could
    be smaller than the original nthreads specified by user).

    The command line option ``--maxFLAG`` is now deprecated (since v1.0.0), 
    please use ``--inclFLAG`` and ``--exclFLAG`` instead, which are more 
    flexible for reads filtering. You could refer to the explain_flags_ page 
    to easily aggregate and convert all flag bits into one integer.
    One example is that the default exclFLAG value (without using UMIs) is 
    1796, which is calculated by adding four flag bits: UNMAP (4), 
    SECONDARY (256), QCFAIL (512) and DUP (1024).


.. _RLIMIT_NOFILE: https://man7.org/linux/man-pages/man2/getrlimit.2.html
.. _explain_flags: https://broadinstitute.github.io/picard/explain-flags.html

