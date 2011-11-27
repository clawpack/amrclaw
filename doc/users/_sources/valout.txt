
.. _valout:

===========
``valout``
===========

Signature:
    ``subroutine valout(lst, lend, time, nvar, naux)``
    

Arguments:
    ``integer, intent(in)``:
        * ``lst``: Coarsest level for detailed output.
        * ``lend``: Finest level for detailed output.
        * ``nvar``: Number of solution variables.
        * ``naux``: Number of aux variables.
    ``double precision, intent(in)``:
        * ``time``: Time of output.
        
        
Description:
    Outputs the results for a general system of conservation laws
    in 2 dimensions.  Writes the results to the file ``fort.q<iframe>``,
    using the format required by the Matlab script ``plotclaw2.m`` or 
    Python plotting tools.

    Setting ``outaux = .true.`` will also output the aux arrays to 
    ``fort.a<iframe>``.