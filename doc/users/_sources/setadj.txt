
.. _setadj:

===========
``setadj``
===========


Signature:
    ``subroutine setadj(mptr, nx, ny, listxlo, listxhi, listylo, listyhi)``


Arguments:
    ``integer, intent(in)``:
        * ``mptr``: Pointer to the grid being modified.
        * ``nx``: Horizontal dimension of grid ``mptr`` (interior only).
        * ``ny``: Vertical dimension of grid ``mptr`` (interior only).
    ``integer, intent(out)``:
        * ``listxlo(ny)``: Left adjacency list of grid ``mptr``.
        * ``listxhi(ny)``: Right adjacency list of grid ``mptr``.
        * ``listylo(ny)``: Bottom adjacency list of grid ``mptr``.
        * ``listyhi(ny)``: Top adjacency list of grid ``mptr``.
        

Description:
    Sets the adjacency lists for the grid ``mptr``.