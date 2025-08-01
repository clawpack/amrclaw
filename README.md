# AMRClaw Package
The AMRClaw version of Clawpack provides Adaptive Mesh Refinement (AMR) capabilities in 2 and 3 space dimensions. (The two-dimensional code can also be used for 1-dimensional problems, see [AMRClaw for 1d problems](https://www.clawpack.org/amrclaw1d.html#amrclaw-1d).)

Block-structured AMR is implemented, in which rectangular patches of the grid at level `L` are refined to level `L+1`.  The general algorithms are described in [BergerLeVeque98](https://doi.org/10.1137/S0036142997315974).

**Links**
 - [Documentation](https://www.clawpack.org/amrclaw.html)
 - [![Tests](https://github.com/clawpack/amrclaw/actions/workflows/testing.yml/badge.svg)](https://github.com/clawpack/amrclaw/actions/workflows/testing.yml)
