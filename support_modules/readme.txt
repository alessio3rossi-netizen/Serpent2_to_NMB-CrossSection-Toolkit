## Support Modules

The `support_modules/` directory contains auxiliary data files and MATLAB functions required during microscopic cross-section extraction and formatting.

### Nuclear Data Files

#### `endfb81.mat`

MATLAB binary file containing the radioactive decay matrix derived from the ENDF/B-VIII.1 nuclear data library and exported from Serpent2.

#### `jeff4.mat`

MATLAB binary file containing the radioactive decay matrix derived from the JEFF-4.0 nuclear data library and exported from Serpent2.

#### `jendl5.mat`

MATLAB binary file containing the radioactive decay matrix derived from the JENDL-5.0 nuclear data library and exported from Serpent2.

These decay matrices are used during the reconstruction and processing of transmutation pathways from Serpent2 depletion matrices.

---

### MATLAB Functions

#### `NMB_Isotopes.m`

Defines the isotopic ordering adopted by the Nuclear Material Balance (NMB) framework.

This ordering is used to:

* Map Serpent2 isotopes to the NMB isotope set.
* Ensure consistent formatting of generated microscopic cross-section libraries.
* Maintain compatibility with NMB input requirements.

#### `nu_ORIGEN.m`

Provides average neutron yields per fission ($\bar{\nu}$) derived from ORIGEN data.

The neutron yields are used during microscopic cross-section processing to reconstruct reaction rates and determine neutron production from fission reactions.

The data are inherited directly from ORIGEN and are applied consistently across all supported reactor and fuel cases.
