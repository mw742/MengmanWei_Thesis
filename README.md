# Mengman Wei’s Ph.D. Thesis
**University of Cambridge**  
**Advisor:** Prof. Jonathan Goodman

## Overview
This repository contains the Ph.D. thesis of Mengman Wei, focused on the development of conformation searching algorithms. The thesis introduces *Diamond Energy*, a systematic conformation searching method designed to enhance our understanding of molecular structures by efficiently identifying low-energy conformations within a diamond lattice.

## Diamond Energy
*Diamond Energy* is a Python program based on the assumption that the structures of local minima in molecular conformations resemble a diamond lattice. This innovative approach allows for systematic searching of conformations using integer arithmetic and simplified energy equations to enhance speed and accuracy.

### Key Features
- **Systematic Conformational Searching**: Capable of exhaustively and accurately performing conformational searches for acyclic alkanes.
- **High Performance**: Demonstrates significant speed in conformational searching, leveraging optimized computations.
- **Extensibility**: Initially developed for acyclic alkanes, the program has been extended to include six-membered rings and saccharides.

### Extensions
The scope of *Diamond Energy* has been broadened through an automatic pipeline that searches for suitable energy evaluation parameters for new atom types. This has facilitated its application to more complex molecular structures.

## Machine Learning Integration
Utilizing the vast amount of data generated by *Diamond Energy*, a machine learning pipeline was developed to analyze the geometry distribution of the diamond lattice framework. This integration aids in refining and improving the conformational searching process.

## Reference
```bibtex
@phdthesis{Wei2024DiamondEnergy,
  title={*Diamond Energy* – a systematic conformation searching method},
  author={Wei, Mengman},
  year={2024},
  school={University of Cambridge}
  doi={https://doi.org/10.17863/CAM.109200}
}
