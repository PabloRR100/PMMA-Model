# PMMA Polymerization Mathematical Model

## Intro

Master's Degree in Chemical Engineering Final Project.  Jun 2017   
This is the work I develop on my visit to the [Technical University of Buenos Aires][itba].    
**Note**: The entire project was develop in Spanish.   

For this project I built a mathematical model for a novel technique for the polymerization of the Poly-Methil-Metacrilate.   

**Tech**: Matlab   

The entire Thesis can be downloaded [here][thesis].   
The raw publication can be found [here][publication].   


## Initiator
The contribution in our work is the use of a novel technicque for the polymerization consisting on the use of symetric cyclic multifunctional initiators:
  - DEKTP: Diethyl ketone triperoxide
  - PDP: Pinacolone diperoxide

The decompositions of both initiators are shown below:
DEKTP:
![][dektp]
PDP:
![][pdp]

## Data
The data was provided by a chemical plant in Mexico after a collaboration with ITBA.
The analysis to calculate the constante of decomposition velocity results in:
![][velocity]
And the analysis to calculate the evolution of the conversion and the molecular weight results in:
Conversion
![][conversion]
Molecular Weight
![][molecular_weight]

## Dynamics
The equations found after the balances of mass, energy and momentum are:
![][equations]

And the equations which govern the cinetics of the reaction are summarized here:
![][dynamics]

## Results
The results for the mathematical model of the evolution of the convesion during the reaction time is:
![][results1]

The results for the mathematica model of the evolution of the molecular weight during the reaction time is:
![][results2]


[itba]: https://www.itba.edu.ar/
[thesis]: https://github.com/PabloRR100/PMMA-Model/blob/master/Tesis%20Final.pdf
[publication]: https://github.com/PabloRR100/PMMA-Model/tree/master/Publication

[dektp]: /imgs/initiator_dektp.png
[pdp]: /imgs/initiator_pdp.png

[velocity]: /imgs/01_determine_velocity.png
[conversion]: /imgs/03_conversions.png
[molecular_weight]: /imgs/02_raw_data.png

[equations]: /imgs/04_equations.png
[dynamics]: /imgs/05_dynamics.png
[results1]: /imgs/06_results.png
[results2]: /imgs/07_results2.png
