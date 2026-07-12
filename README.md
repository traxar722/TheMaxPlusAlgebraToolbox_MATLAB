# TheMaxPlusAlgebraToolbox_MATLAB
This repository presents a toolbox for Max-Plus Algebra developed in the MATLAB environment, including 3 classes that overwrite mathematic operators and MATLAB functions, represent Petri Networks and automatically generate linear Max-Plus systems.

 **The maxplus class**\
 Implements the max-plus algebra as an algebraic structure. It overloads standard arithmetic operators (where addition is mapped to the $\max$ operator and multiplication is mapped to classical addition) and extends these properties to matrix calculus. It also includes advanced algebraic tools essential for system analysis, such as the Kleene Star operator, the maximal eigenvalue. It inludes the possibility to analyze graph propreties such as connectivity and strong connectivity, ciclity/aciclity and also solving linear equations in the Max-Plus sense such as Ax = b and x = Ax + b. 

  **The PetriNetwork class**

  The objective of this class is to represent *timed* and *untimed* Petri Nets using the *initial marking*, *input* and *output* incidence matrixes and the *temporization vector* for timed Nets. The class can analyze the marked graph topology of Petri Nets, classify transitions and positions including the maximal number of tokens for each type of position.

**The MaxPlusSystem class**

  Represents the discrete event input-state-output system from the position timed Petri Nets automatically. Automatically maps it into linear state equations within the max-plus domain starting from P timed Petri Nets. The class can simulate the evolution of the system for a desired number of iterations and inputs, automatically detects deadlock and generates Gantt Diagrams. For the generation of Gantt diagrams, the indices of each operation for each resource are needed , indices wich will be put into a matrix having NxO (where N is the number of resources that we want to represent and O is the number of operations).


 This repository includes a test file including some examples of real life systems and examples of using the Max-Plus Algebra Toolbox.
