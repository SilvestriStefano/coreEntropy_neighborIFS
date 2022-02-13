# coreEntropy_neighborIFS
The up-to-date and documented code for my research on the connection between the Mandelbrot set and the Thurston set.

One endgoal of this project is to create a RESTful API from which one can obtain all the necessary information regarding a rational angle in the context of the Mandelbrot and Thurston set.  


[THEORY SLIDES](https://slides.com/silvas/qsconjcoreentropy)

## TODOs :
- [x] write python code to calculate the Core Entropy 
- [ ] fix eigenvalue problem in the core_entropy
- [x] refactor the KneadToRat function
- [x] create python modules from the code in the jupyter notebooks
- [ ] make some example on how to use this code
- [ ] write down the theory
- [ ] figure out another way to draw graphs other than Gephi
- [ ] write a recursive function that erases vertices with out degree equal to 0.
- [ ] import Mathematica notebooks (?)

### folder's description
 
- `data`: raw data for the project. 

- `docs`: documentation. 

- `results`: results, including checkpoints, as well as figures and tables. 

- `scripts`: scripts - Python and bash alike - as well as .ipynb notebooks.

- `src`: reusable Python modules for the project. 

- `tests`: tests for the code.
