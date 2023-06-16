# coreEntropy_neighborIFS
The up-to-date and documented code for my research on the connection between the Mandelbrot set and the Thurston set.

One endgoal of this project is to create a RESTful API from which one can obtain all the necessary information regarding a rational angle in the context of the Mandelbrot and Thurston set.  


[THEORY SLIDES](https://slides.com/silvas/qsconjcoreentropy)

## TODOs :
- [x] write python code to calculate the Core Entropy 
- [x] fix eigenvalue problem in the core_entropy
- [x] refactor the KneadToRat function
- [x] create python modules from the code in the jupyter notebooks
- [ ] make some example on how to use this code
- [ ] write down the theory
- [ ] figure out another way to draw graphs other than Gephi
- [x] write a recursive function that erases vertices with out degree equal to 0.
- [ ] import Mathematica notebooks (?)
- [ ] write function that obtains the landing parameter in the Mandelbrot set from angle (check Wolf Jung's work)

## Folder's description
 
- `data`: raw data for the project. 

- `docs`: documentation. 

- `results`: results, including checkpoints, as well as figures and tables. 

- `scripts`: scripts - Python and bash alike - as well as .ipynb notebooks.

- `src`: reusable Python modules for the project. 

- `tests`: tests for the code.

## Install requirements

### with pip
Create a virtual environment called `mandel_thurston`  (you can call it whatever you want)
```shell
$ python -m venv mandel_thurston
``` 

activate it
```shell
$ source mandel_thurston/Scripts/activate
``` 

or if you are using PowerShell
```powershell
> mandel_thurston/Scripts/Activate.ps1
``` 

install the required packages
```shell
$(mandel_thurston) pip install -r requirements.txt
$(mandel_thurston) pip install -e . # install this package
```

### with conda
create a virtual environment
```shell
$ conda env create -f mandel_thurston.yml
```

activate it
```shell
$ conda activate mandel_thurston
```

install this package
```shell
$(mandel_thurston) pip install -e .
```


## Run tests
```shell
$(mandel_thurston) pytest -vv --cov=src
``` 
- the flag `-vv` is for verbose
- the flag `--cov` is for the coverage report