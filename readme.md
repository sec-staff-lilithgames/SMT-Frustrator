# WHAT?

 This is a simple SMT-frustrator which aims to protect your function from Symbolic Execution tools like angr or klee or anythings based on SMT solvers. 

# HOW?

 It creates Opaque Predicates based on transformed Pell equation to exploit the weakness of lacking supports on non-linear constraints dealing. And also avoids pattern matching via  Quadratic transformation to add noises. 

# Usage

 1. clone this repository
 2. create a virtual enviroment: python -m venv smt
 3. run the venv: ./smt/Scripts/activate
 4. modify parameters of frustrator or just run: python frustrator.py
