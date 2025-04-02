# Genetic Algorithm

This project implements a genetic algorithm to optimize a quadratic function of the form:

$$
f(x) = a \cdot x^2 + b \cdot x + c
$$

The genetic algorithm performs the following steps:
- **Initialization:** Generates a random population of chromosomes represented as binary strings.
- **Selection:** Uses fitness proportionate selection to pick chromosomes based on their calculated fitness. The selection also implies the elitist method, always keeping the best chromosome for the next generation.
- **Crossover:** Applies a crossover probability to recombine selected chromosomes.
- **Mutation:** Mutates individual genes with a specified mutation probability.
- **Evolution:** Iteratively logs each generation's population, tracking the maximum and mean fitness.

## Project Structure

- **main.py**  
  Contains the core logic for the genetic algorithm, including population initialization, selection, crossover, mutation, and evolution. It logs detailed evolution information into `evolutie.txt`.

- **gui.py**  
  Provides a simple Tkinter-based GUI to input algorithm parameters and start the evolution process.

- **input.txt**  
  (Optional) Contains the input parameters for running the algorithm, such as population size, interval limits, function parameters, precision, crossover/mutation probabilities, and the number of generations.


- **evolutie.txt**  
  The log file produced by the genetic algorithm, which contains information such as the initial population, selection details, crossover events, mutations, and evolution of fitness values across generations.

- **README.md**  
  This file.

- **__pycache__/**  
  Directory containing Python bytecode cache files.

## Requirements

- Python 3.x (tested with Python 3.10+)
- [NumPy](https://numpy.org/) (for numerical computations)

## Installation

1. Clone the repository.
2. Install the required Python packages:
    ```sh
    pip install numpy
    ```

## Usage

### Running via GUI

Execute the following command to start the GUI:

```sh
python main.py