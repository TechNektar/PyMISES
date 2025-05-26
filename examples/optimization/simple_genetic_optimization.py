"""
PyMISES - Simple Genetic Optimization Example

This example demonstrates how to use a simple genetic algorithm to optimize an airfoil shape.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import logging
import time
from scipy.interpolate import interp1d

# Add the PyMISES directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pymises.core.geometry import AirfoilGeometry

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SimpleGeneticOptimizer:
    """
    Simple genetic algorithm for airfoil optimization.
    """

    def __init__(self, base_airfoil='naca0012', n_design_vars=6, population_size=20):
        """
        Initialize the optimizer.

        Parameters
        ----------
        base_airfoil : str
            Name of the baseline airfoil
        n_design_vars : int
            Number of design variables
        population_size : int
            Size of the population
        """
        self.base_airfoil_name = base_airfoil
        self.n_design_vars = n_design_vars
        self.population_size = population_size

        # Load baseline airfoil
        if base_airfoil.lower().startswith('naca'):
            self.base_airfoil = AirfoilGeometry.create_naca(base_airfoil, n_points=101)
        else:
            try:
                self.base_airfoil = AirfoilGeometry.import_from_file(f'{base_airfoil}.dat')
            except:
                print(f"Could not load {base_airfoil}.dat, using NACA 0012 instead")
                self.base_airfoil = AirfoilGeometry.create_naca('0012', n_points=101)

        # Create Hicks-Henne bump functions
        self._create_bump_functions()

        # Initialize population
        self.population = np.random.uniform(-0.05, 0.05, (population_size, 2 * n_design_vars))
        self.fitness = np.zeros(population_size)

        # Best design found
        self.best_design = None
        self.best_fitness = -float('inf')
        self.best_airfoil = None

        # Optimization history
        self.history = {
            'generation': [],
            'best_fitness': [],
            'avg_fitness': [],
            'best_design': []
        }

    def _create_bump_functions(self):
        """
        Create Hicks-Henne bump functions for airfoil parameterization.
        """
        self.bump_locations = np.linspace(0.1, 0.9, self.n_design_vars)
        self.bump_width = 0.2  # Width parameter for bumps

    def _apply_bumps(self, design_vars):
        """
        Apply Hicks-Henne bump functions to modify the airfoil.

        Parameters
        ----------
        design_vars : np.ndarray
            Design variables controlling bump amplitudes

        Returns
        -------
        AirfoilGeometry
            Modified airfoil geometry
        """
        # Split design variables for upper and lower surfaces
        upper_vars = design_vars[:self.n_design_vars]
        lower_vars = design_vars[self.n_design_vars:]

        # Get baseline coordinates
        x = self.base_airfoil.x
        y = self.base_airfoil.y

        # Identify upper and lower surfaces
        upper_mask = y >= 0
        lower_mask = y < 0

        # Create new y-coordinates
        new_y = np.copy(y)

        # Apply bumps to upper surface
        for i, loc in enumerate(self.bump_locations):
            # Hicks-Henne bump function: y = sin(pi * x^t)^2
            # where t controls the location of the bump maximum
            t = np.log(0.5) / np.log(loc)
            bump = np.sin(np.pi * x[upper_mask]**t)**2
            new_y[upper_mask] += upper_vars[i] * bump

        # Apply bumps to lower surface
        for i, loc in enumerate(self.bump_locations):
            t = np.log(0.5) / np.log(loc)
            bump = np.sin(np.pi * x[lower_mask]**t)**2
            new_y[lower_mask] -= lower_vars[i] * bump  # Subtract for lower surface

        # Create new airfoil geometry
        new_coords = np.column_stack((x, new_y))
        new_airfoil = AirfoilGeometry(coords=new_coords, name=f"{self.base_airfoil_name}_optimized")

        return new_airfoil

    def _evaluate_fitness(self, airfoil):
        """
        Evaluate the fitness of an airfoil.

        This is a simplified fitness function that rewards:
        1. Increased thickness for structural strength
        2. Smooth curvature for good aerodynamics
        3. Sharp trailing edge for good aerodynamics

        Parameters
        ----------
        airfoil : AirfoilGeometry
            Airfoil to evaluate

        Returns
        -------
        float
            Fitness value (higher is better)
        """
        # Get coordinates
        x = airfoil.x
        y = airfoil.y

        # Identify upper and lower surfaces
        upper_mask = y >= 0
        lower_mask = y < 0

        # Calculate thickness distribution
        x_upper = x[upper_mask]
        y_upper = y[upper_mask]
        x_lower = x[lower_mask]
        y_lower = y[lower_mask]

        # Sort by x-coordinate
        upper_sort = np.argsort(x_upper)
        lower_sort = np.argsort(x_lower)

        x_upper = x_upper[upper_sort]
        y_upper = y_upper[upper_sort]
        x_lower = x_lower[lower_sort]
        y_lower = y_lower[lower_sort]

        # Interpolate lower surface to match upper surface x-coordinates
        if len(x_lower) > 3:
            f_lower = interp1d(x_lower, y_lower, kind='linear', bounds_error=False, fill_value='extrapolate')
            y_lower_interp = f_lower(x_upper)
            thickness = y_upper - y_lower_interp
        else:
            # Fallback if not enough points
            thickness = np.ones_like(x_upper) * 0.12

        # Calculate maximum thickness
        max_thickness = np.max(thickness)

        # Calculate thickness at 30% chord (structural consideration)
        idx_30 = np.argmin(np.abs(x_upper - 0.3))
        thickness_30 = thickness[idx_30]

        # Calculate curvature (simplified)
        dy_upper = np.diff(y_upper)
        dx_upper = np.diff(x_upper)
        # Avoid division by zero
        dx_upper = np.where(np.abs(dx_upper) < 1e-10, 1e-10, dx_upper)
        slopes = dy_upper / dx_upper
        # Avoid NaN values
        slopes = np.nan_to_num(slopes)
        curvature_upper = np.abs(np.diff(slopes))

        # Penalize high curvature
        curvature_penalty = np.sum(curvature_upper**2)

        # Trailing edge thickness (should be small)
        te_thickness = np.abs(y[-1] - y[0])

        # Leading edge radius (approximation)
        le_radius = thickness[0] / 2

        # Compute fitness
        # Handle potential NaN or inf values
        max_thickness = np.nan_to_num(max_thickness, nan=0.12, posinf=0.2, neginf=0.05)
        thickness_30 = np.nan_to_num(thickness_30, nan=0.12, posinf=0.2, neginf=0.05)
        le_radius = np.nan_to_num(le_radius, nan=0.02, posinf=0.1, neginf=0.01)
        curvature_penalty = np.nan_to_num(curvature_penalty, nan=1.0, posinf=10.0, neginf=0.0)

        fitness = (
            5.0 * max_thickness +  # Reward thickness
            3.0 * thickness_30 +   # Reward thickness at 30% chord
            1.0 * le_radius -      # Reward leading edge radius
            0.5 * curvature_penalty -  # Penalize high curvature
            2.0 * te_thickness     # Penalize thick trailing edge
        )

        # Ensure fitness is finite
        if not np.isfinite(fitness):
            fitness = -10.0  # Penalty for invalid designs

        return fitness

    def evaluate_population(self):
        """
        Evaluate the fitness of the entire population.
        """
        for i in range(self.population_size):
            # Apply design variables to create new airfoil
            airfoil = self._apply_bumps(self.population[i])

            # Evaluate fitness
            self.fitness[i] = self._evaluate_fitness(airfoil)

            # Update best design if needed
            if self.fitness[i] > self.best_fitness:
                self.best_fitness = self.fitness[i]
                self.best_design = self.population[i].copy()
                self.best_airfoil = airfoil

    def select_parents(self):
        """
        Select parents for reproduction using tournament selection.

        Returns
        -------
        np.ndarray
            Selected parents
        """
        # Tournament selection
        tournament_size = 3
        parents = np.zeros((self.population_size, 2 * self.n_design_vars))

        for i in range(self.population_size):
            # Select tournament participants
            participants = np.random.choice(self.population_size, tournament_size, replace=False)

            # Find the best participant
            best_idx = participants[np.argmax(self.fitness[participants])]

            # Select as parent
            parents[i] = self.population[best_idx]

        return parents

    def crossover(self, parents):
        """
        Perform crossover to create new offspring.

        Parameters
        ----------
        parents : np.ndarray
            Selected parents

        Returns
        -------
        np.ndarray
            New offspring
        """
        offspring = np.zeros_like(parents)

        for i in range(0, self.population_size, 2):
            # Select two parents
            parent1 = parents[i]
            parent2 = parents[i+1] if i+1 < self.population_size else parents[0]

            # Crossover point
            crossover_point = np.random.randint(1, 2 * self.n_design_vars - 1)

            # Create offspring
            offspring[i, :crossover_point] = parent1[:crossover_point]
            offspring[i, crossover_point:] = parent2[crossover_point:]

            if i+1 < self.population_size:
                offspring[i+1, :crossover_point] = parent2[:crossover_point]
                offspring[i+1, crossover_point:] = parent1[crossover_point:]

        return offspring

    def mutate(self, offspring, mutation_rate=0.1, mutation_scale=0.02):
        """
        Perform mutation on the offspring.

        Parameters
        ----------
        offspring : np.ndarray
            Offspring to mutate
        mutation_rate : float
            Probability of mutation for each gene
        mutation_scale : float
            Scale of mutation

        Returns
        -------
        np.ndarray
            Mutated offspring
        """
        for i in range(self.population_size):
            for j in range(2 * self.n_design_vars):
                if np.random.random() < mutation_rate:
                    offspring[i, j] += np.random.normal(0, mutation_scale)

        # Clip to bounds
        offspring = np.clip(offspring, -0.1, 0.1)

        return offspring

    def run_optimization(self, n_generations=20):
        """
        Run the genetic optimization.

        Parameters
        ----------
        n_generations : int
            Number of generations to run

        Returns
        -------
        dict
            Optimization results
        """
        print(f"Starting genetic optimization for {self.base_airfoil_name}")
        print(f"Population size: {self.population_size}, Generations: {n_generations}")

        # Start timer
        start_time = time.time()

        # Initial evaluation
        self.evaluate_population()

        # Store initial results
        self.history['generation'].append(0)
        self.history['best_fitness'].append(self.best_fitness)
        self.history['avg_fitness'].append(np.mean(self.fitness))
        self.history['best_design'].append(self.best_design.copy())

        # Run generations
        for gen in range(1, n_generations + 1):
            # Select parents
            parents = self.select_parents()

            # Crossover
            offspring = self.crossover(parents)

            # Mutation
            offspring = self.mutate(offspring)

            # Replace population
            self.population = offspring

            # Evaluate new population
            self.evaluate_population()

            # Store results
            self.history['generation'].append(gen)
            self.history['best_fitness'].append(self.best_fitness if np.isfinite(self.best_fitness) else -10.0)

            # Calculate average fitness, handling NaN values
            valid_fitness = self.fitness[np.isfinite(self.fitness)]
            avg_fitness = np.mean(valid_fitness) if len(valid_fitness) > 0 else -10.0
            self.history['avg_fitness'].append(avg_fitness)

            self.history['best_design'].append(self.best_design.copy())

            # Print progress
            best_fitness_str = f"{self.best_fitness:.4f}" if np.isfinite(self.best_fitness) else "N/A"
            print(f"Generation {gen}: Best fitness = {best_fitness_str}, Avg fitness = {avg_fitness:.4f}")

        # End timer
        end_time = time.time()
        optimization_time = end_time - start_time

        print("\nOptimization completed:")
        print(f"Total time: {optimization_time:.1f} seconds")
        print(f"Best fitness: {self.best_fitness:.4f}")

        # Create result plots
        self._create_result_plots()

        return {
            'best_design': self.best_design,
            'best_fitness': self.best_fitness,
            'best_airfoil': self.best_airfoil,
            'history': self.history,
            'optimization_time': optimization_time
        }

    def _create_result_plots(self):
        """
        Create result plots.
        """
        # Create convergence history plot
        fig1, ax1 = plt.subplots(figsize=(10, 6))

        ax1.plot(self.history['generation'], self.history['best_fitness'], 'b-', label='Best Fitness')
        ax1.plot(self.history['generation'], self.history['avg_fitness'], 'r-', label='Average Fitness')

        ax1.set_xlabel('Generation')
        ax1.set_ylabel('Fitness')
        ax1.set_title('Optimization Convergence History')
        ax1.grid(True)
        ax1.legend()

        fig1.savefig(f'{self.base_airfoil_name}_convergence.png', dpi=300)

        # Create airfoil comparison plot
        fig2, ax2 = plt.subplots(figsize=(10, 6))

        # Plot baseline airfoil
        ax2.plot(self.base_airfoil.x, self.base_airfoil.y, 'b-', label='Baseline')

        # Plot optimized airfoil
        if self.best_airfoil is not None:
            ax2.plot(self.best_airfoil.x, self.best_airfoil.y, 'r-', label='Optimized')

        ax2.set_xlabel('x/c')
        ax2.set_ylabel('y/c')
        ax2.set_title('Airfoil Comparison')
        ax2.grid(True)
        ax2.legend()
        ax2.set_aspect('equal')

        fig2.savefig(f'{self.base_airfoil_name}_comparison.png', dpi=300)

def run_genetic_optimization(airfoil_name='naca0012', n_design_vars=6, population_size=20, n_generations=20):
    """
    Run genetic optimization for an airfoil.

    Parameters
    ----------
    airfoil_name : str
        Name of the baseline airfoil
    n_design_vars : int
        Number of design variables
    population_size : int
        Size of the population
    n_generations : int
        Number of generations to run

    Returns
    -------
    dict
        Optimization results
    """
    # Create optimizer
    optimizer = SimpleGeneticOptimizer(
        base_airfoil=airfoil_name,
        n_design_vars=n_design_vars,
        population_size=population_size
    )

    # Run optimization
    results = optimizer.run_optimization(n_generations=n_generations)

    return results

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="PyMISES Simple Genetic Optimization Example")
    parser.add_argument('--airfoil', type=str, default='naca0012', help="Baseline airfoil name")
    parser.add_argument('--n-vars', type=int, default=6, help="Number of design variables per surface")
    parser.add_argument('--pop-size', type=int, default=20, help="Population size")
    parser.add_argument('--n-gen', type=int, default=20, help="Number of generations")

    args = parser.parse_args()

    # Run optimization
    results = run_genetic_optimization(
        airfoil_name=args.airfoil,
        n_design_vars=args.n_vars,
        population_size=args.pop_size,
        n_generations=args.n_gen
    )

    # Show plots
    plt.show()
