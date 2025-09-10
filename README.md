# Welcome to the FlowsInside Repository!

![flowsinside](https://github.com/user-attachments/assets/f4bebed2-2a9c-4a49-8748-034ff456319b)

FlowsInside: Intracellular Flux Balance Analysis Addon

FlowsInside is an intracellular Flux Balance Analysis (FBA) addon designed to integrate with agent-based modeling platforms such as PhysiCell
.
It provides a biologically grounded way to couple genome-scale metabolic models with dynamic cell simulations, allowing cells to compute metabolic fluxes internally and influence their growth, proliferation, and survival.

✨ Features

🔬 Intracellular Metabolism: Connects FBA models to individual cells in agent-based simulations.

⚡ Efficient Optimization: Built with GLPK
 for linear programming support.

🔄 Dynamic Coupling: Biomass flux can drive cell-cycle progression, volume increase, or custom behaviors.

🧩 Modular Integration: Works as a plugin; easily adaptable to new agent-based frameworks.

📂 File-based Models: Load genome-scale metabolic models (SBML, MATLAB .mat, or custom formats).

🚀 Installation

Clone the repository:

git clone https://github.com/yourusername/FlowsInside.git
cd FlowsInside


Make sure you have:

A C++17 compatible compiler

GLPK library
 installed

(Optional) MATIO
 if you want to load .mat metabolic models

Build with:

make

🛠️ Usage
Example: Linking FBA to an agent
#include "IntracellularFBA.h"

// Create FBA object from a model file
auto fba = IntracellularFBA("models/ecoli.xml");

// Run optimization and get biomass flux
double growth_rate = fba.runFBA();

// Link to cell update
cell->phenotype.cycle.data.transition_rate = growth_rate;

📊 Applications

Cancer metabolism modeling (e.g., glioblastoma microenvironment)

Linking nutrient uptake to cell proliferation

Simulating drug interventions at metabolic level

Testing metabolic network knockouts in agent-based systems

📖 Documentation

Getting Started Guide

API Reference

Examples

🤝 Contributing

Contributions are welcome!
If you’d like to add features or fix bugs:

Fork the repo

Create a feature branch

Open a pull request

📜 License

This project is licensed under the MIT License – see the LICENSE
 file for details.

👨‍🔬 Citation

If you use this addon in scientific work, please cite:

Your Name, Year. FlowsInside: Intracellular Flux Balance Analysis Addon.
GitHub repository, https://github.com/yourusername/FlowsInside
