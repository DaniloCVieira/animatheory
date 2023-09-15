# Simulating Ecological Communities Under Assumptions of Environmental Filtering and Neutrality
This project provides a spatially explicit, individual-based model that simulates ecological communities under different assumptions. Users can explore six different models over varying time steps, visualizing the dynamics through graphical displays.

## Link:
https://danilocvieira.github.io/animatheory/

## Overview
The simulation operates on a grid with the following characteristics:
- Each cell can support only one individual.
- The grid initially contains 49 individuals randomly sampled from a pool of 20 species and placed randomly on the grid.

### Parameters
- Total number of individuals: 49
- Total number of species: 20
- Base birth rate: 0.5
- Base death rate: 0.5
- Individual fitness range: 1 to 20 (unique for each species)

## Simulation Steps
The simulation progresses through a series of steps within 60 (for neutral models) or 100 (for environmental filtering models) time steps. Each step is as follows:

1. **Death**
    - 10 randomly selected individuals die based on a probability (d), which can depend on either:
      - Random selection (neutral model).
      - Environmental filtering (the higher the difference between environmental value (EV) and species fitness, the higher the death probability).

2. **Movement** (if applied)
    - Individuals have a 0.1 probability to move to a different cell.

3. **Immigration** (if applied)
    - New individuals are sampled from the species pool with a 0.05 probability of immigration.

4. **Birth**
    - New births can occur with a probability (b), which can be influenced by either:
      - Random selection.
      - Environmental filtering (the birth probability is higher when EV and fitness are more similar).

5. **Dispersal Limitation** (if applied)
    - New individuals can only occupy neighboring cells.

## Models
Users can explore six different models, each with distinct rules and parameters:
1. Neutral Model: Random death and birth.
2. Neutral Model with Dispersal Limitation: Random death, birth, and dispersal limitation.
3. Neutral Model with Immigration: Random death, birth, immigration, and dispersal limitation.
4. Niche Model: Environmental filtering affecting birth.
5. Niche Model with Dispersal Limitation: Environmental filtering affecting birth with dispersal limitation.
6. Niche Model with Environmental Filtering: Environmental filtering affecting birth, death, and dispersal limitation.

## Graphical Displays
Each model features three graphical displays to visualize the dynamics:
1. **Grid Display**: A 7x7 grid showing individuals represented by animal GIFs (each different species is depicted by a unique GIF) and fitness represented by color gradients.
2. **Species Richness Plot**: A plot depicting time steps versus species richness.
3. **Rank Abundance Plot**: A graph showing species abundance rank versus the number of species.

Background colors (representing EVs) vary depending on the model: white for neutral models and a gradient for environmental filtering models (red, orange, white, light blue, dark blue representing EVs of 1, 5, 10, 15, 20 respectively).

## User Interaction
Users can interact with the simulation in the following ways:
- **Play Button**: Automatically progresses through all time steps sequentially.
- **Step-by-Step Progression**: Allows users to move through time steps one by one, observing how the community evolves over time.


