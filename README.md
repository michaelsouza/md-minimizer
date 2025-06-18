# md-minimizer

### Summary of Noguchi et al. (2024)

Noguchi et al. (2024) conducted a series of computer simulations to investigate how the internal structure of a composite material affects its fracture process. They created a two-dimensional model of a material made of both breakable and unbreakable components. By systematically changing the arrangement of these components, they observed how the material's strength, the way it breaks, and the statistical patterns of fracture change. The key finding is that a regular internal structure, like a reinforcing matrix, can make the material behave in a more ductile (less brittle) manner and that the size of this internal structure dictates the scale of fracture events.

#### The Model: A Virtual Composite Material

The experiment was performed using a **spring network model**. Imagine a 2D sheet of particles arranged in a triangular lattice. These particles are connected to their nearest neighbors by springs. This network of springs represents the atomic or molecular bonds within a material.

* **System Setup:** The model consists of a grid of $N \times N$ particles (where $N=96$ in this study).
* **Two Types of Springs:**
    * **Breakable springs:** These represent the primary, weaker material. Each breakable spring has a randomly assigned breaking point. When the strain (stretch) on the spring exceeds this threshold, it breaks and is permanently removed from the network.
    * **Unbreakable springs:** These represent a stronger, reinforcing material. They are infinitely strong and cannot break.
* **Internal Structure (The Matrix):** The unbreakable springs are not placed randomly. They are arranged in a regular, grid-like pattern to form a matrix within the material. The size of the cells in this matrix is defined by a parameter called $L_{\text{matrix}}$. A smaller $L_{\text{matrix}}$ means a denser, finer matrix, while a larger $L_{\text{matrix}}$ corresponds to a coarser matrix. The study specifically used $L_{\text{matrix}}$ values of 6, 8, and 12. A "normal system" with no unbreakable springs was used as a baseline for comparison.

#### The Simulation Process: Simulating Stretching

The researchers simulated a tensile test on their virtual material to see how it would fracture under stress.

1.  **Applying Strain:** The bottom row of particles in the network was held fixed, while the top row was pulled upwards by a very small, incremental amount ($0.1\%$ of the system's height). This simulates stretching the material.
2.  **Relaxation:** After each small stretch, the system (except the bottom and top rows) was allowed to reach a state of mechanical equilibrium. This means the particles shifted around until all the forces from the springs were balanced.
3.  **Checking for Breaks:** In this new equilibrium state, the strain on every single spring was calculated. If any spring's strain exceeded its predetermined breaking threshold, it was removed.
4.  **Avalanche of Breaks:** The removal of one spring redistributes its load to its neighbors, which can cause them to break. This can trigger a chain reaction, or an **avalanche**, where many springs break in a single step before the system finds a new stable equilibrium.
5.  **Iteration:** This process of stretching, relaxing, and breaking springs was repeated until the entire material fractured into two separate pieces or the total strain reached a maximum limit.
6.  **Data Collection:** Throughout the simulation, the researchers recorded key data, including the overall stress on the system, the number of springs that broke in each avalanche, and the size and location of crack clusters. This entire process was repeated 1000 times with different random configurations of breakable spring thresholds to ensure the results were statistically robust.

#### Key Experiments and Findings

The study focused on three main experimental observations:

##### 1. How the Matrix Changes Material Behavior (Stress-Strain Response)

* **Experiment:** The researchers compared the stress-strain curves of the "normal system" (no matrix) with the "matrix-mixture systems" ($L_{\text{matrix}} = 6, 8, 12$).
* **Finding:** The normal system behaved like a **brittle** material. It stretched, built up stress, and then suddenly failed completely with one catastrophic drop in stress. In contrast, the systems with the unbreakable matrix behaved like a **ductile** material. After an initial elastic phase, they showed a series of smaller, intermittent stress drops, allowing them to withstand much larger overall strain before complete failure. The matrix effectively stops large cracks from propagating through the entire material.

##### 2. The Statistics of Fracture Avalanches

* **Experiment:** They measured the size ($S$) of each avalanche, which is the number of springs that break in a single loading step. They then created a distribution, $P(S)$, showing how often avalanches of different sizes occurred.
* **Finding:** In all systems, small avalanches are much more common than large ones, following a power-law relationship. However, the presence of the matrix introduced a **cut-off size**. The fracture avalanches could not grow larger than a size determined by the matrix cell size ($L_{\text{matrix}}$). A smaller, denser matrix suppressed large avalanches more effectively. The researchers developed a scaling law showing that all the different avalanche distributions would collapse onto a single "master curve" when the avalanche size was scaled by the number of breakable springs within a single matrix cell.

##### 3. The Link Between Crack Growth and Stress Drops

* **Experiment:** The researchers analyzed the geometry of the broken springs, treating connected broken springs as "crack clusters." They introduced a quantity, $\Delta C$, to measure the growth of these clusters in each step and correlated it with the corresponding stress drop, $\Delta\Sigma$.
* **Finding:** A strong relationship was found: when new breaks occurred far from existing cracks, the stress drop was small. However, when new breaks linked up with existing cracks to form a larger, more significant crack, the resulting stress drop was much larger. They showed that the average stress drop increases as a power law with the growth of the crack cluster size ($\Delta C$). Furthermore, they found that the maximum size a crack could attain was limited by the matrix, and this relationship could also be described by a scaling exponent that was consistent with the one found for the stress drop analysis. This confirmed that the matrix directly controls crack growth, which in turn governs the magnitude of stress drops and the overall ductile behavior of the material.

#### Executando os Experimentos

1) Para rodar os experimentos, precisamos primeiro ativar o ambiente virtual Python do Anaconda.

```bash
conda activate lammps
```

2) Com o ambiente ativado, o primeiro passo seria criar uma rede de partículas utilizando o arquivo `create_network.py`. Este script gera um arquivo de entrada para o LAMMPS, que é o simulador utilizado.

```bash
python lammps/create_network.py --N 96 --L_matrix 6 --output networks
```
Aqui voce passa um tamanho de rede $N$ e o tamanho da matriz $L_{\text{matrix}}$. O local de saida `network_6.txt` conterá a configuração inicial da rede.

3) Para executar a simulação de fato, vamos usar o comando

```bash