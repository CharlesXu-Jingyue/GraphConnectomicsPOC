# GraphConnectomicsPOC

By Charles Xu @ Caltech, 20250714

## Project Plan

### 1 Basic Characterization and Preprocessing

#### 1.1 Data Validation and Preprocessing

- **Verify matrix properties**: Check if directed/undirected, presence of self-connections
- **Sparsity analysis**: Calculate connection density $\rho = \frac{\text{nnz}(W)}{n^2}$ and distribution of synaptic weights
- **Weight distribution**: Analyze histogram of non-zero weights, identify outliers
- **Symmetry assessment**: Quantify asymmetry via $\|W - W^T\|_F / \|W\|_F$
- **Initial visualization**: Heatmap with appropriate ordering (e.g., by degree)

#### 1.2 Fundamental Graph Statistics

- **Degree sequences**: 
  - In-degree: $k_i^{\text{in}} = \sum_{j} W_{ji}$
  - Out-degree: $k_i^{\text{out}} = \sum_{j} W_{ij}$
  - Joint degree distribution $P(k^{\text{in}}, k^{\text{out}})$
  - Identify hub neurons ($k > \mu_k + 2\sigma_k$)
- **Weight-degree correlations**: Compute $\text{corr}(k_i, \langle w_i \rangle)$ where $\langle w_i \rangle$ is mean connection weight
- **Reciprocity**: $r = \frac{\sum_{i,j} W_{ij}W_{ji}}{\sum_{i,j} W_{ij}}$

### 2 Local Structure and Motifs

#### 2.1 Local Connectivity Patterns

- **Weighted clustering coefficient** (Onnela et al.):
  $$C_i = \frac{1}{k_i(k_i-1)} \sum_{j,k} (W_{ij}W_{jk}W_{ki})^{1/3}$$
- **Local efficiency**: 
  $$E_{\text{loc}}(i) = \frac{1}{k_i(k_i-1)} \sum_{j,k \in \mathcal{N}_i} \frac{1}{d_{jk}}$$
- **Rich-club coefficient**: 
  $$\phi(k) = \frac{2E_{>k}}{N_{>k}(N_{>k}-1)}$$
  where $E_{>k}$ is edges between nodes with degree $>k$

#### 2.2 Motif Analysis

- **3-node motifs**:
  - Count all 13 directed 3-node motifs
  - Z-score: $z_i = \frac{N_i^{\text{real}} - \langle N_i^{\text{rand}} \rangle}{\sigma_{N_i^{\text{rand}}}}$
  - Focus on:
    - Feedforward loops (motif 9)
    - Bi-fan (motif 5)
    - Bi-parallel (motif 4)
- **4-node motifs** (computational cost permitting):
  - Focus on biologically relevant patterns (e.g., winner-take-all circuits)
- **Motif participation matrix**: $M_{ij}$ = number of motifs where neurons $i$ and $j$ co-participate

### 3 Mesoscale Organization

#### 3.1 Community Detection

- **Multi-resolution modularity optimization**:
  $$Q(\gamma) = \frac{1}{2m} \sum_{ij} \left[W_{ij} - \gamma \frac{k_i^{\text{out}} k_j^{\text{in}}}{2m}\right] \delta(c_i, c_j)$$
  - Scan resolution parameter $\gamma \in [0.5, 2.0]$
  - Identify stable partitions via variation of information
- **Hierarchical clustering**:
  - Distance metric: $d_{ij} = 1 - \frac{W_{ij} + W_{ji}}{\max(W)}$
  - Use Ward's method for agglomeration
- **Overlapping communities**: 
  - Link communities or OSLOM algorithm
  - Quantify overlap: $O_i = \sum_c p_{ic} \log p_{ic}$ where $p_{ic}$ is membership strength

#### 3.2 Module Characterization

- **Participation coefficient**:
  $$P_i = 1 - \sum_{m=1}^{N_M} \left(\frac{k_{im}}{k_i}\right)^2$$
  where $k_{im}$ is connections from node $i$ to module $m$
- **Within-module degree z-score**:
  $$z_i = \frac{k_i^{\text{within}} - \bar{k}_m}{\sigma_{k_m}}$$
- **Inter-module connectivity matrix**: $\mathbf{C}_{mn} = \sum_{i \in m, j \in n} W_{ij}$

#### 3.3 Dimensionality Reduction to Module-Level Connectivity

##### 3.3.1 Community-Based Reduction

- **Module connectivity matrix**:
  $$\tilde{W}_{mn} = \frac{1}{|m||n|} \sum_{i \in m} \sum_{j \in n} W_{ij}$$
  for modules $m \neq n$, and
  $$\tilde{W}_{mm} = \frac{1}{|m|(|m|-1)} \sum_{i,j \in m, i \neq j} W_{ij}$$
- **Module coupling strength**: Normalize by within-module density
  $$S_{mn} = \frac{\tilde{W}_{mn}}{\sqrt{\tilde{W}_{mm} \tilde{W}_{nn}}}$$
- **Effective module dynamics**: Project neuron-level dynamics
  $$\dot{\mathbf{x}}_m = \sum_n \tilde{W}_{mn} f(\mathbf{x}_n) + \mathbf{b}_m$$

##### 3.3.2 Spectral Reduction

- **SVD-based reduction**: Decompose $\mathbf{W} = \mathbf{U}\mathbf{\Sigma}\mathbf{V}^T$
  - Keep top $k$ singular values capturing $\geq 90\%$ variance
  - Reduced connectivity: $\tilde{\mathbf{W}} = \mathbf{U}_k\mathbf{\Sigma}_k\mathbf{V}_k^T$
- **Eigenvector-based modules**:
  - Group neurons by similarity in eigenvector space
  - Define similarity: $s_{ij} = \sum_{k=1}^{r} v_k(i)v_k(j)$ for top $r$ eigenvectors
  - Apply spectral clustering on similarity matrix
- **Schur decomposition for invariant subspaces**:
  $$\mathbf{W} = \mathbf{Q}\mathbf{T}\mathbf{Q}^*$$
  - Identify block structure in upper triangular $\mathbf{T}$
  - Each block represents dynamically invariant subspace

##### 3.3.3 Hybrid Approach

- **Community-constrained SVD**:
  - Compute SVD within each detected community
  - Construct block-diagonal approximation with inter-module connections
- **Participation-weighted reduction**:
  $$\mathbf{W}_{\text{reduced}} = \mathbf{P}^T\mathbf{W}\mathbf{P}$$
  where $\mathbf{P}_{im} = \frac{p_{im}}{\sqrt{\sum_i p_{im}^2}}$ is the normalized participation matrix

### 4 Spectral Analysis

#### 4.1 Eigenspectrum Analysis

- **Full eigendecomposition**: $\mathbf{W}\mathbf{v}_i = \lambda_i \mathbf{v}_i$
  - Plot $\{\lambda_i\}$ in complex plane
  - Identify spectral radius $\rho(\mathbf{W}) = \max_i |\lambda_i|$
- **Spectral gaps**: $\Delta_i = |\lambda_i - \lambda_{i+1}|$
  - Large gaps indicate timescale separation
- **Pseudospectrum**: 
  $$\sigma_\epsilon(\mathbf{W}) = \{z \in \mathbb{C} : \|(z\mathbf{I} - \mathbf{W})^{-1}\| \geq \epsilon^{-1}\}$$
- **Gershgorin circles**: Localize eigenvalues via
  $$\lambda \in \bigcup_i \{z : |z - W_{ii}| \leq \sum_{j \neq i} |W_{ij}|\}$$

#### 4.2 Eigenvector Structure

- **Leading eigenvectors**:
  - Spatial patterns of top 10-20 eigenvectors
  - Correlation with node properties (degree, module membership)
- **Participation ratio**:
  $$\text{PR}_i = \frac{(\sum_j |v_{ij}|^2)^2}{\sum_j |v_{ij}|^4}$$
- **Schur decomposition**: Identify invariant subspaces and their connectivity
- **Eigenvector centrality alignment**: $\text{corr}(\mathbf{v}_1, \mathbf{k})$
- **IPR (Inverse Participation Ratio)**: $\text{IPR}_i = \sum_j |v_{ij}|^4$
  - Low IPR → delocalized mode
  - High IPR → localized mode

### 5 Feedback Architecture Analysis

#### 5.1 Cycle Detection and Analysis

- **Johnson's algorithm** for all simple cycles up to length $L$
- **Cycle weight**: For cycle $c = (i_1, i_2, ..., i_k, i_1)$:
  $$w_c = \prod_{j=1}^k W_{i_j, i_{j+1}}$$
- **Strongly connected components (SCCs)**:
  - Tarjan's algorithm for SCC detection
  - Condensation graph $\mathcal{G}_c$ of SCCs
- **Feedback participation**: $F_i = $ number of cycles containing node $i$

#### 5.2 Signal Flow Analysis

- **Trophic levels** (MacKay et al.):
  $$s_i = 1 + \frac{1}{k_i^{\text{in}}} \sum_{j} W_{ji} s_j$$
- **Flow hierarchy**: $h = \frac{\sum_{ij} W_{ij} \text{sgn}(s_j - s_i)}{\sum_{ij} W_{ij}}$
- **Feedback centrality** (generalized Katz):
  $$\mathbf{x} = \alpha(\mathbf{W} + \mathbf{W}^T)\mathbf{x} + \beta\mathbf{1}$$
- **Loop transfer functions**: For identified feedback loops, compute effective gain

### 6 Dynamical Predictions

#### 6.1 Linear Dynamics Analysis

- **Continuous-time evolution**: $\dot{\mathbf{x}} = \mathbf{W}\mathbf{x}$
  - Solution: $\mathbf{x}(t) = e^{\mathbf{W}t}\mathbf{x}(0)$
- **Timescale hierarchy**: Group eigenvalues by $|\text{Re}(\lambda_i)|$
  - Fast modes: $|\text{Re}(\lambda)| > \lambda_{\text{threshold}}$
  - Slow modes: $|\text{Re}(\lambda)| < \lambda_{\text{threshold}}$
- **Transient amplification**:
  $$G(t) = \|e^{\mathbf{W}t}\|_2^2 = \sigma_{\max}^2(e^{\mathbf{W}t})$$
  - Maximum amplification: $G_{\max} = \sup_{t>0} G(t)$
- **Frequency response**: $\mathbf{H}(\omega) = (i\omega\mathbf{I} - \mathbf{W})^{-1}$

#### 6.2 Nonlinear Dynamics Predictions

- **Fixed point stability** (for $\dot{\mathbf{x}} = -\mathbf{x} + \mathbf{W}f(\mathbf{x})$):
  - Jacobian at fixed point: $\mathbf{J} = -\mathbf{I} + \mathbf{W}\text{diag}(f'(\mathbf{x}^*))$
- **Lyapunov exponents**: For trajectory $\mathbf{x}(t)$:
  $$\lambda_i = \lim_{t \to \infty} \frac{1}{t} \log \|\mathbf{v}_i(t)\|$$
- **Mean-field reduction** for modules:
  $$\dot{m}_\alpha = -m_\alpha + \sum_\beta J_{\alpha\beta} \phi(m_\beta)$$
  where $J_{\alpha\beta} = \tilde{W}_{\alpha\beta} N_\beta$
- **Attractor prediction**:
  - Based on SCC structure and eigenspectrum
  - Estimate attractor dimensionality from participation ratios

### 7 Control-Theoretic Analysis

#### 7.1 Structural Controllability

- **Minimum driver nodes**: Maximum matching in bipartite representation
- **Control matrix**: $\mathbf{B} \in \mathbb{R}^{n \times m}$ where $m$ is number of drivers
- **Controllability Gramian**:
  $\mathbf{W}_c = \int_0^T e^{\mathbf{W}t}\mathbf{B}\mathbf{B}^T e^{\mathbf{W}^T t} dt$ for different input configurations
- **Control energy**: For state transfer $\mathbf{x}_0 \to \mathbf{x}_f$:
  $$E = \min_{\mathbf{u}(t)} \int_0^T \|\mathbf{u}(t)\|^2 dt = (\mathbf{x}_f - e^{\mathbf{W}T}\mathbf{x}_0)^T \mathbf{W}_c^{-1} (\mathbf{x}_f - e^{\mathbf{W}T}\mathbf{x}_0)$$

#### 7.2 Modal Controllability

- **Mode controllability**: For eigenvector $\mathbf{v}_i$:
  $$\phi_i = \mathbf{v}_i^T \mathbf{B}\mathbf{B}^T \mathbf{v}_i^*$$
- **Average controllability**: $\bar{\phi} = \frac{1}{n}\text{Tr}(\mathbf{W}_c)$
- **Modal energy**: Energy to excite mode $i$:
  $$E_i = \frac{1 - e^{2\text{Re}(\lambda_i)T}}{2\text{Re}(\lambda_i)\phi_i}$$

#### 7.3 Target Control

- **Module-specific control**: Identify neurons to control specific modules
- **Output controllability**: If certain neurons are readouts, analyze their controllability
- **Feedback control synthesis**: Design local feedback to stabilize/destabilize specific modes

### 8 Computational Function Predictions

#### 8.1 Integration of Analyses

- **Structure-function mapping**:
  - Map module connectivity $\tilde{W}_{mn}$ to functional interactions
  - Classify modules by motif enrichment:
    - Feedforward processing (high motif 9)
    - Recurrent computation (high motif 6)
    - Convergent integration (high motif 3)
- **Information capacity**:
  $$C = \log_2 \det(\mathbf{I} + \text{SNR} \cdot \mathbf{W}\mathbf{W}^T)$$

#### 8.2 Testable Predictions

##### 8.2.1 Dynamic Signatures

- **Oscillation frequencies**: From imaginary parts of eigenvalues
  $$f_i = \frac{|\text{Im}(\lambda_i)|}{2\pi}$$
- **Response timescales**: $\tau_i = -1/\text{Re}(\lambda_i)$

##### 8.2.2 Perturbation Responses

- **Lesion predictions**: Effect of removing node $i$:
$$\Delta\lambda_{\max} \approx \frac{\mathbf{u}_1(i)\mathbf{v}_1(i)}{\mathbf{u}_1^T\mathbf{v}_1}$$
- Response to stimulation of different modules

##### 8.2.3 Computational Capacity

- **Memory capacity** (for echo state networks):
  $$\text{MC} = \sum_{k=1}^{\infty} \text{corr}^2(u(t-k), y_k(t))$$
- Predict integration vs segregation from modularity

### 9 Visualization and Reporting

#### 9.1 Network Visualizations

- **Module-level graph**: Nodes sized by module size, edges by $\tilde{W}_{mn}$
- **Spectral embedding**: Position nodes using top eigenvectors
  $$\mathbf{X} = [\mathbf{v}_2, \mathbf{v}_3, ..., \mathbf{v}_{k+1}]$$
- **Control node highlighting**: Color by controllability metrics

#### 9.2 Summary Statistics Dashboard

Create comprehensive tables including:
- **Multi-scale network statistics**
- **Module properties and inter-module connectivity**
- **Spectral properties and timescales**
- **Control node rankings**
- **Dynamical predictions summary**

#### 9.3 Computational Notebook

Organize all analyses in executable notebook with:
- Modular functions for each analysis phase
- Intermediate result caching
- Parameter sensitivity analysis
- Interactive visualizations

---

## Implementation

### Project Structure

This project uses a hybrid approach combining modular packages for reusable components and Jupyter notebooks for interactive analysis and visualization.

```
GraphConnectomicsPOC/
├── src/                         # Core analysis modules
│   ├── preprocessing/           # Data validation & preprocessing
│   ├── local_structure/         # Motifs, clustering, local metrics
│   ├── mesoscale/               # Community detection, modules
│   ├── spectral/                # Eigenanalysis, decompositions
│   ├── dynamics/                # Linear/nonlinear predictions
│   ├── control/                 # Controllability analysis
│   └── visualization/           # Plotting utilities
├── notebooks/                   # Analysis workflows
│   ├── 01_preprocessing.ipynb   # Data validation and basic stats
│   ├── 02_local_structure.ipynb # Motifs and local connectivity
│   ├── 03_mesoscale.ipynb       # Community detection
│   ├── 04_spectral.ipynb        # Eigenspectrum analysis
│   ├── 05_dynamics.ipynb        # Dynamical predictions
│   └── 06_control.ipynb         # Control-theoretic analysis
├── data/                        # Raw and processed connectome data
├── results/                     # Output figures, tables, reports
└── tests/                       # Unit tests for analysis functions
```

#### Design Philosophy

- **Modules** provide reusable, testable functions for complex computations (eigendecomposition, community detection, control theory)
- **Notebooks** enable interactive exploration, visualization, and narrative documentation of results
- **Separation of concerns** allows others to use analysis methods independently of full workflows

### Data

**Updated 20250716**

Project data is based on the FlyWire *Drosophila* central brain connectome. The analysis focuses on a subset of 475 neurons from 18 cell types of interest (AOTU046, Delta7, EL, ER2, ER3a, ER3d, ER3m, ER3p, ER3w, ER4d, ER4m, EPG, TuBu01-03, TuBu07-09), generating a 475×475 connectivity matrix with weighted synaptic connections.

### Analysis Approach

This project employs an **unsupervised discovery paradigm**:

1. **Blind Analysis**: Perform all connectivity analyses (clustering, motifs, spectral analysis, etc.) without using known cell type labels
2. **Post-hoc Validation**: Compare discovered patterns with ground truth cell type classifications to validate biological relevance
3. **Novel Discovery**: Identify connectivity patterns that may reveal previously unknown functional relationships

The shuffled neuron ordering ensures unbiased analysis, preventing visual clustering artifacts based on cell type labels.