# PyLLG

Scipy-based Landau-Lifshitz-Gilbert (LLG) equation solver for a magnetic system, described by a classicalmodel including Heisenberg exchange interaction between nearest-neighbor spins, current-driven exchange interaction between lattice spins and nonequilibrium electrons, various anisotropy terms, Dzyaloshinskii-Moriya interaction, and an external field. 

The main equation to be integrated is the Landau-Lifshitz-Gilbert equation (LLG)  
   
<a href="https://www.codecogs.com/eqnedit.php?latex=\partial_t&space;\mathbf{s}_i&space;=&space;-\frac{g}{1&plus;\lambda^2}\left(&space;\mathbf{s}_i\times&space;\mathbf{B}^\mathrm{eff}_i&space;&plus;&space;\lambda&space;\mathbf{s}_i\times&space;[\mathbf{s}_i\times&space;\mathbf{B}_i^\mathrm{eff}]&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\partial_t&space;\mathbf{s}_i&space;=&space;-\frac{g}{1&plus;\lambda^2}\left(&space;\mathbf{s}_i\times&space;\mathbf{B}^\mathrm{eff}_i&space;&plus;&space;\lambda&space;\mathbf{s}_i\times&space;[\mathbf{s}_i\times&space;\mathbf{B}_i^\mathrm{eff}]&space;\right&space;)" title="\partial_t \mathbf{s}_i = -\frac{g}{1+\lambda^2}\left( \mathbf{s}_i\times \mathbf{B}^\mathrm{eff}_i + \lambda \mathbf{s}_i\times [\mathbf{s}_i\times \mathbf{B}_i^\mathrm{eff}] \right )" /></a>  
  
which governs the time-dependent evolution of a spin <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{s}_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{s}_i" title="\mathbf{s}_i" /></a> on site <a href="https://www.codecogs.com/eqnedit.php?latex=i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?i" title="i" /></a>.  
The effective field is computed by  

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{B}^\mathrm{eff}_i&space;=&space;-\mu_M^{-1}&space;\partial_{\mathbf{s}_i}\mathcal{H}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{B}^\mathrm{eff}_i&space;=&space;-\mu_M^{-1}&space;\partial_{\mathbf{s}_i}\mathcal{H}" title="\mathbf{B}^\mathrm{eff}_i = -\mu_M^{-1} \partial_{\mathbf{s}_i}\mathcal{H}," /></a>  
  
with the classical Hamiltonian capturing all dominant interactions  
  
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{H}&space;=&space;-J\sum_{\langle&space;i,j\rangle}&space;\mathbf{s}_i\cdot&space;\mathbf{s}_j&space;-&space;J_{sd}\sum_i&space;\mathbf{s}_i^\mathrm{CD}\cdot&space;\mathbf{s}_i-\sum_{\langle&space;i,j\rangle}\mathbf{A}_{ij}\cdot(\mathbf{s}_i\times&space;\mathbf{s}_j)&space;-K\sum_i&space;(s_{i}^{z})^2&space;&plus;&space;D\sum_i&space;(s_i^y)^2&plus;Q\sum_i&space;[(s_i^x)^4&plus;(s_i^y)^4&plus;(s_i^z)^4]&space;-&space;P\sum_{\langle&space;i,j\rangle}&space;[s_i^x&space;s_j^x&plus;s_i^y&space;s_j^y]&space;-\mu_M&space;\mathbf{H}\sum_i&space;\mathbf{s}_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{H}&space;=&space;-J\sum_{\langle&space;i,j\rangle}&space;\mathbf{s}_i\cdot&space;\mathbf{s}_j&space;-&space;J_{sd}\sum_i&space;\mathbf{s}_i^\mathrm{CD}\cdot&space;\mathbf{s}_i-\sum_{\langle&space;i,j\rangle}\mathbf{A}_{ij}\cdot(\mathbf{s}_i\times&space;\mathbf{s}_j)&space;-K\sum_i&space;(s_{i}^{z})^2&space;&plus;&space;D\sum_i&space;(s_i^y)^2&plus;Q\sum_i&space;[(s_i^x)^4&plus;(s_i^y)^4&plus;(s_i^z)^4]&space;-&space;P\sum_{\langle&space;i,j\rangle}&space;[s_i^x&space;s_j^x&plus;s_i^y&space;s_j^y]&space;-\mu_M&space;\mathbf{H}\sum_i&space;\mathbf{s}_i" title="\mathcal{H} = -J\sum_{\langle i,j\rangle} \mathbf{s}_i\cdot \mathbf{s}_j - J_{sd}\sum_i \mathbf{s}_i^\mathrm{CD}\cdot \mathbf{s}_i-\sum_{\langle i,j\rangle}\mathbf{A}_{ij}\cdot(\mathbf{s}_i\times \mathbf{s}_j) -K\sum_i (s_{i}^{z})^2 + D\sum_i (s_i^y)^2+Q\sum_i [(s_i^x)^4+(s_i^y)^4+(s_i^z)^4] - P\sum_{\langle i,j\rangle} [s_i^x s_j^x+s_i^y s_j^y] -\mu_M \mathbf{H}\sum_i \mathbf{s}_i." /></a>  
  
Here, <a href="https://www.codecogs.com/eqnedit.php?latex=\langle&space;i,&space;j\rangle" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\langle&space;i,&space;j\rangle" title="\langle i, j\rangle" /></a> denotes pairs of nearest neighbors, and the sum over it captures the lattice geometry and the boundary conditions. Within this module, following boundary conditions are implemented: 
- Periodic in both x and y direction
- Periodic in x only, y open
- Periodic in y only, x open
- Open boundaries in x and y direction

on a square lattice of dimensions X,Y.  
    
The parameters of this Hamiltonian are  
- Nearest-neighbor exchange interaction <a href="https://www.codecogs.com/eqnedit.php?latex=J" target="_blank"><img src="https://latex.codecogs.com/gif.latex?J" title="J" /></a> between neighboring spins in the lattice
- Exchange interaction between conducting electrons and localized magnetic moments <a href="https://www.codecogs.com/eqnedit.php?latex=J_{sd}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?J_{sd}" title="J_{sd}" /></a>
- Dzyaloshinskii-Moriya interaction <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}_{ij}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}_{ij}" title="\mathbf{A}_{ij}" /></a>
- magnetic anisotropies <a href="https://www.codecogs.com/eqnedit.php?latex=K,&space;Q" target="_blank"><img src="https://latex.codecogs.com/gif.latex?K,&space;Q" title="K, Q" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=P" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P" title="P" /></a>
- Demagnetization <a href="https://www.codecogs.com/eqnedit.php?latex=D" target="_blank"><img src="https://latex.codecogs.com/gif.latex?D" title="D" /></a>
- External field <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{H}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{H}" title="\mathbf{H}" /></a>
- and <a href="https://www.codecogs.com/eqnedit.php?latex=\mu_M" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_M" title="\mu_M" /></a> is the magnitude of magnetic moment of the localized magnetic moments in the lattice.



 
