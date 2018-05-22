### SU(N) invariant interaction decouplings

The simplest three kinds of form factors maintaining the SU(N) symmetry: (where $$\alpha=1\cdots N$$)

+ $$\frac{g_1}{2}\left(c_{i\alpha}^\dagger c_{i\alpha}-\frac{N}{2}\right)^2$$, with form factor $$1$$.
+ $$\frac{g_{2\pm}}{2}\left(c_{i\alpha}^\dagger c_{i\alpha}\pm c_{j\beta}^\dagger c_{j\beta}\right)^2$$, with form factors $$\sigma_{0,3}$$.
+ $$\frac{g_{3\pm}}{2}\left(c_{i\alpha}^\dagger c_{j\alpha}\pm c_{j\beta}^\dagger c_{i\beta}\right)^2$$, with form factors $$\sigma_{1,2}$$.

The three kinds of relevant interactions are:

+ $$\frac{U}{2} c_{i\alpha}^\dagger c_{i\alpha} c_{i\beta}^\dagger c_{i\beta}$$
+ $$V c_{i\alpha}^\dagger c_{i\alpha} c_{j\beta}^\dagger c_{j\beta}$$
+ $$J c_{i\alpha}^\dagger c_{i\beta} c_{j\beta}^\dagger c_{j\alpha}$$


One solution is

+ $$g_1=\frac{U-V}{2}$$
+ $$g_{2+}=V$$
+ $g_{2-}=0$
+ $$g_{3\pm}=\pm\frac{J}{2}$$

General form of electron-phonon interactions (with onsite phonon):

+ $$\phi_i f_{ijk} c_{j\alpha}^\dagger c_{k\alpha}$$

In special, two well-known examples are given:

+ Holsetin model: $$\phi_ic_{i\alpha}^\dagger c_{i\alpha}$$
+ Peierls model: $$(\phi_j-\phi_i)(c_{i\alpha}^\dagger c_{j\alpha} + h.c.)$$

In general, phonons can also live on bond, then we may expect:

+ $$\phi_{ij}f_{ijkl}c_{k\alpha}^\dagger c_{l\alpha} $$


### phonon

$$Z=\exp\left[\Delta\tau gx_{i\tau} n_i-\Delta\tau\frac{M\Omega^2}{2}x_{i\tau}^2-\Delta\tau\frac{M}{2}\left(\frac{x_{i\tau+\Delta\tau}-x_{i\tau}}{\Delta\tau}\right)^2\right]$$

define $$V=g^2/M\Omega^2$$ and $$\phi_{i\tau}=\Delta\tau gx_{i\tau}$$, we get

$$Z=\exp\left[\phi_{i\tau} n_i-\frac{\phi_{i\tau}^2}{2\Delta\tau V}-\frac{1}{2\Delta\tau V\Omega^2}\left(\frac{\phi_{i\tau+\Delta\tau}-\phi_{i\tau}}{\Delta\tau}\right)^2\right]$$

For each $$\phi_\alpha$$-phonon, we need $$V_\alpha$$ and $$\Omega_\alpha$$. For a local update $$\phi_{i\tau}\rightarrow\phi_{i\tau}'$$, the accept probability is

$$R=(\det{D})^N\exp\left\{(\phi_{i\tau}'-\phi_{i\tau})\left[ -\frac{1}{2\Delta\tau V}(1+\frac{2}{\Omega^2\Delta\tau^2})(\phi_{i\tau}'+\phi_{i\tau})+\frac{1}{\Delta\tau^3 V\Omega^2}(\phi_{i\tau+\Delta\tau}+\phi_{i\tau-\Delta\tau}) \right]\right\}$$

$$R=(\det{D})^N\exp\left\{-(\phi_{i\tau}'-\phi_{i\tau})\left[ \frac{1}{2\Delta\tau V}(\phi_{i\tau}'+\phi_{i\tau})+\frac{1}{\Delta\tau^3 V\Omega^2}(\phi_{i\tau}'+\phi_{i\tau}-\phi_{i\tau+\Delta\tau}-\phi_{i\tau-\Delta\tau}) \right]\right\}$$

#### Local phonon on square lattice (e.g. CuO$_2$)

+ Cu-phonon: $$\phi_i$$
  + Holstein: $$\phi_i n_i$$
  + Frohlich: $$f_\delta\phi_i n_{i+\delta}$$, $$f_\delta$$ can be chosen as $$A_{1g}$$ or $$B_{1g}$$ (not important?)
+ O-phonon: $$\phi_{ix}$$ and $$\phi_{iy}$$
  + density buckling and breathing (Bulut and Scalapino): $$\phi_{i\delta}(n_i\pm n_{i+\delta})$$ 
  + hopping buckling: $$\phi_{i\delta}(c_i^\dagger c_{i+\delta} + h.c.)$$ 
    + $$A_{1g}$$-buckling: $$\phi_{ix}=\phi_{iy}$$ 
    + $$B_{1g}$$-buckling: $$\phi_{ix}=-\phi_{iy}$$ 

### el-ph coupling parameters

$$
H=\frac{M\Omega^2}{2}x^2+\frac{1}{2M}p^2+gxn
$$

Define $$x=\sqrt{\frac{1}{2M\Omega}}(b+b^\dagger)$$ and $$p=i\sqrt{\frac{M\Omega}{2}}(b-b^\dagger)$$, we get
$$
H=\Omega b^\dagger b + g\sqrt{\frac{1}{2M\Omega}}(b+b^\dagger)n
$$
Then, accordingly, the effective el-el interaction is
$$
V_{eff}(\omega)=\frac{g^2}{2M\Omega}\frac{2\Omega}{\omega^2-\Omega^2}=\frac{g^2}{M\Omega^2}\frac{\Omega^2}{\omega^2-\Omega^2}
$$

### benchmark

|                   case                   | 2-components             |                         1-component                          |                 ctqmc                  |  exact   | Johnston13' | Noack91' |   Vekic92'   |
| :--------------------------------------: | ------------------------ | :----------------------------------------------------------: | :------------------------------------: | :------: | :---------: | :------: | :----------: |
|  $1\times4,U=4,V=0.01,\beta=4,\Omega=1$  | $P_{ph}=0.259\pm0.002$   |                       $0.260\pm0.001$                        |                                        | $0.2593$ |   $0.26$    |          |              |
| $1\times4,U=4,V=2,\beta=4,\Omega=\infty$ | $S_{AF}=1.303\pm0.002$   |                                                              |                                        |          |             |          |              |
|        $1\times4,U=2,V=0,\beta=4$        | $S_{AF}=1.303\pm0.001$   |                       $1.303\pm0.001$                        |                                        |          |             |          |              |
| $1\times2,U=4,V=2,\beta=4,\Omega=1,t=0$  | $d=0.0088\pm0.0006$      |                $\langle s \rangle < 10^{-3}$                 |                                        |  0.0090  |             |          |              |
|   $4\times4,U=4,V=4,\Omega=1,\beta=4$    | $S_{AF}=1.022\pm0.005$   |                $\langle s \rangle < 10^{-3}$                 |                                        |          |   $1.37$    |          |              |
|  $4\times4,U=4,V=1.6,\Omega=1,\beta=1$   | $S_{AF}=0.9313\pm0.0002$ | $0.92\pm0.01$, $\langle s \rangle =(1.83\pm0.05)\times10^{-3}$ |                                        |          |    $0.94$     |          |              |
|  $4\times4,U=4,V=1.6,\Omega=1,\beta=2$   | $S_{AF}=1.299\pm0.001$ | $\langle s \rangle < 10^{-3}$ |                                        |          |    $1.34$     |          |              |
|   $1\times2,U=4,V=2,\Omega=1,\beta=1$    | $G_{12}=0.163\pm0.001$   |                                                              | $0.162\pm0.001,\langle s\rangle=0.27$  |          |             |          |              |
|   $1\times2,U=4,V=4,\Omega=1,\beta=1$    | $G_{12}=0.141\pm0.003$   |                                                              | $0.138\pm0.001,\langle s \rangle=0.09$ |          |             |          |              |
|   $8\times8,U=0,V=2,\Omega=1,\beta=12$   | $S_{CDW}=18.2\pm0.8$     |                                                              |                                        |          |             | $23\pm1$ | $25.3\pm0.3$ |
|   $4\times4,U=0,V=2,\Omega=1,\beta=12$   | $S_{CDW}=6.90\pm0.02$    |                        $6.92\pm0.01$                         |                                        |          |             |          | $7.6\pm0.3$  |
|   $4\times4,U=0,V=2,\Omega=1,\beta=8$    |                          |                                                              |                                        |          |             |          |              |
|   $4\times4,U=0,V=2,\Omega=1,\beta=6$    |                          |                    $S_{CDW}=5.06\pm0.05$                     |                                        |          |             |          | $6.5\pm0.5$  |
|   $4\times4,U=0,V=2,\Omega=1,\beta=4$    |                          |                    $S_{CDW}=2.47\pm0.02$                     |                                        |          |             |          | $2.67\pm0.2$ |
|  SSH $4\times4,V=0.5,\Omega=1,\beta=4$   |                          |                    $S_{AF}=0.792\pm0.003$                    | $0.789$ by Beyl(2018) (also with HQMC) |          |             |          |              |

+ The difference between our result and HQMC may be different choice of $\Delta\tau$. 
+ It seems when $\Delta\tau$ becomes small, the results becomes unreliable... 
+ How to fix it: check histogram? global update? reduce autocorrelation by larger ninterval or smaller nsp?
+ By setting nsp=1, nglobal=1, our result is still in disagreement with HQMC. 
+ Is HQMC exact? The authors didn't give any QMC parameters... 

####  More examples to test:
+ 1D Hubbard-Holstein model. In special, $1\times4$ with periodic boundary condition can be viewed as $2\times2$
+ 2-site Hubbard-Holstein model.
+ phonon only couple to 1-site?
+ choose different form?
+ auxiliary imaginary phonon field coupling to fermion spin? overcome the ergodicity problem?

#### Passed tests (parameter space $U,t,V,\Omega$):
+ $U=0$: HQMC(?),CTQMC
+ $t=0$: Lang-Firsov
+ $V=0$: free boson gas
+ $\Omega=\infty$: $U-V$-DQMC
+ $1\times2$ ($1\times4$ or $2\times2$?) in general: CTQMC
+ equivalence between real (spin channel) and imaginary (charge channel) HS decouplings of $U$

` So let's focus on pure Holstein model. Can global update help? Surely, it helps.

+ We can use higher-T auxiliary field as the low-T input, which may strongly reduce the warmup time.
+ The result shows: heat-bath ratio goes to equilibrium more slowly than Metropolis, while the latter may falss to local minima quickly. As a result, people may use corrected Metropolis algorithm instead:
$$
p=\frac{r}{1+\alpha r}, \quad\text{if}\,\, r<1
$$
while
$$
p=\frac{r}{\alpha+r},\quad\text{if}\,\, r>1
$$

### questions:
+ How $U$ affects $T_c$? 
+ Even for pure Holstein model, in fact, the $T_c$ value is still different in different papers, *e.g.* Costa(2017) and Weber(2017).
