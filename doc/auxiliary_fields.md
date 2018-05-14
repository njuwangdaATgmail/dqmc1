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