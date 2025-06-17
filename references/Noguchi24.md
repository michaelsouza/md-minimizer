# Fracture process of composite material in a spring network model 

Haruka Noguchi and Satoshi Yukawa<br>Department of Earth and Space Science,<br>Graduate School of Science, Osaka University,<br>Toyonaka, Osaka 560-0043, Japan

(Dated: May 2, 2024)

---

Abstract

We analyze a two-dimensional spring network model comprising breakable and unbreakable springs. Computer simulations showed this system to exhibit intermittent stress drops in a larger strain regime, and these stress drops resulted in ductile-like behavior. The scaling analysis reveals that the avalanche size distribution demonstrates a cut-off, depending on its internal structure. This study also investigates the relationship between cluster growth and stress drop, and we show that the amount of stress drop increases in terms of power law, corresponding to crack growth. The crack length distribution also demonstrates a cut-off depending on its internal structure. The results show that both the cluster growth-stress drop relationship and the crack size distribution are scaled by the quantity related to the internal structure, and the relevance of the exponent that scales the cluster growth-stress drop relationship to the exponent that scales crack size distribution is verified.

I. INTRODUCTION

Biological materials such as nacre or wood are composed of several materials and constitute a specific structure. It is known that their internal structure serves to improve their mechanical properties, especially fracture toughness. For example, the nacre is a composite material with organic parts and mineral parts arranged in a layered manner; It shows better fracture strength because of this laminated structure [1–3]. The trial to learn from the structure of biomaterial and exploit it to design a product is called “biomimicry” [4], and it is still an important field of engineering.

In material science, the effects of composition and structure on fracture behavior have attracted attention and have been used to create superior materials for failure. A typical man-made structural composite is fiber-reinforced ceramics, a composite of ceramics and fiber. Even though each component is brittle, fiber-reinforced ceramics show ductile-like behavior [5, 6] because fibers prevent crack propagation when cracks propagate in the ceramic phase and meet the fiber. More recently, advancing 3D printing technology has made it possible to create more complex internal structures of composite materials quickly, cheaply, and at a large scale [7], and this advancement has gathered much attention from material science to the relationship between structure and failure [8, 9]. For example, Li et al. [10] created a composite material with a structure consisting of glassy polymer skeletons filled with a highly rubbery thermoplastic elastomer using a 3D printer. By observing the differences in fracture behavior when altering the skeletal structure, they demonstrated that even in composite materials composed of the same type of material, fracture behavior, such as process zone formation, can be controlled by variations in the skeletal structure.

The influence of internal structure can appear not only in fracture behavior but also in scaling behavior. In the layered structures of soft and hard components like nacre, Okumura et al. predicted theoretically [11] and numerically [12, 13] that scaling law with the length of the period between soft and hard layers is valid for the crack tip stress and the crack shape. The other example is about hierarchical structure. Shi et al. [14] derived the scaling law of yield strength between different hierarchy levels and explained the difference in mechanical properties of the nano-scale hierarchical material in the degree of dealloying. Such "structure-based scaling relation" can be a guide to creating composite materials with more complex structures, but it is still not enough for our understanding of how structural properties like the length scale that characterize internal structure appear in fracture behavior.

To bridge this gap, we study the fracture behavior of composite materials, especially scaling behavior for the characteristic internal length scale, with a simple stochastic fracture model by a numerical simulation. There are several types of stochastic fracture models, such as the fiber bundle model [15-17], the random fuse model [18, 19], the spring-network (SN) model [20-23] and so on. They are used to understand disorder-induced statistical aspects of fracture like the power law of released-energy statistics [24-27], the self-affine nature of crack morphology [28-30], intermittent dynamics [26, 31], and pattern formation [32, 33]. The base of our model is the spring-network model. The model comprises two spring types: One spring breaks with the application of a specific amount of load and the other is unbreakable under any load, and these two kinds of springs comprise internal structure. As the previous study about the failure of composite material with the stochastic model, Kun et al. analyzed the fracture of a random mixture of weak and strong fiber composites by using equal-load sharing and local-load sharing fiber bundle models [34, 35]. Tauber et al. [36] considered a spring model that mimics polymer composites. Urabe et al. [23] and Rajesh et al. $[37,38]$ considered bi-material composite by spring model with two kinds of springs that have different young modulus. Compared to their models, we used strong springs in our SN model to form a regular matrix structure.

The contributions of this study in terms of determining the effect of the internal structure of composite materials are threefold. First, the present system demonstrates ductile fractures because of its internal structure. Second, the burst size distribution of the present model shows power-law behavior in the intermediate size scale, and it shows a cut-off in the case of the larger avalanche. The scaling analysis showed that the burst size distribution can be scaled according to the size scale determined by its internal structure. Finally, the crack size distribution is scaled by the internal-structure-based crack length, and the stress drop caused by crack growth is scaled by the crack-opening length depending on the internal structure of the material. Our model boasts simplicity in capturing the fundamental properties of the fracture process in the composite material with a matrix structure. Moreover, the model can be used as a prototype for composite materials with a more complex internal structure.

The organization of this paper is as follows. In Sec. II, the details of the model and simulation method are presented. Results are described in Sec. III, and conclusions are summarized in Sec. IV.

# II. MODEL AND SIMULATION 

In this study, we analyzed a two-dimensional SN model, which represents the composite material as a network of particles connected by Hookean springs. The system comprises $N \times N$ particles on a triangular lattice. We take the $x$-axis as parallel to the edge of the triangular lattice. All particles possess the same mass, which is taken to be a mass unit, and each pair of nearest-neighbor particles is connected by the Hookean spring, the natural length of which is represented by lattice spacing $l_{0}$, which is the unit of length. The periodic boundary condition was imposed in the $x$ direction and the fixed boundary condition was imposed in the $y$ direction for tensile loading. Each spring has a fracture threshold of $l^{*}$, which was randomly selected from a uniform distribution between 0 and 1 . When the strain of the spring, i.e., $\left|l-l_{0}\right| / l_{0}$, becomes larger than the threshold, $l^{*}$, the spring breaks, and it is removed from the system. In this system, the strains caused by the broken springs were distributed among the remaining live springs to reach the mechanical equilibrium. The successive breaking of many springs is possible, a phenomenon, which hereafter, is referred to as burst or avalanche. The potential energy of this model can be formulated as [23]
$$
V=\frac{k}{2} \sum_{\langle i, j\rangle}\left(\left|\boldsymbol{r}_{\boldsymbol{i}}-\boldsymbol{r}_{\boldsymbol{j}}\right|-l_{0}\right)^{2} g_{i j}
$$

where $\boldsymbol{r}_{\boldsymbol{i}}$ is the position vector of the $i$-th particle, $g_{i j}=1$ indicates a live bond, and $g_{i j}=0$ indicates a broken bond. The summing pair $\langle i, j\rangle$ runs the nearest-neighbor pairs on the triangular lattice. Parameter $k$ is a spring constant, taken as unity.

![img-0.jpeg](img-0.jpeg)

FIG. 1: Schematic of the spring network model. $N \times N$ particles are located on a triangular lattice, and nearest-neighbor particles are connected using springs. We take the $x$-axis as parallel to the edge of the triangular lattice. The dotted bonds correspond to breakable springs and the red bond corresponds to unbreakable springs. Unbreakable springs are regularly deployed spatially to constitute a $L_{\text {matrix }} \times L_{\text {matrix }}$ almost square frame. This figure shows the $L_{\text {matrix }}=2$ case. A periodic boundary condition is imposed in the $x$ direction, and a fixed boundary condition is imposed in the $y$ direction for tensile loading.


We replaced some springs with unbreakable springs, $l^{*}=\infty$, for modeling the composite material. These unbreakable springs were regularly deployed spatially to constitute a $L_{\text {matrix }} \times L_{\text {matrix }}$ almost square frame. Here, $L_{\text {matrix }}$ is the number of unbreakable springs for one side of the almost square matrix as shown in Fig. 1. For every $L_{\text {matrix }}$ layer, unbreakable springs are put parallel to the $x$-axis and zigzag for the $y$-direction. We call this system the "matrix-mixture system," and we term the system without unbreakable bonds as the "normal system." In this study, we take $N=96$ and $L_{\text {matrix }}=6,8,12$.

Next, the system was simulated under strain control as follows. The lowest row of particles was fixed, and a small amount of displacement was implemented among the highest row of particles. The system was then allowed to relax to a mechanical equilibrium state. The equilibrium state was explored using the FIRE algorithm [39]. After the system relaxed to the mechanical equilibrium state, we decided on which bonds to break. If a certain spring's strain was over the fracture threshold, that spring was removed from the system. After removing the springs, the mechanical equilibrium configuration was analyzed again without moving the top particles. This loop was continued until the springs stopped breaking in the mechanical equilibrium state. When the system reached this state, we repeated the same procedure. The simulation was finally stopped when the system completely broke into two pieces or the strain of the system reached 1. For one small uniaxial extension step, the strain of the system increased by $0.1 \%$, i.e., $\Delta \epsilon=0.0001$. Here, the statistically independent 1000 configurations were simulated.

# III. RESULT 

## A. Mechanical property

![img-1.jpeg](img-1.jpeg)

FIG. 2: Typical stress-strain response of the system. Every system shows a linear response in a small strain regime. Then, the normal system breaks in a brittle manner. The matrixmixture systems do not demonstrate such an abrupt change but show ductile-like failures with intermittent stress drops. The size distributions of stress drops vary in each system. As the matrix size of the system increases, the larger stress drop also increases, as shown in Fig. 3.

We first discuss the stress-strain curve of the system. In this study, we compute

$$
\Sigma=\frac{k}{H_{0}} \sum_{\langle i, j\rangle}\left(\left|\boldsymbol{r}_{i}-\boldsymbol{r}_{j}\right|-l_{0}\right) \frac{\left|y_{j}-y_{i}\right|}{\left|\boldsymbol{r}_{i}-\boldsymbol{r}_{j}\right|} g_{i j}
$$

as the stress [23], where $H_{0}=N l_{0}$ is the width of the system. Figure 2 shows the typical stress-strain curve of this system. At the beginning of tension application, all systems show elastic-like behavior, though cracks appear in the system. After that, the slope of the stress-strain curve reduces because of the increase in damage. Eventually, the normal and matrix-mixture systems show completely different mechanical responses. The normal system shows an abrupt stress drop. As shown in previous studies [40, 41], this stress drop is compatible with crack propagation from one end to another and breaking the system in two. That is, the normal system shows a brittle fracture. The matrix-mixture systems show similar behavior at the beginning of loading as the normal system. However, these systems behave as a ductile material in the larger strain regime. This ductile behavior could be attributed to the intermittent and instantaneous stress changes. Hereafter, we term this stress change as a stress drop, and its magnitude is denoted as $\Delta \Sigma$. The ductile regime frequently displays small-scale stress drops, which cancel out the increase in stress. The comparison of each $L_{\text {matrix }}$ shows that the system with a small $L_{\text {matrix }}$ value shows less stress drop than the system with a large $L_{\text {matrix }}$ value. To quantify this difference, we investigated the distribution of stress drop, $\Delta \Sigma$, which is denoted as $P(\Delta \Sigma)$, as shown in Fig. 3. All negative instantaneous changes in stress were considered. In the smaller stress drop regime, distribution $P(\Delta \Sigma)$ showed power-law decay, and its exponent is almost the same for the different matrix sizes, $L_{\text {matrix }}$. In the much larger $\Delta \Sigma$ region, the cut-off was observed to depend on $L_{\text {matrix }}$. The result in Fig. 3 suggests that a smaller matrix significantly suppresses large stress drops.

# B. Avalanche 

To understand the effect of suppressing the stress drop by using a matrix structure, we studied the difference between the burst size distribution $P(S)$ between the normal and matrix-mixture SN models, where $S$ is the number of breaking springs during the single loading step with respect to $\Delta \epsilon$. In the normal SN model, the burst size distribution, $P(S)$, behaves as $\sim S^{-\tau}$, and exponent $\tau=2.5$ [40]. We plot the burst size distribution in Fig. 4, wherein all burst events are considered. First, as shown, exponent $\tau$ decreases with the consideration of the internal structure in the matrix-mixture system. This indicates that smaller bursts are more likely to occur in the smaller matrix-mixture system. Second, Figure 4 shows that the burst size distribution demonstrates a cut-off size depending on its matrix size. Based on these observations, the matrix structure increases the burst events until the intermediate scale and suppresses burst events on a larger scale.

![img-2.jpeg](img-2.jpeg)

FIG. 3: Distribution of stress drop, $P(\Delta \Sigma)$. In the smaller $\Delta \Sigma$ regime, $P(\Delta \Sigma)$ shows power-law decay. The larger $\Delta \Sigma$ regime demonstrates a cut-off, depending on $L_{\text {matrix }}$. The distribution is considered by the binning with $1 / N$.

![img-3.jpeg](img-3.jpeg)

FIG. 4: Avalanche size distribution. The black dashed line corresponds to the distribution of the normal SN model, $P(S) \sim S^{-2.5}$.


Next, by using scaling analysis, we discuss the differences among each matrix-mixture system in terms of the burst size distribution. We assume that the burst size distribution has the following scaling form with exponents $a$ and $b$ :

$$
P\left(S, L_{\text {matrix }}\right)=S_{\text {matrix }} f\left(S / S_{\text {matrix }}^{a}\right)
$$

where $S_{\text {matrix }}$ is the number of breakable springs in the matrix calculated as $S_{\text {matrix }}=$ $3 L_{\text {matrix }}\left(L_{\text {matrix }}-1\right)$, and $f(x)$ is the scaling function [35]. The result of the scaling analysis showed that the distribution of burst size, $P(S)$, collapses onto the master curve. We achieved a favorable conformance between the data and master curve $f(x)$, formulated as $f(x)=x^{-\tau} \operatorname{erfc}(x-a)$ with $a=1.0, b=2.38$, and $\tau=2.35$. This result indicates that the burst event follows the power law with the same exponent, $\tau=2.35$, until the burst size was less than $S_{\text {matrix }}$ and showed sigmoidal decay once it increased more than $S_{\text {matrix }}$. The functional form of $f(x)$ indicates that any apparent change in $\tau$ and the decay behavior in Fig. 5 can be attributed to the difference in $S_{\text {matrix }}$. These results show that the stress drops were suppressed by the matrix structure because the fracture events were suppressed by the cut-off size, $S_{\text {matrix }}$.

![img-4.jpeg](img-4.jpeg)

FIG. 5: Scaled avalanche size distribution. The solid line corresponds to the fitting function, $f(x)=x^{-\tau} \operatorname{erfc}(x-1)$, with respect to the scaling hypothesis shown in Eq. (3). Here, we take scaling variables as $a=1.0, b=2.38$, and $\tau=2.35$.

# C. Crack coalescence and stress drop 

A stress drop depends on not only the number of avalanches but also the spatial distribution of crack formation. Figure 6 illustrates the effect of crack coalescence on stress drops. The figure clearly shows the difference between the spatial distributions of appeared cracks of the elastic and ductile regimes. We could clarify the effects of spatial distribution on the stress drop by observing the internal states of the springs. Figure 6 shows the configuration of the SN model at a certain strain rate. The colors represent the changes in the spring length compared with a previous state, i.e., red corresponds to extended springs and blue corresponds to shrunken springs. In addition, green represents the springs broken in the previous step. In both pictures, the number of breaking springs caused by the single tiny displacement is 10. In the elastic regime (Fig. 6(a)), each appeared crack was spatially isolated. As such, the rupture of springs does not significantly affect the entire structure. On the contrary, as the fracturing proceeds to the ductile regime (Fig. 6(b)), the breaking of springs tends to show a more significant effect on the internal structure by coalescing with existing cracks.

![img-5.jpeg](img-5.jpeg)

FIG. 6: Snapshots after breaking events in elastic and ductile regimes for $L_{\text{matrix}}=12$. Springs with reduced (increased) stress compared to those in the previous step are indicated in blue (red), respectively. Broken springs in the previous step are indicated in green. In both snapshots, the number of breaking springs in the previous step is the same, i.e., 10. The left panel shows a smaller stress drop, which is attributed to the burst events occurring at spatially isolated locations. In the right panel, the stress drop is approximately five times larger than in the left panel, and this could be attributed to the fracture events occurring in spatially close locations, accompanied by larger crack growth.


This result suggests that the amount of crack growth is essential for how stress is reduced in certain fracture events [42]. Thus, we quantitatively analyzed the relationship between crack growth and stress drop.

To quantify the effect of crack growth, we introduced two quantities, $C$ and $N_{C}(\epsilon)$, where $C$ corresponds to the size of the crack cluster and $N_{C}(\epsilon)$ is the number of crack clusters with size $C$ at strain $\epsilon$, as shown in Fig. 7. In this study, we assumed that the single crack cluster does not extend over the unbreakable bond. If the cluster seems to extend over the unbreakable bond, we have identified it as two separate clusters.

![img-6.jpeg](img-6.jpeg)

FIG. 7: Schematic of the definition of the cluster size of cracks. The black, green, and red lines correspond to breakable, broken, and unbreakable springs, respectively. The single crack cluster does not extend over the unbreakable bond. Thus, in this case, there are two clusters: $N_{C=3}=1$ and $N_{C=6}=1$.

![img-7.jpeg](img-7.jpeg)

FIG. 8: Relationship between the amount of crack growth, $\Delta C$, and the average amount of stress drop for crack growth, $\overline{\Delta \Sigma}$, in terms of matrix size. Inset: Bare relation. All matrix sizes in the smaller $\Delta C \lesssim 10^1$ regime display the same behavior, which deviates in the larger $\Delta C$ regime. Main: Scaled relation. This deviation can be scaled with respect to $\Delta C$ divided by $L_{matrix}$, where the exponent is $\zeta = 0.55$. For plotting, the horizontal axis, $\Delta C$, has been binned by an integer to reduce fluctuations.


Based on the percolation theory [43], the increment in crack cluster size, $\Delta C$, according to the amount of crack growth is defined as

$$
\Delta C=\sqrt{\sum_{C}\left(C^{2} N_{C}(\epsilon)-C^{2} N_{C}(\epsilon-\Delta \epsilon)\right)}
$$

This quantity is the square root of the difference between the average moments of crack cluster size before and after a tiny loading step in a sample. It is uniquely determined for a tiny loading step and one sample. The stress drop $\Delta \Sigma$ is also uniquely determined for a tiny loading step. Thus, the stress drop $\Delta \Sigma$ during a tiny loading step has one-to-one correspondence with $\Delta C$. The increment of the crack cluster size shows a large value for the large clusters of broken springs and a small value for several small clusters, even with the same number of broken springs. Figure 8 shows the relationship between $\Delta C$, binned by an integer, and the average stress drop $\overline{\Delta \Sigma}$ at each bin over the samples during the whole loading process $\epsilon=0$ to 1 . Its definition is as follows:

$$
\overline{\Delta \Sigma}=\frac{\sum_{\text {sample }} \sum_{\Delta C \in \text { bin }} \Delta \Sigma}{\sum_{\text {sample }} \sum_{\Delta C \in \text { bin }} 1}
$$

The amount of stress drops roughly increases in a power-law manner. The inset of Fig. 8 shows the behavior of the stress drop to deviate near $\Delta C \sim 10^{1}$. This deviation can be scaled with respect to $\Delta C$ over $L_{\text {matrix }}^{\zeta}$, as shown in Fig. 8. Exponent $\zeta \simeq 0.55$ results in the best fit in regime $\Delta C>30$.

Finally, we analyzed the crack size distribution, $P(C)$, for each matrix at the final state, $\epsilon=1$. Here $P(C)$ is defined as

$$
P(C)=\frac{\sum_{\text {sample }} N_{C}(\epsilon=1)}{\sum_{\text {sample }} \sum_{C} N_{C}(\epsilon=1)}
$$

In the inset of Fig. 9, the distribution shows the power-law decay in the small crack length regime. This decay is consistent with that observed in previous studies [44]. In this regime, the matrix structure does not affect the crack size. The effect of matrix structure becomes apparent in the region with large cracks. The distribution shows a cut-off corresponding to the matrix size for a long crack. This result clearly shows that the matrix structure of the system suppresses crack growth. By scaling the crack size by using $L_{\text {matrix }}^{1+\zeta^{\prime}}$, the point at which deviation from the power law begins and the cut-off size can be scaled independently on the matrix size is shown in Fig. 9. We achieved the best fit by using exponent $\zeta^{\prime}=0.56$, and this is almost the same as $\zeta$ in Fig. 8. We ascribe this agreement between $\zeta$ and $\zeta^{\prime}$ to corresponding $L_{\text{matrix}}^{\zeta}$ with an effective crack width. Previous research [45] showed that the roughness exponent of the random fuse model was $2 / 3$, and it was 0.62 for the Born model [46]; These values are close to the value achieved in the current study: $\zeta \simeq 0.55$. Considering $L_{\text{matrix}}^{\zeta^{\prime}}$ as the characteristic length of the crack width in the matrix, $L_{\text{matrix}}^{1+\zeta^{\prime}}$ can be identified according to the typical crack length in the matrix. The crack length distribution follows a power law similar to that of the normal SN model for small crack lengths. However, owing to the matrix structure, the cracks cannot grow larger, and a peak appears in the distribution after the scale of $L_{\text{matrix}}^{1+\zeta^{\prime}}$.

![img-8.jpeg](img-8.jpeg)

FIG. 9: Crack size distribution for each matrix size at the final state, $\epsilon=1$. Inset: Bare distribution. This distribution follows the power law in small-sized cracks. A cut-off size corresponding to its matrix size is observed in the large-sized cracks. Main: Scaled distribution. By scaling the crack size with respect to $L_{\text{matrix}}^{1+\zeta^{\prime}}$, we get the point at which deviation from the power law begins and the cut-off size can be scaled independently on the matrix size. Exponent $\zeta^{\prime}=0.56$ is almost the same as $\zeta$ in Fig. 8 The black dashed line corresponds to $C / L_{\text{matrix}}^{1+\zeta^{\prime}}=1$.


## IV. CONCLUSION

In summary of our results, we analyzed the fracture process of composite materials and their statistical properties according to the internal structure of the material by using the SN model with a mixture of breakable and unbreakable springs. We found that the proposed SN model shows ductile-like fractures because of intermittent stress drops. In addition, we revealed that avalanche size distribution is well scaled by the number of springs in matrix $L_{\text {matrix }}\left(3 L_{\text {matrix }}-1\right)$. The scaling function can be written as $x^{-r} \operatorname{erfc}(x-1)$, suggesting that the fracture event follows the power law, similar to the normal SN model [40], and it decays abruptly when reaching a specific number of springs in a matrix. We also revealed the relationship between crack cluster growth and stress drop. Larger clusters appeared when cracks merged, resulting in a more significant stress drop. On average, the amount of stress drop increased based on the power law, followed by the growth in crack clusters. The cluster size distribution showed a cut-off corresponding to the matrix size, and a typical crack length in the matrix could rescale the size cut-off. The matrix size limits the size of the cluster.

This study showed that the material's internal structure, the regular matrix structure, affects the fracture behavior under quasi-static tensile stress. In particular, the length scale of the internal structure significantly controls the fracture behavior by scaling. Interestingly, these properties appear in the present simple model. Additionally, the simplicity of the present model enables experimental verification of the present results using 3D printing techniques [7]. It is intriguing to consider how the matrix length scale appears in the other type of fracture. For example, how does the effect of the length scale of the internal structure appear in the statistical law of fragmentation [47, 48] or fatigue failure [49]? We can investigate these problems with an extension of the present model because of its simplicity of the present model.

# ACKNOWLEDGMENTS 

The authors thank T. Hatano and N. Sakumichi for the fruitful discussion and useful comments. This work was supported by JSPS KAKENHI Grant Number 19K03652.
[1] J. D. Currey, Proc. R. Soc. Lond. Ser. B. Biol. Sci. 196, 443 (1977).
[2] R. Wang and H. S. Gupta, Annu. Rev. Mater. Res. 41, 41 (2011).
[3] H. D. Espinosa, A. L. Juste, F. J. Latourte, O. Y. Loh, D. Gregoire, and P. D. Zavattieri, Nat. Commun. 2, 173 (2011).
[4] E. Lurie-Luke, Biotechnol. Adv. 32, 1494 (2014).
[5] K. S. Vecchio and F. Jiang, Mater. Sci. Eng., A 649, 407 (2016).
[6] M. Rühle and E. Evans, Prog. Mater. Sci. 33, 85 (1989).
[7] G. Liu, X. Zhang, X. Chen, Y. He, L. Cheng, M. Huo, J. Yin, F. Hao, S. Chen, P. Wang, S. Yi, L. Wan, Z. Mao, Z. Chen, X. Wang, Z. Cao, and J. Lu, Mater. Sci. Eng., R 145, 100596 (2021).
[8] O. Al-Ketan, R. K. A. Al-Rub, and R. Rowshan, Adv. Mater. Technol. 2, 1600235 (2017).
[9] L. S. Dimas, G. H. Bratzel, I. Eylon, and M. J. Buehler, Adv. Funct. Mater. 23, 4629 (2013).
[10] T. Li, Y. Chen, and L. Wang, Compos. Sci. Technol. 167, 251 (2018).
[11] K. Okumura and P.-G. De Gennes, Eur. Phys. J. E 4, 121 (2001).
[12] Y. Aoyanagi and K. Okumura, Phys. Rev. E 79, 066108 (2009).
[13] Y. Hamamoto and K. Okumura, Adv. Eng. Mater. 15, 522 (2013).
[14] S. Shi, Y. Li, B.-N. Ngo-Dinh, J. Markmann, and J. Weissmüller, Science 371, 1026 (2021).
[15] F. T. Peirce, J. Text. Inst. Trans. 17, 355 (1926).
[16] S. Pradhan, A. Hansen, and B. K. Chakrabarti, Rev. Mod. Phys. 82, 499 (2010).
[17] A. Hansen, P. C. Hemmer, and S. Pradhan, The Fiber Bundle Model: Modeling Failure in Materials (Wiley-VCH Verlag GmbH \& Co. KGaA, 2015).
[18] L. de Arcangelis, S. Redner, and H. J. Herrmann, J. Physique Lett. 46, 585 (1985).
[19] A. Shekhawat, S. Zapperi, and J. P. Sethna, Phys. Rev. Lett. 110, 185505 (2013).
[20] A. Hansen, S. Roux, and H. J. Herrmann, J. Phys. France 50, 733 (1989).
[21] S. Zapperi, P. K. V. V. Nukala, and S. Šimunović, Phys. Rev. E 71, 026106 (2005).
[22] S. Arbabi and M. Sahimi, Phys. Rev. B 47, 695 (1993).
[23] C. Urabe and S. Takesue, Phys. Rev. E 82, 016106 (2010).
[24] A. Garcimartin, A. Guarino, L. Bellon, and S. Ciliberto, Phys. Rev. Lett. 79, 3202 (1997).
[25] L. I. Salminen, A. I. Tolvanen, and M. J. Alava, Phys. Rev. Lett. 89, 185503 (2002).
[26] K. J. Måløy, S. Santucci, J. Schmittbuhl, and R. Toussaint, Phys. Rev. Lett. 96, 045501 (2006).
[27] L. I. Salminen, J. M. Pulakka, J. Rosti, M. J. Alava, and K. J. Niskanen, Europhys. Lett. 73, 55 (2005).
[28] S. Santucci, M. Grob, R. Toussaint, J. Schmittbuhl, A. Hansen, and K. Maløy, Europhys. Lett. 92, 44001 (2010).
[29] J. J. Mecholsky, D. E. Passoja, and K. S. Feinberg-Ringel, J. Am. Ceram. Soc. 72, 60 (1989).
[30] P. Daguier, B. Nghiem, E. Bouchaud, and F. Creuzet, Phys. Rev. Lett. 78, 1062 (1997).
[31] D. Bonamy, S. Santucci, and L. Ponson, Phys. Rev. Lett. 101, 045501 (2008).
[32] A. Groisman and E. Kaplan, Europhys. Lett. 25, 415 (1994).
[33] S. Kitsunezaki, Phys. Rev. E 60, 6449 (1999).
[34] R. C. Hidalgo, K. Kovács, I. Pagonabarraga, and F. Kun, Europhys. Lett. 81, 54005 (2008).
[35] K. Kovács, R. C. Hidalgo, I. Pagonabarraga, and F. Kun, Phys. Rev. E 87, 042816 (2013).
[36] J. Tauber, S. Dussi, and J. van der Gucht, Phys. Rev. Mater. 4, 063603 (2020).
[37] A. Mayya, A. Banerjee, and R. Rajesh, J. Mech. Behav. Biomed. 83, 108 (2018).
[38] S. Senapati, A. Banerjee, and R. Rajesh, Phys. Rev. E 107, 055002 (2023).
[39] E. Bitzek, P. Koskinen, F. Gähler, M. Moseler, and P. Gumbsch, Phys. Rev. Lett. 97, 170201 (2006).
[40] S. Zapperi, P. Ray, H. E. Stanley, and A. Vespignani, Phys. Rev. Lett. 78, 1408 (1997).
[41] P. K. V. V. Nukala, S. Zapperi, and S. Šimunović, Phys. Rev. E 71, 066106 (2005).
[42] G. Gagnon, J. Patton, and D. J. Lacks, Phys. Rev. E 64, 051508 (2001).
[43] D. Stauffer and A. Aharony, Introduction to Percolation Theory (Taylor \& Francis, 1992).
[44] C. Spyropoulos, C. H. Scholz, and B. E. Shaw, Phys. Rev. E 65, 056105 (2002).
[45] E. T. Seppälä, V. I. Räisänen, and M. J. Alava, Phys. Rev. E 61, 6312 (2000).
[46] G. Caldarelli, R. Cafiero, and A. Gabrielli, Phys. Rev. E 57, 3878 (1998).
[47] F. Kun, F. K. Wittel, H. J. Herrmann, B. H. Kröplin, and K. J. Måløy, Phys. Rev. Lett. 96, 025504 (2006).
[48] H. A. Carmona, F. K. Wittel, F. Kun, and H. J. Herrmann, Phys. Rev. E 77, 051302 (2008).
[49] F. Kun, H. A. Carmona, J. S. Andrade, and H. J. Herrmann, Phys. Rev. Lett. 100, 094301 (2008).