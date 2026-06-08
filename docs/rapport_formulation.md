# TP Diffusion-Transport-Réaction — Rapport

## Polytech Lyon — MAM 4A

---

## 1) Formulation faible du problème (1)

### Problème fort

Trouver $u$ telle que :

$$
\begin{cases}
-\mu u'' + b\,u' + \sigma\,u = f & \text{sur } (0,1) \\
u(0) = 0 \\
u(1) = 1
\end{cases}
$$

avec $\mu, b, \sigma > 0$ et $f \in L^2(0,1)$.

### Relèvement de la condition de Dirichlet

On pose $u = u_0 + u_D$ avec $u_D(x) = x$ (relèvement affine vérifiant $u_D(0)=0$, $u_D(1)=1$) et $u_0 \in H^1_0(0,1)$.

### Formulation variationnelle

Multiplier par $v \in H^1_0(0,1)$ et intégrer par parties le terme de diffusion :

$$
\int_0^1 \mu\,u_0'\,v'\,dx + \int_0^1 b\,u_0'\,v\,dx + \int_0^1 \sigma\,u_0\,v\,dx = \int_0^1 f\,v\,dx - \int_0^1 \mu\,u_D'\,v'\,dx - \int_0^1 b\,u_D'\,v\,dx - \int_0^1 \sigma\,u_D\,v\,dx
$$

En notant la forme bilinéaire :

$$
a(w, v) = \int_0^1 \mu\,w'\,v'\,dx + \int_0^1 b\,w'\,v\,dx + \int_0^1 \sigma\,w\,v\,dx
$$

et la forme linéaire :

$$
\ell(v) = \int_0^1 f\,v\,dx - a(u_D, v)
$$

La **formulation faible** s'écrit :

> Trouver $u_0 \in V = H^1_0(0,1)$ tel que $\forall v \in V$ :
> $$a(u_0, v) = \ell(v)$$
> puis $u = u_0 + u_D$.

**Remarque** : $a(\cdot,\cdot)$ est coercive sur $H^1_0(0,1)$ (grâce au terme $\mu \int u'v'$ et au lemme de Poincaré). L'existence et l'unicité sont assurées par le théorème de Lax-Milgram.

---

## 2) Discrétisation par éléments finis P1 de Lagrange — Identification de A et b

### Maillage

Partition uniforme de $(0,1)$ en $n$ sous-intervalles $[x_j, x_{j+1}]$ avec $h = 1/n$ et $x_j = jh$ pour $j = 0, \ldots, n$.

### Espace d'approximation

$V_h = \{v_h \in C^0([0,1]) : v_h|_{[x_j,x_{j+1}]} \in \mathbb{P}_1, \; v_h(0)=v_h(1)=0 \}$

Base : fonctions chapeau $\varphi_i$ ($i=1,\ldots,n-1$) telles que $\varphi_i(x_j) = \delta_{ij}$.

### Matrices élémentaires

Sur un élément $[x_j, x_{j+1}]$ de longueur $h$, les matrices élémentaires sont :

**Raideur (diffusion)** :
$$K^e = \frac{\mu}{h} \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$$

**Convection (transport)** :
$$C^e = \frac{b}{2} \begin{pmatrix} -1 & 1 \\ -1 & 1 \end{pmatrix}$$

**Masse (réaction)** :
$$M^e = \frac{\sigma h}{6} \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$$

**Vecteur charge élémentaire** (pour $f$ constant) :
$$F^e = \frac{fh}{2} \begin{pmatrix} 1 \\ 1 \end{pmatrix}$$

### Système global $Ax = b$ (équation (2))

Après assemblage et application des conditions de Dirichlet ($u_0 = 0$, $u_n = 1$) :

$$A = K + C + M \in \mathbb{R}^{(n-1)\times(n-1)}$$

où les nœuds intérieurs $i = 1, \ldots, n-1$ portent les inconnues, et :

$$A_{i,i} = \frac{2\mu}{h} + \frac{\sigma h}{3}, \qquad A_{i,i-1} = -\frac{\mu}{h} - \frac{b}{2} + \frac{\sigma h}{6}, \qquad A_{i,i+1} = -\frac{\mu}{h} + \frac{b}{2} + \frac{\sigma h}{6}$$

Le second membre $b_i$ incorpore la source $f$ et les conditions de Dirichlet :

$$b_i = \int_0^1 f\,\varphi_i\,dx + (\text{termes de bord issus de } u(1)=1)$$

**$A$ est une matrice tridiagonale.**

---

## 5) Cas transport dominant ($\sigma=0$, $b \gg \mu$) — Méthodes stabilisées

### Nombre de Péclet élémentaire

$$\text{Pe}_h = \frac{|b| h}{2\mu}$$

Lorsque $\text{Pe}_h > 1$, le Galerkin standard produit des **oscillations parasites**.

### Viscosité artificielle (upwind)

On ajoute une diffusion numérique $\mu_{\text{art}}$ :
- **Upwind plein** : $\mu_{\text{art}} = |b|h/2$ (stable mais trop diffusif)
- **Upwind optimal** : $\mu_{\text{art}} = \xi |b|h/2$ avec $\xi = \coth(\text{Pe}_h) - 1/\text{Pe}_h$

### SUPG (Streamline Upwind Petrov-Galerkin)

Test modifié $w_i = \varphi_i + \tau b \varphi_i'$ :

$$\sum_e \int_e \tau b \varphi_i' \left(b u_h' + \sigma u_h - f\right) dx = 0$$

Avec $\tau_{\text{opt}} = \frac{h}{2|b|}\left(\coth \text{Pe}_h - \frac{1}{\text{Pe}_h}\right)$.

### GALS (Galerkin Least Squares)

Test modifié $w_i = \varphi_i + \tau \mathcal{L}(\varphi_i)$ où $\mathcal{L}(w) = -\mu w'' + bw' + \sigma w$ :

Pour P1 ($\varphi_i'' = 0$) : $\mathcal{L}(\varphi_i) = b\varphi_i' + \sigma\varphi_i$.

Avec $\tau_j = h/(2b)$, GALS équivaut à l'upwind plein lorsque $\sigma = 0$.

---

## 6) Cas réaction dominant ($b=0$, $\sigma \gg \mu$)

La solution exacte présente des **couches limites** en $x=0$ et $x=1$ d'épaisseur $\sim \sqrt{\mu/\sigma}$.

Le critère de stabilité analogue est $\sigma h^2/(6\mu) < 1$ soit $h < \sqrt{6\mu/\sigma}$.

Pour $\mu = 1/6000$ et $\sigma = 1$ : $\alpha = \sqrt{6000} \approx 77.5$, couche limite d'épaisseur $\sim 1/77.5 \approx 0.013$.
