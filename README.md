# TP Diffusion-Transport-Réaction

**Polytech Lyon — MAM 4A**

Résolution par éléments finis P1 de Lagrange du problème 1D :

```
-μ u'' + b u' + σ u = f   sur (0,1)
u(0) = 0,  u(1) = 1
```

## Structure du projet

```
src/
  ├── feasmbl2.m          # Assemblage matrice/vecteur système (fourni)
  ├── feeldof1.m           # Connectivité dofs par élément 1D (fourni)
  ├── fef1l.m              # Vecteur charge pour f=1, élément linéaire (fourni)
  ├── feode2l.m            # Matrice élémentaire Galerkin (fourni)
  ├── feaplyc2.m           # Application conditions de Dirichlet (fourni)
  ├── Diff_Trans_Reac1D.m  # Script pilote d'origine fourni (template, verbatim)
  ├── fef_load.m           # Vecteur charge pour f(x) quelconque (Gauss 2 pts)
  ├── compute_error.m      # Erreurs L2, H1-semi, max nodale (Gauss 3 pts)
  ├── solve_dtr.m          # Solveur EF flexible (Galerkin/lumped/upwind/SUPG/GALS)
  ├── q4a_diffusion_transport.m   # Q4a : μ=b=f=1, σ=0
  ├── q4b_diffusion_reaction.m    # Q4b : μ=σ=1, b=0, f(x)=x
  ├── q5a_transport_dominant_galerkin.m  # Q5a : Galerkin, n=10,40,60
  ├── q5b_upwind_supg.m           # Q5b : upwind + SUPG
  ├── q5c_gals_vs_upwind.m        # Q5c : GALS vs upwind
  ├── q5d_convergence.m           # Q5d : tableau/plot de convergence (transport)
  ├── q6a_reaction_dominant.m     # Q6a : réaction dominante
  ├── q6b_reaction_stabilized.m   # Q6b : Galerkin/lumping/GALS/SUPG (réaction)
  ├── q6c_convergence.m           # Q6c : tableau/plot de convergence (réaction)
  └── run_all.m                    # Exécuter toutes les questions
figures/                            # Figures générées (.png)
docs/
  └── rapport_formulation.md       # Formulation faible + identification A, b
```

## Exécution

Depuis MATLAB ou GNU Octave :

```matlab
cd src
run_all
```

Les figures sont sauvegardées dans `figures/`.

## Méthodes implémentées

| Méthode | Description |
|---------|-------------|
| `galerkin` | Galerkin standard P1 |
| `galerkin_lumped` | Galerkin + mass-lumping (stabilisation réaction) |
| `upwind_full` | Viscosité artificielle pleine μ_art = \|b\|h/2 |
| `upwind_opt` | Viscosité artificielle optimale (Scharfetter-Gummel) |
| `supg` | Streamline Upwind Petrov-Galerkin |
| `gals` | Galerkin Least Squares |
