%% Run all TP exercises
%  Execute from the src/ directory

% Setup for Octave off-screen rendering (use fltk + xvfb-run)
if exist('OCTAVE_VERSION','builtin')
  graphics_toolkit('fltk');
  warning('off','all');
end

fprintf('====== TP Diffusion-Transport-Reaction ======\n\n');

fprintf('--- Q4a: Diffusion-Transport ---\n');
q4a_diffusion_transport;

fprintf('\n--- Q4b: Diffusion-Reaction ---\n');
q4b_diffusion_reaction;

fprintf('\n--- Q5a: Transport dominant - Galerkin ---\n');
q5a_transport_dominant_galerkin;

fprintf('\n--- Q5b: Upwind et SUPG ---\n');
q5b_upwind_supg;

fprintf('\n--- Q5c: GALS vs Upwind ---\n');
q5c_gals_vs_upwind;

fprintf('\n--- Q6a: Reaction dominant ---\n');
q6a_reaction_dominant;

fprintf('\n--- Q6b: Reaction dominant - methodes stabilisees ---\n');
q6b_reaction_stabilized;

fprintf('\n--- Q5d: Convergence (transport) ---\n');
q5d_convergence;

fprintf('\n--- Q6c: Convergence (reaction) ---\n');
q6c_convergence;

fprintf('\n====== Toutes les figures sont dans ../figures/ ======\n');
