using ModelingToolkit
using DifferentialEquations: solve
using CairoMakie
using Latexify

# Configuration generale des figures
CairoMakie.activate!(; px_per_unit=2)

### Modele phyto remediation ###

# Parametres de la simulation
time_max = 500.0

# Definition des variables d'interet
@variables t C(t) C_P(t) C_H(t) P(t) H(t)

# Parametres du modele
@parameters θ=0.0 μ=0.0 a=0.2 ϵ=0.6 b=1e-1 γ_H=5e-2 m_P=0.08 m_H=0.08 r=0.2 K=10

# Compartiments au temps initial
u0 = [1, 0, 0, 1, 1]

# Creation d'un operateur de derivee
D = Differential(t)

# Processus modelises
contamination = θ - μ*C
absorption_plante = b*P*C
croissance_plante = r*P*(1 - P/K)
predation = a*P*H
conversion = ϵ*predation
mortalite_plante = m_P*P
mortalite_herbivore = m_H*H
mortalite_herbivore_cont = γ_H * (C_H/H) * H

# Definition du systeme d'equations differentielles
sir_equations = [
    D(C) ~ contamination - absorption_plante + mortalite_plante * (C_P/P) + mortalite_herbivore_cont * (C_H/H) + mortalite_herbivore * (C_H/H),
    D(C_P) ~ absorption_plante - mortalite_plante * (C_P/P) - predation * (C_P/P),
    D(C_H) ~ predation *(C_P/P) - mortalite_herbivore * (C_H/H) - mortalite_herbivore_cont * (C_H/H),
    D(P) ~ croissance_plante - predation - mortalite_plante,
    D(H) ~ conversion - mortalite_herbivore - mortalite_herbivore_cont
]

# Afficher le systeme d'equations differentielles en LaTeX
Latexify.latexify(sir_equations) |> render

# Construire une representation formelle du systeme ODE
@named sir_system = ODESystem(sir_equations)

# Creer un probleme ODE a resoudre numeriquement
# (systeme, valeurs initiales, temps)
prob = ODEProblem(sir_system, u0, (0.0, time_max))

# Resoudre le systeme ODE
sol = solve(prob, saveat=0:1:time_max)

## Afficher les resultats

# Creation d'une figure
fig = Figure(; resolution=(1000,500))

# Figure 1: axes
timecourse_population = Axis(fig[1,1]; xlabel="Time", ylabel="Population")
timecourse_contaminant = Axis(fig[1,2]; xlabel="Time", ylabel="Contaminant")
phaseportrait = Axis(fig[1,3]; xlabel="Plantes", ylabel="Herbivores")

# Figure 1a: courbes
lines!(timecourse_population, sol[t], sol[P], label="Plantes", color=:green)
lines!(timecourse_population, sol[t], sol[H], label="Herbivores", color=:red)

# Figure 1b : courbes
lines!(timecourse_contaminant, sol[t], sol[C], label="Environnement", color=:blue)
lines!(timecourse_contaminant, sol[t], sol[C_P], label="Plantes", color=:green)
lines!(timecourse_contaminant, sol[t], sol[C_H], label="Herbivores", color=:red)

# Figure 1c : courbes
lines!(phaseportrait, sol[P], sol[H], color=:black)

# Ajouter une legende
Legend(fig[2, 1], timecourse_population; tellheight=true, tellwidth=false, framevisible=false, orientation=:horizontal)
Legend(fig[2, 2], timecourse_contaminant; tellheight=true, tellwidth=false, framevisible=false, orientation=:horizontal)

# Definir les limites des axes
ylims!(timecourse_population; low=0.0, high=maximum(sol[P]+sol[H]))
ylims!(timecourse_contaminant; low=0.0, high=maximum(sol[C]+sol[C_P]+sol[C_H]))
xlims!(phaseportrait; low=0.0, high=maximum(sol[P]))
ylims!(phaseportrait; low=0.0, high=maximum(sol[H]))

# Embellir
tightlimits!(timecourse_population)
tightlimits!(timecourse_contaminant)
tightlimits!(phaseportrait)

# Afficher la figure
current_figure()
save("phyto.png", current_figure())

#Hello
