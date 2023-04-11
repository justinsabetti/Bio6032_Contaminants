####################################################################
# Modèle de contamination du sol
####################################################################

##################################
# Configuration
##################################
using ModelingToolkit
using DifferentialEquations: solve
using CairoMakie
using Latexify
import JLD
import CSV

# Configuration des figures
CairoMakie.activate!(; px_per_unit=2)

##################################
# Définition
##################################

# Définition des compartiments du modèle
# C(t): quantité de contaminant dans le sol (au cours du temps)
# C_P(t): quantité de contaminant dans la biomasse des plantes (au cours du temps)
# C_H(t): quantité de contaminant dans la biomasse des herbivores (au cours du temps)
# P(t): biomasse des plantes (au cours du temps)
# H(t): biomasse des herbivores (au cours du temps)
@variables t C(t) C_P(t) C_H(t) P(t) H(t)

# Paramètres du modèle
## θ: taux d'apport du contaminant [θ]=g/s
## μ: taux d'export du contaminant [μ]=(g de contaminant)/s/(g de plante)
## a: taux d'attaque [a]=
## ϵ: taux d'efficacité de la conversion de biomasse [ϵ]=
## b: taux d'absorption de contaminant par les plantes [b]=
## γ_H: taux de mortalité des herbivore due au contaminant [γ_H]=
## m_P: taux de mortalité naturelle des plantes [m_P]=
## m_H: taux de mortalité naturelle des herbivores [m_H]=
## r: taux de croissance des plantes [r]=
## K: capacité de support pour les plantes [K]=
@parameters θ=0.0 μ=0.0 a=0.2 ϵ=0.6 b=1e-1 γ_H=5e-2 m_P=0.08 m_H=0.08 r=0.2 K=10

# Paramètres de la simulation

## Compartiments au temps initial: [C(0), C_P(0), C_H(0), P(0), H(0)]
sim_init_conditions = [1, 0, 0, 1, 1]

## Etendue temporelle
sim_time_max = 500.0

# Création d'un opérateur de derivée
D = Differential(t)

# Processus modélisés

## Processus de contamination du sol
## Apport - Export
contamination = θ - μ*C
## Absorption du contaminant par les plantes
## Choix de modélisation: L'absorption dépend du métabolisme des plantes
absorption_plante = b*P*C
## Croissance de la biomasse des plantes
## Choix de modélisation: croissance logistique
croissance_plante = r*P*(1 - P/K)
## Prédation des plantes par les herbivores
## Choix de modélisation: action de masse
predation = a*P*H
## Conversion de la biomasse de plantes ingérée en biomasse d'herbivores
conversion = ϵ*predation
## Mortalité naturelle des plantes
mortalite_plante = m_P*P
## Mortalité naturelle des herbivores
mortalite_herbivore = m_H*H
## Mortalité des herbivore due a la présence de contaminant dans leur biomasse
mortalite_herbivore_cont = γ_H * (C_H/H) * H

# Système d'équations differentielles
model_equations = [
    # Variation du contaminant dans le sol
    D(C) ~ contamination - absorption_plante + mortalite_plante * (C_P/P) + mortalite_herbivore_cont * (C_H/H) + mortalite_herbivore * (C_H/H),
    # Variation du contaminant dans la biomasse des plantes
    D(C_P) ~ absorption_plante - mortalite_plante * (C_P/P) - predation * (C_P/P),
    # Variation du contaminant dans la biomasse des herbivores
    D(C_H) ~ predation *(C_P/P) - mortalite_herbivore * (C_H/H) - mortalite_herbivore_cont * (C_H/H),
    # Variation de la biomasse des plantes
    D(P) ~ croissance_plante - predation - mortalite_plante,
    # Variation de la biomasse des herbivores
    D(H) ~ conversion - mortalite_herbivore - mortalite_herbivore_cont
]

##################################
# Resolution
##################################

# Construire une représentation formelle du système ODE
@named model_system = ODESystem(model_equations)

# Créer un problème ODE a résoudre numériquement
# (système, valeurs initiales, temps)
prob = ODEProblem(model_system, sim_init_conditions, (0.0, sim_time_max))

# Résoudre le système ODE
sol = solve(prob, saveat=0:1:sim_time_max)

##################################
# Affichage des résultats
##################################

# Création d'une figure
fig = Figure(; resolution=(1000,500))

# Figure 1: axes
timecourse_population = Axis(fig[1,1]; xlabel="Time", ylabel="Population")
timecourse_contaminant = Axis(fig[1,2]; xlabel="Time", ylabel="Contaminant")
phaseportrait = Axis(fig[1,3]; xlabel="Plants", ylabel="Herbivores")

# Figure 1a: Évolution dans le temps
# de la biomasse de Plantes et d'Herbivores
lines!(timecourse_population, sol[t], sol[P], label="Plants", color=:green)
lines!(timecourse_population, sol[t], sol[H], label="Herbivores", color=:red)

# Figure 1b : Évolution dans le temps 
# de la quantité de contaminant 
# dans l'environnement, les plantes et les herbivores 
lines!(timecourse_contaminant, sol[t], sol[C], label="Environment", color=:blue)
lines!(timecourse_contaminant, sol[t], sol[C_P], label="Plants", color=:green)
lines!(timecourse_contaminant, sol[t], sol[C_H], label="Herbivores", color=:red)

# Figure 1c : 
lines!(phaseportrait, sol[P], sol[H], color=:black)

# Ajouter une légende
Legend(fig[2, 1], timecourse_population; tellheight=true, tellwidth=false, framevisible=false, orientation=:horizontal)
Legend(fig[2, 2], timecourse_contaminant; tellheight=true, tellwidth=false, framevisible=false, orientation=:horizontal)

# Définir les limites des axes
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

##################################
# Sauvegarde des résultats
##################################

output_folder="output"

# Création du dossier de sortie des résultats
if ~isdir(output_folder)
    mkdir(output_folder)
end

# Identifiant de la simulation courante
sim_identifier = "$(hash(model_equations))"

# Emplacement de la sauvegarde
 result_path = joinpath(output_folder, sim_identifier)

# Création du dossier au besoin
if !ispath(result_path)
    mkpath(result_path)
end


# Sauvegarde du modèle
model_path = joinpath(result_path, "model.png")
Latexify.latexify(model_equations) |> render
CairoMakie.save(model_path, current_figure())

# Sauvegarde des paramètres
simulation_parameters = (theta=θ, mu=μ, a=a, epsilon=ϵ, b=b, gamma_H=γ_H, m_P=m_P, m_H=m_H, r=r, K=K)
param_path = joinpath(result_path, "parameters.csv")
CSV.write(param_path, simulation_parameters)

# Sauvegarde des resultats
fig_path = joinpath(result_path, "results_overview.png")
CairoMakie.save(fig_path, current_figure())
