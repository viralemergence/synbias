using Distributions
using Plots
using Plots.PlotMeasures
theme(:solarized_light)

x = range(0, 5, length = 100)
λa = 0.5
dist = Exponential(λa)
y = cdf(dist, x)
λb = 1
dist2 = Exponential(λb)
y2 = cdf(dist2, x)
λc = 1.5
dist3 = Exponential(λc)
y3 = cdf(dist3, x)
p1 = plot(x, [y y2 y3], xlabel = "", ylabel = "Number of Species Detected",
	labels = ["high facilitation" "neutral" "high competition"],
	title = "A. Species Accumulation Curves",
	xticks = false,
	yticks = false,
	left_margin = 5mm,
	top_margin = 5mm,
	bottom_margin = 5mm)
p2 = plot(x, [y y2 y3], xlabel = "Number of Sampling Units", ylabel = "Two-Way Species Interactions Detected",
	labels = ["high facilitation" "neutral" "high competition"],
	title = "B. Interaction Accumulation Curves",
	xticks = false,
	yticks = false,
	left_margin = 5mm,
	top_margin = 5mm,
	bottom_margin = 5mm)

p = plot(p1, p2, layout = (2, 1), size = (600, 800), aspect_ratio = :auto)
savefig("Figures/conceptual_figure.png")
