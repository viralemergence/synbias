using Distributions
using PlotlyJS

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

p1 = plot(
	[
		scatter(x = x, y = y, name = "high facilitation", yaxis = "y", line = attr(color = "#004488FF")),
		scatter(x = x, y = y2, name = "neutral", yaxis = "y", line = attr(color = "#DDAA33FF")),
		scatter(x = x, y = y3, name = "high competition", yaxis = "y2", line = attr(color = "#BB5566FF")),
	],
	Layout(
		title = "Species & Interaction Accumulation Curves",
		xaxis = attr(title = "Number of Sampling Units", showticklabels = false, showgrid = false),
		yaxis = attr(title = "Number of Species Detected", showticklabels = false, showgrid = false),
		yaxis2 = attr(
			title = "Number of Interactions Detected",
			overlaying = "y",
			side = "right",
			showticklabels = false,
			showgrid = false,
		),
	))

savefig(p1, "Figures/conceptual_figure.png")
