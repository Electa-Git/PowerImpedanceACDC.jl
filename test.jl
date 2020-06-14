using PGFPlotsX
x = range(0; stop=2, length = 100)
exp_plot = PlotInc(Table(x, exp.(x)))
exp_legend = LegendEntry(raw"$\exp(x)$")
log_plot = PlotInc(Table(x, log.(x)))
log_legend = LegendEntry(raw"$\log(x)$")

axs1 = @pgf Axis(exp_plot, exp_legend, log_plot, log_legend)
axs2 = @pgf SemiLogYAxis(exp_plot, exp_legend, log_plot, log_legend)
axs3 = @pgf SemiLogXAxis(exp_plot, exp_legend, log_plot, log_legend)
axs4 = @pgf LogLogAxis(exp_plot, exp_legend, log_plot, log_legend)

@pgf GroupPlot(
    { group_style = { group_size="2 by 2" },
      no_markers,
      legend_pos="north west",
      xlabel=raw"$x$",
    },
    axs1, axs2, axs3, axs4)
