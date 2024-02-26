import polars as pl
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

def plot_and_whisker(
        incoming_data,
        comparison_column,
        response_column,
        ylabel,
        xticklabels,
        custom_colours = None,
        subsetting_column = None,
        subsetting_value = None,
        print_stats_table = True,
        facet = False,
        facet_on = None
    ):
    """
    Create a violin-and-whisker plot with significance bars.
    """

    sns.set_theme(style = "ticks")

    def _annotate(**kwargs):
        data = pl.DataFrame(kwargs.pop('data'))
        # Assume the order of groups is the same as labels on x axis
        df_list = xticklabels
        # Stat test for difference
        # Default uses mannwhitneyu
        # Initialise a list of combinations of groups
        significant_combinations = []
        # Check from the outside pairs of boxes inwards
        ls = list(range(1, len(df_list) + 1))
        combinations = [
            (ls[x], ls[x + y]) for y in reversed(ls) for x in range((len(ls) - y))
        ]
        for combination in combinations:
            data1 = data.filter(
                pl.col(comparison_column) == df_list[combination[0] - 1]
            ).to_series()
            data2 = data.filter(
                pl.col(comparison_column) == df_list[combination[1] - 1]
            ).to_series()
            # Significance
            U, p = stats.mannwhitneyu(data1, data2, alternative='two-sided')
            if p < 0.05:
                significant_combinations.append([combination, p])

        # Add significance labels
        # Get the y-axis limits
        bottom, top = 0, 1.1
        y_range = top - bottom

        # Significance bars
        for i, significant_combination in enumerate(significant_combinations):
            # Columns corresponding to the datasets of interest
            x1 = significant_combination[0][0] - 1
            x2 = significant_combination[0][1] - 1
            # What level is this bar among the bars above the plot?
            level = len(significant_combinations) - i
            # Plot the bar
            bar_height = (y_range * 0.07 * level) + top
            bar_tips = bar_height - (y_range * 0.02)
            plt.plot(
                [x1, x1, x2, x2],
                [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k'
            )
            # Significance level
            p = significant_combination[1]
            if p < 0.001:
                sig_symbol = '***'
            elif p < 0.01:
                sig_symbol = '**'
            elif p < 0.05:
                sig_symbol = '*'
            text_height = bar_height + (y_range * 0.01)
            plt.text(
                (x1 + x2) * 0.5,
                text_height,
                sig_symbol,
                ha='center',
                va='bottom',
                c='k'
            )


    # Subset data if needed
    if subsetting_column is not None:
        plot_title = subsetting_value
        data = incoming_data.filter(
            pl.col(subsetting_column) == subsetting_value
            ).select(
                pl.col(response_column)
            ).insert_column(
                1,
                incoming_data.filter(
                    pl.col(subsetting_column) == subsetting_value
                ).select(
                    pl.col(comparison_column)
                ).to_series()
            )
    else:
        data = incoming_data.select([response_column, comparison_column])

    if facet:
        if subsetting_column is not None:
            data = data.insert_column(
                2,
                incoming_data.filter(
                    pl.col(subsetting_column) == subsetting_value
                ).select(
                    pl.col(facet_on)
                ).to_series()
            )
        else:
            data = data.insert_column(
                2,
                incoming_data.select(
                    pl.col(facet_on)
                ).to_series()
            )

        # Generate plot
        ax = sns.FacetGrid(data.to_pandas(), col=facet_on, margin_titles=True)
        ax.map_dataframe(
            sns.violinplot,
            x = comparison_column,
            y = response_column,
            hue = comparison_column,
            palette = custom_colours if custom_colours else "vlag",
            inner = "quart",
            fill = False
        )
        ax.map_dataframe(
            sns.stripplot,
            x = comparison_column,
            y = response_column,
            hue = comparison_column,
            size = 4,
            color = ".3",
            alpha = 0.2
        )
        ax.set(ylabel = ylabel)
        ax.fig.subplots_adjust(top=0.75)
        ax.fig.suptitle(
            subsetting_value if subsetting_value else comparison_column
        )
        sns.despine(trim = True)
        ax.set_axis_labels("")
        ax.map_dataframe(_annotate)
        if len(xticklabels) > 4:
            ax.set_xticklabels(rotation = 90)
        plt.show()

    else:
        # Assume the order of groups is the same as labels on x axis
        df_list = xticklabels
        # Stat test for difference
        # Default uses mannwhitneyu
        # Initialise a list of combinations of groups
        significant_combinations = []
        # Check from the outside pairs of boxes inwards
        ls = list(range(1, len(df_list) + 1))
        combinations = [
            (ls[x], ls[x + y]) for y in reversed(ls) for x in range((len(ls) - y))
        ]
        for combination in combinations:
            data1 = data.filter(
                pl.col(comparison_column) == df_list[combination[0] - 1]
            ).to_series()
            data2 = data.filter(
                pl.col(comparison_column) == df_list[combination[1] - 1]
            ).to_series()
            # Significance
            U, p = stats.mannwhitneyu(data1, data2, alternative='two-sided')
            if p < 0.05:
                significant_combinations.append([combination, p])

        # Generate plot
        f, ax = plt.subplots(figsize = (7, 6))
        sns.violinplot(
            data.to_pandas(),
            x = comparison_column,
            y = response_column,
            hue = comparison_column,
            inner = "quart",
            fill = False,
            palette = custom_colours if custom_colours else "vlag"
        )
        sns.stripplot(
            data.to_pandas(),
            x = comparison_column,
            y = response_column,
            hue = comparison_column,
            size = 4,
            color = ".3",
            alpha = 0.2
        )
        ax.set(ylabel = ylabel)
        plt.title(
            label = plot_title if plot_title else comparison_column,
            loc = 'left',
            fontweight = 'bold'
        )
        sns.despine(trim = True)

        # Add significance labels
        # Get the y-axis limits
        bottom, top = ax.get_ylim()
        y_range = top - bottom

        # Significance bars
        for i, significant_combination in enumerate(significant_combinations):
            # Columns corresponding to the datasets of interest
            x1 = significant_combination[0][0] - 1
            x2 = significant_combination[0][1] - 1
            # What level is this bar among the bars above the plot?
            level = len(significant_combinations) - i
            # Plot the bar
            bar_height = (y_range * 0.07 * level) + top
            bar_tips = bar_height - (y_range * 0.02)
            plt.plot(
                [x1, x1, x2, x2],
                [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k'
            )
            # Significance level
            p = significant_combination[1]
            if p < 0.001:
                sig_symbol = '***'
            elif p < 0.01:
                sig_symbol = '**'
            elif p < 0.05:
                sig_symbol = '*'
            text_height = bar_height + (y_range * 0.01)
            plt.text((x1 + x2) * 0.5, text_height, sig_symbol, ha='center', va='bottom', c='k')
        
        if len(xticklabels) > 4:
            ax.xticks(rotation = 90)
        plt.show()

    if print_stats_table:
        first = [item[0][0] - 1 for item in significant_combinations]
        last = [item[0][1] - 1 for item in significant_combinations]
        stats_table = pl.DataFrame(
                significant_combinations,
                schema = ['index', "p.value"]
            ).select(
                'p.value'
            ).insert_column(
                0,
                pl.DataFrame(
                    [df_list[i] for i in first],
                    schema = ["group_1"]
                ).to_series()
            ).insert_column(
                1,
                pl.DataFrame(
                    [df_list[i] for i in last],
                    schema = ["group_2"]
                ).to_series()
            ).with_columns(
                signature = pl.lit(plot_title)
            )
        print(stats_table)


def plot_stacked(
        data,
        sample_id,
        custom_colours,
        sample_column='Tumor_Sample_Barcode',
        comparison_column='method',
        response_column='exposure',
        color_on='signature',
        ylabel='Relative SBS exposure'
):
    small_maf = data.filter(
        pl.col(sample_column) == sample_id
    )

    # Generate plot
    sns.set_theme(style = "ticks")
    ax = sns.histplot(
        small_maf.filter(pl.col(response_column) > 0).to_pandas(),
        x=comparison_column,
        hue=color_on,
        weights=response_column,
        multiple='stack',
        palette=custom_colours
    )
    sns.move_legend(
        ax,
        "lower left",
        bbox_to_anchor=(-0.1, -0.5),
        ncol=5,
        title=None,
        frameon=False
    )
    ax.set(ylabel = ylabel)
    plt.title(
            label = sample_id,
            loc = 'left',
            fontweight = 'bold'
        )
    sns.despine(trim=True)
    plt.show()
