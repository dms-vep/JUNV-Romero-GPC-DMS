---
aside: false
---

# Effects of mutations on pseudovirus entry in Vesper mice TfR1 expressing 293T cells 

## Interactive plot of mutation effects
The plot below shows how mutations affect cell entry, and allows you to zoom and mouseover sites and mutations. Positive values indicate better entry than the parental GPC while negative values indicate impaired entry. 

Here is an explanation of the key plot elements:
 - The zoom bar at the top of the plot shows different regions of GPC, and can be used to zoom in on specific sites.
 - The line plot summarizes the effects of mutations at each site, using the site smmary statistic specified by the interactive option at the bottom of the plot (eg, sum or mean effect of mutations at a site).
  - The heatmap shows how each individual mutation affects cell entry. The `x`'s indicate the amino-acid identity in the parental GPC strain, light gray indicates mutations that were not measured, and dark gray indicates mutations that are filtered out by one of the interactive sliders at the bottom of the plot. If you mouseover points, you can also see the measurement in each individual library, which provides a sense of the reproducibility between experiments. 
  - The options at the bottom of the plot allow you to adjust what is displayed. For instance, the *times_seen* slider filters for only mutations seen in at least that many unique variants per-library; larger values will generally give more accurate measurements. Similarly, you can filter by how many different experimental selections the mutation was measured in.

## Effects of mutations on entry in *C. musculinus* TfR1 expressing cells

Click on the expansion box in the upper right of the plot to enlarge it for easier viewing.

<Figure caption="Interactive plot showing effects of mutations on entry into Vesper mice TfR1 expressing 293T cells">
    <Altair :showShadow="true" :spec-url="'htmls/293T-TFRC-KO-C-musculinus-TfR1_func_effects.html'"></Altair>
</Figure>

## Effects of mutations on entry in *C. callosus* TfR1 expressing cells

Click on the expansion box in the upper right of the plot to enlarge it for easier viewing.

<Figure caption="Interactive plot showing effects of mutations on entry into Vesper mice TfR1 expressing 293T cells">
    <Altair :showShadow="true" :spec-url="'htmls/293T-TFRC-KO-C-callosus-TfR1_func_effects.html'"></Altair>
</Figure>

## Comparison of effects of GPC mutations on entry in cells expressing *C. musculinus* versus *C. callosus* TfR1

The plot below compares the effects of GPC mutations on entry into cells expressing *C. musculinus* versus *C. callosus* TfR1.

Here is an explanation of the key plot elements:
 - The line plot shows the difference in effects of mutations at each site between a *comparator* and *reference* cell line (use the dropdowns below the plot to choose which cell line is the *comparator* and *reference*). If you click on sites on the line plot, then the heatmaps below the line plot will zoom to show that site in the center, and you can mouseover the heatmaps to see the effects of mutations in each cell. The heatmaps show how each individual mutation affects cell entry. The `x`'s indicate the amino-acid identity in the parental GPC strain and light gray indicates mutations that were not measured. The options below the chart also allow you to show different metrics for how the difference between mutation effects is calculated for the line plot, as well as the floor assigned to the mutation effects in the plots. 
 - The zoom bar below the site level summary shows different regions of GPC, and can be used to zoom in on specific sites.
 - The scatter plot to the right of the heatmaps shows the effect of each mutation in each of the two cell lines being compared.

Click on the expansion box in the upper right of the plot to enlarge it for easier viewing.

<Figure caption="Interactive plot comparing the effects of mutations on entry into Vesper mice TfR1 expressing 293T cells">
    <Altair :showShadow="true" :spec-url="'htmls/compare_cell_entry_site_zoom.html'"></Altair>
</Figure>

## Numerical values of mutation effects on cell entry
The pre-filtered numerical data plotted on this page can be found [here](https://github.com/dms-vep/JUNV-Romero-GPC-DMS/blob/main/results/summaries/entry_293T-various-TfR1s.csv).