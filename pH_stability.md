---
aside: false
---

# Effects of mutations on low pH stability 

## Interactive plot of mutation effects
The plot below shows how mutations affect stability, and allows you to zoom and mouseover sites and mutations. Positive values indicate better stability than the parental GPC while negative values indicate impaired stability. 

Here is an explanation of the key plot elements:
 - The zoom bar at the top of the plot shows different regions of GPC, and can be used to zoom in on specific sites.
 - The line plot summarizes the effects of mutations at each site, using the site smmary statistic specified by the interactive option at the bottom of the plot (eg, sum or mean effect of mutations at a site).
  - The heatmap shows how each individual mutation affects stability. The `x`'s indicate the amino-acid identity in the parental GPC strain, light gray indicates mutations that were not measured, and dark gray indicates mutations that are filtered out by one of the interactive sliders at the bottom of the plot. If you mouseover points, you can also see the measurement in each individual library, which provides a sense of the reproducibility between experiments. 
  - The options at the bottom of the plot allow you to adjust what is displayed. For instance, the *times_seen* slider filters for only mutations seen in at least that many unique variants per-library; larger values will generally give more accurate measurements. Similarly, you can filter by how many different experimental selections the mutation was measured in.

Click on the expansion box in the upper right of the plot to enlarge it for easier viewing.

<Figure caption="Interactive plot showing effects of mutations on low pH stability">
    <Altair :showShadow="true" :spec-url="'htmls/stability_mut_effect.html'"></Altair>
</Figure>

## Interactive plot showing both cell entry and stability measurements

Click on the expansion box in the upper right of the plot to enlarge it for easier viewing.

<Figure caption="Interactive plot showing effects of mutations on all phenotypes">
    <Altair :showShadow="true" :spec-url="'htmls/cell_entry_and_pH_stability_overlaid.html'"></Altair>
</Figure>

## Numerical values of mutation effects on stability
The pre-filtered numerical data plotted on this page can be found [here](https://github.com/dms-vep/JUNV-Romero-GPC-DMS/blob/main/results/summaries/cell_entry_and_pH_stability.csv). Note that this file also includes the cell entry data, which can be useful for interpreting the stability measurements. 