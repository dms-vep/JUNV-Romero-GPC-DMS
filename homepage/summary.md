---
aside: false
---

# Summary of how mutations affect both cell entry and low pH stability

## Interactive plot of both cell entry and stability
The plot below shows how mutations affect stability and cell entry, and allows you to zoom and mouseover sites and mutations. Positive cell entry values indicate better entry than the parental GPC while negative values indicate impaired entry. Positive escape values indicate increased stability while negative values indicate decreased stability. Note that the two different shades of gray in the heatmaps have differing meanings: light gray means a mutation was *missing (not measured)* in the library, whereas dark gray means a mutation *was measured but was so deleterious for cell entry* it is not possible to reliably estimate its effect on stability (the threshold for how deleterious a mutation must be for cell entry to be shown in dark gray is controlled by the cell entry slider at the bottom of the plot).

Here is an explanation of the key plot elements:
 - The zoom bar at the top of the plot shows different regions of GPC, and can be used to zoom in on specific sites.
  - The heatmaps show how each individual mutation affects each phenotype. The `x`'s indicate the amino-acid identity in the parental GPC strain, light gray indicates mutations that were not measured, and dark gray indicates mutations that are filtered out by one of the interactive sliders at the bottom of the plot (eg, they have a very negative effect on cell entry). You can mouse over any of the points for details.
  - The options at the bottom of the plot allow you to interactively adjust what is displayed. For instance, the *minimum mutation entry in C musculinus TfR1 cells* only shows mutations with at least some minimal cell entry function (and grays out mutations that are more deleterious). You can also select to floor the escape at zero to show / not show "negative" escape values.

Click on the expansion box in the upper right of the plot to enlarge it for easier viewing, or [click here](/htmls/cell_entry_and_pH_stability_overlaid.html){target="_self"} to open the plot in a stand-alone window.

<Figure caption="Interactive plot showing effects of mutations on all phenotypes">
    <Altair :showShadow="true" :spec-url="'htmls/cell_entry_and_pH_stability_overlaid.html'"></Altair>
</Figure>


