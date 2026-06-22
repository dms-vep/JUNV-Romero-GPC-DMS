"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""


# Configure dms-viz JSONs ---------------------------------------------------------------

# read configuration for `configure_dms_viz`
with open("data/inputs_for_dms-viz/dms_viz_config.yml") as f:
    dms_viz_config = yaml.YAML(typ="safe", pure=True).load(f)

rule configure_dms_viz:
    """Configure a JSON for `dms-viz`."""
    input:
        data_csv=lambda wc: dms_viz_config[wc.viz_name]["data_csv"],
        sitemap_csv=lambda wc: dms_viz_config[wc.viz_name]["sitemap_csv"],
        nb="analysis_notebooks/configure_dms_viz.ipynb",
    output:
        dms_viz_json="results/dms-viz/{viz_name}/{viz_name}.json",
        pdb_file="results/dms-viz/{viz_name}/{viz_name}.pdb",
        input_data_csv="results/dms-viz/{viz_name}/{viz_name}_data.csv",
        input_sitemap_csv="results/dms-viz/{viz_name}/{viz_name}_sitemap.csv",
        nb="results/notebooks/configure_dms_viz_{viz_name}.ipynb",
    params:
        params_yaml=lambda wc: yaml_str(
            {
                key: dms_viz_config[wc.viz_name][key]
                for key in [
                    "pdb_id",
                    "pdb_type",
                    "name",
                    "melt_condition_metric_cols",
                    "metric",
                    "opt_params",
                ]
            }
        ),
    conda:
        "conda_envs/dms-viz.yml"
    log:
        "results/logs/configure_dms_viz_{viz_name}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p data_csv {input.data_csv} \
            -p sitemap_csv {input.sitemap_csv} \
            -p dms_viz_json {output.dms_viz_json} \
            -p pdb_file {output.pdb_file} \
            -p input_data_csv {output.input_data_csv} \
            -p input_sitemap_csv {output.input_sitemap_csv} \
            -y "{params.params_yaml}" \
            &> {log}
        """

docs["dms-viz visualizations"] = {
    "dms-viz JSON files": {
        viz_name: rules.configure_dms_viz.output.dms_viz_json.format(viz_name=viz_name)
        for viz_name in dms_viz_config
    },
    "Notebooks prepping dms-viz JSONs": {
        viz_name: rules.configure_dms_viz.output.nb.format(viz_name=viz_name)
        for viz_name in dms_viz_config
    },
}

rule visualize_mutation_distributions:
    """
    Create altair version of amino-acid count
    distribution plots.
    """
    input:
        variant_data = "results/variants/codon_variants.csv",
        nb="analysis_notebooks/visualize_mutation_distributions.ipynb",
    params:
        out_dir="results/summary_of_libraries/",
    output:
        out_file="results/summary_of_libraries/present_and_absent_mutations_in_DMS_libraries.csv",
        nb="results/notebooks/visualize_mutation_distributions.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/visualize_mutation_distributions.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p variant_data {input.variant_data} \
            -p out_dir {params.out_dir} \
            -p out_file {output.out_file} \
            &> {log}
        """

rule compare_cell_entry:
    """Compare entry across cell lines."""
    input:
        nb="analysis_notebooks/compare_cell_entry.ipynb",
        mut_effects_csv="results/summaries/entry_293T-various-TfR1s.csv",
    output:
        nb="results/notebooks/compare_cell_entry.ipynb",
        site_diffs_csv="results/compare_cell_entry/site_diffs.csv",
        mut_scatter_chart="results/compare_cell_entry/compare_cell_entry_scatter.html",
        site_zoom_chart="results/compare_cell_entry/compare_cell_entry_site_zoom.html",
    params:
        params_yaml=lambda wc: yaml_str(
            {
                # cells and their names in input CSV file
                "cells": {
                    "C callosus TfR1": "C callosus TfR1",
                    "C musculinus TfR1": "C musculinus TfR1",
                },
                # for calculating differences and display, floor mutation effects at this
                "floor_mut_effects": -6,
            }
        ),
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/compare_cell_entry.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p mut_effects_csv {input.mut_effects_csv} \
            -p site_diffs_csv {output.site_diffs_csv} \
            -p mut_scatter_chart {output.mut_scatter_chart} \
            -p site_zoom_chart {output.site_zoom_chart} \
            -y "{params.params_yaml}" \
            &> {log}
        """

# read configuration for wrapped heatmaps
with open("data/wrapped_heatmap_config.yaml") as f:
    wrapped_heatmap_config = yaml.YAML(typ="safe", pure=True).load(f)

rule wrapped_heatmap:
    """Make row-wrapped heatmaps."""
    input:
        data_csv=lambda wc: wrapped_heatmap_config[wc.wrapped_hm]["data_csv"],
    output:
        chart_html="results/wrapped_heatmaps/{wrapped_hm}_wrapped_heatmap.html",
    params:
        params_dict=lambda wc: wrapped_heatmap_config[wc.wrapped_hm],
    log:
        notebook="results/notebooks/wrapped_heatmap_{wrapped_hm}.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    notebook:
        "analysis_notebooks/wrapped_heatmap.ipynb"


rule map_scores_onto_pdb:
    """Map scores onto PDB structure."""
    input:
        # input files
        dms_scores="results/summaries/cell_entry_and_pH_stability.csv",
        pdb_file="data/9MT6.pdb",
        nb="analysis_notebooks/pdb_mapping.ipynb",
    output:
        nb="results/notebooks/pdb_mapping.ipynb",
    params:
        pdb_dir="results/mapped_scores_onto_pdb/",
    log:
        "results/logs/map_scores_onto_pdb.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p dms_scores {input.dms_scores} \
            -p pdb_file {input.pdb_file} \
            -p pdb_dir {params.pdb_dir} \
            &> {log}
        """
    


docs["Additional analyses and data files"] = {
    "Variant mutation analysis" : {
        "Notebook creating altair version of mutation distribution" : rules.visualize_mutation_distributions.output.nb,
    },
    "Compare cell entry": {
        "Scatter plots of mutation effects in different cells": rules.compare_cell_entry.output.mut_scatter_chart,
        "Zoomable site chart of differences": rules.compare_cell_entry.output.site_zoom_chart,
    },
    "Row-wrapped heatmaps" : {
        wrapped_hm: rules.wrapped_heatmap.output.chart_html.format(wrapped_hm=wrapped_hm)
        for wrapped_hm in wrapped_heatmap_config
    },
    "Map DMS scores onto PDB structure":{
        "Notebook mapping scores onto pdb structure" : rules.map_scores_onto_pdb.output.nb,
    },
}
