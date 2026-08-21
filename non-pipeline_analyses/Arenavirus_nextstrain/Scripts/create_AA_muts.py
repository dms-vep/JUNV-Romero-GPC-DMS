import json
import argparse


def create_AA_muts(input_file, output_file):
    GENE_NAME = "GPC"  # must match what you want shown in Auspice

    with open(input_file) as f:
        data = json.load(f)

    new_nodes = {}
    for node_name, node_data in data["nodes"].items():
        aa_mut_list = []
        for mut in node_data.get("muts", []):
            # muts look like "A42V" — from, position, to
            from_aa = mut[0]
            to_aa = mut[-1]
            pos = int(mut[1:-1])
            aa_mut_list.append(f"{from_aa}{pos}{to_aa}")
        
        new_nodes[node_name] = {
            "aa_muts": {
                GENE_NAME: aa_mut_list
            }
        }

    output = {"nodes": new_nodes}
    with open(output_file, "w") as f:
        json.dump(output, f, indent=2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str, required=True, help="input")
    parser.add_argument('--output', type=str, required=True, help="output")
    args = parser.parse_args()

    create_AA_muts(args.input, args.output)