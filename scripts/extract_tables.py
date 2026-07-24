"""Extract MMFFANG, MMFFDef, MMFFStbn, MMFFDfsb tables from RDKit source."""

import re

# Read the Params.cpp file we already have
with open(
    "/Users/rliu/.local/share/opencode/tool-output/tool_de0f1f0d0001of23VMexxz9iwE", "r"
) as f:
    content = f.read()

# The file is a single line JSON, so search within it

# Extract defaultMMFFAngle table
angle_match = re.search(r'defaultMMFFAngle\s*=\s*"(.*?)";', content, re.DOTALL)
if angle_match:
    table = angle_match.group(1)
    # Unescape
    table = table.replace("\\n", "\n").replace("\\t", "\t").replace('\\"', '"')
    with open(
        "/Users/rliu/Projects/WebMM/mmff_params_extracted/defaultMMFFAngle.tsv", "w"
    ) as f:
        f.write(table)
    lines = [l for l in table.split("\n") if l.strip() and not l.startswith("*")]
    print(f"MMFFAngle: {len(lines)} data lines")
else:
    print("Could not find defaultMMFFAngle")

# Extract defaultMMFFDef table
def_match = re.search(r'defaultMMFFDef\s*=\s*"(.*?)";', content, re.DOTALL)
if def_match:
    table = def_match.group(1)
    table = table.replace("\\n", "\n").replace("\\t", "\t").replace('\\"', '"')
    with open(
        "/Users/rliu/Projects/WebMM/mmff_params_extracted/defaultMMFFDef.tsv", "w"
    ) as f:
        f.write(table)
    lines = [l for l in table.split("\n") if l.strip() and not l.startswith("*")]
    print(f"MMFFDef: {len(lines)} data lines")
else:
    print("Could not find defaultMMFFDef")

# Extract defaultMMFFStbn table (stretch-bend)
stbn_match = re.search(r'defaultMMFFStbn\s*=\s*"(.*?)";', content, re.DOTALL)
if stbn_match:
    table = stbn_match.group(1)
    table = table.replace("\\n", "\n").replace("\\t", "\t").replace('\\"', '"')
    with open(
        "/Users/rliu/Projects/WebMM/mmff_params_extracted/defaultMMFFStbn.tsv", "w"
    ) as f:
        f.write(table)
    lines = [l for l in table.split("\n") if l.strip() and not l.startswith("*")]
    print(f"MMFFStbn: {len(lines)} data lines")
else:
    print("Could not find defaultMMFFStbn")

# Extract defaultMMFFDfsb table (default stretch-bend)
dfsb_match = re.search(r'defaultMMFFDfsb\s*=\s*"(.*?)";', content, re.DOTALL)
if dfsb_match:
    table = dfsb_match.group(1)
    table = table.replace("\\n", "\n").replace("\\t", "\t").replace('\\"', '"')
    with open(
        "/Users/rliu/Projects/WebMM/mmff_params_extracted/defaultMMFFDfsb.tsv", "w"
    ) as f:
        f.write(table)
    lines = [l for l in table.split("\n") if l.strip() and not l.startswith("*")]
    print(f"MMFFDfsb: {len(lines)} data lines")
else:
    print("Could not find defaultMMFFDfsb")

# Extract defaultMMFFProp table (has crd, val, mltb, arom, sbmb etc.)
prop_match = re.search(r'defaultMMFFProp\s*=\s*"(.*?)";', content, re.DOTALL)
if prop_match:
    table = prop_match.group(1)
    table = table.replace("\\n", "\n").replace("\\t", "\t").replace('\\"', '"')
    with open(
        "/Users/rliu/Projects/WebMM/mmff_params_extracted/defaultMMFFProp.tsv", "w"
    ) as f:
        f.write(table)
    lines = [l for l in table.split("\n") if l.strip() and not l.startswith("*")]
    print(f"MMFFProp: {len(lines)} data lines")
else:
    print("Could not find defaultMMFFProp")
