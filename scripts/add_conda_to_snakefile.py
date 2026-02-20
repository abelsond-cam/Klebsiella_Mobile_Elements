#!/usr/bin/env python3
"""Add conda environment to all MGEfinder rules in the Snakefile."""

import re
from pathlib import Path

def main():
    snakefile = Path("workflow/mgefinder.end2end.snakefile")
    
    with open(snakefile) as f:
        content = f.read()
    
    # Find all rules that use mgefinder
    mgefinder_rules = [
        "formatbam", "find", "pair", "inferseq_assembly", "inferseq_reference", 
        "inferseq_overlap", "make_database", "inferseq_database", "clusterseq", 
        "genotype", "summarize", "makefasta"
    ]
    
    # Add conda directive to each MGEfinder rule
    for rule_name in mgefinder_rules:
        # Pattern to match rule definition
        pattern = f"rule {rule_name}:(.*?)shell:"
        match = re.search(pattern, content, re.DOTALL)
        
        if match:
            rule_content = match.group(1)
            # Check if conda directive already exists
            if "conda:" not in rule_content:
                # Add conda directive before shell
                replacement = f"rule {rule_name}:{rule_content}    conda:\n        \"../envs/mgefinder.yaml\"\n    shell:"
                content = content.replace(match.group(0), replacement)
                print(f"Added conda environment to rule {rule_name}")
            else:
                print(f"Rule {rule_name} already has conda environment")
        else:
            print(f"Could not find rule {rule_name}")
    
    # Write the modified content back
    with open(snakefile, 'w') as f:
        f.write(content)
    
    print("Updated Snakefile with conda environments")

if __name__ == "__main__":
    main()