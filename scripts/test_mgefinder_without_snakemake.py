#!/usr/bin/env python3
"""Test if we can use MGEfinder tools without the internal snakemake import."""

import subprocess
import sys

def test_individual_tools():
    """Test MGEfinder individual commands that might not need snakemake import."""
    
    tools_to_test = [
        ["mgefinder", "find", "--help"],
        ["mgefinder", "pair", "--help"], 
        ["mgefinder", "inferseq-assembly", "--help"],
        ["mgefinder", "makefasta", "--help"],
        ["mgefinder", "formatbam", "--help"]
    ]
    
    print("Testing individual MGEfinder tools...")
    
    working_tools = []
    broken_tools = []
    
    for tool_cmd in tools_to_test:
        try:
            result = subprocess.run(tool_cmd, capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                working_tools.append(" ".join(tool_cmd[:2]))
                print(f"✓ {' '.join(tool_cmd[:2])}: WORKS")
            else:
                broken_tools.append(" ".join(tool_cmd[:2]))
                print(f"✗ {' '.join(tool_cmd[:2])}: FAILED - {result.stderr.strip()[:100]}")
        except Exception as e:
            broken_tools.append(" ".join(tool_cmd[:2]))
            print(f"✗ {' '.join(tool_cmd[:2])}: ERROR - {str(e)[:100]}")
    
    print(f"\nSummary:")
    print(f"Working tools: {len(working_tools)} - {working_tools}")
    print(f"Broken tools: {len(broken_tools)} - {broken_tools}")
    
    return len(working_tools) > 0

if __name__ == "__main__":
    if test_individual_tools():
        print("\n>>> Some MGEfinder tools work! We might be able to use external Snakemake.")
        sys.exit(0)
    else:
        print("\n>>> All MGEfinder tools fail. Need to fix Snakemake compatibility.")
        sys.exit(1)