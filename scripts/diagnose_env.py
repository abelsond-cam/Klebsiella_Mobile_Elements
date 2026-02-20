#!/usr/bin/env python3
"""Diagnose environment issues."""

import sys
import subprocess

def check_package(name):
    try:
        # Python 3.9 compatible version - avoid nested quotes in f-strings
        cmd = f"import {name}; print('{name}: ' + str(getattr({name}, '__version__', 'unknown version')))"
        result = subprocess.run([sys.executable, "-c", cmd], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print(result.stdout.strip())
        else:
            print(f"{name}: FAILED - {result.stderr.strip()}")
    except Exception as e:
        print(f"{name}: ERROR - {e}")

def check_pulp_detailed():
    try:
        import pulp
        version = str(getattr(pulp, '__version__', 'unknown'))
        print("pulp version: " + version)
        
        solver_attrs = [attr for attr in dir(pulp) if 'solver' in attr.lower()]
        print("pulp attributes: " + str(solver_attrs))
        
        # Try the failing call
        try:
            if hasattr(pulp, 'list_solvers'):
                solvers = pulp.list_solvers(onlyAvailable=True)
                print("pulp.list_solvers() works: " + str(solvers))
            elif hasattr(pulp, 'listSolvers'):
                solvers = pulp.listSolvers(onlyAvailable=True)  
                print("pulp.listSolvers() works: " + str(solvers))
            else:
                print("No list_solvers or listSolvers function found")
        except Exception as e:
            print("pulp solver listing FAILED: " + str(e))
            
    except ImportError as e:
        print("Cannot import pulp: " + str(e))

def main():
    print("Python: " + sys.version)
    print("Environment: " + sys.executable)
    print("")
    
    packages = ["snakemake", "pulp", "yaml", "pandas"]
    for pkg in packages:
        check_package(pkg)
    
    print("\nDetailed pulp check:")
    check_pulp_detailed()
    
    print("\nTrying mgefinder:")
    try:
        result = subprocess.run(["mgefinder", "--version"], capture_output=True, text=True)
        if result.returncode == 0:
            print("mgefinder: " + result.stdout.strip())
        else:
            print("mgefinder FAILED: " + result.stderr.strip())
    except Exception as e:
        print("mgefinder ERROR: " + str(e))

if __name__ == "__main__":
    main()