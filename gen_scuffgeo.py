"""
Generate .scuffgeo files for all group/displacement combinations.
Usage: python gen_scuffgeo.py
"""
import os

# Parameters
W = 1.0

# Group definitions
groups = {
    'a': {'H': 5, 'delta': 0.6, 'displacements': [1.0, 2.0, 4.0, 5.0, 10.0]},
    'b': {'H': 5, 'delta': 0.2, 'displacements': [1.0, 2.0, 4.0, 5.0, 10.0]},
    'c': {'H': 5, 'delta': 0.1, 'displacements': [1.0, 2.0, 4.0, 5.0, 10.0]},
    'd': {'H': 2, 'delta': 0.6, 'displacements': [0.5, 1.0, 2.0, 5.0]},
    'e': {'H': 2, 'delta': 0.2, 'displacements': [0.5, 1.0, 2.0, 5.0]},
    'f': {'H': 2, 'delta': 0.1, 'displacements': [0.5, 1.0, 2.0, 5.0]},
}

# Create output directory
os.makedirs("scuffgeo", exist_ok=True)

count = 0
for group, params in groups.items():
    for d in params['displacements']:
        filename = f"scuffgeo/{group}_d{d:.1f}.scuffgeo"
        with open(filename, 'w') as f:
            f.write(f"# scuff-em geometry file for group {group}, displaced = {d}\n")
            f.write(f"# W = {W}, H = {params['H']}, delta = {params['delta']}\n")
            f.write(f"# p = 2*(W+delta) = {2*(W + params['delta'])}\n\n")

            # Lower object (not displaced)
            f.write("OBJECT Lower\n")
            f.write(f"  MESHFILE ./mshs/{group}3-low.msh\n")
            f.write("  MATERIAL ./materials/silica-27p.txt\n")
            f.write("ENDOBJECT\n\n")

            # Upper object (displaced in +y direction)
            f.write("OBJECT Upper\n")
            f.write(f"  MESHFILE ./mshs/{group}3-up.msh\n")
            f.write("  MATERIAL ./materials/silica-27p.txt\n")
            f.write(f"  DISPLACED 0 {d} 0\n")
            f.write("ENDOBJECT\n")

        count += 1
        print(f"Generated: {filename}")

print(f"\nTotal: {count} scuffgeo files generated")
