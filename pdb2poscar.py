#!/usr/bin/env python3
"""
PDB to VASP POSCAR Converter
Converts Protein Data Bank (PDB) files to VASP POSCAR format
Atoms are grouped by element type as required by VASP
"""

import numpy as np
import sys
import os
from collections import defaultdict


def read_pdb(pdb_file):
    """
    Read PDB file and extract atomic coordinates and elements
    Returns atoms grouped by element type
    """
    # Use dictionary to group atoms by element
    atoms_by_element = defaultdict(list)
    box_vectors = None

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Extract coordinates
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())

                    # Extract element symbol
                    element = line[76:78].strip()
                    if not element:
                        # If element column is empty, try to get from atom name
                        element = line[12:16].strip()[0:2].strip()
                        # Clean element name (remove numbers and spaces)
                        element = ''.join([c for c in element if not c.isdigit()]).strip()

                    # Standardize element symbols (first letter uppercase, second lowercase)
                    if len(element) > 1:
                        element = element[0].upper() + element[1].lower()
                    else:
                        element = element.upper()

                    atoms_by_element[element].append([x, y, z])

                except (ValueError, IndexError):
                    continue

            # Try to read unit cell parameters from CRYST1 record
            elif line.startswith('CRYST1'):
                try:
                    a = float(line[6:15].strip())
                    b = float(line[15:24].strip())
                    c = float(line[24:33].strip())
                    alpha = float(line[33:40].strip())
                    beta = float(line[40:47].strip())
                    gamma = float(line[47:54].strip())

                    # Convert to Cartesian vectors
                    if alpha == 90 and beta == 90 and gamma == 90:
                        box_vectors = np.array([
                            [a, 0, 0],
                            [0, b, 0],
                            [0, 0, c]
                        ])
                    else:
                        print("Warning: Non-orthogonal unit cell detected. Using approximate box.")
                        # Convert angles to radians
                        alpha_rad = np.radians(alpha)
                        beta_rad = np.radians(beta)
                        gamma_rad = np.radians(gamma)

                        # Calculate box vectors for triclinic cell
                        box_vectors = np.array([
                            [a, 0, 0],
                            [b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0],
                            [c * np.cos(beta_rad),
                             c * (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad),
                             c * np.sqrt(1 - np.cos(beta_rad) ** 2 - (
                                         (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(
                                     gamma_rad)) ** 2)]
                        ])

                except (ValueError, IndexError):
                    print("Warning: Could not read unit cell parameters from CRYST1 record")

    return atoms_by_element, box_vectors


def convert_fractional(coords, box_vectors):
    """
    Convert Cartesian coordinates to fractional coordinates
    """
    if box_vectors is None:
        return coords

    # Calculate inverse of box matrix
    inv_box = np.linalg.inv(box_vectors)

    # Convert to fractional coordinates
    fractional_coords = np.dot(coords, inv_box.T)

    return fractional_coords


def get_element_order(elements):
    """
    Define a standard order for elements
    You can customize this order based on your needs
    """
    # Common order in VASP calculations
    preferred_order = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                       'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
                       'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                       'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
                       'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
                       'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
                       'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
                       'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
                       'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']

    # Sort elements according to preferred order
    sorted_elements = []
    for element in preferred_order:
        if element in elements:
            sorted_elements.append(element)

    # Add any remaining elements not in preferred order (alphabetically)
    remaining = [e for e in elements if e not in sorted_elements]
    sorted_elements.extend(sorted(remaining))

    return sorted_elements


def write_poscar(atoms_by_element, box_vectors, output_file, system_name="Generated from PDB"):
    """
    Write VASP POSCAR file with atoms grouped by element type
    """
    # Get all elements and sort them
    all_elements = list(atoms_by_element.keys())
    sorted_elements = get_element_order(all_elements)

    # Prepare data for writing
    element_counts = []
    all_coords = []

    for element in sorted_elements:
        coords = atoms_by_element[element]
        element_counts.append(len(coords))
        all_coords.extend(coords)

    # Convert to numpy arrays
    all_coords = np.array(all_coords)

    with open(output_file, 'w') as f:
        # System name
        f.write(f"{system_name}\n")

        # Scaling factor
        f.write("1.0\n")

        # Lattice vectors
        if box_vectors is not None:
            for vector in box_vectors:
                f.write(f"  {vector[0]:16.10f}  {vector[1]:16.10f}  {vector[2]:16.10f}\n")
        else:
            # If no box info, create a bounding box with some padding
            min_coords = np.min(all_coords, axis=0)
            max_coords = np.max(all_coords, axis=0)
            padding = 10.0  # Angstrom padding

            box_vectors = np.array([
                [max_coords[0] - min_coords[0] + padding, 0, 0],
                [0, max_coords[1] - min_coords[1] + padding, 0],
                [0, 0, max_coords[2] - min_coords[2] + padding]
            ])

            for vector in box_vectors:
                f.write(f"  {vector[0]:16.10f}  {vector[1]:16.10f}  {vector[2]:16.10f}\n")

            # Center atoms in the box
            center_shift = (box_vectors[0, 0] + box_vectors[1, 1] + box_vectors[2, 2]) / 2 - (
                        max_coords + min_coords) / 2
            all_coords += center_shift

        # Element names
        f.write("  " + "  ".join(sorted_elements) + "\n")

        # Element counts
        f.write("  " + "  ".join(map(str, element_counts)) + "\n")

        # Coordinate type
        if box_vectors is not None:
            f.write("Direct\n")
            # Convert to fractional coordinates
            fractional_coords = convert_fractional(all_coords, box_vectors)

            # Write coordinates grouped by element
            start_idx = 0
            for count in element_counts:
                for i in range(start_idx, start_idx + count):
                    coord = fractional_coords[i]
                    f.write(f"  {coord[0]:16.10f}  {coord[1]:16.10f}  {coord[2]:16.10f}\n")
                start_idx += count
        else:
            f.write("Cartesian\n")
            # Write coordinates grouped by element
            start_idx = 0
            for count in element_counts:
                for i in range(start_idx, start_idx + count):
                    coord = all_coords[i]
                    f.write(f"  {coord[0]:16.10f}  {coord[1]:16.10f}  {coord[2]:16.10f}\n")
                start_idx += count


def main():
    if len(sys.argv) < 2:
        print("Usage: python pdb_to_poscar.py <pdb_file> [output_file]")
        print("Example: python pdb_to_poscar.py protein.pdb POSCAR")
        sys.exit(1)

    pdb_file = sys.argv[1]

    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        base_name = os.path.splitext(pdb_file)[0]
        output_file = f"POSCAR_{base_name}"

    if not os.path.exists(pdb_file):
        print(f"Error: PDB file '{pdb_file}' not found!")
        sys.exit(1)

    try:
        print(f"Converting {pdb_file} to {output_file}...")

        # Read PDB file
        atoms_by_element, box_vectors = read_pdb(pdb_file)

        total_atoms = sum(len(coords) for coords in atoms_by_element.values())
        if total_atoms == 0:
            print("Error: No atoms found in PDB file!")
            sys.exit(1)

        print(f"Found {total_atoms} atoms")
        print(f"Elements: {list(atoms_by_element.keys())}")
        for element, coords in atoms_by_element.items():
            print(f"  {element}: {len(coords)} atoms")

        if box_vectors is not None:
            print("Unit cell information found in PDB file")
        else:
            print("No unit cell information found - creating bounding box")

        # Write POSCAR file
        system_name = f"System from {pdb_file}"
        write_poscar(atoms_by_element, box_vectors, output_file, system_name)

        print(f"Successfully created {output_file}")

    except Exception as e:
        print(f"Error during conversion: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()


