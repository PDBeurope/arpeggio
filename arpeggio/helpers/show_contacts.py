#!/usr/bin/env python

import argparse
import itertools
import json
import os
from collections import namedtuple

pymol_config = {
    "dashcolor": {
        "normal": "cyan",
        "disordered_region": "deepsalmon",
        "xbond": {"vdw_clash": "blue", "vdw": "tv_blue", "proximal": "slate"},
        "hbond": {"vdw_clash": "red", "vdw": "deepsalmon", "proximal": "salmon"},
        "polar": {"vdw_clash": "red", "vdw": "deepsalmon", "proximal": "salmon"},
        "weak_hbond": {"vdw_clash": "orange", "vdw": "amber", "proximal": "peach"},
        "weak_polar": {"vdw_clash": "orange", "vdw": "amber", "proximal": "peach"},
        "metalcomplex": {
            "covalent": "purple",
            "vdw_clash": "purple",
            "vdw": "medium_purple",
        },
        "ionic": {
            "covalent": "yellow",
            "vdw_clash": "yellow",
            "vdw": "maize",
            "proximal": "cream",
        },
        "aromatic": {
            "vdw_clash": "cyan",
            "vdw": "electric_blue",
            "proximal": "baby_blue",
        },
        "hydrophobic": {
            "vdw_clash": "green",
            "vdw": "emerald",
            "proximal": "moss_green",
        },
        "carbonyl": {
            "vdw_clash": "rose",
            "vdw": "brilliant_rose",
            "proximal": "thulian_pink",
        },
        "undefined": {
            "covalent": "white",
            "vdw_clash": "grey90",
            "vdw": "grey60",
            "proximal": "grey30",
        },
    },
    "dashradius": {
        "normal": 0.02,
        "xbond": {"vdw_clash": 0.12, "vdw": 0.08, "proximal": 0.04},
        "hbond": {"vdw_clash": 0.12, "vdw": 0.08, "proximal": 0.04},
        "weak_hbond": {"vdw_clash": 0.12, "vdw": 0.08, "proximal": 0.04},
        "polar": {"vdw_clash": 0.12, "vdw": 0.08, "proximal": 0.04},
        "weak_polar": {"vdw_clash": 0.12, "vdw": 0.08, "proximal": 0.04},
        "metalcomplex": {"covalent": 0.16, "vdw_clash": 0.12, "vdw": 0.08},
        "ionic": {"covalent": 0.16, "vdw_clash": 0.12, "vdw": 0.08, "proximal": 0.04},
        "aromatic": {"vdw_clash": 0.12, "vdw": 0.08, "proximal": 0.04},
        "hydrophobic": {"vdw_clash": 0.12, "vdw": 0.08, "proximal": 0.04},
        "carbonyl": {"vdw_clash": 0.12, "vdw": 0.08, "proximal": 0.04},
        "undefined": {
            "covalent": 0.16,
            "vdw_clash": 0.12,
            "vdw": 0.08,
            "proximal": 0.04,
        },
    },
    "dashgap": {
        "normal": 0.0,
        "xbond": {"vdw_clash": 0.25, "vdw": 0.35, "proximal": 0.45},
        "hbond": {"vdw_clash": 0.25, "vdw": 0.35, "proximal": 0.45},
        "weak_hbond": {"vdw_clash": 0.25, "vdw": 0.35, "proximal": 0.45},
        "polar": {"vdw_clash": 0.25, "vdw": 0.35, "proximal": 0.45},
        "weak_polar": {"vdw_clash": 0.25, "vdw": 0.35, "proximal": 0.45},
        "metalcomplex": {"covalent": 0.20, "vdw_clash": 0.25, "vdw": 0.35},
        "ionic": {"covalent": 0.20, "vdw_clash": 0.25, "vdw": 0.35, "proximal": 0.45},
        "aromatic": {"vdw_clash": 0.25, "vdw": 0.35, "proximal": 0.45},
        "hydrophobic": {"vdw_clash": 0.25, "vdw": 0.35, "proximal": 0.45},
        "carbonyl": {"vdw_clash": 0.25, "vdw": 0.35, "proximal": 0.45},
        "undefined": {
            "covalent": 0.0,
            "vdw_clash": 0.25,
            "vdw": 0.35,
            "proximal": 0.45,
        },
    },
    "dashlength": {
        "normal": 1.0,
        "xbond": {"vdw_clash": 0.08, "vdw": 0.06, "proximal": 0.04},
        "hbond": {"vdw_clash": 0.08, "vdw": 0.06, "proximal": 0.04},
        "weak_hbond": {"vdw_clash": 0.08, "vdw": 0.06, "proximal": 0.04},
        "polar": {"vdw_clash": 0.08, "vdw": 0.06, "proximal": 0.04},
        "weak_polar": {"vdw_clash": 0.08, "vdw": 0.06, "proximal": 0.04},
        "metalcomplex": {"covalent": 0.10, "vdw_clash": 0.08, "vdw": 0.06},
        "ionic": {"covalent": 0.10, "vdw_clash": 0.08, "vdw": 0.06, "proximal": 0.04},
        "aromatic": {"vdw_clash": 0.08, "vdw": 0.06, "proximal": 0.04},
        "hydrophobic": {"vdw_clash": 0.08, "vdw": 0.06, "proximal": 0.04},
        "carbonyl": {"vdw_clash": 0.08, "vdw": 0.06, "proximal": 0.04},
        "undefined": {
            "covalent": 1.0,
            "vdw_clash": 0.08,
            "vdw": 0.06,
            "proximal": 0.04,
        },
    },
    "colors": {"carbon": "grey", "non-std-res": "green", "variations": "magenta"},
}


def build_parser():
    """Parse the command-line arguments."""
    description = (
        "Script for showing the inter-atomic interactions calculated "
        "by Arpeggio in PyMOL. \nPlease note: **You must already "
        "have run Arpeggio on your structure to generate the "
        "required input files**."
    )

    parser = argparse.ArgumentParser(description=description, prog="show_contacts")

    parser.add_argument(
        "filename",
        type=str,
        help="Path to the structure file (PDB or mmCIF) to be analyzed.",
    )

    parser.add_argument(
        "-op",
        type=str,
        default="",
        help="The output postfix you used with Arpeggio (if any).",
    )

    parser.add_argument(
        "-s",
        "--selection",
        action="append",
        nargs="+",
        type=str,
        default=None,
        help=(
            "Only report contacts that involve residues in the selection. "
            "This flag can be used multiple times and the selection syntax "
            "is '/<chain_id>/<res_num>/<atom_name>' for a single selection or "
            "'[/<chain_id>/<res_num>/<atom_name>, /<chain_id>/<res_num>/<atom_name>]' "
            "for a pairwise selection."
        ),
    )

    parser.add_argument(
        "-v", "--verbose", help="increase output verbosity", action="store_true"
    )

    return parser


def default_pymol_header():
    """Define the default PyMOL settings."""
    pymol_header = []
    add = pymol_header.append

    add("reinitialize")

    add("\n# set PyMOL variables for prettier molecules")
    add("set valence, 1")
    add("set stick_rad, 0.15")
    add("set line_width, 1")
    add("set mesh_width, 0.3")

    add("\n# decrease font size")
    add("set label_size, 10")

    add("\n# make spheres look prettier (i.e., water, ions...')")
    add("set sphere_scale, 0.15")

    add("\n# set dash properties")
    add("set dash_round_ends, 0")
    add("set dash_gap, 0.25")
    add("set dash_length, 0.05")

    add("\n# color definitions")
    add("set_color maize, (251, 236, 93)  # IONIC VDW")
    add("set_color cream, (255, 253, 208)  # IONIC PROXIMAL")
    add("set_color electric_purple, (191, 0, 255)  # METAL COMPLEX VDW CLASH")
    add("set_color medium_purple, (147, 112, 219)  METAL COMPLEX VDW")
    add("set_color amber, (255, 191, 0)")
    add("set_color peach, (255, 229, 180)")
    add("set_color electric_blue, (125, 249, 255)")
    add("set_color baby_blue, (224, 255, 255)")
    add("set_color emerald, (80, 200, 120)")
    add("set_color moss_green, (173, 200, 173)")
    add("set_color rose, (255, 0, 127)")
    add("set_color brilliant_rose, (246, 83, 166)")
    add("set_color thulian_pink, (222, 111, 161)")
    add("set_color neon_red, (255, 0, 102)")

    add("\n# settings for DNA")
    add("set cartoon_ladder_mode, 1")
    add("set cartoon_ring_finder, 1")
    add("set cartoon_ring_mode, 3")
    add("set cartoon_nucleic_acid_mode, 4")
    add("set cartoon_ring_transparency, 0.5")

    add("\n# settings for figure quality")
    add("set line_smooth, 1")
    add("set antialias, 2")
    add("set cartoon_fancy_helices, 1")
    add("set depth_cue, 1")
    add("set specular,1")
    add("set surface_quality, 1")
    add("set stick_quality, 15")
    add("set sphere_quality, 2")
    add("set cartoon_sampling, 14")
    add("set ribbon_sampling, 10")
    add("set transparency_mode, 2")
    add("set stick_ball, 1")
    add("set stick_ball_ratio, 1.5")
    add("rebuild")

    add("\n# defer updating/redrawing until all interactions are added")
    add("set defer_updates, 1")

    return "\n".join(pymol_header)


def get_interaction_information(interaction):
    """Extract basic information about the interactions."""
    bgn = interaction["bgn"]
    atom_bgn = f"{bgn['auth_asym_id']}/{bgn['auth_seq_id']}/{bgn['auth_atom_id']}"
    end = interaction["end"]
    atom_end = f"{end['auth_asym_id']}/{end['auth_seq_id']}/{end['auth_atom_id']}"
    contact = interaction["contact"]
    interaction_type = interaction["type"]

    Interaction = namedtuple(
        "Interaction", ["atom_bgn", "atom_end", "contact", "interaction_type"]
    )

    return Interaction(atom_bgn, atom_end, contact, interaction_type)


def atom_atom_interaction(interaction):
    """Interaction between two atoms."""
    dist_flag = "undefined"
    interaction_type = "undefined"
    used_interaction_types = set()

    interactions = []
    mutual_exclusive = ("clash", "covalent", "vdw_clash", "vdw", "proximal")
    feature_contacts = [c for c in interaction.contact if c not in mutual_exclusive]

    # 1: mutually exclusive interactions
    dist_flag = [me for me in mutual_exclusive if me in interaction.contact][0]

    if dist_flag == "clash" or (
        dist_flag == "covalent" and "metal_complex" not in feature_contacts
    ):
        # skip clashed and covalent bonds other than metal_complexes for now
        return None, None

    # 2: feature contacts
    if feature_contacts:
        for interaction_contact in feature_contacts:
            interactions.append((interaction_contact, dist_flag))
    else:
        interactions.append(("undefined", dist_flag))

    # generate PyMOL script for the interaction(s)
    pymol_interaction = []
    add = pymol_interaction.append

    for interaction_type, flag in interactions:
        add(
            f"distance {interaction_type}-{flag}, {interaction.atom_bgn}, {interaction.atom_end}"
        )
        add(f"show sticks, byres {interaction.atom_bgn}")
        add(f"show sticks, byres {interaction.atom_end}")
        add(f"select binding_site, binding_site + byres {interaction.atom_bgn}")
        add(f"select binding_site, binding_site + byres {interaction.atom_end}")

        used_interaction_types.add((interaction_type, flag))
    add("")

    return pymol_interaction, used_interaction_types


def plane_plane_interaction(interaction_info, interaction):
    """Interaction between two aromatic rings."""
    pymol_interaction = []
    add = pymol_interaction.append

    ring_centers = [interaction[res].get("center") for res in ("bgn", "end")]
    add(
        f"# ring-ring interaction: {interaction_info.atom_bgn} - {interaction_info.atom_end}"
    )
    for ring_center in ring_centers:
        add(f"pseudoatom ring_centers, pos={ring_center}")

    if interaction["interacting_entities"] == "INTRA_NON_SELECTION":
        return None, None

    add(f"pseudoatom pt1, pos={ring_centers[0]}")
    add(f"pseudoatom pt2, pos={ring_centers[1]}")
    add(f"distance ring-{interaction_info.contact[0]}, pt1, pt2")
    add("delete pt1")
    add("delete pt2")
    add("")

    return pymol_interaction, set({("ring-ring", None)})


def atom_plane_interaction(interaction_info, interaction):
    """Interaction between atom and aromatic ring."""
    pymol_interaction = []
    add = pymol_interaction.append

    add(
        f"# atom-ring interaction: {interaction_info.atom_bgn} - {interaction_info.atom_end}"
    )
    add(f'pseudoatom pt1, pos={interaction["end"]["center"]}')

    if interaction["interacting_entities"] == "INTRA_NON_SELECTION":
        return None, None

    for itype in interaction_info.contact:
        add(f"distance {itype}, pt1, {interaction_info.atom_bgn}")
    add("delete pt1")
    add("")

    return pymol_interaction, set({("atom-ring", None)})


def group_group_interaction(interaction_info, interaction):
    """Interaction between amide/amide."""
    pymol_interaction = []
    add = pymol_interaction.append

    if interaction["interacting_entities"] == "INTRA_NON_SELECTION":
        return None, None

    pseudo_atoms = [
        f'pseudoatom pt1, pos={interaction["bgn"]["center"]}',
        f'pseudoatom pt2, pos={interaction["end"]["center"]}',
    ]

    add(
        f"# amide-amide interaction: {interaction_info.atom_bgn} - {interaction_info.atom_end}"
    )
    for pa in pseudo_atoms:
        add(f"{pa}")
    add("distance amide-amide, pt1, pt2")
    interaction_type = set({("amide-amide", None)})

    add("delete pt1")
    add("delete pt2")
    add("")

    return pymol_interaction, interaction_type


def group_plane_interaction(interaction_info, interaction):
    """Interaction between amide/aromatic ring."""
    pymol_interaction = []
    add = pymol_interaction.append

    if interaction["interacting_entities"] == "INTRA_NON_SELECTION":
        return None, None

    pseudo_atoms = [
        f'pseudoatom pt1, pos={interaction["bgn"]["center"]}',
        f'pseudoatom pt2, pos={interaction["end"]["center"]}',
    ]

    add(
        f"# amide-ring interaction: {interaction_info.atom_bgn} - {interaction_info.atom_end}"
    )
    for pa in pseudo_atoms:
        add(f"{pa}")
        add("distance amide-ring, pt1, pt2")
        interaction_type = set({("amide-ring", None)})

    add("delete pt1")
    add("delete pt2")
    add("")

    return pymol_interaction, interaction_type


def default_pymol_footer(used_interaction_types):
    """Finish up with general PyMOL settings."""
    pymol_footer = []
    add = pymol_footer.append

    add("\n# *** specify general settings for each interaction type *** #")
    for interaction_type, flag in used_interaction_types:
        if interaction_type == "ring-ring":
            add("# ring-ring interaction")
            add("hide everything, ring_centers")
            add("show spheres, ring_centers")
            add("set dash_radius, 0.25, ring-*")
            add("set dash_gap, 0.1, ring-*")
            add("color white, ring-*")
        elif interaction_type == "atom-ring":
            add("# atom-ring interaction")
            add("set dash_radius, 0.15, *PI")
            add("set dash_gap, 0.1, *PI")
            add("color white, CARBONPI")
            add("color blue, DONORPI")
            add("color green, HALOGENPI")
            add("color red, CATIONPI")
            add("color yellow, METSULPHURPI")
        elif interaction_type == "amide-ring":
            add("# amide-ring interaction")
            add("set dash_radius, 0.25, amide-ring")
            add("set dash_gap, 0.2, amide-ring")
            add("set dash_length, 0.5, amide-ring")
            add("color white, amide-ring")
        elif interaction_type == "amide-amide":
            add("# amide-amide interaction")
            add("set dash_radius, 0.25, amide-amide")
            add("set dash_gap, 0.2, amide-amide")
            add("set dash_length, 0.5, amide-amide")
            add("color blue, amide-amide")
        else:  # atom-atom
            add("# atom-atom interaction")
            label = f"{interaction_type}-{flag}"

            add(f"color {pymol_config['dashcolor'][interaction_type][flag]}, {label}")
            add(
                f"set dash_radius, {pymol_config['dashradius'][interaction_type][flag]}, {label}"
            )
            add(
                f"set dash_gap, {pymol_config['dashgap'][interaction_type][flag]}, {label}"
            )
            add(
                f"set dash_length, {pymol_config['dashlength'][interaction_type][flag]}, {label}"
            )
        add("")

    add("delete pt1")

    add("\nhide labels")
    add("util.cbaw")
    add("bg_color white")
    add("show cartoon")
    add("set cartoon_side_chain_helper, 1")
    add("hide lines")
    add("hide everything, het")
    add("show sticks, het")
    add("show spheres, het")
    add("disable undefined-proximal")

    add("\nset defer_updates, 0")

    return "\n".join(pymol_footer)


def parse_selection(selection):
    """Parse the user-provided selection string."""
    selected_residues = []

    for sel in selection:
        # single selection
        if len(sel) == 1:
            selected_residues.append(sel[0].strip("/"))
        # pairwise selection
        elif len(sel) == 2:
            selected_residues.append(
                (sel[0].strip("[").strip("/"), sel[1].strip("]").strip("/"))
            )

    return selected_residues


def interaction_involves_selected_residues(interaction_info, selected_residues):
    """Output the interaction only if it involves residues in selection."""
    show_interaction = False

    for sel in selected_residues:
        if isinstance(sel, str):  # single selection
            if sel in interaction_info.atom_bgn or sel in interaction_info.atom_end:
                show_interaction = True
        elif isinstance(sel, tuple):  # pairwise selection
            if (
                sel[0] in interaction_info.atom_bgn
                and sel[1] in interaction_info.atom_end
            ):
                show_interaction = True

    return show_interaction


def main():
    """Do all the magic."""
    parser = build_parser()
    args = parser.parse_args()

    # handle input/output files
    filename = os.path.abspath(args.filename)
    basename, ext = os.path.splitext(args.filename)

    output_postfix = args.op
    json_filename = filename.replace(f"{ext}", output_postfix + ".json")
    script_filename = filename.replace(f"{ext}", output_postfix + ".pml")

    # residue selection (optional)
    if args.selection:
        selected_residues = parse_selection(args.selection)
        if args.verbose:
            print(f"selected residues = {selected_residues}")

    # load the JSON output from Arpeggio
    with open(json_filename, "rb") as contacts_file:
        all_interactions = json.load(contacts_file)

    # go through all interactions and output PyMOL syntax to show them
    pymol_interaction_script = []
    used_interaction_types = set()

    for interaction in all_interactions:
        interaction_info = get_interaction_information(interaction)

        if args.selection:
            show_interaction = interaction_involves_selected_residues(
                interaction_info, selected_residues
            )
            if args.verbose:
                print(
                    f"interaction {interaction_info} included because it is contained in the residue selection ({selected_residues})"
                )
            if not show_interaction:
                if args.verbose:
                    print(
                        f"interaction {interaction_info} NOT included because it is not contained in the residue selection ({selected_residues})"
                    )
                continue

        if interaction_info.interaction_type == "atom-atom":
            script, int_type = atom_atom_interaction(interaction_info)
            if script and int_type:
                pymol_interaction_script.append(script)
                used_interaction_types.update((int_type))

        elif interaction_info.interaction_type == "plane-plane":
            script, int_type = plane_plane_interaction(interaction_info, interaction)
            if script and int_type:
                pymol_interaction_script.append(script)
                used_interaction_types.update((int_type))

        elif interaction_info.interaction_type == "atom-plane":
            script, int_type = atom_plane_interaction(interaction_info, interaction)
            if script and int_type:
                pymol_interaction_script.append(script)
                used_interaction_types.update((int_type))

        elif interaction_info.interaction_type == "group-group":
            script, int_type = group_group_interaction(interaction_info, interaction)
            if script and int_type:
                pymol_interaction_script.append(script)
                used_interaction_types.update((int_type))

        elif interaction_info.interaction_type == "group-plane":
            script, int_type = group_group_interaction(interaction_info, interaction)
            if script and int_type:
                pymol_interaction_script.append(script)
                used_interaction_types.update((int_type))

        else:
            print(f"WARNING: interaction {interaction} was not processed!")

    # write the PyMOL script to a file
    with open(script_filename, "w") as output_file:
        output_file.write(default_pymol_header())

        output_file.write("\n\n# load the structure\n")
        output_file.write(f"load {filename}\n")
        output_file.write("select binding_site, None\n\n")

        flattened_pymol_interaction_script = list(
            itertools.chain.from_iterable(pymol_interaction_script)
        )
        output_file.write("\n".join(flattened_pymol_interaction_script))

        output_file.write(default_pymol_footer(used_interaction_types))


if __name__ == "__main__":
    main()
