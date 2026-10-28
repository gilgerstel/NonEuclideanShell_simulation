from generate_inputs import create_input_files, make_in_file
import os

# ------------------------------------------------------------
# Simple flat-sheet baseline test
# ------------------------------------------------------------

sim_name = "flat_baseline"
output_dir = "baseline_output"

os.makedirs(output_dir, exist_ok=True)

# Material-coordinate domain:
# 0 <= u <= 1
# 0 <= v <= 1
domain_params = (0.0, 1.0, 0.0, 1.0)

# Create triangular mesh.
#
# max_area controls mesh resolution.
# 0.02 gives a fairly small mesh, good for a first test.
#
# tols=(0.05, 0) fixes thin strips near u=0 and u=1.
create_input_files(
    sim_name=sim_name,
    output_dir=output_dir,
    domain_params=domain_params,
    max_area=0.02,
    tols=(0.05, 0.0)
)

vertices_file = os.path.join(
    output_dir,
    f"{sim_name}_Vertices"
)

faces_file = os.path.join(
    output_dir,
    f"{sim_name}_Faces"
)

# ------------------------------------------------------------
# Parameters passed to the C++ shell simulator
# ------------------------------------------------------------

params = {
    "vertices_file": vertices_file,
    "faces_file": faces_file,

    # Keep this small for our first test.
    "loops": 50,
    "save_every": 10,

    # Reference first fundamental form:
    # abar = identity
    "abar": [
        "1", "0",
        "0", "1"
    ],

    # Reference second fundamental form:
    # bbar = 0
    "bbar": [
        "0", "0",
        "0", "0"
    ],

    "cbar": [
    "0", "0",
    "0", "0"
    ],

    "c_prefactor": "0",
    "K1111": "0",
    "K1112": "0",
    "K1122": "0",
    "K1212": "0",
    "K1222": "0",
    "K2222": "0",

    # Material parameters
    "thickness": "0.01",
    "E": "1.0",
    "nu": "0.3",

    # Initial embedding:
    # r(u,v) = (u,v,0)
    "pos0": [
        "u",
        "v",
        "0"
    ],

    # Turn connection energy completely off
    "eta": "0",
    "lambdaG": "0",
    "muG": "0",

    # Reference connection.
    # Irrelevant because eta = 0, but the input format
    # still expects all eight components.
    "gamma": [
        [
            ["0", "0"],
            ["0", "0"]
        ],
        [
            ["0", "0"],
            ["0", "0"]
        ]
    ],

    # No continuation tricks for the baseline test
    "thickness_adjust": "1",
    "metric_adjust": "1",
    "connection_adjust": "1",

    # Start fresh, not from a saved state
    "restart": 0,
}

input_file = make_in_file(
    params=params,
    output_dir=output_dir,
    sim_name=sim_name
)

print()
print("Created baseline simulation.")
print("Vertices:", vertices_file)
print("Faces:   ", faces_file)
print("Input:   ", input_file)
print()
print("Run with:")
print(f"./RunShell < {input_file}")