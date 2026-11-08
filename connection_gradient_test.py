from generate_inputs import create_input_files, make_in_file
import os

# ------------------------------------------------------------
# Exact cylinder geometry test
# ------------------------------------------------------------

sim_name = "connection_gradient_test"
output_dir = "optimizer_output"

os.makedirs(output_dir, exist_ok=True)

# We use R = 1.
#
# Exact embedding:
#
#   r(u,v) = ( sin(u), v, cos(u) )
#
# so theoretically:
#
#   a = [[1,0],
#        [0,1]]
#
# and the nonzero principal curvature has magnitude 1.

domain_params = (-0.75, 0.75, -0.5, 0.5)

create_input_files(
    sim_name=sim_name,
    output_dir=output_dir,
    domain_params=domain_params,
    max_area=0.05,
    tols=(0.0, 0.0)
)

vertices_file = os.path.join(
    output_dir,
    f"{sim_name}_Vertices"
)

faces_file = os.path.join(
    output_dir,
    f"{sim_name}_Faces"
)

params = {
    "vertices_file": vertices_file,
    "faces_file": faces_file,

    # IMPORTANT:
    # zero optimization iterations.
    # We want to measure the supplied cylinder,
    # not let the program deform it yet.
    "loops": 0,
    "save_every": 1,

   "abar": [
        "(1+0.3*u)^2", "0",
        "0", "1"
    ],

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

    "thickness": "0.01",
    "E": "1.0",
    "nu": "0.3",

    "pos0": [
        "u+0.15*u*u",
        "v",
        "0"
    ],

    "eta": "1",
    "lambdaG": "1",
    "muG": "1",

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

    "thickness_adjust": "1",
    "metric_adjust": "1",
    "connection_adjust": "1",

    "restart": 0,
}

input_file = make_in_file(
    params=params,
    output_dir=output_dir,
    sim_name=sim_name
)

print()
print("Created cylinder geometry test.")
print("Vertices:", vertices_file)
print("Faces:   ", faces_file)
print("Input:   ", input_file)
print()
print("Run with:")
print(f"./RunShell < {input_file}")