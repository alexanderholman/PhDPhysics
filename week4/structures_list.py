from classes.Lattice import Lattice

# Obj tpl
# {
#     "comment": "",
#     "lattice": Lattice(
#         a1=[0.0, 0.0, 0.0],
#         a2=[0.0, 0.0, 0.0],
#         a3=[0.0, 0.0, 0.0]
#     ),
#     "species": [
#         {
#             "name": "Si",
#             "positions": [
#                 [0.0, 0.0, 0.0],
#             ]
#         }
#     ]
# },

structures = [
    {
        "comment": "diy-bulk-si", # not really diy
        "species": [
            {
                "name": "Si",
                "positions": [
                    [0.0, 0.0, 0.0],
                    [0.5, 0.5, 0.0],
                    [0.5, 0.0, 0.5],
                    [0.0, 0.5, 0.5],
                    [0.25, 0.25, 0.25],
                    [0.75, 0.75, 0.25],
                    [0.75, 0.25, 0.75],
                    [0.25, 0.75, 0.75]
                ]
            }
        ]
    },
    # # Structures from Materials Project
    # # source https://next-gen.materialsproject.org/materials?_limit=75&chemsys=Si
    # {
    #     "comment": "mp-149",
    #     "lattice": Lattice(
    #         a1=[5.4437023729394527, 0.0, 0.0000000000000003],
    #         a2=[0.0000000000000009, 5.4437023729394527, 0.0000000000000003],
    #         a3=[0.0, 0.0, 5.4437023729394527]
    #     ),
    #     "species": [
    #         {
    #             "name": "Si",
    #             "positions": [
    #                 [0.75, 0.75, 0.25],
    #                 [0.0, 0.5, 0.5],
    #                 [0.75, 0.25, 0.75],
    #                 [0.0, 0.0, 0.0],
    #                 [0.25, 0.75, 0.75],
    #                 [0.5, 0.5, 0.0],
    #                 [0.25, 0.25, 0.25],
    #                 [0.5, 0.0, 0.5],
    #             ]
    #         }
    #     ]
    # },
]