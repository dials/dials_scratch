from dials_scratch.jbe.sys_abs.screw_axes import (
    screw_axis_21a, screw_axis_21b, screw_axis_21c,
    screw_axis_41a, screw_axis_41b, screw_axis_41c, screw_axis_42c,
    screw_axis_61c, screw_axis_62c, screw_axis_63c, screw_axis_31c,
    ScrewAxisObserver
)
from libtbx.table_utils import simple_table


def score_screw_axes(laue_group_info, reflection_table):
    """Get the relevant screw axes and score them. Print pretty table."""
    print('scoring individual screw axes')
    scores = []
    rows = []
    for axis in laue_group_info["unique_axes"]:
        a = axis()
        a.register_observer(event="selected data for scoring", observer=ScrewAxisObserver())
        if "equivalent_axes" in laue_group_info:
            if axis in laue_group_info["equivalent_axes"]:
                for equivalent in laue_group_info["equivalent_axes"][axis]:
                    a.add_equivalent_axis(equivalent())
        score = a.score_axis(reflection_table)
        rows.append([a.name, "%.5f" % score, str(a.n_refl_used[0]), str(a.n_refl_used[1])])
        scores.append(score)

    st = simple_table(rows, column_headers=["Screw axis", "Score", "No. present", "No. absent"])
    print(st.format())
    return scores

def score_space_groups(screw_axis_scores, laue_group_info):
    # Allowed sgs is a tuple of space group and axis conditions
    # Screw axis scores is a list of the scores.
    print('scoring space groups')
    allowed_sgs = laue_group_info["space_groups"]
    rows = []
    max_score = 0.0
    best_sg = None
    for space_g, conditions in allowed_sgs:
        sg_score = 1.0
        for condition, score in zip(conditions, screw_axis_scores):
            if condition:
                sg_score *= score
            else:
                sg_score *= (1.0 - score)
        rows.append([space_g, "%.4f" % sg_score])
        if sg_score > max_score:
            max_score = sg_score
            best_sg = space_g

    st = simple_table(rows, column_headers=['Space group', 'score'])
    print(st.format())
    print("Recommended space group: %s" % best_sg)
    if "enantiomorphic pairs" in laue_group_info:
        if best_sg in laue_group_info["enantiomorphic pairs"]:
            print("Space group with equivalent score (enantiomorphic pair): %s" %
                laue_group_info["enantiomorphic pairs"][best_sg])

    return best_sg

axial_zones = {
    'P 1 2/m 1': ["2(1) b"], #P2, P21
    'P m m m': ["2(1) a", "2(1) b", "2(1) c"], #P222, P2221, P21212, P212121
    'C m m m': ["2(1) c"], #C222 C2221
    'P 4/m': ["4(1) c"], # P4, P41, P42, P43
    'I 4/m': ["4(1) c"],# I4, I41, not I42
    'P 4/m m m': ["4(1) c", "2(1) a", "2(1) b"],#symmetry equivalent
        #P422, #P4212, #P4122, P4222, P4322, #P41212, #P42212, P43212
    'I 4/m m m': ["4(1) c"], #I422, I4122
    'P -3': ["3(1) c"], #P3 P31 P32
    'P -3 1 m': ["3(1) c"], #P312, P3112, P3212
    'P -3 m 1': ["3(1) c"], #P321 P3121, P3231
    'P 6/m': ["6(1) c"], #P6, P61, P62, P63, P64, P65
    'P 6/m m m': ["6(1) c"], #P622, P6122, P6222, P6322, P6422, P6522
    'P m -3': ["2(1) a", "2(1) b", "2(1) c"],#symmetry equivalent #P23, P213
    'P m -3 m': ["4(1) a", "4(1) b", "4(1) c"],#symmetry equivalent #P423, P4132, P4232, P4332
    'I m -3 m': ["4(1) a", "4(1) b", "4(1) c"],#symmetry equivalent # I432, I4132
    'F m -3 m': ["4(1) a", "4(1) b", "4(1) c"],#symmetry equivalent # F432, F4132
}

p2m = {
    "unique_axes" : [screw_axis_21b],
    "space_groups" : [
        ("P 2", [False]),
        ("P 21", [True]),
    ],
}

Cmmm = {
    "unique_axes" : [screw_axis_21c],
    "space_groups" : [
        ("C 2 2 2", [False]),
        ("C 2 2 21", [True]),#What about different settings?
    ],
}

Immm = {
    "unique_axes" : [screw_axis_21c, screw_axis_21b, screw_axis_21a],
    "space_groups" : [
        ("I 2 2 2", [False, False, False]),
        ("I 21 21 21", [True, True, True]), #is this correct?
    ],
}

pmmm = {
    "unique_axes" : [screw_axis_21a, screw_axis_21b, screw_axis_21c],
    "space_groups" : [
        ("P 2 2 2", [False, False, False]),
        ("P 2 2 21", [False, False, True]), #these three will have to be combined/reindexed?
        ("P 2 21 2", [False, True, False]),
        ("P 21 2 2", [True, False, False]),
        ("P 21 21 2", [True, True, False]),#these three will have to be combined/reindexed?
        ("P 21 2 21", [True, False, True]),
        ("P 2 21 21", [False, True, True]),
        ("P 21 21 21", [True, True, True]),
    ],
}

p6m = {
    "unique_axes" : [screw_axis_61c, screw_axis_62c, screw_axis_63c],
    "space_groups": [
        ('P 6', [False, False, False]),
        ('P 61', [True, False, False]),
        ('P 62', [False, True, False]),
        ('P 63', [False, False, True]),
    ],
    "enantiomorphic pairs" : {
        'P 61' : 'P 65',
        'P 62' : 'P 64',
    }
}

p6mmm = {
    "unique_axes" : [screw_axis_61c, screw_axis_62c, screw_axis_63c],
    "space_groups": [
        ('P 6 2 2', [False, False, False]),
        ('P 61 2 2', [True, False, False]),
        ('P 62 2 2', [False, True, False]),
        ('P 63 2 2', [False, False, True]),
    ],
    "enantiomorphic pairs" : {
        'P 61 2 2' : 'P 65 2 2',
        'P 62 2 2' : 'P 64 2 2',
    }
}

p3 = {
    "unique_axes" : [screw_axis_31c],
    "space_groups":[
        ('P 3', [False]),
        ('P 31', [True]),
    ],
    "enantiomorphic pairs" : {
        'P 31' : 'P 32',
    },
}

p31m = {
    "unique_axes" : [screw_axis_31c],
    "space_groups":[
        ('P 3 1 2', [False]),
        ('P 31 1 2', [True]),
    ],
    "enantiomorphic pairs" : {
        'P 31 1 2' : 'P 32 1 2',
    },
}

p3m1 = {
    "unique_axes" : [screw_axis_31c],
    "space_groups":[
        ('P 3 2 1', [False]),
        ('P 31 2 1', [True]),
    ],
    "enantiomorphic pairs" : {
        'P 31 2 1' : 'P 32 2 1',
    },
}

p4mmm = {
    "unique_axes" : [screw_axis_41c, screw_axis_21a, screw_axis_42c],
    "equivalent_axes" : {screw_axis_21a : [screw_axis_21b]},
    "space_groups": [# name and which unique axes expect.
        ('P 4 2 2', [False, False, False]),
        ('P 4 21 2', [False, True, False]),
        ('P 41 2 2', [True, False, False]),
        ('P 42 2 2', [False, False, True]),
        ('P 41 21 2', [True, True, False]),
        ('P 42 21 2', [False, True, True]),
    ],
    "enantiomorphic pairs" : {
        'P 41 2 2' : 'P 43 2 2',
        'P 41 21 2' : 'P 43 21 2',
    }
}

p4m = {
    "unique_axes" : [screw_axis_41c, screw_axis_42c],
    "space_groups": [
        ('P 4', [False, False]),
        ('P 41', [True, False]),
        ('P 42', [False, True]),
    ],
    "enantiomorphic pairs" : {
        'P 41' : 'P 43',
    }
}

I4m = {
    "unique_axes" : [screw_axis_41c],
    "space_groups": [
        ('I 4', [False]),
        ('I 41', [True]),
    ],
}

I4mmm = {
    "unique_axes" : [screw_axis_41c],
    "space_groups": [
        ('I 4 2 2', [False]),
        ('I 41 2 2', [True]),
    ],
}

pm3 = {
    "unique_axes" : [screw_axis_21a],
    "equivalent_axes" : {screw_axis_21a : [screw_axis_21b, screw_axis_21c]},
    "space_groups" : [
        ("P 2 3", [False]),
        ("P 21 3", [True]),
    ]
}

Im3 = {
    "unique_axes" : [screw_axis_21a],
    "equivalent_axes" : {screw_axis_21a : [screw_axis_21b, screw_axis_21c]},
    "space_groups" : [
        ("I 2 3", [False]),
        ("I 21 3", [True]),
    ]
}

pm3m = {
    "unique_axes" : [screw_axis_41a, screw_axis_42c],
    "equivalent_axes" : {screw_axis_41a : [screw_axis_41b, screw_axis_41c]},
    "space_groups" : [
        ("P 4 3 2", [False, False]),
        ("P 41 3 2", [True, False]),
        ("P 42 3 2", [False, True]),
    ],
    "enantiomorphic pairs" : {
        'P 41 3 2' : 'P 43 3 2',
    }
}

Im3m = {
    "unique_axes" : [screw_axis_41a],
    "equivalent_axes" : {screw_axis_41a : [screw_axis_41b, screw_axis_41c]},
    "space_groups" : [
        ("I 4 3 2", [False]),
        ("I 41 3 2", [True]),
    ]
}

Fm3m = {
    "unique_axes" : [screw_axis_41a],
    "equivalent_axes" : {screw_axis_41a : [screw_axis_41b, screw_axis_41c]},
    "space_groups" : [
        ("F 4 3 2", [False]),
        ("F 41 3 2", [True]),
    ]
}


laue_groups = {
    'P 1 2/m 1' : p2m,
    'P m m m ' : pmmm,
    'C m m m' : Cmmm,
    'I m m m' : Immm,
    'P 4/m' : p4m,
    'I 4/m' : I4m,
    'P 4/m m m' : p4mmm,
    'I 4/m m m' : I4mmm,
    'P 6/m' : p6m,
    'P 6/m m m' : p6mmm,
    'P -3': p3,
    'P -3 1 m': p31m,
    'P -3 m 1': p3m1,
    'P m -3': pm3,
    'I m -3': Im3,
    'P m -3 m': pm3m,
    'I m -3 m': Im3m,
    'F m -3 m': Fm3m,
}
