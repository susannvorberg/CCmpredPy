import numpy as np
import operator


def select_CB_or_CA(res):
    if 'CB' in res:
        return res['CB']
    else:
        return res['CA']

select_CA = operator.itemgetter('CA')


def get_coordinates(resis, atom_selection=select_CB_or_CA):
    # transform into 75x3 array of atom coordinates
    return np.array([atom_selection(res).get_coord() for res in resis])


def distance_map(resis, atom_selection=select_CB_or_CA, squared=False):
    """Calculate a distance map from a selection of residues"""

    coords = get_coordinates(resis, atom_selection)

    # delta is 75x75x3 array of coordinate differences
    delta = coords[np.newaxis, :] - coords[:, np.newaxis]
    dist = np.sum(delta * delta, axis=2)

    if not squared:
        dist = np.sqrt(dist)

    return dist
