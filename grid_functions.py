# -*- coding: utf-8 -*-
"""
Created on Jul 13 16:30 2017

@author: wpreimes
"""
import ast
import numpy as np
import os

def cells_for_continent(continents):
    # type: (str) -> list
    # continents in file: "Australia", "North_America"
    continents_cells_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'continents_cells.txt')
    cells = []
    if isinstance(continents, str):
        continents = [continents]
    for continent in continents:
        with open(continents_cells_file, 'r') as f:
            s = f.read()
            cells += ast.literal_eval(s)[continent]

    return cells

def grid_points_for_cells(grid, cells):

    if type(cells) == str:
        cells = cells_for_continent(cells)

    if type(cells) == list:
        grid_points = []
        for cell in cells:
            grid_points+=np.ndarray.tolist(grid.grid_points_for_cell(cell)[0])
        return grid_points


if __name__ == '__main__':
    cells = cells_for_continent('Australia')
    print cells