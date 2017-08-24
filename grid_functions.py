# -*- coding: utf-8 -*-
"""
Created on Jul 13 16:30 2017

@author: wpreimes
"""
import ast
import numpy as np

def cells_for_continent(continents):
    # type: (str) -> list
    # continents in file: "Australia", "North_America"
    cells = []
    if isinstance(continents, str):
        continents = [continents]
    for continent in continents:
        with open(r"H:\HomogeneityTesting_data\continents_cells.txt", 'r') as f:
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